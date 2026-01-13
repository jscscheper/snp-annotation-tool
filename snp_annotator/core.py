import csv
from pathlib import Path
from typing import Dict, List, Union, Tuple, Optional
from Bio import AlignIO, SeqIO
from Bio.Align import substitution_matrices
from Bio.Seq import MutableSeq
from Bio.SeqRecord import SeqRecord


CODON_TABLE: Dict[str, str] = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
    'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W'
}

class SNPAnnotationError(Exception):
    pass

class InvalidPositionError(SNPAnnotationError):
    pass

class SNPAnnotator:
    """
    Handles SNP insertion and checks it against an MSA.
    """
    
    def __init__(self, msa_path: Union[str, Path], sequence_path: Union[str, Path]):
        self.msa_path = Path(msa_path)
        self.sequence_path = Path(sequence_path)
        self._gene_record: Union[SeqRecord, None] = None
        self._msa: Union[AlignIO.MultipleSeqAlignment, None] = None

    def load_data(self) -> None:
        """Grabs the FASTA and MSA files."""
        if not self.msa_path.exists():
            raise FileNotFoundError(f"MSA file not found: {self.msa_path}")
        if not self.sequence_path.exists():
            raise FileNotFoundError(f"Sequence file not found: {self.sequence_path}")
            
        self._gene_record = SeqIO.read(self.sequence_path, "fasta")
        self._msa = AlignIO.read(self.msa_path, format="clustal")

    def insert_snp(self, position: int, nucleotide: str) -> MutableSeq:
        """
        Puts a SNP at the given position (1-based).
        """
        if self._gene_record is None:
            self.load_data()
            
        if position < 1:
             raise InvalidPositionError(f"Position must be greater than 0, got {position}")
             
        sequence_len = len(self._gene_record.seq)
        if position > sequence_len:
            raise InvalidPositionError(
                f"Position {position} is out of range. Sequence length is {sequence_len}."
            )

        gene_seq = MutableSeq(self._gene_record.seq)
        
        index = position - 1
        original_nuc = gene_seq[index]
        

        gene_seq[index] = nucleotide
        
        return gene_seq

    def translate_sequence(self, dna_sequence: MutableSeq) -> List[str]:
        """
        Translates DNA to amino acids.
        """
        seq_str = str(dna_sequence)
        

        codon_list = [seq_str[i:i+3] for i in range(0, len(seq_str), 3)]
        

        amino_sequence = [
            CODON_TABLE[codon] if codon in CODON_TABLE else "-" 
            for codon in codon_list
        ]
        return amino_sequence

    def determine_outcome(self, amino_sequence: List[str], position: int) -> float:
        """
        Figures out the conservation percentage.
        """
        if self._msa is None:
            self.load_data()


        aa_index = (position - 1) // 3
        
        if aa_index >= len(amino_sequence):
             raise InvalidPositionError(f"Amino acid position {aa_index} out of range.")

        target_aa = amino_sequence[aa_index]
        
        match_count = 0
        total_sequences = len(self._msa)
        
        for alignment in self._msa:

            if aa_index < len(alignment.seq):
                if alignment.seq[aa_index] == target_aa:
                    match_count += 1
                    
        return round((match_count / total_sequences) * 100, 2)

    def calculate_blosum_score(self, amino_sequence: List[str], position: int, nucleotide: str) -> Dict[str, Union[float, str]]:
        """
        Gets the biological impact score (BLOSUM62).
        """
        if self._msa is None:
            self.load_data()

        aa_index = (position - 1) // 3
        if aa_index >= len(amino_sequence):
            raise InvalidPositionError(f"Amino acid position {aa_index} out of range.")
            
        mutated_aa = amino_sequence[aa_index]
        if mutated_aa == "-":
             return {"score_original": 0.0, "score_mutated": 0.0, "delta": 0.0, "verdict": "unknown"}

        orig_seq_str = str(self._gene_record.seq)
        orig_codon = orig_seq_str[(position-1)//3*3 : (position-1)//3*3+3]
        original_aa = CODON_TABLE.get(orig_codon, "-")

        if original_aa == "-":
             return {"score_original": 0.0, "score_mutated": 0.0, "delta": 0.0, "verdict": "unknown"}
        try:
            matrix = substitution_matrices.load("BLOSUM62")
        except Exception:
             
             raise ImportError("Could not load BLOSUM62 matrix.")

        score_original_sum = 0.0
        score_mutated_sum = 0.0
        count = 0
        
        for alignment in self._msa:
            if aa_index < len(alignment.seq):
                msa_aa = alignment.seq[aa_index]
                if msa_aa in matrix.alphabet and original_aa in matrix.alphabet and mutated_aa in matrix.alphabet:
                    score_original_sum += matrix[original_aa, msa_aa]
                    score_mutated_sum += matrix[mutated_aa, msa_aa]
                    count += 1
        
        if count == 0:
            return {"score_original": 0.0, "score_mutated": 0.0, "delta": 0.0, "verdict": "unknown"}
            
        avg_original = score_original_sum / count
        avg_mutated = score_mutated_sum / count
        delta = avg_original - avg_mutated
        
        verdict = "deleterious" if delta > 1.0 else "neutral"
        
        return {
            "score_original": round(avg_original, 2),
            "score_mutated": round(avg_mutated, 2),
            "delta": round(delta, 2),
            "verdict": verdict
        }

    def process_batch(self, csv_path: Union[str, Path]) -> List[Dict]:
        """
        Runs through a CSV file.
        """
        results = []
        with open(csv_path, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            reader.fieldnames = [name.strip().lower() for name in reader.fieldnames]
            
            if 'position' not in reader.fieldnames or 'nucleotide' not in reader.fieldnames:
                 raise ValueError("CSV must contain 'position' and 'nucleotide' columns.")
                 
            for row in reader:
                try:
                    pos = int(row['position'])
                    nuc = row['nucleotide'].strip().upper()
                    

                    mutated = self.insert_snp(pos, nuc)
                    amino_seq = self.translate_sequence(mutated)
                    
                    cons_score = self.determine_outcome(amino_seq, pos)
                    cons_verdict = self.get_verdict(cons_score)
                    
                    blosum_stats = self.calculate_blosum_score(amino_seq, pos, nuc)
                    
                    results.append({
                        "position": pos,
                        "nucleotide": nuc,
                        "conservation_score": cons_score,
                        "conservation_verdict": cons_verdict,
                        "blosum_delta": blosum_stats.get("delta"),
                        "blosum_verdict": blosum_stats.get("verdict")
                    })
                    
                except Exception as e:
                    results.append({"position": row.get('position'), "nucleotide": row.get('nucleotide'),  "error": str(e)})
        return results

    @staticmethod
    def get_verdict(conservation_score: float) -> str:
        """Returns 'deleterious' or 'neutral'."""
        if conservation_score < 90:
            return "deleterious"
        return "neutral"
