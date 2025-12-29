
import unittest
from unittest.mock import MagicMock, patch
from Bio.Seq import Seq, MutableSeq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from snp_annotator.core import SNPAnnotator, InvalidPositionError

class TestSNPAnnotator(unittest.TestCase):

    def setUp(self):
        self.annotator = SNPAnnotator("dummy_msa.clustal", "dummy_seq.fasta")

    def test_insert_snp_valid(self):
        # Mock loaded data
        self.annotator._gene_record = SeqRecord(Seq("ATGC"))
        
        # Insert 'T' at position 1 (A -> T)
        mutated = self.annotator.insert_snp(1, "T")
        self.assertEqual(str(mutated), "TTGC")
        
        # Insert 'G' at position 4 (C -> G)
        mutated = self.annotator.insert_snp(4, "G")
        self.assertEqual(str(mutated), "ATGG")

    def test_insert_snp_invalid_position(self):
        self.annotator._gene_record = SeqRecord(Seq("ATGC"))
        
        with self.assertRaises(InvalidPositionError):
            self.annotator.insert_snp(0, "A")
            
        with self.assertRaises(InvalidPositionError):
            self.annotator.insert_snp(5, "A")

    def test_translate_sequence(self):
        # ATG -> M, GCC -> A
        seq = MutableSeq("ATGGCC")
        aa = self.annotator.translate_sequence(seq)
        self.assertEqual(aa, ['M', 'A'])
        
        # Partial codon
        seq = MutableSeq("ATGAT") 
        aa = self.annotator.translate_sequence(seq)
        self.assertEqual(aa, ['M', '-']) 

    def test_determine_outcome(self):
        # MSA has 3 sequences
        seq1 = SeqRecord(Seq("MA"), id="seq1")
        seq2 = SeqRecord(Seq("MA"), id="seq2")
        seq3 = SeqRecord(Seq("MG"), id="seq3")
        self.annotator._msa = MultipleSeqAlignment([seq1, seq2, seq3])
        
        # Test match
        outcome = self.annotator.determine_outcome(['M', 'A'], 1)
        self.assertEqual(outcome, 100.0)
        
        # Test mismatch
        outcome = self.annotator.determine_outcome(['M', 'A'], 4)
        self.assertEqual(outcome, 66.67)


    def test_get_verdict(self):
        self.assertEqual(self.annotator.get_verdict(90.0), "neutral")
        self.assertEqual(self.annotator.get_verdict(95.5), "neutral")
        self.assertEqual(self.annotator.get_verdict(89.9), "deleterious")
        self.assertEqual(self.annotator.get_verdict(0.0), "deleterious")

    def test_calculate_blosum_score(self):
        # Setup: MSA with conserved 'M' at pos 1
        seq1 = SeqRecord(Seq("MA"), id="seq1")
        seq2 = SeqRecord(Seq("MA"), id="seq2")
        self.annotator._msa = MultipleSeqAlignment([seq1, seq2])
        self.annotator._gene_record = SeqRecord(Seq("ATGC")) # M, ...
        
        # Test 1: M -> T (Deleterious)
        # We need translated seq of mutation. M -> T
        # Original: M
        # Mutated: T
        # MSA: M, M
        # Matrix(M, M) is usually high (e.g. 5)
        # Matrix(T, M) is usually low 
        
        # Determine amino sequence input. 
        # If we mutated pos 1 (ATG) to T... wait, A->T is T...
        # Let's just simulate the amino_sequence passed in.
        amino_seq = ['T', 'A']
        
        # We need real biopython matrix for this test or mock it.
        # Let's hope biopython is installed in test env (it is).
        try:
            stats = self.annotator.calculate_blosum_score(amino_seq, 1, "T")
            # We expect Delta > 0 (Original score > Mutated score)
            self.assertGreater(stats['delta'], 0)
            self.assertEqual(stats['verdict'], 'deleterious')
            
        except ImportError:
            self.skipTest("Biopython substitution matrices not available")

    @patch("builtins.open", new_callable=MagicMock)
    @patch("csv.DictReader")
    def test_process_batch(self, mock_reader_class, mock_open):
        # Mock CSV content instance
        mock_reader_instance = MagicMock()
        mock_reader_instance.__iter__.return_value = iter([
            {'position': '1', 'nucleotide': 'A'}, 
            {'position': '4', 'nucleotide': 'G'}
        ])
        # Mock fieldnames (initial read)
        mock_reader_instance.fieldnames = ['position', 'nucleotide']
        
        # When csv.DictReader() is called, return our mock instance
        mock_reader_class.return_value = mock_reader_instance
        
        # Mock internal methods to avoid complex setup
        self.annotator.insert_snp = MagicMock(return_value=MutableSeq("ATGC"))
        self.annotator.translate_sequence = MagicMock(return_value=['M', 'A'])
        self.annotator.determine_outcome = MagicMock(return_value=100.0)
        self.annotator.get_verdict = MagicMock(return_value="neutral")
        self.annotator.calculate_blosum_score = MagicMock(return_value={'delta': 0.0, 'verdict': 'neutral'})
        
        results = self.annotator.process_batch("dummy.csv")
        
        self.assertEqual(len(results), 2)
        self.assertEqual(results[0]['position'], 1)
        self.assertEqual(results[0]['nucleotide'], 'A')
        self.assertEqual(results[0]['conservation_score'], 100.0)

if __name__ == '__main__':
    unittest.main()
