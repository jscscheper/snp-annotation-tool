import argparse
import sys
import csv
import logging
from pathlib import Path
from snp_annotator.core import SNPAnnotator, InvalidPositionError, SNPAnnotationError

def parse_arguments() -> argparse.Namespace:
    """
    Command-line arguments.
    """
    parser = argparse.ArgumentParser(description="SNP Annotation Tool")

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-n", "--nucleotide",
                        help="The Single Nucleotide Polymorphism (SNP) to insert (A, C, G, T).",
                        choices=['A', 'C', 'G', 'T'])
    group.add_argument("--batch",
                       help="Path to CSV file for batch processing columns: position, nucleotide).",
                       type=str)
                       
    parser.add_argument("-p", "--position",
                        help="1-based position of the SNP in the nucleotide sequence (required if -n is used).", 
                        type=int)
                        
    parser.add_argument("-m", "--msa",
                        help="Path to the Multiple Sequence Alignment (MSA) CLUSTAL file.",
                        required=True, type=str)
    parser.add_argument("-s", "--dna_sequence",
                        help="Path to the FASTA file containing the DNA sequence.", 
                        required=True, type=str)
                        
    parser.add_argument("--method",
                        help="Scoring method to use.",
                        choices=['conservation', 'blosum'],
                        default='conservation')
                        
    return parser.parse_args()

def main() -> None:
    """
    Main.
    """

    logging.basicConfig(
        level=logging.INFO,
        format='%(levelname)s: %(message)s',
        stream=sys.stdout
    )
    
    args = parse_arguments()
    
    if args.nucleotide and not args.position:
        logging.error("--position is required when --nucleotide is specified.")
        sys.exit(1)
    
    try:
        annotator = SNPAnnotator(msa_path=args.msa, sequence_path=args.dna_sequence)
        annotator.load_data()
        
        if args.batch:
            logging.info(f"Processing batch file: {args.batch}")
            results = annotator.process_batch(args.batch)
            
            writer = csv.DictWriter(sys.stdout, fieldnames=results[0].keys())
            writer.writeheader()
            writer.writerows(results)
            
        else:
            logging.info(f"Processing SNP: {args.nucleotide} at position {args.position}")
            
            mutated_seq = annotator.insert_snp(args.position, args.nucleotide)
            logging.info(f"EVENT: Nucleotide at position {args.position} replaced by {args.nucleotide}.")
            
            amino_seq = annotator.translate_sequence(mutated_seq)
            
            if args.method == 'conservation':
                outcome = annotator.determine_outcome(amino_seq, args.position)
                verdict = annotator.get_verdict(outcome)
                print("OUTCOME (Conservation):")
                print(f"The inserted SNP {args.nucleotide} on position {args.position} has {outcome}% similarity.")
                print(f"Verdict: {verdict.upper()}")
                
            elif args.method == 'blosum':
                stats = annotator.calculate_blosum_score(amino_seq, args.position, args.nucleotide)
                print("OUTCOME (BLOSUM62):")
                print(f"Original Score: {stats['score_original']}")
                print(f"Mutated Score: {stats['score_mutated']}")
                print(f"Delta: {stats['delta']}")
                print(f"Verdict: {stats['verdict'].upper()}")

    except (SNPAnnotationError, FileNotFoundError) as e:
        logging.error(f"{e}")
        sys.exit(1)
    except Exception as e:
        logging.critical(f"An unexpected error occurred: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
