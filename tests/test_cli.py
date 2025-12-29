
import unittest
from unittest.mock import patch, MagicMock
from io import StringIO
import sys
from snp_annotator.cli import main

class TestCLI(unittest.TestCase):
    
    @patch('sys.stdout', new_callable=StringIO)
    @patch('snp_annotator.cli.SNPAnnotator')
    def test_main_cli_conservation(self, mock_annotator_class, mock_stdout):
        # Setup Mock
        mock_instance = MagicMock()
        mock_annotator_class.return_value = mock_instance
        
        mock_instance.insert_snp.return_value = "ATGC"
        mock_instance.translate_sequence.return_value = ['M', 'A']
        mock_instance.determine_outcome.return_value = 100.0
        mock_instance.get_verdict.return_value = "neutral"
        
        test_args = ['cli.py', '-n', 'A', '-p', '1', '-m', 'dummy_msa', '-s', 'dummy_seq']
        with patch.object(sys, 'argv', test_args):
            main()
            
        output = mock_stdout.getvalue()
        self.assertIn("Processing SNP: A at position 1", output)
        self.assertIn("100.0% similarity", output)

    @patch('sys.stdout', new_callable=StringIO)
    @patch('snp_annotator.cli.SNPAnnotator')
    def test_main_cli_blosum(self, mock_annotator_class, mock_stdout):
        mock_instance = MagicMock()
        mock_annotator_class.return_value = mock_instance
        
        mock_instance.calculate_blosum_score.return_value = {
            "score_original": 5.0,
            "score_mutated": 2.0,
            "delta": 3.0,
            "verdict": "deleterious"
        }
        
        test_args = ['cli.py', '-n', 'T', '-p', '1', '-m', 'dummy_msa', '-s', 'dummy_seq', '--method', 'blosum']
        with patch.object(sys, 'argv', test_args):
            main()
            
        output = mock_stdout.getvalue()
        self.assertIn("OUTCOME (BLOSUM62)", output)
        self.assertIn("DELETERIOUS", output)

    @patch('sys.stdout', new_callable=StringIO)
    @patch('snp_annotator.cli.SNPAnnotator')
    def test_main_cli_batch(self, mock_annotator_class, mock_stdout):
        mock_instance = MagicMock()
        mock_annotator_class.return_value = mock_instance
        
        mock_instance.process_batch.return_value = [
            {'position': 1, 'nucleotide': 'A', 'verdict': 'neutral'}
        ]
        
        test_args = ['cli.py', '--batch', 'batch.csv', '-m', 'dummy_msa', '-s', 'dummy_seq']
        with patch.object(sys, 'argv', test_args):
            main()
            
        output = mock_stdout.getvalue()
        self.assertIn("Processing batch file: batch.csv", output)
        # CSV output check
        self.assertIn("position,nucleotide,verdict", output.replace('\r\n', '\n'))
        self.assertIn("1,A,neutral", output.replace('\r\n', '\n'))

if __name__ == '__main__':
    unittest.main()
