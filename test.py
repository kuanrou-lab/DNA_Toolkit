import unittest  # Import the unittest module for unit testing
from unittest.mock import patch  # Import patch for mocking dependencies in tests
import tkinter as tk
from Bio import SeqIO
from io import StringIO
from Bio.Seq import Seq
from dnatoolkit import DNAToolkitApp, SequenceViewer  # Import the main application classes from the dnatoolkit
import random
import os
from tkinter import Tk
from tkinter import messagebox

TEST_FILE = "test_folder/test_sequence.fasta"  # Define the test file path
EXPORT_TEST_FILE = "test_folder/exported_sequence.fa"  # Define exported file path
TEST_FILE_VALID = "test_folder/test_sequence.fasta"
TEST_FILE_INVALID = "test_folder/invalid_file.pdf"  # Example invalid file


class TestDNAToolkit(unittest.TestCase):
    """This class tests various functionalities of the DNA Toolkit GUI application."""

    @classmethod
    def setUpClass(cls):
        cls.root = tk.Tk()
        cls.root.withdraw()
        cls.app = DNAToolkitApp(cls.root)
        with open(TEST_FILE, "r") as file:
            cls.test_sequence = file.read().strip()  # Read FASTA sequence
        fasta_io = StringIO(cls.test_sequence)
        cls.seq_record = SeqIO.read(fasta_io, "fasta")
        cls.sequence_str = str(cls.seq_record.seq)  # Extract sequence only

    @patch("tkinter.filedialog.askopenfilename", return_value="test_folder/valid_sequence.fasta")
    def test_valid_file_selection(self, mock_file_dialog):
        """Test the selection and correct opening of a valid FASTA file."""

        with patch("builtins.open", unittest.mock.mock_open(read_data=">Valid FASTA Sequence\nATCG")):
            result = self.app.process_input()
            self.assertIsNone(result, "Valid file selection should not raise errors")
            print("Testing valid file format, pass ")

    @patch("tkinter.filedialog.askopenfilename", return_value="test_folder/invalid_file.pdf")
    def test_invalid_file_selection(self, mock_file_dialog):
        """Test the selection of an invalid file format and ensure an error message is displayed."""

        with patch.object(messagebox, "showwarning") as mock_warning:
            self.app.process_input()

            # Verify that the warning messagebox is called with the correct message
            mock_warning.assert_called_with("Input Error", "Invalid file format! Only .txt and .fasta files are allowed.")
            print("Testing invalid file format, pass ")

    @patch("tkinter.filedialog.askopenfilename", return_value="")  # Simulate no file selected
    @patch("tkinter.messagebox.showwarning")
    def test_no_file_selected(self, mock_warning, mock_file_dialog):
        """Test if a warning pops up when the user does not select a file."""

        self.app.process_input()

        # Verify that the warning messagebox is called with the correct message
        mock_warning.assert_called_with("Input Error", "Please select a sequence file.")
        print("Testing no file selected warning, pass ")

    @patch("tkinter.filedialog.asksaveasfilename", return_value=EXPORT_TEST_FILE)
    def test_export_function(self, mock_save_dialog):
        """Test the export function to ensure that the DNA sequence can be correctly saved to a file."""

        sequence_viewer = SequenceViewer(None, self.test_sequence)
        sequence_viewer.save_to_file()
        self.assertTrue(os.path.exists(EXPORT_TEST_FILE), "Exported file was not created!")  # Verify that the exported file was created successfully
        with open(EXPORT_TEST_FILE, "r") as file:
            saved_content = file.read().strip()

        # Ensure the exported file content matches the original sequence
        self.assertEqual(saved_content, sequence_viewer.seq_input.strip(), "Exported file content mismatch!")
        print("Testing export function, pass")
        os.remove(EXPORT_TEST_FILE)

    @patch("tkinter.messagebox.showwarning")
    def test_invalid_sequence_content(self, mock_warning):

        incorrect_file_path = "test_folder/incorrect_content.txt"
        with open(incorrect_file_path, "r") as file:
            wrong_sequence_content = file.read().strip()
        self.app.validate_and_open(wrong_sequence_content)

        # Assert that messagebox.showwarning was called with the expected warning message
        mock_warning.assert_called_with("Input Error", "Invalid sequence content. Only A, T, C, and G are allowed.")
        print("Testing invalid sequence content, pass ")

    @patch("tkinter.messagebox.showwarning")
    def test_invalid_sequence_format(self, mock_warning):

        incorrect_file_path = "test_folder/incorrect_format.txt"
        with open(incorrect_file_path, "r") as file:
            wrong_sequence_format= file.read().strip()
        self.app.validate_and_open(wrong_sequence_format)

        # Assert that messagebox.showwarning was called with the expected warning message
        mock_warning.assert_called_with("Input Error", "Invalid sequence format. The sequence must start with '>' (FASTA format).")
        print("Testing invalid sequence format, pass ")

    @patch("tkinter.messagebox.showwarning")  # Mock messagebox.showwarning
    def test_sequence_length_exceeds_limit(self, mock_warning):

        bases = ["A", "T", "G", "C"]
        random_sequence = "".join(random.choices(bases, k=2001))  # The test sequence is 3000 bp
        extended_sequence = self.test_sequence + random_sequence  # Form a sequence with 5001 bp (>5000bp)

        self.app.sequence_entry.insert("1.0", extended_sequence)  # Simulate user input in the text box
        self.app.submit_input()  # Call the submit_input function to trigger validation

        mock_warning.assert_called_with("Input Error", "Sequence length exceeds 5000 bp.")
        print("Testing sequence length exceeds limit, pass ")

    def test_generate_reverse_complement(self):

        ground_truth_file = "test_folder/test_complement_strand.fa"
        with open(ground_truth_file, "r") as file:
            ground_truth_sequence = file.read().strip()
        ground_truth_seq = "".join(ground_truth_sequence.split("\n")[1:])
        sequence_viewer = SequenceViewer(None, self.test_sequence)
        sequence_viewer.generate_reverse_complement()
        result_sequence = "".join(sequence_viewer.seq_input.split("\n")[1:])
        self.assertEqual(result_sequence, ground_truth_seq, "Reverse complement mismatch!")
        print("Testing generate reverse complement, pass")

    def test_generate_rna(self):

        ground_truth_file = "test_folder/test_rna_strand.fa"
        with open(ground_truth_file, "r") as file:
            ground_truth_sequence = file.read().strip()
        ground_truth_seq = "".join(ground_truth_sequence.split("\n")[1:])
        sequence_viewer = SequenceViewer(None, self.test_sequence)
        sequence_viewer.generate_rna()   # Call the generate_rna function
        result_sequence = "".join(sequence_viewer.seq_input.split("\n")[1:])
        self.assertEqual(result_sequence, ground_truth_seq, "RNA transcription mismatch!")
        print("Testing generate RNA, pass")

    def test_generate_protein(self):

        ground_truth_file = "test_folder/test_translate_to_protein.fa"
        with open(ground_truth_file, "r") as file:
            ground_truth_sequence = file.read().strip()
        ground_truth_seq = "".join(ground_truth_sequence.split("\n")[1:])
        sequence_viewer = SequenceViewer(None, self.test_sequence)
        sequence_viewer.generate_protein()
        result_sequence = "".join(sequence_viewer.seq_input.split("\n")[1:])
        self.assertEqual(result_sequence, ground_truth_seq, "Protein translation mismatch!")
        print("Testing generate protein, pass")

    def test_base_counting(self):

        sequence_viewer = SequenceViewer(None, self.test_sequence)
        actual_base_counts = sequence_viewer.analyze_sequence()
        expected_base_counts = {base: self.sequence_str.count(base) for base in "ATGC"}
        self.assertEqual(actual_base_counts, expected_base_counts, "Mismatch in nucleotide base counts!")
        print("Testing base count, pass")

    def test_complement_strand_then_base_count(self):

        # Step 1: Create a sequence viewer and click "Complement Strand" first
        sequence_viewer_complement = SequenceViewer(self.root, self.test_sequence)
        sequence_viewer_complement.handle_action("complement")  # Generate complement strand

        # Store the modified sequence after complement transformation
        expected_complement = str(SeqIO.read(StringIO(self.test_sequence), "fasta").seq.reverse_complement())
        modified_sequence = "".join(sequence_viewer_complement.seq_input.split("\n")[1:])  # Remove FASTA header

        # Ensure sequence is modified after clicking complement strand
        self.assertNotEqual(self.test_sequence, sequence_viewer_complement.seq_input, "Sequence did not change after complement!")
        self.assertEqual(modified_sequence, expected_complement, "Complement sequence is incorrect!")

        # Step 2: Click "Base Count" after complementing
        base_counts_complement = sequence_viewer_complement.analyze_sequence()

        # Step 3: Create a new sequence viewer and perform direct "Base Count" without complement
        sequence_viewer_original = SequenceViewer(self.root, self.test_sequence)
        base_counts_original = sequence_viewer_original.analyze_sequence()

        # Step 4: Verify Complement Rules (A ↔ T, G ↔ C)
        self.assertEqual(base_counts_complement["A"], base_counts_original["T"], "Mismatch in A/T complement count!")
        self.assertEqual(base_counts_complement["T"], base_counts_original["A"], "Mismatch in T/A complement count!")
        self.assertEqual(base_counts_complement["G"], base_counts_original["C"], "Mismatch in G/C complement count!")
        self.assertEqual(base_counts_complement["C"], base_counts_original["G"], "Mismatch in C/G complement count!")
        print("Testing click complement strand then click base count, pass")


if __name__ == '__main__':
    unittest.main()