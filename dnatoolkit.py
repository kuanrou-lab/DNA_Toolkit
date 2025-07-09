# Import necessary libraries
import tkinter as tk  # Import tkinter for GUI development
from tkinter import messagebox, Toplevel, filedialog  # Import tkinter components for dialogs and file selection
import matplotlib.pyplot as plt  # # Import Matplotlib for data visualization
from collections import Counter  # Import Counter for nucleotide counting
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from Bio import SeqIO  # Import Biopython for processing biological sequences
from io import StringIO  # Used for processing string input
from Bio.Seq import Seq  # Import Seq object from Biopython, supports DNA/RNA operations
import textwrap  # Used for formatting strings
import re  # Import regular expressions for sequence format validation

# Color mapping, used to highlight different nucleotides in the GUI
BASE_COLORS = {'A': 'red', 'T': 'green', 'G': 'blue', 'C': 'orange', 'U': 'purple'}
# Limit the maximum length of the DNA sequence to avoid processing overly large data
MAX_SEQUENCE_LENGTH = 5000


class DNAToolkitApp:
    """Main application class, used to create the GUI interface for the DNA sequence processing tool."""

    def __init__(self, root):
        """
        Initialize the main interface.
        :param root: Tkinter main window
        """
        self.root = root  # Bind Tkinter main window
        self.root.title("DNA Toolkit")  # window title
        self.create_widgets()  # Call the method to create interface widgets

    def create_widgets(self):
        """Create interface components, including input boxes and buttons."""

        # Create a top frame for placing buttons
        top_frame = tk.Frame(self.root)
        top_frame.pack(anchor="w", padx=10, pady=10, fill="x")

        # "Choose File" button, clicking it calls the process_input method to select a FASTA file
        tk.Button(top_frame, text="Choose File", command=self.process_input).pack(side="left", padx=10)

        # Label to indicate that users can manually enter a FASTA format DNA sequence
        tk.Label(top_frame, text="or enter the FASTA format sequence:").pack(side="left")

        # Create a text input box for manually entering a FASTA format DNA sequence
        self.sequence_entry = tk.Text(self.root, font=("Courier", 12), width=80, height=10)
        self.sequence_entry.pack(pady=10)

        # Create a button frame to hold the "Clear" and "Submit" buttons
        button_frame = tk.Frame(self.root)
        button_frame.pack(pady=10)

        # "Clear" button, clicking it calls the clear_input method to clear the text input box
        tk.Button(button_frame, text="Clear", command=self.clear_input).grid(row=0, column=1, padx=5, pady=10)
        # "Submit" button, clicking it calls the submit_input method to process the entered sequence
        tk.Button(button_frame, text="Submit", command=self.submit_input).grid(row=0, column=2, padx=5, pady=10)

    def process_input(self):
        """Allows users to select any file, but only processes .txt and .fasta files."""

        # Open a file selection dialog, allowing users to choose a file
        file_path = filedialog.askopenfilename(
            title="Select a Sequence File",  # Set the dialog title
            filetypes=[("All Files", "*.*")]  # Allow selection of all file types
        )

        # If the user does not select any file, show a warning and return None
        if not file_path:
            messagebox.showwarning("Input Error", "Please select a sequence file.")  # Show warning to prompt user to choose a file
            return None  # Return directly, do not execute subsequent code

        # Check if the file extension meets the requirements (only accepts .txt, .fa, .fasta)
        if not file_path.endswith((".txt", ".fa", ".fasta")):
            messagebox.showwarning("Input Error",
                                   "Invalid file format! Only .txt and .fasta files are allowed.")  # Show format error
            return None  # Return directly, do not execute subsequent code

        try:
            # Attempt to open the file and read its content
            with open(file_path, "r") as f:
                seq_input = f.read().strip()  # Read file content and remove leading/trailing spaces

            # Call the validate_and_open method to validate and process the sequence
            self.validate_and_open(seq_input)
        except Exception as e:

            # If an error occurs while reading the file, show an error message
            messagebox.showerror("File Error", f"Error reading file: {e}")

    def submit_input(self):
        """Process the DNA sequence entered manually by the user."""

        # Retrieve the content entered by the user from the text box and remove leading/trailing spaces
        seq_input = self.sequence_entry.get("1.0", tk.END).strip()

        # **Optimization: Show a warning if the user does not enter any content**
        if not seq_input:
            messagebox.showwarning("Input Error", "Please choose a file or enter the .fasta sequence")
            return  # Terminate execution to prevent passing an empty value to validate_and_open
        self.validate_and_open(seq_input)

    def validate_and_open(self, seq_input):
        """
        Validate whether the user-input DNA sequence follows the FASTA format and open the display window.
        :param seq_input: The DNA sequence input by the user (FASTA format)
        """

        # Check if the input starts with ">", which is a basic requirement of the FASTA format
        if not seq_input.startswith(">"):
            messagebox.showwarning("Input Error",
                                   "Invalid sequence format. The sequence must start with '>' (FASTA format).")
            return  # Terminate function execution

        # Extract the FASTA sequence part (remove the first-line description)
        sequence_content = "".join(seq_input.split("\n")[1:])

        # Check if the sequence length exceeds the maximum allowed length
        if len(sequence_content) > MAX_SEQUENCE_LENGTH:
            messagebox.showwarning("Input Error", f"Sequence length exceeds {MAX_SEQUENCE_LENGTH} bp.")
            return  # Terminate function execution

        # Use a regular expression to check if the sequence contains only DNA nucleotides A, T, C, and G
        if not re.fullmatch(r"[ATCG]+", sequence_content):
            messagebox.showwarning("Input Error", "Invalid sequence content. Only A, T, C, and G are allowed.")
            return  # Terminate function execution

        # If all validations pass, open a new window to display the sequence
        SequenceViewer(self.root, seq_input)

    def clear_input(self):
        """Clear the DNA sequence from the text input box."""
        self.sequence_entry.delete("1.0", tk.END)  # Delete all content from the text box


class SequenceViewer:
    """Sequence viewer class, used to display and manipulate DNA sequences."""

    def __init__(self, master, seq_input):
        """
        Initialize the sequence viewer window.
        :param master: The main window object
        :param seq_input: Validated DNA sequence in FASTA format
        """
        self.original_sequence = seq_input  # Store the original DNA sequence
        self.seq_input = seq_input  # Record the current DNA sequence being manipulated

        # Create a new window for displaying the DNA sequence
        self.win = Toplevel(master)
        self.win.title("DNA Toolkit - Sequence Viewer")  # Set the window title

        # Create a text widget to display the DNA sequence
        self.text_widget = tk.Text(self.win, wrap="word", font=("Courier", 12))
        self.text_widget.pack(side="left", fill="both", expand=True)
        self.text_widget.tag_configure("bold", font=("Courier", 12, "bold"))  # Set bold format

        # Set colors based on nucleotide type
        for base, color in BASE_COLORS.items():
            self.text_widget.tag_configure(base, foreground=color)

        self.colorize_sequence()  # Apply coloring to the DNA sequence

        # Create a button frame
        self.button_frame = tk.Frame(self.win)
        self.button_frame.pack(side="right", fill="y")

        self.create_buttons()  # Create functional buttons

    def create_buttons(self):
        """Create buttons for manipulating the DNA sequence."""

        button_params = {"width": 20, "height": 2, "font": ("Arial", 12)}

        # Create "Base Count" button
        self.btn_base_count = tk.Button(self.button_frame, text="Base Count",
                                        command=lambda: self.handle_action("count"), **button_params)
        self.btn_base_count.grid(row=0, column=0, pady=5)

        # Create "Complement Strand" button
        self.btn_complement = tk.Button(self.button_frame, text="Complement Strand",
                                        command=lambda: self.handle_action("complement"), **button_params)
        self.btn_complement.grid(row=1, column=0, pady=5)

        # Create "RNA Strand" button
        self.btn_rna = tk.Button(self.button_frame, text="RNA Strand", command=lambda: self.handle_action("rna"),
                                 **button_params)
        self.btn_rna.grid(row=2, column=0, pady=5)

        # Create "Translate to Protein" button
        self.btn_protein = tk.Button(self.button_frame, text="Translate to Protein",
                                     command=lambda: self.handle_action("protein"), **button_params)
        self.btn_protein.grid(row=3, column=0, pady=5)

        # Create "Export" button
        self.btn_export = tk.Button(self.button_frame, text="Export", command=self.save_to_file, **button_params)
        self.btn_export.grid(row=4, column=0, pady=5)

        # Create "Restore Sequence" button
        self.btn_return = tk.Button(self.button_frame, text="Restore", command=self.reset_sequence, **button_params)
        self.btn_return.grid(row=5, column=0, pady=5)
        self.btn_return.grid_remove()  # Hide the button by default

    def handle_action(self, action):
        """Handle the button clicks and execute the corresponding DNA processing functions."""

        if action == "count":
            self.analyze_sequence()  # Calculate the nucleotide composition of the DNA sequence
            self.update_button_states("count")
        elif action == "complement":
            self.generate_reverse_complement()  # Generate the complementary strand
            self.update_button_states("complement")
        elif action == "rna":
            self.generate_rna()  # Transcribe DNA into RNA
            self.update_button_states("rna")
        elif action == "protein":
            self.generate_protein()  # Translate DNA into a protein sequence
            self.update_button_states("protein")

    def update_button_states(self, mode):
        """
        Update button states based on the current mode of DNA processing.
        :param mode: The current DNA processing mode
        """

        # Enable all buttons by default
        self.btn_base_count["state"] = tk.NORMAL
        self.btn_complement["state"] = tk.NORMAL
        self.btn_rna["state"] = tk.NORMAL
        self.btn_protein["state"] = tk.NORMAL
        self.btn_export["state"] = tk.NORMAL
        self.btn_return.grid_remove()  # Hide "Restore" button initially

        if not hasattr(self, "disabled_actions"):
            self.disabled_actions = set()

        if mode == "restore":
            # Reset all buttons to their original state
            self.disabled_actions.clear()  # Clear tracked disabled actions
            self.btn_base_count["state"] = tk.NORMAL
            self.btn_complement["state"] = tk.NORMAL  # Re-enable "Complement Strand"
            self.btn_rna["state"] = tk.NORMAL
            self.btn_protein["state"] = tk.NORMAL
            self.btn_export["state"] = tk.NORMAL
            self.btn_return.grid_remove()  # Hide "Restore" since we're restored
            return  # Stop further processing

        elif mode == "complement":
            # Disable "Complement Strand" button and show "Restore"
            self.btn_complement["state"] = tk.DISABLED
            self.btn_return.grid()
            self.disabled_actions.add("complement")  # Track that "complement" was clicked

        elif mode == "count":
            # When "Base Count" is clicked, if complement was applied, keep it disabled
            if "complement" in self.disabled_actions:
                self.btn_complement["state"] = tk.DISABLED
            self.btn_return.grid()

        elif mode == "rna":
            self.btn_base_count["state"] = tk.DISABLED
            self.btn_complement["state"] = tk.DISABLED
            self.btn_rna["state"] = tk.DISABLED
            self.btn_return.grid()

        elif mode == "protein":
            self.btn_base_count["state"] = tk.DISABLED
            self.btn_complement["state"] = tk.DISABLED
            self.btn_rna["state"] = tk.DISABLED
            self.btn_protein["state"] = tk.DISABLED
            self.btn_return.grid()

        # Ensure "Restore" button appears only if the sequence has changed
        if self.seq_input != self.original_sequence:
            self.btn_return.grid()

    def reset_sequence(self):
        """Restore the DNA sequence to its original input state."""

        self.seq_input = self.original_sequence  # Reset to original sequence
        self.colorize_sequence()  # Reapply sequence coloring
        self.update_button_states("restore")  # Use "restore" mode to re-enable all buttons

    def colorize_sequence(self):
        """Highlight the DNA sequence with colors and display it in the GUI."""

        self.text_widget.delete("1.0", tk.END)  # Clear the text box content
        desc, seq = self.seq_input.split("\n", 1)  # Split the FASTA description and nucleotide sequence
        self.text_widget.insert(tk.END, desc + "\n", "bold")  # Insert the description and apply bold formatting

        # Iterate through the DNA sequence and apply different colors based on nucleotide type
        for base in seq:
            tag = base if base in BASE_COLORS else None  # Get the corresponding color tag
            self.text_widget.insert(tk.END, base, tag)  # Insert nucleotide into the text box and apply color formatting

    def analyze_sequence(self):
        """Analyze the nucleotide composition of the DNA sequence and visualize the results."""

        fasta_io = StringIO(self.seq_input)  # Convert the string into a file-like object for SeqIO parsing
        seq_record = SeqIO.read(fasta_io, "fasta")  # Parse the DNA sequence in FASTA format
        counts = Counter(str(seq_record.seq))  # Count the number of nucleotides

        # Create a new window to display the analysis results
        result_win = Toplevel(self.win)
        result_win.title("Base Count")  # Set the window title

        # Display the sequence length
        tk.Label(result_win, text=f"Sequence Length: {len(seq_record.seq)} bp").pack()

        # Display the count and percentage of each nucleotide
        for base in "ATGC":
            tk.Label(result_win,
                     text=f"{base}: {counts[base]} ({(counts[base] / len(seq_record.seq)) * 100:.2f}%)").pack()

        # Create a bar chart
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.barh(["C", "G", "T", "A"], [counts[base] for base in "CGTA"], color=[BASE_COLORS[base] for base in "CGTA"])
        ax.set_xlabel("Count")  # Set X-axis label
        ax.set_ylabel("Base")  # Set Y-axis label
        ax.set_title("Nucleotide Counts")  # Set the chart title

        # Embed the Matplotlib figure into the Tkinter window
        canvas = FigureCanvasTkAgg(fig, master=result_win)
        canvas.draw()
        canvas.get_tk_widget().pack()
        return dict(counts)

    def generate_reverse_complement(self):
        """Generate the reverse complement of the DNA sequence."""

        seq = str(SeqIO.read(StringIO(self.seq_input), "fasta").seq.reverse_complement())  # Compute the reverse complement strand
        self.update_sequence(seq, " (Reverse Complement)")  # Update GUI display

    def generate_rna(self):
        """Transcribe the DNA sequence into an RNA sequence."""

        seq = str(SeqIO.read(StringIO(self.seq_input), "fasta").seq.transcribe())  # Perform RNA transcription
        self.update_sequence(seq, " (RNA)")  # Update GUI display

    def generate_protein(self):
        """Translate the DNA sequence into a protein sequence."""

        seq_record = SeqIO.read(StringIO(self.seq_input), "fasta")  # Read the DNA sequence in FASTA format
        seq = str(seq_record.seq)  # Extract the DNA string

        # Since DNA translation requires groups of 3 nucleotides per amino acid, ensure the sequence length is a multiple of 3
        trimmed_length = len(seq) - (len(seq) % 3)  # Compute the translatable DNA sequence length
        trimmed_seq = seq[:trimmed_length]  # Extract the translatable part

        protein_seq = str(Seq(trimmed_seq).translate())  # Perform protein translation

        self.update_sequence(protein_seq, " (Protein)", colorize=False)  # Update GUI display

    def update_sequence(self, new_seq, suffix, colorize=True):
        """
        Update the current DNA sequence and display it in the GUI.
        :param new_seq: The new DNA/RNA/protein sequence
        :param suffix: A label for the new sequence (e.g., RNA, Protein, Reverse Complement)
        :param colorize: Whether to apply color highlighting to the sequence; enabled by default
        """

        self.seq_input = f">{SeqIO.read(StringIO(self.seq_input), 'fasta').description}{suffix}\n{textwrap.fill(new_seq, 70)}"
        self.text_widget.delete("1.0", tk.END)  # Clear the current text box content
        self.text_widget.insert(tk.END, self.seq_input)  # Insert the new sequence

        if colorize:  # If color highlighting is needed
            self.colorize_sequence()  # Apply color highlighting

    def save_to_file(self):
        """Allow the user to save the current DNA sequence to a file."""

        file_path = filedialog.asksaveasfilename(defaultextension=".txt",
                                                 filetypes=[("Text file", "*.txt"), ("FASTA file", "*.fa")])
        if file_path:
            with open(file_path, "w") as file:
                file.write(self.seq_input)  # Write the current DNA sequence to the file


if __name__ == "__main__":
    # Create a Tkinter window
    root = tk.Tk()
    app = DNAToolkitApp(root)  # Instantiate the application
    root.mainloop()  # Run the Tkinter main loop