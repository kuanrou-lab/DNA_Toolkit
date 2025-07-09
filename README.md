# DNA Toolkit

Hi! Welcome to DNA Toolkit! 😊🎉 If this is your first time hearing about DNA Toolkit, let me first tell you what it is! 😁🔬

DNA Toolkit is a versatile GUI application for manipulating and analyzing DNA sequences. 🧬 Whether you need to transcribe a DNA strand into an RNA strand, translate an RNA strand into a protein chain, or generate a complementary sequence for a DNA strand, this toolkit enables you to accomplish these tasks quickly, supporting your bioinformatics analysis.

## Features

Currently, DNA Toolkit includes the following features:
 - Calculate the nucleotide composition of a DNA sequence and generate a statistical chart.
 - Obtain the complementary strand of a DNA sequence (reverse complement).
 - Transcribe a DNA sequence into an RNA sequence.
 - Translate a DNA sequence into a protein sequence.
 - Save transcribed or translated sequences in FASTA and TXT formats.

To avoid lagging, the program currently limits the sequence length to 5000 bp. Future updates will work on optimizing this restriction. There are also plans to add ORFs (Open Reading Frames) prediction and a search function for specific DNA fragments to make the tool even more useful.

Actually, in the second edition, the tool introduced ORFs prediction and supported translating ORFs into proteins. However, in the third edition, this feature was removed and replaced with a direct DNA-to-protein translation function. The main reasons for removing ORFs prediction were limited development time, code length constraints, and the author's programming skills at the time, which made it difficult to get the expected results in a short period. So, it was taken out temporarily. In future updates, ORFs prediction will be reintroduced after optimization, allowing start and stop codons to be highlighted in different colors, as well as enabling ORFs to be underlined for better visualization. ORF sequences will also be extractable and displayed in FASTA format in the Sequence Viewer window, with an option to export files.

Finally, we would like to invite you to follow the future development of DNA Toolkit! In the near future, we will continue to explore and implement more features for processing DNA sequences, striving to provide you with an even better user experience. 🌟🔍

## Requirements

To run this project, your computer or laptop must meet the following requirements:
 - Python version 3.0 and above
 - **Python libraries:** matplotlib, Biopython

## Running the Application

This project requires no installation and is easy to run. Simply download or clone this project, then enter the following command in the terminal or command prompt:
```bash
python dnatoolkit.py
```
After running the command, a GUI interface will pop up, allowing you to input a DNA sequence or select a file for analysis.

## Running Unit Tests

This project includes unit tests. If you need to run them, please enter the following command in the terminal or command prompt:
```bash
python -m unittest test.py
```
The test file `test.py` relies on the data in the test folder for testing. Therefore, before running the tests, please ensure that the test folder and `test.py` are located in the same directory. This folder can be found at the project's download location.

## How to Use

1. **Enter DNA Sequence:** Input a DNA sequence in the text box, or click `Choose File` to select a FASTA or TXT file. The sequence must be in FASTA format (starting with `>`). Then the entered or selected sequence will be displayed in a new window called Sequence Viewer.
2. **Clear Input:** Click Clear to remove the current input sequence.
3. **Perform DNA Operations:**
   - Click `Base Count` to analyze nucleotide composition and generate a histogram.
   - Click `Complement Strand` to calculate the complementary strand.
   - Click `RNA Strand` to transcribe the sequence into RNA.
   - Click `Translate to Protein` to translate the sequence into a protein sequence.
4. **Restore Original Sequence:** After generating a complementary strand, RNA, or protein sequence, click `Restore` to revert to the original DNA sequence.
5. **Export Sequence:** Click `Export` to export the sequence currently displayed in the Sequence Viewer window in FASTA or TXT format.

At any address where you can access this project, we also provide you with a FASTA file named `Try_Using_DNAToolkit`, which contains the DNA sequence of human chromosome 11 from position 121,101,243 bp to 121,104,242 bp (a total of 3,000 bp). You can try using the DNA Toolkit to process this sequence to quickly familiarize yourself with the tool. Of course, you may also choose any DNA sequence of interest when using the DNA Toolkit for the first time. This application is incredibly easy to use—you can even try it out right now in just three minutes! ⏳💡 However, please note that this application currently only supports sequences of up to 5,000 bp!

## License

This project is licensed under the MIT License, allowing free use and modification.

## Author

Developed by Kuanrou Fan
