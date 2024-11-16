## ORF Finder

This tool is designed to identify Open Reading Frames (ORFs) in a nucleotide genome sequence provided in the FASTA format. ORFs are continuous stretches of DNA that potentially encode proteins. Identifying ORFs is crucial in bioinformatics for understanding the protein-coding potential of a genome.

## How to Use

Prerequisites
Python 3.x
Running the Script
To run the script, use the following command:


python3 ORF_finder.py genome.fasta -m 100

genome.fasta: Input filename of the genome sequence in FASTA format.

-m or --minOrfSize: (Optional) The minimum length of an ORF in amino acids. Default is set to 50.

Output: The script generates a FASTA file containing the identified ORFs and base pair location.

## Features

Finds ORFs in six reading frames: three forward frames and three reverse frames.

Excludes ORFs containing ambiguous nucleotides ('N') in the codons.

Identifies the longest ORF among the identified ones.

Provides warnings regarding amino acids with 'N' present in the codon and certain exceptions related to codon usage.
Notes

This tool is intended for research purposes and may require further validation or customization based on specific use cases.
For any questions or issues, feel free to contact the developer.
