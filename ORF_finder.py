# HOW TO RUN SCRIPT - python3 ORF_finder.py genome.fasta -m 100
import argparse
parser = argparse.ArgumentParser(description='find ORFs in a nucleotide genome sequence')
parser.add_argument("fileName", help="Input filename of the genome sequence (.fasta)")     
parser.add_argument("-m", "--minOrfSize", help="the minimum length of an ORF in amino acids", type=int, default=50)
args = parser.parse_args()




inputfile = args.fileName
min_orf_length = args.minOrfSize

print()
print("########################-WARNING-##############################\nAll amino acids with 'N' present in the codon have been excluded.\nHowever as only GGN can be Glycine for example, amino acids where this is the case such as Leucine, Valine, Serine, Proline, Threonie, Alanine and Argenine have been included for that specific codon only.")


order = []                                           # Empty list
seqs = {}                                            # Empty dictionary
print()
def fastaread(filename):
    fastafile = open(filename, "r")                  # Open file
    lines = fastafile.readlines()                     # Read lines in one go into a list, 'lines'
    for line in lines:                                # Iterate over fasta file, line by line
        if line.startswith('>'):                      # If statement for if line starts with ">"
                words = line.split()                  # Separates the header words so the if len(words) >= 1: works and can print 1st word
                if len(words) >= 1:                   # Print only first word of species (I.claudius)         
                    print(words[0])
    current_key = None                                # Initialize a variable to keep track of the current key
    for line in lines: 
            line = line.strip().upper()               # Remove whitespace and make upper case to work with amino dictionary
            if line.startswith('>'):
                current_key = line                    # If the line starts with ">" Set that line to the variable "current key"
                seqs[current_key] = []                # Initialise an empty list for the key in the seqs dictionary
            elif current_key is not None:             # If the "current key" key in dictionary is empty(None)
                    seqs[current_key].append(line)     # If line doesn't start ">" and current_key has been assigned, append line to the key list                                                                      associated with the sequence header in the seqs dictionary 

    
    first_seq_key = list(seqs.keys())[0]                    # Putting the first sequence in the dictionary into variable first_seq_key
    return "".join(seqs[first_seq_key])                             # Returning the concatenated sequence with.join
    
seq = fastaread(inputfile)                                          # Storing the contents of fastaread function into variable "seq"

def reverse(s):                                                     # Creating a function to reverse the sequence
    rev_s = s[::-1]                                                 # s is the sequence provided and [::-1] slices it backwards
    return rev_s

def rev_complement(s):
    complement_dictionary = {"a": "T", "c": "G", "g": "C", "t": "A", "n": "N"} # Nucleotide complement dict (In caps to work with amino dict)
    mytable = reverse(s).maketrans(complement_dictionary)           # .maketrans produces a mapping table that I can use with .translate
    new_sequence_reverse_complement = reverse(s).translate(mytable) # Nucleotides reversed and replaced with complement in dict using .translate
    return new_sequence_reverse_complement                          # Returns variable made from last line of code      

#seq -------------------------------------- forward sequence
rev_complement_seq = rev_complement(seq) # - reverse sequence

gencode = { "AAA": "K", "AAC": "N", "AAG": "K", "AAT": "N", "ACA": "T", "ACC": "T",
    "ACG": "T", "ACT": "T", "AGA": "R", "AGC": "S", "AGG": "R", "AGT": "S",
    "ATA": "I", "ATC": "I", "ATG": "M", "ATT": "I", "CAA": "Q", "CAC": "H",
    "CAG": "Q", "CAT": "H", "CCA": "P", "CCC": "P", "CCG": "P", "CCT": "P",
    "CGA": "R", "CGC": "R", "CGG": "R", "CGT": "R", "CTA": "L", "CTC": "L",
    "CTG": "L", "CTT": "L", "GAA": "E", "GAC": "D", "GAG": "E", "GAT": "D",
    "GCA": "A", "GCC": "A", "GCG": "A", "GCT": "A", "GGA": "G", "GGC": "G",
    "GGG": "G", "GGT": "G", "GTA": "V", "GTC": "V", "GTG": "V", "GTT": "V",
    "TAA": "", "TAC": "Y", "TAG": "", "TAT": "Y", "TCA": "S", "TCC": "S",    # Stop codons set to "" so dont appear in final orf seq and dont
    "TCG": "S", "TCT": "S", "TGA": "", "TGC": "C", "TGG": "W", "TGT": "C",   # affect amino length
    "TTA": "L", "TTC": "F", "TTG": "L", "TTT": "F", "CTN": "L", "GTN": "V",  # Added XXN to dictionary as only GGN can be Glycine for example
    "TCN": "S", "CCN": "P", "ACN": "T", "GCN": "A", "CGN": "R", "GGN": "G" } # Same for Leucine, Valine, Serine, Proline, Threonie, Alanine, Arg
                                                                             
    
def find_orfs(dna, frameoffset, min_orf_length, start_codon="ATG", stop_codons=["TAA", "TGA", "TAG"]):
    orfs = []                                                  # Empty list to store ORFs
    in_orf = False                                             # To track if currently inside an ORF, we are not in orf yet so set to False 
    start_coordinate = 0                                       # Set start coordinate to 0
    start_codon_position = None                                # Start codon set to none as there is no start codon yet

    for frame in range(frameoffset, len(dna), 3):              # Go through 3 chars of DNA, starting at specified frameoffset
        codon = dna[frame:frame + 3]                           # Slice from start frame to frame + 3, put in the variable codon

        if codon == start_codon:                               # Check if variable codon is a start codon(ATG)
            if not in_orf:                                     # If start codon and not in orf 
                current_orf = codon                            # Current orf variable = codon variable from above 
                in_orf = True                                  # Set in orf to true as we are now reading the first ATG which is in orf
                start_coordinate = (frame - frameoffset) + 1   # Frame-frameoffset is orig start point in DNA. +1 to get it back to frame offset
                start_codon_position = frame                   # Setting start codon position as frame
            else:
                current_orf += amino_acid                      # If start codon but already inside an ORF, append the codon to the current ORF being made.

        amino_acid = gencode.get(codon, None)                  # Find amino acid in dictionary, or none if not a valid amino acid

        if 'N' in codon:
            in_orf = False                                     # Exclude the ORF if it contains 'N'

        if in_orf and amino_acid is not None:                  # If inside an ORF and codon is a valid amino acid, it adds the amino acid to the current ORF
            current_orf += amino_acid

        if codon in stop_codons and in_orf:                    # Check if the codon is a stop codon
            in_orf = False                                     # Not in ORF/ends ORF
            if len(current_orf) - 3 >= min_orf_length:         # -3 as ATG was appearing before first Methionine affecting len calculation
                orfs.append((start_coordinate, start_codon_position, current_orf[3:])) #If current ORF's length is greater than or equal to min_orf_length, append a tuple of start_coordinate, start_codon_position, and ORF sequence to the orfs list.

    return orfs


def fasta(orf_dict, species_name):       
    fasta_strings = {}                                         # Empty dict
    
    for key, orfs in orf_dict.items():                                 # Extract each key (Frame1,2,3 etc) and its list of orfs
        fasta_strings[key] = []                                        # Each key(Frame) added to fasta_strings dict as keys and the values are an empty list
        for i, (start_coord, start_codon_pos, orf) in enumerate(orfs, 1): #Look through list from pos 1. Look at the 3 variables for each orf
            amino_acid_length = len(orf)                               # Calculate the length of the amino acid sequence
            header = f">{species_name}_{key}_{i:04}_len{amino_acid_length}_start{start_coord}"#-------------#Create header with species_name, key(Frames),
            fasta_entry = f"{header}\n{orf}"                           # Combine header above orf            #frame number, amino len, nucleotide position
            fasta_strings[key].append(fasta_entry)                     # Append fasta_entry to each key(Frame)
    
    return fasta_strings

species_name = "CLAUD"                                                 # Create new dict for all 6 frames and their corresponding ORFs from find_orfs
ORF_DICT = {
    "Frame_1": find_orfs(seq, frameoffset=0, min_orf_length=args.minOrfSize),
    "Frame_2": find_orfs(seq, frameoffset=1, min_orf_length=args.minOrfSize),
    "Frame_3": find_orfs(seq, frameoffset=2, min_orf_length=args.minOrfSize),
    "Rev_Frame_1": find_orfs(rev_complement_seq, frameoffset=0, min_orf_length=args.minOrfSize),
    "Rev_Frame_2": find_orfs(rev_complement_seq, frameoffset=1, min_orf_length=args.minOrfSize),
    "Rev_Frame_3": find_orfs(rev_complement_seq, frameoffset=2, min_orf_length=args.minOrfSize)
}

#fasta_format = fasta(ORF_DICT, species_name)                         # Using fasta function with ORF_DICT and species name
#for key, value in fasta_format.items():                              # For loop the keys and values in fast_format, .items access key/value for printing 
#    print(key)
#    print("\n".join(value))                                          # Join values list but separate values with new lines
    

    
s_dict = str(ORF_DICT)                                         # Convert dict to string
delimiter = ","
string = s_dict.split(delimiter)                               #Converts the string into list
longest_orf = max(string, key = len)                           #Finding the Longest ORF from the list using max and key is len of each string

total_orfs = sum(len(orfs) for orfs in ORF_DICT.values())      # Sum of orf entries for each frame in dict values
print(f"Total number of ORFs: {total_orfs}\n")
print("Length of longest ORF:", len(longest_orf),)
print()


output_filename = "I.claudias_ORFs.fasta"

with open(output_filename, "w") as output_file:                     # Writing Orfs to FASTA file
    fasta_format = fasta(ORF_DICT, species_name)
    for key, value in fasta_format.items():
        output_file.write(key + "\n")
        output_file.write("\n".join(value) + "\n")    
    
    
    
print("File has been saved as: I.claudias_ORFs.fasta")    
    
    
if __name__ == "__main__":
	pass
else:
	print("run as module\n")
