import sys
from Bio import SeqIO

def sequence_cleaner(fasta_file, min_length=0, por_n=100):
    # Create our hash table to add the sequences
    sequences={}

    # Using the Biopython fasta parse we can read our fasta input
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        # Take the current sequence
        sequence = str(seq_record.seq).upper()
        # Check if the current sequence is according to the user parameters
        if (len(sequence) >= min_length and
            (float(sequence.count("N"))/float(len(sequence)))*100 <= por_n):
        # If the sequence passed in the test "is it clean?" and it isn't in the
        # hash table, the sequence and its id are going to be in the hash
            if sequence not in sequences:
                sequences[sequence] = seq_record.id
       # If it is already in the hash table, we're just gonna concatenate the ID
       # of the current sequence to another one that is already in the hash table
            else:
                sequences[sequence] += "_" + seq_record.id


    # Write the clean sequences

    for sequence in sequences:
        print(">" + sequences[sequence] + "\n" + sequence)


userParameters = sys.argv[1:]

try:
    if len(userParameters) == 1:
        sequence_cleaner(userParameters[0])
    elif len(userParameters) == 2:
        sequence_cleaner(userParameters[0], float(userParameters[1]))
    elif len(userParameters) == 3:
        sequence_cleaner(userParameters[0], float(userParameters[1]),
                         float(userParameters[2]))
    else:
        print("There is a problem!")
except:
    print("There is a problem!")
