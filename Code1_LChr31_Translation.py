# Import the SeqIO module 
from Bio import SeqIO

# Function to translate the nucleotide sequences to protein sequences
def translate_sequences(input_fasta, output_fasta):
    with open(output_fasta, 'w') as out_file: # Open the output FASTA file
        for record in SeqIO.parse(input_fasta, 'fasta'): # Parse the input FASTA file
            for strand, name in [(+1, 'forward'), (-1, 'reverse')]: # Iterate over both the forward and reverse strands
                seq = record.seq if strand == 1 else record.seq.reverse_complement() # Reverse sequence if on the reverse strand
                for frame in range(3): # Translate in all three reading frames
                    subseq = seq[frame:] # Adjust the sequence for the current reading frame
                    if len(subseq) % 3 != 0: # Ensure the sequence length is a multiple of 3
                        subseq = subseq[:-(len(subseq) % 3)]
                    translated_seq = subseq.translate() # Translate to a protein sequence
                    out_file.write(f'>{record.id}_{name}_frame{frame+1}\n{translated_seq}\n') # Write the translated sequence to the output file

input_fasta = 'Chr31_LmajorFriedlin2021.fasta'
output_fasta = 'Translated_Chr31_LmajorFriedlin2021.fasta'
translate_sequences(input_fasta, output_fasta)

# Function to remove stop codons from the translated sequences 
def remove_stop_codons(input_fasta, output_fasta):
    with open(output_fasta, 'w') as out_file: # Open the translated FASTA file
        for record in SeqIO.parse(input_fasta, 'fasta'): # Parse the input FASTA file
            protein_seq = str(record.seq).replace('*', '') # Convert the sequence to a string and remove all stop codons (*)
            out_file.write(f'>{record.id}\n{protein_seq}\n') # Write the modified sequence to the output file

input_fasta = 'Translated_Chr31_LmajorFriedlin2021.fasta'
output_fasta = 'NoStop_Translated_Chr31_LmajorFriedlin2021.fasta'
remove_stop_codons(input_fasta, output_fasta)
