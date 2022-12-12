# Load modules and libraries
from Bio import SeqIO
from Bio import AlignIO
import pandas as pd
import numpy as np

# Load sequence data and parse FASTA

for seq_record in SeqIO.parse("../data/dog_breeds.fa", "fasta"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))
    


    
    
    
    
    
    
    
#fasta_seqs = SeqIO.parse(open("../data/dog_breeds.fa", "fasta")

#genomic_data = "../data/..fasta"

#sequences = []
#for seq_record in SeqIO.parse("../data/dog_breeds.fa", "fasta")
#    sequences.append(

#for trg_record in SeqIO.parse("../data/mystery.fa", "fasta"