# Import necessary libraries and load modules

import Bio as Bio

from Bio import SeqIO
from Bio import AlignIO
from Bio import Phylo
from Bio import Align
from Bio.Seq import Seq
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor

import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

def load_alignment_get_distance(file:str):
    
    ### Function for loading given multiple sequence alignment data and outputting DNA distance in a matrix format
    ###
    ### Input:
    ###       File - string path to a DNA (MAFFT) profile-multiple sequence allignment file (clustal, .aln)
    ###
    ### Output:
    ###        A Biopython DistanceMatrix (Bio.Phylo.TreeConstruction.DistanceMatrix) and distance calculator
        
    # Load the multiple sequence alignment file
    read_aln = AlignIO.read(file, "clustal")
    
    # Initialise DNA distance calculator
    cal_dist = DistanceCalculator("identity")
    
    # Generate a distnace matrix from your alignment
    dm = cal_dist.get_distance(read_aln)
    
    return(dm, cal_dist)
    
def distance_matrix_to_df(dm):   
    
    ### Function for formating the BioPython distance matrix into a pandas dataframe (df)
    ###
    ### Input:
    ###        Biopython's DistanceMatrix
    ###
    ### Output:
    ###        Pandas Dataframe with alignment distance matrix
    
    # Create a pandas dataframe by appending all distances into a list 
    distances = []

    for dislst in dm.matrix:
          distances.append(dislst)

    dstmtxDF = pd.DataFrame(distances)
    
    # Replace all 0.0000 in the df with NaN to avoid identical matches with same sequences
    dstmtxDF.replace(0.000000, np.nan, inplace=True)
    
    # Set all the dog ID's as the column headers
    dstmtxDF.columns = dm.names
    
    # Create a new column with dog ID's for each row
    dstmtxDF['seq_id'] = dm.names
    
    # Set the row ID's as our index column in a new variable
    orgDF = dstmtxDF.set_index("seq_id")
    
    # Create a new column that returns the dog ID of the closest distance for each of the rows respectively
    orgDF['closest_dist'] = orgDF.idxmin(axis=1)
    
    return(orgDF)

def seq_match_in_df(orgDF: pd.DataFrame,
                    match_seq:str):
    
    ### Function for matching and displaying the closest sequence to the mystery sequence using our pandas df
    ###
    ### Input:
    ###        Pandas df with alignment distance matrix and a matching dog sequence ID to our mystery sequence
    ###
    ### Output:
    ###        The seqence that is closest related to the input sequence string (our mystery sequence ID)
    
    # Identifies the closest sequence match to our mystery sequence and stores it in a variable to return (Output Goal)
    mystery_match = orgDF.loc[match_seq, "closest_dist"]
        
    return(mystery_match)

def prepare_seqs(msa):
    
    ### Function for preparing the indentified seqences of interes
    ###
    ### Input:
    ###        Our profile-multiple sequence alignment file (note: dog IDs of interest were pre-determined from our df) 
    ###
    ### Output:
    ###        Individual mystery and target sequence fasta files saved in output folder
      
    # Read the multiple sequence alignement file and store in variable
    testing_seqs = [i for i in AlignIO.read(msa, "clustal")]
    
    # Initialise a list for both sequences of interest
    seq_1 = []
    seq_2 = []
    
    # Append the mystery and target sequence/fasta details into the new lists
    for target in testing_seqs:
        if target.id=="New|gb|KM061522":
            seq_1.append(target)

    for target in testing_seqs:
        if target.id=="gb|AY656744.1|":
            seq_2.append(target)
    
    # Test to see if the lists have been populated with our intended records
    assert len(seq_1) > 0, "Assertion error: Sequence 1 not loaded in list"
    assert len(seq_2) > 0, "Assertion error: Sequence 2 not loaded in list"
    
    # Writing the new sequence list pair to new fasta files in our output folder, because we cannot feed lists into a pairwise aligner
    SeqIO.write(seq_1, "../output/KM061522.fasta", "fasta") 
    SeqIO.write(seq_2, "../output/AY656744.fasta", "fasta")
                
    
def pairwise_alignment():
    
    ### Function for a pairwise seqence alignment of selected matching sequences of interest
    ###
    ### Input:
    ###        Our mystery sequence and its closest match target determined by the distance matrix df
    ###
    ### Output:
    ###        A pairwise sequence alignment, text file written into output folder to observe differences. Number of matched nucleotides in the alignment are displayed
    
    # Pairwise aligner does not accept lists, so a new fasta file was generated to feed biopython's .seq object directly into the aligner
    myst_seq = r"../output/KM061522.fasta"
    targ_seq = r"../output/AY656744.fasta"
    
    # Read fasta sequences and store as input variables
    input_1 = SeqIO.read(myst_seq, "fasta")
    input_2 = SeqIO.read(targ_seq, "fasta")
    
    # Test to ensure sequences are the same length for alignment
    assert len(input_1.seq) == len(input_2.seq), "Assertion error: Sequence lengths are not equal"
    seq_len = len(input_1.seq)
    print(f"\nTotal nucleotide sequence length to align is: {seq_len}\n") 
    
    # Initialise our aligner parameters for global alignment and scoring (using as a nucleotide match counter)
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 1
    aligner.mismatch_score = 0
    
    # Start global pairwise alignment of the sequences. Prints the number of nucleotide matches and alignment display on terminal/console (Output Goal)
    for compare in aligner.align(input_1.seq, input_2.seq):
        print("Nucleotides matched in sequence = %.1f" % compare.score)
        print("\n\n\nA global pairwise sequence alignment has been successfully generated. Please observe the representation below for any differences:\n\n\n")
        print(compare)
        break

    # Set variable outside of the for loop (nucleotide match count to be returned)
    seq_match = compare.score
    
    # Returning input_1 for a later function, we will use its length in a percentage calculation along with seq_match
    return(input_1, seq_match)

def alignment_match_score(input_1, seq_match):  
    
    ### Function for calculating the percentage match of the two sequences aligned
    ###
    ### Input:
    ###        input_1 sequence variable and the seq_match score 
    ###
    ### Output:
    ###        Percentage of total sequence aligned and found matching
    
    # Nucleotide sequence length input into the aligner
    seq_len = len(input_1.seq)
   
    # Calculate the percentage match in the alignment
    pct_match = (seq_match/seq_len)*100
    
    # Print out the percentage match and the details of the IDs fed into the aligner (Output Goal & Stretch Goal 1)
    return(print(f"\n\n\nMystery sequence (KM061522) has a {pct_match:.2f}% match with AY656744.1 (English Springer Spaniel).\n\n\n"))

def construct_phylo_tree(cal_dist, dm):  
    
    ### Function for constructing a phylogenetic tree for the multiple seqence alignment distance matrix
    ###
    ### Input:
    ###        Distance calculator and distance matrix 
    ###
    ### Output:
    ###        Generates a neighbour-joining phylogenetic tree figure and saves it to the output folder 
    
    # Variable that combines our distance calculator and tree constructor
    constructor = DistanceTreeConstructor(cal_dist)
    
    # Generates our neighbour-joining phylogeny tree from our msa distance matrix
    nj_tree = constructor.nj(dm)
    
    # Plot a neighbour joining tree from constructor and distance matrix, constrained layout minimises overcrowding
    fig = plt.figure(figsize=(100, 100), dpi=300, layout="constrained")

    # Set plot font and label sizes
    matplotlib.rc("font", size=2)
    matplotlib.rc("xtick", labelsize=2)
    matplotlib.rc("ytick", labelsize=2)
    
    # Sort clades according to terminal nodes
    nj_tree.ladderize()
    
    # Set plot axes parameters
    axes = fig.add_subplot(1, 1, 1)
    
    # Get current figure so it can save to file
    fig1 = plt.gcf()
    
    # Draw the phylogeny tree
    Phylo.draw(nj_tree, axes=axes)
    
    # Set our output folder path to save file
    save_results_to = "../output/"
    
    # Save phylogeny plot as .png in the output folder (Stretch Goal 2)
    fig1.savefig(save_results_to + 'phylogenetic_tree.png', dpi=300, bbox_inches='tight')

    # Print message to prompt the checking of a phylo tree generated in output folder
    print("\nA phylogenetic tree has been successfully generated. A .png file has now been saved in the output folder.\n\n\n")
    
    return(nj_tree)

def draw_ascii_tree(nj_tree):  

    ### Function for presenting basic ascii representation of our tree on terminal/console
    ###
    ### Input:
    ###        nj_tree phylogeny tree
    ###
    ### Output:
    ###        Phylogenetic tree drawn in ascii characters which helps clarify branches and nodes
    
    print("Below is a basic ascii tree representation for visual clarity. Please don't forget to scroll up and review all output from the beginning. Thank you!\n\n")
    
    # Basic ascii tree representation displayed on terminal (Stretch Goal 2)
    return(Phylo.draw_ascii(nj_tree))

if __name__ == '__main__':
    
    # Define file path
    msa = r"../data/mafft_alignments/second_MSA.aln"
    # Mystery sequence to match
    match_seq = "New|gb|KM061522"
    # Get distance matrix and calculator
    dm, cal_dist = load_alignment_get_distance(msa)
    # Generate pandas df from distance matrix
    orgDF = distance_matrix_to_df(dm)
    # Closest match in distance matrix df depnding on ID input
    mystery_match = seq_match_in_df(orgDF, match_seq)
    # Prints our desired outputs onto terminal
    print(f"\n\n\n\nThe closest related dog to our mystery sequence is: {mystery_match} (English Springer Spaniel)\n\n\n")
    print("Loading sequence data... \n\nAppending sequence closest to the mystery sequence... \n")
    # To fetch target sequences needed for aligning
    prepare_seqs(msa)
    print("Sequences successfully appended and saved into the output folder.\n")
    # Needed for match calculations
    input_1, seq_match = pairwise_alignment()
    alignment_match_score(input_1, seq_match)
    # Phylogeny tree construction
    nj_tree = construct_phylo_tree(cal_dist, dm)
    # Tree representation on terminal
    draw_ascii_tree(nj_tree)