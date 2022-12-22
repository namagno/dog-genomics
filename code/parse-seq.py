# Load modules and import libraries

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
    ### Function for loading alignment data and outputting DNA distance in a matrix format
    ###
    ### Input:
    ###       file - string path to DNA allignment file
    ###
    ### Output:
    ###        Biopython DistanceMatrix
        
    # Load alignment file
    read_aln = AlignIO.read(file, "clustal")
    
    # Initialise DNA distance calculator
    cal_dist = DistanceCalculator("identity")
    
    # Get distnace matrix 
    dm = cal_dist.get_distance(read_aln)
    
    return(dm, cal_dist)
    

def distance_matrix_to_df(dm):   
    ### Function for formating the BioPython distance matrix into pandas dataframe
    ###
    ### Input:
    ###        Biopython DistanceMatrix
    ###
    ### Output:
    ###        Pandas Dataframe with alignment distance matrix
    
    # Create a pandas dataframe by appending all distances into a list 
    distances = []

    for dislst in dm.matrix:
          distances.append(dislst)

    dstmtxDF = pd.DataFrame(distances)
    
    # Replace all 0.0000 in df with NaN to avoid identical matches
    dstmtxDF.replace(0.000000, np.nan, inplace=True)
    
    # Set the dog ID's as column headers
    dstmtxDF.columns = dm.names
    
    # Create a new column with dog ID's for each row
    dstmtxDF['seq_id'] = dm.names
    
    # Set the ID's as our index column in a new variable
    orgDF = dstmtxDF.set_index("seq_id")
    
    # Create a new column containing the closest distances for each of the rows respectively
    orgDF['closest_dist'] = orgDF.idxmin(axis=1)
    
    return(orgDF)

def seq_match_in_df(orgDF: pd.DataFrame,
                    match_seq:str):
    ### Function for matching the closest sequence to the mystery sequence using our pandas dataframe
    ###
    ### Input:
    ###        Pandas Dataframe with alignment distance matrix and a sequence ID to match to
    ###
    ### Output:
    ###        The seqence that is closest related to input sequence string, in this case our mystery sequence ID
    
    #Identify the closest sequence match to our mystery sequence
    mystery_match = orgDF.loc[match_seq, "closest_dist"]
        
    return(mystery_match)

def prepare_seqs(msa):
    
    ### Function for preparing seqs of interest
    ###
    ### Input:
    ###        Target sequence and its closest match determined by the distance matrix DF
    ###
    ### Output:
    ###        Sequecnce files saved in output
      
    testing = [i for i in AlignIO.read(msa, "clustal")]
    
    # Initialise list for sequences of interest

    seq_1 = []
    seq_2 = []
    
    # Attempting to append the sequence closest to mystery sequence into new file
    for target in testing:
        if target.id=="New|gb|KM061522":
            seq_1.append(target)

    for target in testing:
        if target.id=="gb|AY656744.1|":
            seq_2.append(target)
    
    # Check if list has been populated
    assert len(seq_1) > 0, "Assertion error: Sequence 1 not loaded in list"
    assert len(seq_2) > 0, "Assertion error: Sequence 2 not loaded in list"
    
    
    # Writing the new sequence pair to a new output file
    SeqIO.write(seq_1, "../output/KM061522.fasta", "fasta") 
    SeqIO.write(seq_2, "../output/AY656744.fasta", "fasta")
                
    
def pairwise_alignment():
    
    ### Function for pairwise seqence alignment of selected matching seqences of interest
    ###
    ### Input:
    ###        Target sequence and its closest match determined by the distance matrix DF
    ###
    ### Output:
    ###        A pairwise sequence alignment, its differences and percentage match
    
    # Pairwise aligner does not accept lists, so a new file needs to be generated to feed into the aligner
    myst_seq = r"../output/KM061522.fasta"
    targ_seq = r"../output/AY656744.fasta"
    
    # Read fasta sequences
    input_1 = SeqIO.read(myst_seq, "fasta")
    input_2 = SeqIO.read(targ_seq, "fasta")
    
    # Ensure sequences are the same length for alignment
    assert len(input_1.seq) == len(input_2.seq), "Assertion error: Sequence lengths are not equal"
    
    # Initialise aligner parameters for scoring
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 1
    aligner.mismatch_score = 0
    
    for compare in aligner.align(input_1.seq, input_2.seq):
        print("Nucleotides matched in sequence = %.1f" % compare.score)
        print("\nA pairwise sequence alignment has been generated in the output folder.\n") # compare, maybe save this to output
        break

    seq_match = compare.score
    
    # Returning input_1 for a later function to to use its length in a calculation and seq_match

    return(input_1, seq_match)
   

def alignment_match_score(input_1, seq_match):  
    
    ### Function for calculating the percentage match of the two sequences provided
    ###
    ### Input:
    ###        input_1 sequence variable and compare.score 
    ###
    ### Output:
    ###        Percentage of sequence match
    
    seq_len = len(input_1.seq)
   
    pct_match = (seq_match/seq_len)*100
    
    return(print(f"Mystery sequence (KM061522) has a {pct_match:.2f}% match with AY656744.1 (English Springer Spaniel)\n\n"))

def construct_phylo_tree(cal_dist, dm):  
    
    ### Function for calculating the percentage match of the two sequences provided
    ###
    ### Input:
    ###        input_1 sequence variable and compare.score 
    ###
    ### Output:
    ###        Percentage of sequence match
    
    constructor = DistanceTreeConstructor(cal_dist)
    
    nj_tree = constructor.nj(dm)
    
    # Plot a neighbour joining tree from constructor and distance matrix

    fig = plt.figure(figsize=(100, 100), dpi=300, layout="constrained")

    #plt.subplots(layout="constrained")
    matplotlib.rc("font", size=2)
    matplotlib.rc("xtick", labelsize=2)
    matplotlib.rc("ytick", labelsize=2)
    
    # Sort clades according to terminal nodes
    nj_tree.ladderize()
    
    # Set axes parameters
    axes = fig.add_subplot(1, 1, 1)
    
    # Get current figure so it can save
    fig1 = plt.gcf()
    
    # Draw the phylogeny tree
    Phylo.draw(nj_tree, axes=axes)
    
    # Set path to save file
    save_results_to = "../output/"
    
    # Save phylogeny plot as png in the output folder
    fig1.savefig(save_results_to + 'phylogenetic_tree.png', dpi=300, bbox_inches='tight')

    return(nj_tree)

def draw_ascii_tree(nj_tree):  

    ### Function for drawing basic ascii representation of our tree
    ###
    ### Input:
    ###        nj_tree
    ###
    ### Output:
    ###        Phylogenetic tree drawn in ascii characters for node clarity
    
    # Basic representation printed on terminal for tree label clarity 
    return(Phylo.draw_ascii(nj_tree))
    

if __name__ == '__main__':
    
    # Define file path
    msa = r"../data/mafft_alignments/second_MSA.aln"
    
    # Seq to match
    match_seq = "New|gb|KM061522"

    # Get distance matrix and calculator
    dm, cal_dist = load_alignment_get_distance(msa)
    
    # Generate pandas df from distance matrix
    orgDF = distance_matrix_to_df(dm)
    
    # Closest match in distance matrix df depnding on ID input
    mystery_match = seq_match_in_df(orgDF, match_seq)

    print(f"\nThe closest related dog to our mystery sequence is: {mystery_match} \n")
    print("Loading sequence data... \n\nAppending sequence closest to the mystery sequence... \n")
    #
    prepare_seqs(msa)
    print("Sequences successfully appended and saved into the output folder. \n")
    
    # 
    input_1, seq_match = pairwise_alignment()
    
    #
    alignment_match_score(input_1, seq_match)
    
    nj_tree = construct_phylo_tree(cal_dist, dm)

    draw_ascii_tree(nj_tree)