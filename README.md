# Dog genomics :dog:

## BBK UoL, MSc Bioinformatics - Biocomputing module - Dog genomics project

## Overview

This project is intended to identify the closest DNA sequence in a given FASTA file (containing multiple sequences), matched to an individual "mystery sequence" contained in a separate FASTA file.

I was provided with the "dog_breeds.fa" and "mystery.fa" files which have been stored as a reference in the "data" folder. 

Initially, I had installed locally the ClustalW, MUSCLE & MAFFT sequence alignment software tools. However, it was a challenge to get the Biopython command-line wrappers to work with them properly. However, when I did, ClustalW struggled to generate the multiple sequence alignment in a reasonable amount of time, maybe because of the large fasta file I was working with (1.7mb), limited computing power or perhaps just my dodgy code. Instead of using code with a Biopython wrapper for this type of software, I instead opted to use an online server as it was much faster in generating the alignments that I needed to work with.

I had chosen to use the MAFFT multiple sequence alignment tool (https://mafft.cbrc.jp/alignment/server/) to first align the "dog_breeds.fa" file using default settings (first_MSA.aln/fasta). After this was generated, I saved the aligned file and ran it again paired together with the "mystery.fa" added in a profile-msa (https://mafft.cbrc.jp/alignment/server/add_sequences.html). This was useful because you have the option to insert "New|" at the head of title for each new sequence introduced, in this case it was just our single sequence in "mystery.fa" and it helped me easily identify the mystery sequence whilst working with it. Default settings were run again, and I saved the aligned files (second_MSA.aln/fasta). Both the first and second MSAs are saved for reference under "mafft_alignments" in the "data" folder. I only needed to use the second (profile-MSA) alignment file (.aln). It's worth mentioning that some extended information is lost during the alignment, such as the dogs' breed, sub-species etc.., so I have used the FASTA IDs as a core way to identify all the sequences.

I needed to use pandas to generate a dataframe that could effectively organise and isolate the shortest distances along column and row IDs and find their closest related matches accordingly. Many thanks to Dr Tristan Cragnolini for suggesting the use of ".idxmin()" in the distance matrix dataframe. After I made a few tweaks to the dataframe values and assigned the appropriate IDs to the headers and index column, it was a simple and effective function to find the closest related sequence of any given sequence (ID).  

Numpy is basic a requirement for Biopython to operate and I have also used it to manipulate data in the pandas dataframe to avoid identical sequence matching using the distance matrix values.

Matplotlib is used to generate the phylogenetic tree along with several modules that come with Biopython. 

## How to Install and Run the Project

I used Anaconda Navigator to manage my environments when working on this project. I have exported my dog_genomics package dependencies into a YAML file (env.yml) listing what will need to be configured to successfully run my python code "seq_matcher.py" on your system. It's worth noting that you will need to have the latest version of Biopython (1.80) installed for my code to run. If you are working in an anaconda environment, unfortunately the conda installation route only fetches up to version 1.78. You will need to install Biopython 1.80 directly from your command-line using pip instead (https://biopython.org/wiki/Download). 

I have run my python code "seq_matcher.py" using the command prompt within my anaconda environment (as per the YAML specs) and it is fully functional.   

## How to Use the Project

This project contains a standalone piece of python code intended to be run within the confines of this project folder, to meet the objectives set out by the Biocomputing module, which are:

##### Output: Closest sequence to "Mystery" sequence, and the difference

##### Stretch Goal 1: Probabilities

##### Stretch Goal 2: Reconstructed phylogeny

I have included extensive documentation for each line of my code, and in brackets, (i.e. Output, Stretch Goal 1 or Stretch Goal 2) which objective I am trying to solve at that given instance. 

## Issues and Future Improvements

If I had more time to dedicate to the project I would introduce more tests, perhaps more focused unit, or integration testing in a separate python file. I would probably also attempt to apply object-oriented programming and create some classes for the more distinct aspects of the project, rather than sticking solely to functional programming as seen throughout this project. 

As mentioned already, my code has been generated for a niche task using curated data, but with some tweaking it might be more generally applicable to output sequence matching and phylogenetic tree building for more data input types. Please be aware that I have not tested any other file formats not already present in the code, such as "embl, genbank, seqxml etc.." for any sequence file inputs, or "stockholm, phylip, nexus etc.." for phylogenetic tree building. It would be nice to enable the user to have some input in a flexible piece of code, where you could specify more parameters if prompted, such as "local" instead of "global" alignments, or an "UPMGA" tree instead of "NJ". In hindsight, I would have also started by converting the nucleotide sequences into amino acids as they may have been a better source for comparison when calculating alignment probabilities and statistics as stated in (https://www.ncbi.nlm.nih.gov/BLAST/tutorial/Altschul-1.html). If I had more time, I may have also tried working with a PAML command-line wrapper to calulate the P-values from a Chi-squared distribution of my phylogenetic tree data (https://biopython.org/docs/dev/api/Bio.Phylo.PAML.chi2.html).

