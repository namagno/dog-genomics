# Dog genomics 

## BBK UoL, MSc Bioinformatics - Biocomputing module - Dog genomics project

## Overview

This project is intended to identify the closest DNA sequence in a given FASTA file, matched to a separate "mystery sequence" contained in a separate FASTA file.

I was provided with the "dog_breeds.fa" and "mystery.fa" files which have been stored as a reference in the "data" folder. 

To begin, I had installed locally the ClustalW, MUSCLE & MAFFT sequence alignment tools as software. However, it was a challange to get the Biopython commandline wrappers to work properly. When I did, because of the large fasta file I was working with (1.7mb), it struggled to generate the multiple sequence alignment in a reasonable amount of time. Instead of using code with a Biopython wrapper for this software, I opted to use an online server as it was much faster in generating the alignments that I needed to work with.

I had chosen to use the MAFFT multiple sequence alignment tool (https://mafft.cbrc.jp/alignment/server/) to first align the "dog_breeds.fa" using default settings (first_MSA.aln/fasta). After this was generated, I saved the aligned file and ran it again together with the "mystery.fa" added in a profile-msa (https://mafft.cbrc.jp/alignment/server/add_sequences.html). This was useful because you have the option to insert "New|" at the head of title of each new sequence introduced, in this case it was our single sequence in "mystery.fa" and it helped me identify the sequence when working with it. Default settings were run and I saved the aligned files (second_MSA.aln/fasta). Both first and second MSAs were saved for reference under "mafft_alignments" in the "data" folder. I only needed to use the second (profile-MSA) alignment file (.aln). It's worth mentioning that some extended information is lost in the alignment, such as the dogs' breed, sub-species etc.. so I have used the FASTA IDs as a core way to identify the sequences.

I needed to use pandas to generate a dtaframe that could effectively organise and isolate the shortest distances along column and row IDs and find closest related matches accordingly. Many thanks to Tristan Cragnolini for suggesting the use of ".idxmin()" in the distance matrix dataframe, it was an effective method after I made a few tweaks to the dataframe values and assigned appropriate IDs to the headers and index column.

Numpy is a requirement for Biopython to operate and I have used it to manipulate data in the pandas dataframe to avoid identical sequence matching using the distance matrix values.

## How to Install and Run the Project

I have used Anaconda Navigator as an easy way to manage my environments and I have exported my dog_genomics package dependencies into a YAML file (env.yml) listing what will need to be configured to sucessfully run my python code "seq_matcher.py". It's worth noting that you will need to have the latest version of Biopython (1.80) installed for my code to run. If you are working in an anaconda environment, the conda installation route only fetches up to version 1.78. You will need to install Biopython 1.80 directly from your commandline using pip (https://biopython.org/wiki/Download). 

I have run the python code using the command prompt within my annaconda environment (as per the YAML specs) and it is fully functional.   

## How to Use the Project

This project contains a standalone piece of python code intended to be run within the confines of this project folder to meet the objectives set out by the module, which are:

• Output: the closest sequence, and the difference
• Stretch Goal 1: Probabilities
• Stretch Goal 2: Reconstructed phylogeny

I have included documentation in my code and in brackets (i.e. Output, Stretch Goal 1 or Stretch Goal 2) which objective I am trying to solve at that instance. 

## Issues and future improvements

If I had more time to dedicate to the project I would introduce more tests, perhaps more focused unit or integration testing in a seperate python file. I would probably also attempt to move more towards object oriented programming and create some classes for more distinct aspects of the project, rather than sticking solely to functional programming as seen in this project. 

As mentioned alread, my code has been generated for a single use basis as per the niche objectives and provided data, but with some tweaking it might be more generally applicable to output sequence matching and phylogenetic tree building for more data input types. Please be aware I have not tested any other file formats not already present in the code, such as "embl, genbank, seqxml etc.." for sequence file inputs, or "stockholm, phylip, nexus etc.." for phylogenetic tree building. It would be nice to have the user have more input in a more modular piece of code where you can specify more parameters if prompted, such as 

