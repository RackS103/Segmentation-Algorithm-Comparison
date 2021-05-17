#  A Comparison of two Compositional Segmentation Algorithms for Genomic Sequences
Rac Mukkamala, May 2021

## Abstract
Of the various segmentation algorithms created to predict the locations of compositionally homogeneous domains within genomic sequences, two of the most widely used algorithms are IsoPlotter (Elhaik et al. 2010b) and IsoSegmenter (Cozzi et al. 2015). However, these two algorithms yield significantly different predictions, and no study to date has thoroughly examined their differences. Here, I present a detailed comparison of the IsoPlotter and IsoSegmenter algorithms, using a library of simulated random genomic sequences as a benchmark to test algorithm performance and accuracy. Each simulated genomic sequence consisted of multiple simulated compositional domains which were assigned distinct guanine-cytosine (GC) percentages based on the isochore families model (Bernardi 2000). Of the 2,000 simulated sequences generated in this study, 1,100 consisted of domains assigned equal lengths, and the other 900 sequences contained domains assigned variable lengths based on a power-law distribution. My results show that IsoPlotter significantly outperforms IsoSegmenter under a variety of test scenarios, and that IsoSegmenter consistently predicts the existence of large (>200,000bp) domains regardless of underlying genomic architecture. However, there is room for both algorithms to be improved upon, such as IsoPlotter’s tendency to underpredict compositional domain sizes.

## Repository Contents
This repository contains all supplementary data, figures, and scripts used as a part of this research project. Below is a summary of the directories and contents of this repo:
- **src**: Contains all the Python source code files used to generate simulated sequences and score the performance of isoSegmenter/isoPlotter on these sequences.
- **R**: Contains the R code used to create plots and conduct data analysis
- **data**: Contains CSV files which list the performance and number of correct predictions of isoPlotter/isoSegmenter on each of the 2,000 simulated sequences created.
- **figures**: Contains all figures created from the data and added to the final manuscript
- **docs**: Contains a final copy of the manuscript as well as a list of cited references
