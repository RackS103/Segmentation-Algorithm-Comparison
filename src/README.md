# Python Source Code

## Information about each script
- `Generate_Seq.py`: generates simulated DNA sequences with domains of homogeneous GC Content. Simulated Sequences can either have equal-length domains or variable-length domains.
- `Testing_Library.py`: 
- `ViewGC.py`: Shows a graph of window-based GC content using the FASTA sequence and JSON sequence data files inputted.
- `score_eqlengths.py`: Compiles data for correlational analysis between ground truth domain length for EQUAL-LENGTH SEQUENCES and average predicted length by both isoPlotter and isoSegmenter.
- `score_varlengths.py`:Compiles data for correlational analysis between ground truth number of domains for VARIABLE-LENGTH SEQUENCES and average predictions by both isoPlotter and isoSegmenter.
- `score_isoplotter.py`: Scores all predictions made by isoPlotter (equal length and variable length) and outputs summary data
- `score_isosegmenter.py`: Scores all predictions made by isoSegmenter (equal length and variable length) and outputs summary data

## Python Package Requirements
- biopython (version 1.78)
- numpy (version 1.20.2)
- matplotlib (version 3.4.1)
