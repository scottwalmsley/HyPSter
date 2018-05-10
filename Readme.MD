# HyPster
### *Hybrid peptide identification from mass spectrometry data using word counts.*

Developed in conjunction with Thomas Delong for the analysis of hybrid peptide detection in the Type I diabetes model.
Hybrid peptides can form from a theoretical mechanism within cells.  These hybrid peptides are not evaluable using traditional mass spectrometry identification routines.  As such we are developing a custom workflow and algorithims to determine the identity of un-discovered hybrid peptides. This component reads the k-mer potential from a proteome fasta file.

# Dependencies
C11.
Build with your favorite build tool (VC, gcc).  The executable has worked on both Linux and Windows.

# Usage
This module is a c++ program written for command line usage.

```
USAGE: getNmerFreqs.exe fasta_file Kmer_size <output_file>
``` 





