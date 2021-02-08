# 2021_TrmD
Analysis of ribosome profiling and RNA-seq data from E. coli cells lacking the tRNA methylase TrmD

Developed by Allen Buskirk using code originally by Fuad Mohammad and Nick Guydosh

Dependencies:

Jupyter notebook
python 2.7
skewer v0.2.2
bowtie v1.2.2
BCBio
Biopython

The raw sequencing data are available as FASTQ files and processed WIG files at the GEO, accession number GSE165592.

Processing of the FASTQ files, mapping with bowtie, and storing ribosome density are all described in the iPython notebook. Bowtie indexes are given in a separate file, and the GFF file for annotation of version 2 of the MG1655 genome is given as coli2.gff.

The iPython notebook also describes how lists of genes and genecounts are generated and how pauses are calculated depending on the amino acid encoded in the A site codon, or all 61 sense codons. 

Analyses in the paper using DESEQ or XTAIL to determine genes whose expression were statistically significantly altered by TrmD depletion were carried out in R. The script and input files are given as separate folders. 



