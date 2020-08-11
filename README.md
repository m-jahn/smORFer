# smORFer
smORFer is a small ORF (smORF) detection algorithm that integtates genome, ribosome profiling and translation inition stalling data. Ribosome profiling sequencing data (Ribo-Seq) generates ribosome protected fragments (RPFs) that can exactly located ribosome on mRNA. Translation inition stalling sequencing (TIS-seq) is generated by blocking ribosomes at the start codon and subsequent sequencing of the protected fragments. Both Ribo-Seq and TIS-seq use special antibiotics to stall the ribosomes.

# Introduction to the algorithm

Our algorithm cotains 3 moduls:

* module A: Genome search module to detect all possible ORFs from a gives genome (FASTA)
* module B: Detect smORF candidates using ribosome profiling data
* module C: Detect smORF candidates using translation inition sequencing data.

# Code and example data

We provide our code together with an truncated ***E.coli*** genome input (~100000 first bases) and outputs as examples.

