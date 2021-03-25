# smORFer - small ORF detection algorithm
smORFer is a small ORF (smORF) detection algorithm that integtates genome, ribosome profiling and translation inition stalling data. Ribosome profiling sequencing data (Ribo-Seq) generates ribosome protected fragments (RPFs) that can exactly located ribosome on mRNA. Translation initiation stalling sequencing (TIS-seq) is generated by blocking ribosomes at the start codon and subsequent sequencing of the protected fragments. Both Ribo-Seq and TIS-seq use special antibiotics to stall the ribosomes.

In the following you can find information regarding:

* [General introduction](#introduction-to-the-algorithm)
* [Example data and replicates](#example-data-and-replicates)
* [Code and data](#code-and-example-data) including:
  * [Module A](#module-a)
  * [Module B](#module-b)
  * [Module C](#module-c)
  * [Helper scripts](#helper-scripts)
* [Installation and requirements](#installation-and-requirements)
* [File formats](#file-formats)
* [Contact](#contact)


## Introduction to the algorithm

Our algorithm contains scripts written in bash, Perl and R. The algorithm contains 3 moduls:

* module A: Genome search module to detect all possible ORFs from a gives genome (FASTA)
* module B: Detect smORF candidates using ribosome profiling data
* module C: Detect smORF candidates using translation inition sequencing data.

*An schematic view on the workflow will be link upon publication*

smORFer is a modular algorithm and can be used according to your available data. For further details see the detailed module sections [module A](#module-a), [module B](#module-b), [module C](#module-c).

smORFer is developed and tested for prokaryotes however it will work also for eukaryotes. The main restriction for eukaryotes it the ORF detection procedure of modul A developed by [Baek et al.](https://academic.oup.com/g3journal/article/7/3/983/6027630) (FileS1). This procedure searches for ORFs in continues stretches of sequences and it does not include exon-intron structures.     

Thoughout the different step different file formats are used: **BED** format in the 6-column version (**BED6**) for ORF/gene locations, **BAM** format to store sequencing reads and some **CSV** files. For further information on **BED** and **BAM** format see the section [file formats](#file-formats)

## Example data and replicates

We provide our code together with an truncated input files: ***E.coli*** genome input (~100000 first bases) in **FASTA** format(`ecoli_100k_nt.fa`) and all protein coding genes in **BED** format (`ecoli_genes.bed`). Additionally, we provide examples of Ribo-Seq, calibrated Ribo-Seq and TIS-Seq sequencing reads in the **BAM** format that are already mapped/aligned such that the workflow can be executed using the examples. Each of the **BAM** files contain ~500,000 reads. To execute the examples call the bash script in each subfolder. The exact example calls are give in the  [Code](#code-and-example-data) section.

The algortihm does not use replicate sequencing data directly. A ususal approach is first to check that the replicates show that same behavior for high expressed genes regarding the read counts. Second, the replicates can be merged into one file to increase the detection limit for low expressed genes.  

## Code and example data

In the following we will give a detailed explanation of the different steps to run and execute the modules including code example of each step. In addtion, we will explain how to execute a full module in one-step execution using standard parameters. 

The complete workflow includes command-line tools (e.g. bedtools), a perl script and several R scripts. To simplify the excution we provide bash/shell scripts such that all scripts can be called from the command-line.  

### Module A

Module A is the genome search module and detects all possible ORFs from a gives genome /transcriptome (FASTA). The main step is the detection of putative ORFs within a given genome (FASTA file). There two optional step: the detection of sequences periodicity patter using Fourier Transformation and a filtering for selected regions.

The aim of the first module is to generate a set of putative (sm)ORFs that is used by the other modules. However, you are free to generate any set of regions in **BED** format on you own (see )  

##### Inputs module A

* [required] **FASTA** file containing a genome or transcriptome. Only the upper letters should be used A/T/G/C/N.
* [optional] **BED** file for filtering specific regions of interest 

##### Output module A

* **BED** file containing putative (sm)ORFs or interesting regions to analyse with the other modules
* [optional] **BED** file of putative (sm)ORFs that show a 3-nt sequence periodicity that might be promising candidates.

##### Details on module A

**Step 1**: putative (sm)ORF search
This step searches for all possible ORFs in all reading frame and for a given ORF length. Please note that for large genomes this search can take some time (usually mins to hours). This search assumes only one chromosome in the **FASTA** file.

```
# call script 
bash putative_orfs.sh fasta.fa genome_name out_name min_ORF_length max_ORF_length
# call using example data
bash modulA_pORFs_genome_search/1_pORF/putative_orfs.sh example_data/ecoli_100k_nt.fa U00096.3 pORFs 9 150
```

The output will be written into a folder called `output` inside the location of the folder of the script. Here, this would be `modulA_pORFs_genome_search/1_pORF/output`. The example command create two file: `pORFs.bed` a **BED** file with all putative ORFs and `pORF.txt` a list of the same ORFs including detected start and stop codons. 

If you want to change the used start and stop codons you have to modify lines 16 and 17 of the `porf_bedformat.pl`:

```
@STA=("ATG","GTG","TTG","CTG"); #start codons
@STO=("TGA","TAG","TAA"); #stop codons
```

If you want to perform the ORF search for eukaryotes please use the dedicated script for eukaryotes, which first seperates the chromosomes in the **FASTA** file. Because of the exon-intro structure of eukaryotes you should use transcripts. However, you are free to search on any stretch of nucleotide sequence. Please, note that the search for eukaryotes can take very long (up to days) when the sequences is large. 

```
# call script 
bash putative_orfs_eukaryotes.sh fasta.fa genome_name out_name min_ORF_length max_ORF_length
```

There are few more files written to the output, that are used to handle the different transcripts/chromosomes.

**Step 2 [optional]**: select or unselect regions of interest

This step is used to filter the putative ORFs for specific regions of interest. For the smORFer publication we searched in non-annotated regions, thus we *un-selected* known ORFs. 

```
# call script 
bash unselect_regions.sh regions_to_unselect.bed putative_orfs.bed out_name
# call using example data
bash modulA_pORFs_genome_search/2_region_selection/unselect_regions.sh example_data/ecoli_genes.bed modulA_pORFs_genome_search/1_pORF/output/pORFs.bed pORFs_filtered
```
Again an output folder at the location of the script is created that contains a **BED** file with only putative ORFs that do not intersect with the given regions. 

If you want the oposite filter and select putative ORFs that intersect with your given regions use the `select_regions.sh` in the same way. 

```
# call script 
bash select_regions.sh regions_to_select.bed putative_orfs.bed
```

**Step 3 [optional]**: Fourier transfor of the sequence periodicity 

This step uses the input regions a checks for 3-nt sequence periodicity to prefilter ORFs that might be translated. 

```
# call:
bash FT_GCcontent.sh fasta.fa orfs.bed
# call using example data: 
bash modulA_pORFs_genome_search/3_FT_GCcontent/FT_GCcontent.sh example_data/ecoli_100k_nt.fa modulA_pORFs_genome_search/2_region_selection/output/pORFs_filtered.bed
```

In the output folder you will find the ORFs that show a sequence (GC content) 3-nt periodic pattern. The output file is call `FT_passed.bed`. 

If you wish to modify the cutoff for the FT ratio signal detection, open `FT_GCcontent.R` in `modulA_pORFs_genome_search/3_FT_GCcontent/` and modify line 

```
cutoff <- 3
```

##### One-step module A 

We do not recommend to run the full module A in a one-step procedure. However, you can do this by calling:

```
# call:
bash run_moduleA.sh fasta.fa genome_name out_name min_ORF_length max_ORF_length regions_to_select.bed out_name2
# call using example data: 
bash modulA_pORFs_genome_search/run_moduleA.sh example_data/ecoli_100k_nt.fa U00096.3 pORFs 9 150 example_data/ecoli_genes.bed pORFs_filtered
```

The output files can be found in the according subfolders (see details above step 1 to step 3). Please note, that step 2 runs in selection (not un-selection) mode. 


### Module B

Module B is working on Ribo-Seq data that Ribosome Protected Fragment (RPFs). It thus required aligned/mapped Ribo-Seq data in **BAM** file format. This mapping can be performed using various tool (e.g. [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) or [BWA](http://bio-bwa.sourceforge.net/)). In addition, a set of ORF location is need as it generate by [module A](#module-a). If you want to make use of the Fourier Transform signal detection for RPFs you also need to provide a **BAM** file of the calibrated reads reduced to the one nucleotide that results from the calibration. This calibration is a multi-step procedure. An external manual can be found [here](https://github.com/AlexanderBartholomaeus/MiMB_ribosome_profiling) and detailed explanation will be available soon in [Methods in Molecular Biology](https://www.springer.com/gp/book/9781071611494)  

##### Inputs module B

* [required] **BAM** file mapped Ribo-Seq sequencing reads. 
* [required] **BED** file with ORFs, as generated by [module A](#module-a) 
* [optional] **BAM** file with calibrated Ribo-Seq reads.

Please not the for the mapping the genome name has to be used as in the resulting


##### Output module B

* **BED** file containing translated ORFs
* [optional] **BED** file containing 3-nt translated ORFs

##### Details on module B

**Step 4**: count Ribo-Seq / RPF reads

Here we count mapped RPF reads within each ORF using bedtools and generate a **BED** file with the highest expressed ORFs. 

```
# call script 
bash count_RPF.sh pORFs.bed RiboSeq.bam 
# call using example data
bash modulB_RPF_analysis/4_count_RPF/count_RPF.sh modulA_pORFs_genome_search/2_region_selection/output/pORFs_filtered.bed example_data/RPF.bam 
```

Within the folder of count_RPF.sh an `output` folder will be created. It will contain 3 files: 

* **RPF_counts.txt**: contains results from bedtool's coverageBed command.
* **RPF_high.bed**: **BED** file containing on the ORF with expression of 100 RPK (reads per kilobase of length). 
* 

You may change the cutoff for detecting ORFs as translated by modifying `count_RPF.sh` file line 13. Change the >= 5 to any number you like.

```
# get candidates > 5 RPF (validated)
awk '$7 >= 5' "${script_path}/output/RPF_counts.txt" > "${script_path}/output/RPF_translated.txt"

```

You my change the cutoff of detecting ORF as high expressed that can be subject to the Fourier Transform in the next steps. To do so, modifying `count_RPF.sh` file line 16. Change the `>= 0.1` to any number you like. Please not 100 RPK = 1 RPF per 10 nucleotides ORF length = 0.1 RPF per 1 nt ORF length. 

```
awk '$7/$9 >= 0.1 {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' "${script_path}/output/RPF_counts.txt" > "${script_path}/output/RPF_high.bed"
```

**Step 5[optional]**: Fourier Transform to find 3-nt translated ORF

To find ORF that clearly show a 3-nt periodic signal Fourier Transform is performed on calibrated Ribo-Seq sequencing. 


```
# call script 
bash FT_RPF.sh high_expressed_ORFs.bed calibrated_RiboSeq.bam
# call using example data
bash modulB_RPF_analysis/5_FT_RPF/FT_RPF.sh modulB_RPF_analysis/4_count_RPF/output/RPF_high.bed example_data/RPF_calibrated.bam
```

The output is a **BED** like file with all ORF passing the cutoff with an additional column showing the FT ratio signal. Please note that our example call did not reveal any ORF because the example **BAM** file contain too few reads and the complete file would be too large to be stored here.

**Step 6[optional]**: find most probable start

The Step is experimentally. You may find it useful to detect the most probable start codon for a sets of ORF with the same stop codon if you have only Ribo-Seq and no Tis-Seq data available. 

```
# get most probable start codon from RPF reads
# bash find_best_start.sh RPF_validated.bed RFP_calibrated.bed
# call using example data
# bash modulB_RPF_analysis/6_find_most_probable_start/find_best_start.sh modulB_RPF_analysis/4_count_RPF/output/RPF_translated.txt example_data/RPF_calibrated.bam

```

##### One-step module B

We do not recommend to run the full module B in a one-step procedure. However, you can do this by calling:

```
# call:
bash run_moduleB.sh pORFs.bed RPF.bam RPF_calibrated.bam
# call with example data:
bash modulB_RPF_analysis/run_moduleB.sh modulA_pORFs_genome_search/2_region_selection/output/pORFs_filtered.bed example_data/RPF.bam example_data/RPF_calibrated.bam

```

Please note, only step 4 and step 5 are performed because step 6 is experimental. The output files can be found in the according subfolders (see details above step 4 and step 5).

### Module C

Module C uses mapped Tis-Seq reads to obtain true translated start codons. First, start codons locations are identified and offset (+- one codon) is added. Second, Tis-Seq reads are counted and ORF with more reads than the cutoff are selected. 

##### Inputs module C

* [required] **BAM** file mapped TIS-Seq sequencing reads. 
* [required] **BED** file with ORFs, as generated by [module A](#module-a)  

Please not the for the mapping the genome name has to be used as in the resulting. If you want extract the middle nucleotide from each TIS read to get a more precise signal you can use the [parse_readMiddle.pl](#helper-scripts) perl script.

##### Output module C

* **BED**/**TXT** files containing ORFs with detect start codon TIS-Seq coverage. 

##### Details on module C


**Step 7**: obtain start codon location to prepare for counting

First, the start codon location has to be identified and a new **BED** file with the start codon +- and offset one codon will be created for the counting of TIS-Seq reads in the next step. For our analysis the TIS-Seq reads show nice enrichment for the offset of +- one codon when considering the middle nucleotide of each TIS-Seq read.


```
# call script 
bash start_codon.sh pORFs.bed 
# call using example data
bash modulC_TIS_analysis/7_get_start_codon/start_codon.sh modulA_pORFs_genome_search/2_region_selection/output/pORFs_filtered.bed
```

The **BED** output file can be found in the output folder. If you wish to change the offset, open the `start_codon.R` file and modify the offset in line 5 and 6. The offsets are given in nucleotides:

```
# upstream offset in nucleotides
up_offset <- 3
down_offset <- 3

```

**Step 8**: count Tis-Seq reads and identify translated ORF start codons

Next, TIS-Seq reads will counted and ORF over the cutoff will be returned. If you want extract the middle nucleotide from each TIS read to get a more precise signal you can use the [parse_readMiddle.pl](#helper-scripts) perl script.  

```
# call script 
bash count_TIS.sh start_codons.bed TisSeq.bam 
# call using example data
bash modulC_TIS_analysis/8_count_TIS/count_TIS.sh modulC_TIS_analysis/7_get_start_codon/output/start_codons.bed example_data/TIS.bam 
```

The output folder will contain two **BED** lke files. `TIS_counts.txt` contains all ORFs including a new output columns from bedtool's `coverageBed` commands. Column 7 contains the obtained read counts. The other file `TIS_candidates.txt` contains only the ORF candidates that passed the cutoff. If you wish to change the cutoff you can open the `count_TIS.sh` and change `>= 5` in line 11 to any number you like:

```
awk '$7 >= 5' output/TIS_counts.txt > output/TIS_candidates.txt
```

##### One-step module C

You can run complete module C by calling:

```
# call:
bash run_moduleC.sh pORFs.bed TisSeq.bam
# call using example data: 
bash modulC_TIS_analysis/run_moduleC.sh modulA_pORFs_genome_search/2_region_selection/output/pORFs_filtered.bed example_data/TIS.bam 
```

The output files can be found in the according subfolders (see details above step 7 and step 8). 

### Helper scripts

To manipulate and modify input file we provide few scripts in the `helper_scripts` folder. 

**GFF to BED parser**: `parse_gff.R` This script can parse GFF ot BED files. To use the script open the script and change the input location. In addtion you might need to change the pattern for the gene name detection in line 21.

**Read middle nucleotide**: `parse_readMiddle.pl` extracts the middle nucleotid from each mapped sequencing read. This script can be but used to extract the middle nucleotide for the mapped TIS reads. Below you find an example on few full length TIS reads:

```
# call script 
samtools view -h example_data/TIS_full_read.bam | perl helper_scripts/parse_readMiddle.pl | samtools view -b > out.bam

```

## Installation & Requirements

### Requirements

* Linux operating system 
* bedtools
* Perl
* R (including the following packages):
  * Biostrings (Bioconductor)
  * seqinr (CRAN)
 
#### Installation

##### bedtools

Bedtools can be installed on Fedora/Centos and Debian/Ubuntu using package managers. For more information see [bedtools installation](https://bedtools.readthedocs.io/en/latest/content/installation.html). 

##### Perl 

Perl is usally installed on Linux. If this is not the case you can install it using package managers and search for `perl`. 

##### Installing R

R can be install on Linux by installing it using the package managers. The package is usally called `r-base` (Debian/Ubuntu) or `R` (Fedora/Centos). 

##### Installing R packages

The installation of the required R packages can be done using our installation script (`R` is required first, see above). Move to the downloaded github directory and install using:

```
Rscript install_dependencies.R
```

If the script finishes without error the two dependencies (`seqinr` and `Biostrings` are installed)

## File formats

### BED format

The **BED** format is a tab seperated text file that is usually stored with the extension `.bed` (e.g. `filename.bed`) but the extension is not necessary. We use the **BED6** format with 6 columns, which is the minimal requirement to include strand information. For more information on the **BED** format see [bedtools webpage](https://bedtools.readthedocs.io/en/latest/content/general-usage.html).

### BAM format

The **BAM** format is a common an efficient format to store aligned/mapped sequencing reads. **BAM** file are binary files and not human readable. You can convert the **BAM** file to human readable **SAM** file. For more information see [samtools](https://samtools.github.io/hts-specs/SAMv1.pdf).
 
## More information / contact
For more information about smORFer contact:
bartholomaeus.alexander@gmail.com or 
baban.kolte@chemie.uni-hamburg.de
