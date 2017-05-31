# Assembly of thermophilic culture/enrichment
In Fall 2014 we aerobically cultured thermophiles on TSA50 at 60C for 3 days. At the end of incubation, the "overgrown" plates were scraped to remove the colonies on the plate. DNA was extracted from these scrapings to represent an enrichment of thermophiles that we shotgun metagenome sequenced and 16S rRNA gene amplicon sequenced. The following work represents efforts to quality control the data and assemble it using megahit. The resulting assembled contigs will be pooled into metagenome assembled genomes (MAGs) by using the relative abundance data of the reads from JGI's sequencing of 12 metagenomes across the chronosequence.
## Software used in this workflow
1. khmer
1. fastqc
1. megahit


## Prepping data for metagenomic analysis

### Interleaving reads
Data starts off as two separate fastq.gz files. The following script interleaves the reads that have mates in each file.
** The following script is called `Interleave_Reads_Culture.qsub` and is found in the current directory**
```
#! /bin/bash
#PBS -l walltime=8:00:00
#PBS -l mem=100Gb
#PBS -l nodes=1:ppn=8
#PBS -e /mnt/scratch/sorens75/Cen13_Culture/
#PBS -o /mnt/scratch/sorens75/Cen13_Culture/
#PBS -N Interleave_Reads_Culture
#PBS -A bicep
#PBS -M jacksonwsorensen@gmail.com
#PBS -m abe
module load GNU/4.8.2
module load khmer/2.0
cd /mnt/scratch/sorens75/Cen13_Culture/
interleave-reads.py /mnt/research/ShadeLab/Shade/20150123_metaG_Centralia/Cen13_05102014Pooled_ctDNA_GCCAAT_L002_R* > combined.fastq
qsub Quality_Filter_Culture.qsub
```
Ran quickly with no issues. 51 GB final size of the combined.fastq file.
### Quality Filtering reads
For my next trick, I make use of fastqc to filter out reads with an average Q score less than 30.
** The following script is called `Quality_Filter_Culture.qsub` and is found in the current directory**
```
#! /bin/bash
#PBS -l walltime=12:00:00
#PBS -l mem=100Gb
#PBS -l nodes=1:ppn=8
#PBS -e /mnt/scratch/sorens75/Cen13_Culture/
#PBS -o /mnt/scratch/sorens75/Cen13_Culture/
#PBS -N Quality_Filtering_Culture
#PBS -A bicep
#PBS -M jacksonwsorensen@gmail.com
#PBS -m abe
module load fastx
cd /mnt/scratch/sorens75/Cen13_Culture/
fastq_quality_filter -Q33 -q 30 -p 50 -i combined.fastq -o combined_filtered.fastq
qsub Match_Reads_Culture.qsub
```
This removed about 9 gigabytes of data. The resulting file `combined_filtered.fastq` is approximately 42 GB.

### Extracting Paired Reads from the Quality Filtered Data
Megahit will throw errors if any of your reads are missing their mate. Some of our reads have likely lost their mate in the previous step due to low quality score.
** The following script is called `Match_Reads_Culture.qsub` and is found in the current directory**
```
#! /bin/bash
#PBS -l walltime=8:00:00
#PBS -l mem=100Gb
#PBS -l nodes=1:ppn=8
#PBS -e /mnt/scratch/sorens75/Cen13_Culture/
#PBS -o /mnt/scratch/sorens75/Cen13_Culture/
#PBS -N Match_Reads_Culture
#PBS -A bicep
#PBS -M jacksonwsorensen@gmail.com
#PBS -m abe
module load GNU/4.8.2
module load khmer/2.0
cd /mnt/scratch/sorens75/Cen13_Culture/
extract-paired-reads.py combined_filtered.fastq
```
This results in two files `combined_filtered.fastq.pe` and `combined_filtered.fastq.se`. The paired end file contains ~34GB of reads and the single end file contains ~ 7.4GB of data. Not terrible but not great either. If only megahit could make use of paired end and single end data at the same time.

### Assembling quality filtered reads
We've finally come to the assembly part of things. Assembly traditionally takes a huge amount of memory. Megahit has drastically cut down on the amount of memory needed for assembly.
** The following script is called `megahit.qsub` and is found in the current directory**
```
#! /bin/bash
#PBS -l walltime=07:00:00:00
#PBS -l mem=400Gb
#PBS -l nodes=1:ppn=12
#PBS -e /mnt/scratch/sorens75
#PBS -o /mnt/scratch/sorens75
#PBS -N Megahit_QC_Culture
#PBS -A bicep
#PBS -M jacksonwsorensen@gmail.com
#PBS -m abe
module load megahit
cd /mnt/scratch/sorens75
megahit --12 combined_filtered.fastq.pe --k-list 27,37,47,57,67,77,87,97,107 -o Megahit_QC_Assembly/ -t $PBS_NUM_PPN
```
