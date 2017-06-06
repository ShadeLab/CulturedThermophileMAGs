# Jane's Notes

## Interleaving Reads
Last week I submitted my first job, which was to interleave the two data files with forward and reverse reads. Each of the data files, from the cultured thermophiles are 26 Gb.

Software used to interleave Cen13_05102014Pooled_ctDNA_GCCAAT_L002_R1_001.fastq and Cen13_05102014Pooled_ctDNA_GCCAAT_L002_R2_001.fastq:
khmer

Job script for interleave:
```
#! /bin/bash

#PBS -l walltime=8:00:00
#PBS -l nodes=1:ppn=8
#PBS -l mem=100 Gb
#PBS -e /mnt/ls15/scratch/users/f0002184/
#PBS -o /mnt/ls15/scratch/users/f0002184/
#PBS -N Interleave_Reads
#PBS -M jlee4946@gmail.com
#PBS -m abe

module load GNU/4.8.2
module load khmer/2.0
cd /mnt/ls15/scratch/users/f0002184/
interleave-reads.py Cen13_05102014Pooled_ctDNA_GCCAAT_L002_R1_001.fastq Cen13_05102014Pooled_ctDNA_GCCAAT_L002_R2_001.fastq -o combined_reads.fastq
```
This job script, called interleave.qsub can be found in /mnt/ls15/scratch/users/f0002184/

The results were one file, named combined_reads.fastq, is in the same directory with 51 Gb.

## Quality Control of the Reads
Once I got my combined reads file, I used software fastx in order to quality control and eliminate reads with 51% or greater q-scores less than 30 (meaning 1 in 1000 errors).

Job script for fastx:
```
#! /bin/bash    

#PBS -l walltime=12:00:00
#PBS -l nodes=1:ppn=8
#PBS -l mem=100Gb
#PBS -N Quality_Filtering_Culture
#PBS -e /mnt/ls15/scratch/users/f0002184/
#PBS -o /mnt/ls15/scratch/users/f0002184/
#PBS -M jlee4946@gmail.com
#PBS -m abe

module load fastx
cd /mnt/ls15/scratch/users/f0002184/
fastq_quality_filter -Q33 -q 30 -p 50 -i combined_reads.fastq -o combined_filtered.fastq
```
The job script can be found in the same directory. Its name is quality_filter.qsub. The results from this job is a FASTQ file, called combined_filtered.fastq. This file size is 42 Gb, so about 9 Gb of data were eliminated.

## Extracting Match-Paired Reads
Continuing from the filtered data, I had to filter the data more so that I have a file with only paired reads, since filtering the reads by q-scores may have gotten rid of a reverse read or a forward read without eliminating its pair. I did this by submitting the job script below, called Match_Paired_Culture.qsub in the same directory.

Job script:
```
#! /bin/bash

#PBS -l walltime=8:00:00
#PBS -l mem=100Gb
#PBS -l nodes=1:ppn=8
#PBS -e /mnt/ls15/scratch/users/f0002184/
#PBS -o /mnt/ls15/scratch/users/f0002184/
#PBS -N Match_Reads_Culture
#PBS -M jlee4946@gmail.com
#PBS -m abe
module load GNU/4.8.2
module load khmer/2.0
cd /mnt/ls15/scratch/users/f0002184/
extract-paired-reads.py combined_filtered.fastq
```

The result of this job is two files, combined_filtered.fastq.pe (pe for paired end) of 34 Gb and combined_filtered.fastq.se (se for single end) of 7.4 Gb. I will only use the paired end file.

## Assembling the Paired Reads
Lastly, I used the paired end data file in order to assemble the reads into contigs by submitting another job. I wrote a script for this, called megahit.qsub in the same directory, found below.

```
#! /bin/bash

#PBS -l walltime=07:00:00:00
#PBS -l mem=400Gb
#PBS -l nodes=1:ppn=12
#PBS -e /mnt/ls15/scratch/users/f0002184/
#PBS -o /mnt/ls15/scratch/users/f0002184/
#PBS -N Megahit_QC_Culture
#PBS -M jlee4946@gmail.com
#PBS -m abe
module load megahit
cd /mnt/ls15/scratch/users/f0002184/
megahit --12 combined_filtered.fastq.pe --k-list 27,37,47,57,67,77,87,97,107 -o Megahit_QC_Assembly/ -t $PBS_NUM_PPN
```
The job ran overnight and it was done around noon. The result is a directory, called Megahit_QC_Assembly. Within that directory, there is a .fa file called final.contigs.fa with the assembled contigs. The size of this file is 1.1 Gb. I made a copy of this file just in case and copied it into my original directory.
___
## Mapping
There are a couple steps to mapping the reads from the other sites to the assembled contigs. Primarily I need to index the contigs (the Mega-Assembly or MA) and then map the reads onto them.

Since I have the final.contigs.fa file in my original directory, I made a new directory for the sole purpose of mapping, called MAPPING_MEGA_ASSEMBLY which can be found in /mnt/ls15/scratch/users/f0002184/ and copied the final contigs file into it.
```
cp final.contigs.fa /mnt/ls15/scratch/users/f0002184/MAPPING_MEGA_ASSEMBLY
```

### Indexing the Contigs
I used a job to do this. It was very speedy! Only 3 minutes or so.
```
#! /bin/bash

#PBS -l walltime=2:00:00
#PBS -l mem=250 Gb
#PBS -l nodes=1:ppn=12
#PBS -e /mnt/scratch/ls15/users/f0002184/MAPPING_MEGA_ASSEMBLY
#PBS -o /mnt/scratch/ls15/users/f0002184/MAPPING_MEGA_ASSEMBLY
#PBS -N index_mega_assembly
#PBS -M jlee4946@gmail.com
#PBS -m abe

module load bbmap
cd /mnt/ls15/scratch/users/f0002184/MAPPING_MEGA_ASSEMBLY
bbmap.sh ref=final.contigs.fa build=1 -Xmx215g
```
This job is called index_assembled_contigs.qsub in the MAPPING_MEGA_ASSEMBLY directory.
It created a directory called ref, in which there are two directories, genome and index, in which there is a directory called 1 in each of them, where there are files for chrom1-3, and in /genome/1/ there is an info.txt file and a summary.txt file. The entire ref directory is 33 Gb.

## Mapping Reads from Cen01, Cen03, Cen04, Cen05, Cen06, Cen07, Cen10, Cen12, Cen14, Cen15, Cen16, and Cen17
I copied the reads from all the other Centralia sites from /mnt/research/ShadeLab/Sorensen/JGI_Metagenomes/ into the MAPPING_MEGA_ASSEMBLY directory. They all vary in size.

I then wrote job scripts for each of those sites to map the reads, the one for Cen01 named map_Cen01_MA.qsub accordingly, the one for Cen02 named map_Cen02_MA.qsub etc. I gave 2 days for map_Cen01_MA.qsub and 7 days for every other one, submitted between 13:22 and 17:00 on 5 June 2017. The job script for Cen01 looks like this:
```
#! /bin/bash

#PBS -l walltime=48:00:00
#PBS -l mem=250Gb
#PBS -l nodes=1:ppn=8
#PBS -e /mnt/ls15/scratch/users/f0002184/MAPPING_MEGA_ASSEMBLY
#PBS -o /mnt/ls15/scratch/users/f0002184/MAPPING_MEGA_ASSEMBLY
#PBS -N map_Cen01_MA
#PBS -M jlee4946@gmail.com
#PBS -m abe

module load bbmap
cd /mnt/ls15/scratch/users/f0002184/MAPPING_MEGA_ASSEMBLY
bbmap.sh in=/mnt/ls15/scratch/users/f0002184/MAPPING_MEGA_ASSEMBLY/Cen01.anqdp.fastq build=1 -Xmx215g out=Cen01_MA.sam
```
