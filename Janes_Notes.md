# Jane's Notes
This project takes quality controlled sequenced DNA from a cultured pool and assembles the reads using the khmer package and metahit respectively. Reads from 12 other sites are mapped against the assembled contigs from the cultured pool and respective abundances are estimated. The contigs are then sorted and binned based on abundance patterns as well as tetranucleotide sequences using MetaBAT. Using CheckM, the bins are quality controlled and taxonomy determined. Response patterns are noted and gene content annotated in order to determine any gene loss and/or gain.

## Table of Contents
* [Interleaving Reads](https://github.com/ShadeLab/CulturedThermophileMAGs/blob/master/Janes_Notes.md#interleaving-reads)
* [Quality Control](https://github.com/ShadeLab/CulturedThermophileMAGs/blob/master/Janes_Notes.md#quality-control-of-the-reads)
* [Extracting Match Pair Reads](https://github.com/ShadeLab/CulturedThermophileMAGs/blob/master/Janes_Notes.md#extracting-match-paired-reads)
* [Assembling Paired Reads](https://github.com/ShadeLab/CulturedThermophileMAGs/blob/master/Janes_Notes.md#assembling-the-paired-reads)
* [MetaQUAST](https://github.com/ShadeLab/CulturedThermophileMAGs/blob/master/Janes_Notes.md#metaquast)
* [Mapping](https://github.com/ShadeLab/CulturedThermophileMAGs/blob/master/Janes_Notes.md#mapping)
* [Indexing the Contigs](https://github.com/ShadeLab/CulturedThermophileMAGs/blob/master/Janes_Notes.md#indexing-the-contigs)
* [Mapping Reads](https://github.com/ShadeLab/CulturedThermophileMAGs/blob/master/Janes_Notes.md#mapping-reads-from-cen01-cen03-cen04-cen05-cen06-cen07-cen10-cen12-cen14-cen15-cen16-and-cen17)
* [Converting .sam to .bam](https://github.com/ShadeLab/CulturedThermophileMAGs/blob/master/Janes_Notes.md#converting-sam-to-bam)
* [Binning](https://github.com/ShadeLab/CulturedThermophileMAGs/blob/master/Janes_Notes.md#binning)
* [CheckM](https://github.com/ShadeLab/CulturedThermophileMAGs/blob/master/Janes_Notes.md#checkm)

#### 30 May 2017 - 6 June 2017
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

The results were one file, named combined_reads.fastq, and are in the same directory with 51 Gb.

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

## MetaQUAST
I used MetaQuast to get summary statistics on the assembled contigs using MetaQUAST version 2.3. I created a new directory called quast-2.3, where I loaded and executed the program interactively.
```
wget https://downloads.sourceforge.net/project/quast/quast-2.3.tar.gz
tar xzvf quast-2.3.tar.gz
cd quast-2.3
```

After this, I ran the program with the final.contigs.fa file.
```
python metaquast.py -o /mnt/ls15/scratch/users/f0002184/quast-2.3 final.contigs.fa
```
I don't need a job for this, but it does take a couple minutes. This resulted in several directories, of which one is quast_output. Once I cd into that, there is a text file named report.txt. When I use more to see the text file, there is a nice summary of my contigs!

|Assembly | final.contigs|
| --------|:------------|
|# contigs (>= 0 bp) | 508410|
|# contigs (>= 1000 bp) | 222766|
|Total length (>= 0 bp) | 898424779|
|Total length (>= 1000 bp) | 698652595|
|# contigs | 508410|
|Largest contig | 413537|       
|Total length | 898424779|    
|GC (%) | 63.90|   
|N50 | 2651|       
|N75 | 1087|        
|L50 | 59145|        
|L75 | 198925|       
|# N's per 100 kbp | 0.00|
___
## Mapping
There are a couple steps to mapping the reads from the other sites to the assembled contigs. Primarily I need to index the contigs (the Mega-Assembly or MA) and then map the reads onto them.

Since I have the final.contigs.fa file in my original directory, I made a new directory for the sole purpose of mapping, called MAPPING_MEGA_ASSEMBLY which can be found in /mnt/ls15/scratch/users/f0002184/ and copied the final contigs file into it.
```
cp final.contigs.fa /mnt/ls15/scratch/users/f0002184/MAPPING_MEGA_ASSEMBLY
```

## Indexing the Contigs
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

## Converting .sam to .bam
#### 6 June 2017
Now that the jobs are submitted and are in queue, I think it's going to take a few days. I've taken the liberty to go ahead and write some more job scripts in preparation. Once these jobs are done, I'm going to have .sam files that will need to be converted to .bam files in order for MetaBAT to bin these genomes. I will use SAMTools/1.3 to do this. Here is an example of one job script, but I will have to do this for every sample site.
```
#! /bin/bash

#PBS -l walltime=10:00:00
#PBS -l mem=100Gb
#PBS -l nodes=1:ppn=8
#PBS -e /mnt/ls15/scratch/users/f0002184
#PBS -o /mnt/ls15/scratch/users/f0002184
#PBS -N Sam_to_Bam_Cen01
#PBS -M jlee4946@gmail.com
#PBS -m abe

module load GNU/4.8.3
module load SAMTools/1.3
cd /mnt/ls15/scratch/users/f0002184
samtools view -bS Cen01_MA.sam > Cen01_MA.bam
samtools sort -o Cen01_MA.bam -T Cen01_Sorted -@ 8 -m 8G Cen01_MA.bam
```
This job is called sam_to_bam_Cen01.qsub in /mnt/ls15/scratch/users/f0002184.

I ran into my first small and rather silly problem here. I originally had put 100 Gb in my script, but clearly that didn't work, so keep in mind, there are no spaces between them!!!!!!! As seen above, it's 100Gb, NOT 100 Gb! Also, converting these files only took around 3 hours each.

#### 7 June 2017
map_Cen01_MA begun execution and finished! It took about 5 hours, which is much shorter than what I was expecting. It created a new file, Cen01_MA.sam in /mnt/ls15/scratch/users/f0002184/MAPPING_MEGA_ASSEMBLY. I submitted the job to convert it to a bam file (s/o to Jackson for updating SAMTools yay!). I created a new directory for only .bam files in the MAPPING_MEGA_ASSEMBLY directory called BAM_Files.

Still waiting on all the other .sam files though.

#### 8 June 2017
map_Cen03_MA, map_Cen05_MA, map_Cen06_MA, map_Cen07_MA, map_Cen12_MA, map_Cen15_MA are done! Submitted the jobs to convert them into .bam files. Had some issues regarding memory and changed the walltime from 1 hour to 1 day to 10 hours to 5 hours but other than that, everything is good! Other mapping jobs began to run! Yay!

The .sam files seem to be around 50-60Gb.

Cen01_MA.bam (13 Gb), Cen03_MA.bam (15 Gb), Cen05_MA.bam (15 Gb) finished converting from .sam to .bam! I copied them into my BAM_Files directory.

## Binning
#### 9 June 2017
Everything has been converted to BAM files and copied to my BAM_Files directory!

I indexed the files in BAM_Files in order for MetaBAT to work using
```
module load GNU/4.8.3
module load SAMTools/1.3
samtools index -b BAM_files/*
```
in my bash terminal under the MAPPING_MEGA_ASSEMBLY directory.

Now I will create a depth file.

This is where I ran into my first BIG PROBLEM! After trying to figure out this ominous error message:
```
[E::hts_open] fail to open file
```
for ETERNITY, turns out my Cen03_MA.bam file was only something like 15M instead of like 15 Gb like the other .bam files (s/o to Jackson for coming in clutch!!). I re-ran my job script to convert Cen03_MA.sam to Cen03_MA.bam but then I looked at my original Cen03_MA.bam which was 15 Gb so I guess something went wrong copying the original to the BAM_Files directory. Either way, I am now waiting for the depth file! I did everything in the command line, although as of now I'm kind of regretting my decision and wondering if I should have submitted a job. Anyway, I entered:
```
cd /mnt/ls15/scratch/users/f0002184/MAPPING_MEGA_ASSEMBLY/BAM_Files
tmux new -s METABAT
module load GNU/4.8.3
module load MetaBAT/20160622
jgi_summarize_bam_contig_depths --outputDepth depth.txt *.bam
```
This took a couple hours so maybe I should've submitted a job but it's done now!
I made a Genome_Binning directory under BAM_Files and put the result of this, depth.txt in it. depth.txt is 87 Mb.

I submitted my binning job!!! I hope it works smoothly, but I'm glad it has the weekend to run!
Here is the script, called binning.qsub in /mnt/ls15/scratch/users/f0002184/MAPPING_MEGA_ASSEMBLY/BAM_Files/Genome_Binning.
```
#! /bin/bash

#PBS -l walltime=48:00:00
#PBS -l mem=200Gb
#PBS -l nodes=2:ppn=20
#PBS -l feature='intel14|intel16'
#PBS -e /mnt/ls15/scratch/users/f0002184/MAPPING_MEGA_ASSEMBLY/BAM_Files/Genome_Binning
#PBS -o /mnt/ls15/scratch/users/f0002184/MAPPING_MEGA_ASSEMBLY/BAM_Files/Genome_Binning
#PBS -N Binning_Genomes
#PBS -M jlee4946@gmail.com
#PBS -m abe

module load GNU/4.8.3
module load MetaBAT/20160622
cd /mnt/ls15/scratch/users/f0002184/MAPPING_MEGA_ASSEMBLY/Genome_Binning
metabat -i cd /mnt/ls15/scratch/users/f0002184/MAPPING_MEGA_ASSEMBLY/final.contigs.fa -v -a depth.txt -o METABAT_VerySpecific --saveTNF saved.tnf --saveDistance saved.dist -t 40 --veryspecific
```

As you can see above, it says #PBS -l nodes=2! 2 nodes!!! #PBS -l feature='intel14|intel16' specifies which nodes should work with the job. Could it be...COULD IT BE...PARALLEL COMPUTING!?!

I submitted an additional binning job, deservingly called binning_2.qsub using MetaBAT/20170607. The script is the same, other than the name changes and
```
module load MetaBAT/20170607
```
#### 12 June 2017
My binning job aborted! I'm thinking it's because of the 'cd' in metabat -i cd /mnt/ls15/scratch/users/f0002184/MAPPING_MEGA_ASSEMBLY/final.contigs.fa so I deleted that and re-submitted the job.

binning_2 still hadn't executed so I changed it in that job script as well.

Also changed ppn=40 from ppn=20.

Also check out this cool [link](https://wiki.gacrc.uga.edu/wiki/MetaBAT-Sapelo)
with some tips!

So I tried running the commands from the script separately just to see if they worked, and they did, and only gave me 4 "bins". Seemed weird, but apparently it's okay! So now I have 4 bins, all in a separate directory called Genome_Bins, called METABAT_VerySpecific_Trial.1.fa - METABAT_VerySpecific_Trial.4.fa. They are between 2.0 Mb and 4.8 Mb. I used the 2016 MetaBAT version for this.

## CheckM
#### 13 June 2017
Lots of progress today!
Spent a lot time yesterday just trying to install CheckM so that I could quality control and determine taxonomy, but apparently installing programs on the HPCC is much more complicated than I would have ever imagined.
This morning, I struggled a bit with the CheckM installation again. I tried updating it and it didn't work, but maybe it finished updating eventually or something because it worked even though I didn't change anything some time later (s/o to Pat at the HPCC center & Jackson for trying to troubleshoot what was wrong even though apparently nothing was).
One thing to note is remembering to load HMMER, pplacer and prodigal using module load, something like this:
```
module load HMMER
module load pplacer
module load prodigal
```

So now that I got CheckM to work, I did everything in the command line. It took a few minutes but nothing too bad. Here are the commands I used for the very specific bins for CheckM:
```
checkm lineage_wf -t 8 -x fa Genome_Bins/ CheckM/
```
I had to make a new directory for the program to write into, CheckM. This directory must be completely empty in order for it to work.
I then put the summary into a file so I could keep it:
```
checkm qa CheckM/lineage.ms CheckM/ > CheckM_Results.txt
```
The result text file has the bin ID, the marker lineage, # genomes, # markers, # marker sets, 0, 1, 2, 3, 4, 5+, Completeness, Contamination and Strain heterogeneity. For my 4 bins, it looks like completeness is between 54.53 and 94.94, which I think (?) is pretty good? All contamination values are below 10, which is good too.

However, I ran the MetaBAT job script again with the 2017 version, using the very sensitive tag instead of very specific, because it's better for more simple communities.
Here are the commands:
```
module load GNU/4.8.3
module load MetaBAT/20170607
cd /mnt/ls15/scratch/users/f0002184/MAPPING_MEGA_ASSEMBLY/BAM_Files/Genome_Binning
metabat -i /mnt/ls15/scratch/users/f0002184/MAPPING_MEGA_ASSEMBLY/final.contigs.fa -v -a depth.txt -o METABAT_VerySensitive --saveTNF saved.tnf --saveDistance saved.dist -t 40 --verysensitive
```
This gave me 12 bins! I ran CheckM on this as well, and my results are in /mnt/ls15/scratch/users/f0002184/MAPPING_MEGA_ASSEMBLY/BAM_Files/Genome_Binning/CheckM_Sensitive.
These bins have much more varying completenesses, but other than 2 contamination values of 13.56 and 29.95, they're below 6.36! Although to be fair, the completeness level is lower on this own than on the very specific bins.

So I thought I would try using the 2017 version of MetaBAT with the very specific bins. I got 5 bins this time! The only different thing between the CheckM results from the 2016 MetaBAT and the 2017 MetaBAT seems to be well 1) there is one more bin from the 2017 one and 2) that bin is root(U1D1) (?) and its completeness is only 4.17 with 0 contamination.

#### 14 June 2017
Today I downloaded the very specific 2016 CheckM result file and after much struggling to get it formatted correctly, I had it in a .csv format in Excel before I could read it into R. Although since for this one there is only 4 bins I could just look at it to sort the bins based on the thresholds but it was a good exercise to get familiar with R. Here is what I put in R:
```
setwd("/Users/janelee/Documents/MSU_REU/CheckM")
CheckM_Specific <- read.csv("CheckM_Results.csv")
CheckM_Specific <- as.data.frame(CheckM_Specific)
CheckM_Specific
Clean_Contam <- CheckM_Specific[CheckM_Specific$Contamination<5,]
Clean_Contam
Clean_Complete <- CheckM_Specific[CheckM_Specific$Completeness>90,]
Clean_Complete
Clean <- Clean_Contam[Clean_Contam$Completeness>90,]
Clean
```
I get one bin each for when contamination < 5 and for when completeness > 90, but no bins for both contamination < 5 and completeness > 90. Since there are only 4 bins I could have easily just looked at it as well. The most complete bin has a contamination value of 6.37 while the least contaminated bin has a completeness value of 67.18.

I read my very sensitive CheckM file into R and did the formatting a little differently which made it much easier.
```
setwd("/Users/janelee/Documents/MSU_REU/CheckM")
CheckM_Sensitive <- read.table("CheckM_VerySensitive_Results.txt", header=FALSE, blank.lines.skip = TRUE)
CheckM_Sensitive <- as.data.frame(CheckM_Sensitive)
colnames(CheckM_Sensitive)<- c("BinID", "MarkerLineage", "#Genomes", "#Markers", "#MarkerSets", "0","1", "2", "3", "4", "5+", "Completeness", "Contamination", "StrainHeterogeneity")
CheckM_Sensitive
Clean_Contamination <- CheckM_Sensitive[CheckM_Sensitive$Contamination<5,]
Clean_Contamination
Clean_Complete <- CheckM_Sensitive[CheckM_Sensitive$Completeness>90,]
Clean_Complete
```
For these bins, 9 of them have a contamination level less that 5, but the highest completeness is 69.03, and the second highest completeness is 20.61. There was only one bin that had a completeness level greater than 90, at 95.4. The contamination is 6.36. The Marker Lineage, k__Bacteria(UID3187), is the same for this bin as the bin above from the very specific bins. k__Bacteria(UID3187), under very specific and very sensitive flags, has a > 90 completeness value and ~6.36 contamination value in both bins.

## MiGA
#### 14 June 2017
I downloaded the most complete genome of the very specific genomes (k__Bacteria(UID3187)) and uploaded it to MiGA's NCBI Prok to "query the NCBI Genome database (Prokaryotes) with your own genomes to determine the most likely taxonomic classification and novelty rank" ([MiGA](http://enve-omics.ce.gatech.edu:3000/)). MiGA is still currently running and I'm not sure how long it's going to take.

#### 15 June 2017
Last night I uploaded the remaining 3 bins with the very specific tag, and they were all finished by this morning. MiGA gives you the closest relative in the MiGA database, the taxonomy of the genome, starting from domain down to subspecies, with p-values. It additionally gives you a completeness %, contamination %, a quality score # of predicted proteins, average length (amino acids), total number of contigs, and total length (in base pairs). 
