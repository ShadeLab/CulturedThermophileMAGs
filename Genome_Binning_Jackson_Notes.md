# Binning Genomes from assembled metagenomes.

## Software
#### bbmap
Creates the abundance data
#### samTools
Converts abundance data into useable form for metaBAT
#### metaBAT
Bins the contigs into MAGs
#### CheckM
Checks the quality of the MAGs

##### May 4th, 2016
Previously I worked to create a large MEGA-assembly that was made from the metagenome reads from all 12 sites JGI sequenced for the CSP (See other workflow for details). Here I am attempting to bin out metagenome assembled genomes (MAGs) from this MEGA-assembly. To do this I will use the software metaBAT. metaBAT works using abundance data and tetranucleotide frequency of the assembled contigs to group the contigs into these MAGs. Before I can use metaBAT, I need to determine the abundance of each contig in each site using the tool bbmap.

### Creating abundance data for contigs using bbmap.
Determining the abundance of each contigs across the chronosequence requires mapping the reads from each site to the contigs from the MEGA-assembly. This process has 2 steps. 1) The reference(in the case the MEGA-assembly) has to be indexed. 2) Reads from each site are mapped to the indexed contigs.

Copied the final.contigs.fa from the MEGA-assembly into a new directory for the mapping analysis
```
mkdir MAPPING_MEGA_ASSEMBLY
cp ../MEGA_ASSEMBLY/Megahit_QC_Assembly/final.contigs.fa MA_Contigs.fa
```

**The following is a qsub script for indexing the assembled contigs**


```
#! /bin/bash
#PBS -l walltime=04:00:00
#PBS -l mem=250Gb
#PBS -l nodes=1:ppn=12
#PBS -e /mnt/scratch/sorens75/MAPPING_MEGA_ASSEMBLY
#PBS -o /mnt/scratch/sorens75/MAPPING_MEGA_ASSEMBLY
#PBS -N index_MA
#PBS -A bicep
#PBS -M jacksonwsorensen@gmail.com
#PBS -m abe
module load bbmap
cd /mnt/scratch/sorens75/MAPPING_MEGA_ASSEMBLY
bbmap.sh ref=MA_Contigs.fa build=1 -Xmx215g
```
This script only took 50 minutes to run and worked perfectly!

The next step is to map each site's reads to the contigs to create the abundance data.

**The following is a qsub script for mapping the reads from Cen01 to the MEGA-assembly**

```
#! /bin/bash
#PBS -l walltime=08:00:00
#PBS -l mem=250Gb
#PBS -l nodes=1:ppn=12
#PBS -e /mnt/scratch/sorens75/MAPPING_MEGA_ASSEMBLY
#PBS -o /mnt/scratch/sorens75/MAPPING_MEGA_ASSEMBLY
#PBS -N map_Cen01_MA
#PBS -A bicep
#PBS -M jacksonwsorensen@gmail.com
#PBS -m abe
module load bbmap
cd /mnt/scratch/sorens75/MAPPING_MEGA_ASSEMBLY
bbmap.sh in=/mnt/research/ShadeLab/Sorensen/Cen01.anqdp.fastq.gz build=1 -Xmx215g out=Cen01_MA.sam
```

##### May 9th, 2016
This job failed because it exceeded the walltime. Resubmitted with the following job script. The only differences are a longer walltime and less ppn.
```
#! /bin/bash
#PBS -l walltime=48:00:00
#PBS -l mem=250Gb
#PBS -l nodes=1:ppn=8
#PBS -e /mnt/scratch/sorens75/MAPPING_MEGA_ASSEMBLY
#PBS -o /mnt/scratch/sorens75/MAPPING_MEGA_ASSEMBLY
#PBS -N map_Cen01_MA
#PBS -A bicep
#PBS -M jacksonwsorensen@gmail.com
#PBS -m abe
module load bbmap
cd /mnt/scratch/sorens75/MAPPING_MEGA_ASSEMBLY
bbmap.sh in=/mnt/research/ShadeLab/Sorensen/Cen01.anqdp.fastq.gz build=1 -Xmx215g out=Cen01_MA.sam
```
##### May 12th, 2016
**And it worked!** The output sam file `Cen01_MA.sam` is 65GB in size(Pretty effin big). The mapping only took ~16 hours. After this succeeded I submitted jobs to map the rest of the JGI mgDNA reads to our MEGA-assembly. **Of note here is that the default minimum identity of mapped reads is 0.76 by default** It maybe be worthwhile to change this in the future.
##### Update May 16th, 2016
Cen03 and Cen04 finished mapping without issue.

Cen05 exceeded walltime, hrmmm. Clearly needs more time. Resubmitted with two days walltime.

##### Update May 18th, 2016
Cen06 completed without issue.

Cen07, 10, 12, 14, 15, 16, 17 all exceeded the requested walltime. Upped the walltime to 3 days and resubmitted each job.  

##### Update May 24th, 2016
Mapping job for each sample finished in the allotted amount of walltime.   

### Converting sam to bam
Alright, now that we have sam files (well one at the moment but we should have the rest within a week or two), let's get started converting these to a useful format for metaBAT to use. metaBAT needs ordered bam files in order to bin contigs. So the next couple steps are going to be converting these sam files into bam files.

Working with the `Cen01_MA.sam` file, let's try to get it into a bam format. Maybe we can get lucky and not have to submit a job. (I have a bad feeling about this)
```
module load samTools
samtools view -bS Cen01_MA.sam > Cen01_MA.bam
```
After about 15 minutes I'm thinking doing this interactively isn't the best way. Let's get a job going. Also it looks like you can convert sam straight to an order bam file so we can make use of the job script to do this as well.
```
#! /bin/bash
#PBS -l walltime=12:00:00
#PBS -l mem=100Gb
#PBS -l nodes=1:ppn=8
#PBS -e /mnt/scratch/sorens75/MAPPING_MEGA_ASSEMBLY
#PBS -o /mnt/scratch/sorens75/MAPPING_MEGA_ASSEMBLY
#PBS -N Cen01_Sam_to_Bam
#PBS -A bicep
#PBS -M jacksonwsorensen@gmail.com
#PBS -m abe
module load samTools
cd /mnt/scratch/sorens75/MAPPING_MEGA_ASSEMBLY
samtools view -bS Cen01_MA.sam | samtools sort - -o Cen01_MA.bam
```
##### May 15th 2016
Converting SAM to BAM was unsuccessfuly. Likely just not enough walltime. Error output follows
```
Lmod Warning: samTools not found, loading: SAMTools/1.2
[bam_sort_core] merging from 114 files...
=>> PBS: job killed: walltime 43216 exceeded limit 43200
Unknown option: 1
Unknown option: 1
```
Resubmitted the job with greater walltime requested.

```
#! /bin/bash
#PBS -l walltime=48:00:00
#PBS -l mem=100Gb
#PBS -l nodes=1:ppn=8
#PBS -e /mnt/scratch/sorens75/MAPPING_MEGA_ASSEMBLY
#PBS -o /mnt/scratch/sorens75/MAPPING_MEGA_ASSEMBLY
#PBS -N Cen01_Sam_to_Bam
#PBS -A bicep
#PBS -M jacksonwsorensen@gmail.com
#PBS -m abe
module load samTools
cd /mnt/scratch/sorens75/MAPPING_MEGA_ASSEMBLY
samtools view -bS Cen01_MA.sam | samtools sort - -o Cen01_MA.bam
```

##### May 23rd, 2016
Converting Cen01_MA.sam to a bam file still did not work. Eror output is below.
```
Lmod Warning: samTools not found, loading: SAMTools/1.2
[bam_sort_core] merging from 114 files...
=>> PBS: job killed: walltime 172831 exceeded limit 172800
Unknown option: 1
Unknown option: 1
```
Looks like the walltime again was not enough. Also, looking at my job script I think their may have been an error in the command. `samtools view -bS Cen01_MA.sam | samtools sort - -o Cen01_MA.bam`. Resubmitted two jobs, one will older command and one with a modified command. Also increased the memory for both of them.  

Same command more walltime
```
#! /bin/bash
#PBS -l walltime=96:00:00
#PBS -l mem=200Gb
#PBS -l nodes=1:ppn=8
#PBS -e /mnt/scratch/sorens75/MAPPING_MEGA_ASSEMBLY
#PBS -o /mnt/scratch/sorens75/MAPPING_MEGA_ASSEMBLY
#PBS -N Cen01_Sam_to_Bam
#PBS -A bicep
#PBS -M jacksonwsorensen@gmail.com
#PBS -m abe
module load samTools
cd /mnt/scratch/sorens75/MAPPING_MEGA_ASSEMBLY
samtools view -bS Cen01_MA.sam | samtools sort -o Cen01_MA_sorted.bam
```

New command more walltime
```
#! /bin/bash
#PBS -l walltime=96:00:00
#PBS -l mem=200Gb
#PBS -l nodes=1:ppn=8
#PBS -e /mnt/scratch/sorens75/MAPPING_MEGA_ASSEMBLY
#PBS -o /mnt/scratch/sorens75/MAPPING_MEGA_ASSEMBLY
#PBS -N Cen01_Sam_to_Bam_Exact
#PBS -A bicep
#PBS -M jacksonwsorensen@gmail.com
#PBS -m abe
module load samTools
cd /mnt/scratch/sorens75/MAPPING_MEGA_ASSEMBLY/Cen01_Exact
samtools view -bS ../Cen01_MA.sam | samtools sort - Cen01_MA_sorted
```
##### May 24th, 2016
Okay so I'm getting tired of this nonsense with converting sam to bam and making a sorted bam file... After resubmitting with 4 days of walltime yesterday and seeing that it would start for another week I looked into seeing if there is anyway to speed up the way the program is running. Looks like the program defaults to only using 716 MB per thread and only 1 thread. I have been requesting a large amount of resources but they haven't been making use of them. Okay let's try to actually make use of these resources this time. Requested less walltime for this job because it should theorectically be running at least 8 times faster than before. This has the added bonus of getting my job through the queue quicker. Specified the number of threads with `-@` and the memory per thread with `-m`.  

```
#! /bin/bash
#PBS -l walltime=06:00:00
#PBS -l mem=100Gb
#PBS -l nodes=1:ppn=8
#PBS -e /mnt/scratch/sorens75/MAPPING_MEGA_ASSEMBLY
#PBS -o /mnt/scratch/sorens75/MAPPING_MEGA_ASSEMBLY
#PBS -N Cen03_Sam_to_Bam
#PBS -A bicep
#PBS -M jacksonwsorensen@gmail.com
#PBS -m abe
module load samTools
cd /mnt/scratch/sorens75/MAPPING_MEGA_ASSEMBLY
samtools view -bS Cen03_MA.sam | samtools sort -o Cen03_MA.bam -@ 8 -m 8G
```
##### June 30th, 2016

So after much reflection and fighting with HPCC and SAMTools I HAVE CREATED SORTED BAM FILES FROM SAM FILES! It turns out SAMTools/1.2 has a serious coding inefficiency that makes sorting and merging bam files take way too long (would not finish even with 7 days of walltime on hpcc) [See Here](https://github.com/samtools/samtools/pull/337). Asked HPCC to update to SAMTools 1.3 and now everything works. Sample qsub script listed below.

```
#! /bin/bash
#PBS -l walltime=24:00:00
#PBS -l mem=100Gb
#PBS -l nodes=1:ppn=8
#PBS -e /mnt/scratch/sorens75/MAPPING_MEGA_ASSEMBLY
#PBS -o /mnt/scratch/sorens75/MAPPING_MEGA_ASSEMBLY
#PBS -N Cen01_Sam_to_Bam
#PBS -A bicep
#PBS -M jacksonwsorensen@gmail.com
#PBS -m abe
module load GNU/4.8.3
module load SAMTools/1.3
cd /mnt/scratch/sorens75/MAPPING_MEGA_ASSEMBLY
samtools view -bS Cen01_MA.sam > Cen01_MA.bam
samtools sort -o Cen01_MA.bam -T Cen01_Sorted -@ 8 -m 8G Cen01_MA.bam
```

Once done, copied all sorted BAM files to research space `/mnt/research/ShadeLab/Sorensen/MEGA_ASSEMBLY/Mapping_Mega_Assembly/Sorted_BAM_files`

##### July 2nd, 2016
It also turns out that in order for metabat to use these files they must be indexed as well.
```
module load GNU/4.8.3
module load SAMTools/1.3
samtools index -b Sorted_BAM_files/*
```

##### July 5th, 2016
I finally have all the pieces to perform genome binning using metabat. The version of metabat on hpcc is a little bit old(Oct-10-2014, current release is June-22-2016), but I have asked them to update through a request form.

First, I need to create a depth file so that metabat can use the relative abundance informations of the different contigs.

```
tmux new -s METABAT
module load GNU/5.2
module load metabat
cd /mnt/research/ShadeLab/Sorensen/MEGA_ASSEMBLY/
mkdir Genome_Binning
jgi_summarize_bam_contig_depths --outputDepth Genome_Binning/depth_20141007.txt  Mapping_Mega_Assembly/Sorted_BAM_files/*.bam
```

Cool, it worked using tmux. I tried to do the genome binning without submitting a job but it did not finish in the 2 hours hpcc gives on development nodes.
```
cd Genome_Binning
metabat -i ../Megahit_QC_Assembly/final.contigs.fa -a depth_20141007.txt -o 20141007_bin
```

##### July 6th, 2016

Since it doesn't look like the job can finish in time on a development node I decided to submit a job to do the actual binning.
```
#! /bin/bash
#PBS -l walltime=72:00:00
#PBS -l mem=500Gb
#PBS -l nodes=1:ppn=8
#PBS -e /mnt/research/ShadeLab/WorkingSpace/Sorensen/MEGA_ASSEMBLY/Genome_Binning
#PBS -o /mnt/research/ShadeLab/WorkingSpace/Sorensen/MEGA_ASSEMBLY/Genome_Binning
#PBS -N Binning_Genomes
#PBS -A bicep
#PBS -M jacksonwsorensen@gmail.com
#PBS -m abe
module load GNU/5.2
module load metabat
cd /mnt/research/ShadeLab/WorkingSpace/Sorensen/MEGA_ASSEMBLY/Genome_Binning
metabat -i ../Megahit_QC_Assembly/final.contigs.fa -a depth_20141007.txt -o 20141007_bin
```

And the job started today as well. I'm so lucky.

##### July 9th, 2016
Genome Binning didn't finish in time, and the logfile was pretty uninformative. Let's up it to a week of wall time.

```
#! /bin/bash
#PBS -l walltime=168:00:00
#PBS -l mem=500Gb
#PBS -l nodes=1:ppn=8
#PBS -e /mnt/research/ShadeLab/WorkingSpace/Sorensen/MEGA_ASSEMBLY/Genome_Binning
#PBS -o /mnt/research/ShadeLab/WorkingSpace/Sorensen/MEGA_ASSEMBLY/Genome_Binning
#PBS -N Binning_Genomes
#PBS -A bicep
#PBS -M jacksonwsorensen@gmail.com
#PBS -m abe
module load GNU/5.2
module load metabat
cd /mnt/research/ShadeLab/WorkingSpace/Sorensen/MEGA_ASSEMBLY/Genome_Binning
metabat -i ../Megahit_QC_Assembly/final.contigs.fa -a depth_20141007.txt -o 20141007_bin
```
Job started on July 11th.

##### July 18th, 2016
Binning did not finish after a week of walltime. Job was killed for exceeding the maximum allowed walltime. Meet with HPCC Thursday to see if there are anyways to improve the script.

##### July 21st, 2016
Talked with folks at HPCC they suggested trying to increase the number of processors and also suggested that the amount of RAM we are using is probably a bit excessive. In addition, they also taught me how to print out the resources actually used by the job while it ran. You can figure out on average how many processors your job actually used by dividing the cpu hours by the walltime. Below is updated job script.

```
#! /bin/bash
#PBS -l walltime=48:00:00
#PBS -l mem=200Gb
#PBS -l nodes=2:ppn=20
#PBS -e /mnt/research/ShadeLab/WorkingSpace/Sorensen/MEGA_ASSEMBLY/Genome_Binning
#PBS -o /mnt/research/ShadeLab/WorkingSpace/Sorensen/MEGA_ASSEMBLY/Genome_Binning
#PBS -N Binning_Genomes_2016_2N_20PPN_200GB
#PBS -l feature='intel14|intel16'
#PBS -M jacksonwsorensen@gmail.com
#PBS -m abe
module load GNU/4.8.3
module load MetaBAT/20160622
cd /mnt/research/ShadeLab/WorkingSpace/Sorensen/MEGA_ASSEMBLY/Genome_Binning
(sleep $((72*60*60-5*60)); qstat -f $PBS_JOBID )&
metabat -i ../Megahit_QC_Assembly/final.contigs.fa -v -a depth_20141007.txt -o 20141007_bin --saveTNF saved.tnf --saveDistance saved.dist -t 40
qstat -f $PBS_JOBID
```
Cool this worked! By increasing the number of processors for 8 to 40 and specify the number in the script with the `-t` flag I was able to get the binning to complete in less than a day. METABAT identified 1213 separate bins from our MEGA_ASSEMBLY. This is pretty cool. Let's go forward and try to figure out how complete/clean these bins are by using [CheckM](http://ecogenomics.github.io/CheckM/).

#### Let's check the quality of these bins using CheckM

Alright I think we can get away with checking the quality of these bins without submitting a job.

```
mkdir Genome_Bins
mv 20141007_bin.*.fa Genome_Bins
tmux new -s CheckM
source /mnt/research/ShadeLab/Sorensen/software/load_CheckM.sh
checkm lineage_wf -t 8 -x fa Genome_Bins/ CheckM_Output
```

In the above command I started a tmux session and then loaded an anaconda version of CheckM because HPCC's isn't working very well. Unfortunately, because CheckM is weird in its output I needed to re-run the summary of the analysis so that I could have the output in a file.

```
checkm qa CheckM_Output/lineage.ms CheckM_Output/ > CheckM_Results.txt
```
Downloaded the Results file and read it into RStudio to quickly parse through the results. Before importing in R though, there is some small manual alterations that need to be made. Delete the header section and the last line of the table. The result should be a table of genome bins and their different stats. Use this alterred table in R.

```
setwd("~/GitHub_Repos/ShadeLab/2014_MetaG_Assembly/JGI_MEGA_ASSEMBLY/")
MAGS <- read.table("~/GitHub_Repos/ShadeLab/2014_MetaG_Assembly/JGI_MEGA_ASSEMBLY/CheckM_Results.txt", quote="\"", comment.char="")
MAGS <- MAGS[,-3]
colnames(MAGS)<- c("BinID", "MarkerLineage", "#Genomes", "#Markers", "#MarkerSets", "0","1", "2", "3", "4", "5+", "Completeness", "Contamination", "StrainHeterogeneity")

Clean <- MAGS[MAGS$Contamination<5,]
Clean_Complete <- Clean[Clean$Completeness>90,]
```
This binning effort led to 29 Bins that are >90% complete and with <5% contamination. Not a bad effort, but I think we can still do a bit better. Looking at [Bendal & Stevens](http://www.nature.com/ismej/journal/vaop/ncurrent/full/ismej2015241a.html) they use MAGs that are as low as 50% complete as estimated by CheckM. Doesn't look like they use the contamination parameter from CheckM at all. Instead they make sure they use contigs whose abundances throughout their time series have correlation coefficients >.995 with the median contig abundance in the bin.

##### July 23rd, 2016
Let's re-do the binning but this time with the `--veryspecific` flag to potentially clean up some of our bins.
```
#! /bin/bash
#PBS -l walltime=48:00:00
#PBS -l mem=200Gb
#PBS -l nodes=2:ppn=20
#PBS -e /mnt/research/ShadeLab/WorkingSpace/Sorensen/MEGA_ASSEMBLY/Genome_Binning
#PBS -o /mnt/research/ShadeLab/WorkingSpace/Sorensen/MEGA_ASSEMBLY/Genome_Binning
#PBS -N Binning_Genomes_2016_2N_20PPN_200GB
#PBS -l feature='intel14|intel16'
#PBS -M jacksonwsorensen@gmail.com
#PBS -m abe
module load GNU/4.8.3
module load MetaBAT/20160622
cd /mnt/research/ShadeLab/WorkingSpace/Sorensen/MEGA_ASSEMBLY/Genome_Binning
(sleep $((72*60*60-5*60)); qstat -f $PBS_JOBID )&
metabat -i ../Megahit_QC_Assembly/final.contigs.fa -v -a depth_20141007.txt -o METABAT_VerySpecific --saveTNF saved.tnf --saveDistance saved.dist -t 40 --veryspecific
qstat -f $PBS_JOBID
```
##### July 25th, 2016
Cool this worked. Let's go ahead and make use CheckM to assess completeness and contamination.
```
tmux attach -t CheckM
checkm lineage_wf -t 8 -x fa VerySpecific_Bins/ CheckM_VerySpecific_Output
checkm qa CheckM_VerySpecific_Output/lineage.ms CheckM_VerySpecific_Output/ > CheckM_VerySpecific_Results.txt
```
