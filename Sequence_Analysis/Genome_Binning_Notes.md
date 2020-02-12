# Russell Ranch Metagenome Assembled Genomes Workflow
## Rationale/Setup
Russell Ranch is a long term research site set up to investigate the longterm influences of different management practices. Our goal was to understand how these different management practices may influence the local viral communities, and their respective hosts. In July and October of 2018, Laura Zinke, Joanne Emerson, and Sara Geonczy sampled 6 different Tomato/Corn rotation plots, 3 representing that were managed with an Organic practice (Poultry manure for fertlizer + Winter Cover Crop) and 3 that were managed with a Traditional practice (Mineral fertilizer + Fallow winter season). At each plot 2 replicate cores were taken at different spots. From each of these samples (2 management practices * 3 plots per practice * 2 reps per plot * 2 time points == 24 total samples), BULK metagenomes were produced via direct DNA extraction with the Powersoil kit + Illumina PE 150bp sequencing, and Viromes were produced by first physically separating viral particles from soil and microbes, then extracting DNA from these viral particles and sequencings using Illumina PE 150bp sequencing. 6 Viromes failed to produce enough sequences for analysis. All six of these samples were from the July time point and were evenly split between the conventional and organic management practices. The following is a workflow notes for performing metagenome assembly and Genome binning. Below is a list of filenames for all the samples.
```
RR131Jul_Bulk_S47_L001_R1_001.fastq.gz	RR641Jul_Bulk_S58_L001_R1_001.fastq.gz
RR131Jul_Bulk_S47_L001_R2_001.fastq.gz	RR641Jul_Bulk_S58_L001_R2_001.fastq.gz
RR131Oct_Bulk_S12_L001_R1_001.fastq.gz	RR641Oct_Bulk_S50_L001_R1_001.fastq.gz
RR131Oct_Bulk_S12_L001_R2_001.fastq.gz	RR641Oct_Bulk_S50_L001_R2_001.fastq.gz
RR132Jul_Bulk_S49_L001_R1_001.fastq.gz	RR642Jul_Bulk_S45_L001_R1_001.fastq.gz
RR132Jul_Bulk_S49_L001_R2_001.fastq.gz	RR642Jul_Bulk_S45_L001_R2_001.fastq.gz
RR132Oct_Bulk_S20_L001_R1_001.fastq.gz	RR642Oct_Bulk_S51_L001_R1_001.fastq.gz
RR132Oct_Bulk_S20_L001_R2_001.fastq.gz	RR642Oct_Bulk_S51_L001_R2_001.fastq.gz
RR231Jul_Bulk_S41_L001_R1_001.fastq.gz	RR681Jul_Bulk_S43_L001_R1_001.fastq.gz
RR231Jul_Bulk_S41_L001_R2_001.fastq.gz	RR681Jul_Bulk_S43_L001_R2_001.fastq.gz
RR231Oct_Bulk_S40_L001_R1_001.fastq.gz	RR681Oct_Bulk_S57_L001_R1_001.fastq.gz
RR231Oct_Bulk_S40_L001_R2_001.fastq.gz	RR681Oct_Bulk_S57_L001_R2_001.fastq.gz
RR232Jul_Bulk_S38_L001_R1_001.fastq.gz	RR682Jul_Bulk_S36_L001_R1_001.fastq.gz
RR232Jul_Bulk_S38_L001_R2_001.fastq.gz	RR682Jul_Bulk_S36_L001_R2_001.fastq.gz
RR232Oct_Bulk_S54_L001_R1_001.fastq.gz	RR682Oct_Bulk_S22_L001_R1_001.fastq.gz
RR232Oct_Bulk_S54_L001_R2_001.fastq.gz	RR682Oct_Bulk_S22_L001_R2_001.fastq.gz
RR451Jul_Bulk_S30_L001_R1_001.fastq.gz	RR891Jul_Bulk_S29_L001_R1_001.fastq.gz
RR451Jul_Bulk_S30_L001_R2_001.fastq.gz	RR891Jul_Bulk_S29_L001_R2_001.fastq.gz
RR451Oct_Bulk_S25_L001_R1_001.fastq.gz	RR891Oct_Bulk_S44_L001_R1_001.fastq.gz
RR451Oct_Bulk_S25_L001_R2_001.fastq.gz	RR891Oct_Bulk_S44_L001_R2_001.fastq.gz
RR452Jul_Bulk_S31_L001_R1_001.fastq.gz	RR892Jul_Bulk_S33_L001_R1_001.fastq.gz
RR452Jul_Bulk_S31_L001_R2_001.fastq.gz	RR892Jul_Bulk_S33_L001_R2_001.fastq.gz
RR452Oct_Bulk_S26_L001_R1_001.fastq.gz	RR892Oct_Bulk_S56_L001_R1_001.fastq.gz
RR452Oct_Bulk_S26_L001_R2_001.fastq.gz	RR892Oct_Bulk_S56_L001_R2_001.fastq.gz
```

### Read Trimming
Wanted to start by trimming reads to have decent quality Q scores.
```
for file in $(<All_Samples.txt)
do
        trimmomatic PE -threads 8 -phred33 /share/emersonlab/zinke/RR_190502/${file}L001_R1_001.fastq.gz /share/emersonlab/zinke/RR_190502/${file}L001_R2_001.fastq.gz \
        ${file}R1_paired.fastq.gz ${file}R2_paired.fastq.gz \
        ${file}R1_unpaired.fastq.gz ${file}R2_unpaired.fastq.gz \
        ILLUMINACLIP:/share/emersonlab/annie/programs/Trimmomatic-0.38/adapters/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:30 MINLEN:50
done
```
The file `All_Samples.txt` is a list of the prefixes/sample names for each of the BULK metagenomes. We wrap through this list to perform all of the trimming.

Unfortunately, the output arguments for this command above were listed in the wrong order, so we ended up with the PAIRED R2 samples being names R1_Unpaired and the UNPAIRED R1 Samples being named R2_paired. I renamed the files in order to have the appropriate name for each file.

### Removing PhiX
Performed filtering in order to remove any reads that may be from PhiX that the sequencing facility did not detect.
```
for file in $(<All_Samples.txt)
do
	bbduk.sh in1=../Trimmed_Reads/${file}R1_paired.fastq.gz in2=../Trimmed_Reads/${file}R2_paired.fastq.gz threads=8 \
	out1=${file}R1_filtered.fastq.gz out2=${file}R2_filtered.fastq.gz stats=${file}stats.txt \
	ref=/share/emersonlab/annie/programs/bbduk/phix_adapters.fa k=31 \
	hdist=1
done
```
### Metagenome assembly
Because of the number and size of samples, decided to assemble each of the metagenomes separately using MEGAHIT. Below is an example slurm job script used for assembly.
```
#!/bin/bash --login
########## SBATCH Lines for Resource Request ##########

#SBATCH --time=08:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=1                 # number of different nodes - could be an exact number or a range of nodes (same as -N)
#SBATCH --ntasks=1                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=12           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem-per-cpu=8G            # memory required per allocated CPU (or core) - amount of memory (in bytes)
#SBATCH --job-name Name_of_Job      # you can give your job a name for easier identification (same as -J)

########## Command Lines to Run ##########


module load megahit
cd /share/emersonlab/sorensen/RussellRanch/FreshAnalysis/megahit
  megahit -1 ../PhiX_Removed/RR131Jul_Bulk_S47_R1_filtered.fastq.gz \
   -2 ../PhiX_Removed/RR131Jul_Bulk_S47_R2_filtered.fastq.gz \
   --continue -t 12 --k-min 27 	--min-contig-len 2000 --presets meta-large \
   --out-dir RR131Jul_Bulk_S47_assembly --out-prefix RR131Jul_Bulk_S47_
scontrol show job
```
After assembly, I concatenated all the assembled contigs from each sample togeter using the code below.
```
for file in $(<All_Samples.txt)
do
	cat ${file}assembly/${file}.contigs.fa >> Total_Assembly.fa
done
```

### DETOUR: Attempting to find CRISPRS in the assembled contigs.

One of the goals of binning MAG's is to be able to tie specific viral populations with specific hosts. Several tools currently exist to search for and or assemble CRISPR arrays in with assembled genomes, contigs, or unassembled metagenome reads. The approach used in [(Emerson et al 2018)](https://doi.org/10.1038/s41564-018-0190-y) was to first assemble CRISPR arrays from unassembled reads using CRASS ([Paper](https://doi.org/10.1093/nar/gkt183), [GitHub](https://github.com/ctSkennerton/crass)) and then look for the repeats in assembled Bacterial/Archaeal contigs, and look for the Spacers in Viral contigs identified using VirSorter ([Paper](https://peerj.com/articles/985/), [GitHub](https://github.com/simroux/VirSorter)). Annie and I have both been unsuccessful in making CRASS run without errors. Used CRT ([Paper](https://doi.org/10.1186/1471-2105-8-209), [Software](http://www.room220.com/crt/)) and MinCED ([GitHub](https://github.com/ctSkennerton/minced))to look for CRISPR arrays in my assembled contigs from the Bulk soil samples.
#### 1. CRT

Below is the code used to download and run CRT on the assembled contigs from the bulk soil samples.

```
wget http://www.room220.com/crt/CRT1.2-CLI.jar.zip
unzip CRT1.2-CLI.jar.zip
java -cp CRT1.2-CLI.jar crt ../megahit/Total_Assembly.fa Total_Assembly.out
```

CRT Detected 63 CRISPRS arrays. Unfortunately, the output of CRT is rather hard to parse, and I still need to figure out how to handle the output for optimal efficiency.

#### 2. MinCED

Below is the code used to download and run MinCED on the assembled contigs.

```
git clone https://github.com/ctSkennerton/minced.git
../../../minced-master/minced ../megahit/Total_Assembly.fa Total_Assembly.crisprs
```

MinCED only found 28 CRISPR arrays in the assembled contigs and unfortunately the output is almost identical to CRT so I still need to detemine the best way to parse this file.

### Clustering Metagenome Assemblies

Because I chose to assemble all of the metagenomes separately from each other, I need to perform some clustering of the resulting assembled contigs in order to reduce any redundancy in the dataset. Currently am trying to use psi-cd-hit ([GitHub](https://github.com/weizhongli/cdhit/wiki/3.-User's-Guide#PSICDHIT_clustering)) for this clustering as it is meant to deal with particularly long sequences. This step has not finished running as of yet, but I have posted the code used to install and submit below.
```
# First to install
conda create --name CDHIT
conda activate cdhit
conda install -c bioconda cd-hit
git clone https://github.com/weizhongli/cdhit.git

# Now to submit Job

#!/bin/bash --login
########## SBATCH Lines for Resource Request ##########

#SBATCH --time=08:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=1                 # number of different nodes - could be an exact number or a range of nodes (same as -N)
#SBATCH --ntasks=1                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=12           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem-per-cpu=8G            # memory required per allocated CPU (or core) - amount of memory (in bytes)
#SBATCH --job-name Name_of_Job      # you can give your job a name for easier identification (same as -J)

########## Command Lines to Run ##########

source ~/.bashrc
conda activate CDHIT
cd /share/emersonlab/sorensen/RussellRanch/FreshAnalysis/Cluster_Assemblies/cdhit-master/psi-cd-hit
./psi-cd-hit.pl -i ../../Total_Assembly_Sorted.fa -o ../../Cluster_Total_Assembly.fa -c 0.9 -G 1 -g 1 -prog blastn -circle 1 -exec local -para 3 -blp 4
scontrol show job $SLURM_JOB_ID
```
The above took to long to complete so the job got killed. Recently discovered a way for CD-HIT to save the progress of a job using the `-rs` flag which saves checkpoints every so many steps. Edited job script to increase resources as well has adding the `-rs` flag in case it does not finish in the allotted time.

```
#!/bin/bash --login
########## SBATCH Lines for Resource Request ##########

#SBATCH --time=24:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=1                 # number of different nodes - could be an exact number or a range of nodes (same as -N)
#SBATCH --ntasks=1                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=24           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem-per-cpu=8G            # memory required per allocated CPU (or core) - amount of memory (in bytes)
#SBATCH --job-name Name_of_Job      # you can give your job a name for easier identification (same as -J)

########## Command Lines to Run ##########


source ~/.bashrc
conda activate CDHIT
cd /share/emersonlab/sorensen/RussellRanch/FreshAnalysis/Cluster_Assemblies/cdhit-master/psi-cd-hit
./psi-cd-hit.pl -i ../../Total_Assembly_Sorted.fa -o ../../Cluster_Total_Assembly_MoreResources.fa -c 0.9 -G 1 -g 1 -prog blastn -circle 1 -exec local -para 6 -blp 4 -rs 5000
#./psi-cd-hit.pl -i ../../Total_Assembly_Sorted.fa -o ../../Clustered_Total_Assembly.fa -c 0.90 -aS 0.85
scontrol show job $SLURM_JOB_ID
```
