# Russell Ranch Metagenome Assembled Genomes Workflow
## Rationale/Setup
Russell Ranch is a long term research site set up to investigate the longterm influences of different management practices. Our goal was to understand how these different management practices may influence the local viral communities, and their respective hosts. In July and October of 2018, Laura Zinke, Joanne Emerson, and Sara Geonczy sampled 6 different Tomato/Corn rotation plots, 3 representing that were managed with an Organic practice (Poultry manure for fertlizer + Winter Cover Crop) and 3 that were managed with a Traditional practice (Mineral fertilizer + Fallow winter season). At each plot 2 replicate cores were taken at different spots. From each of these samples (2 management practices * 3 plots per practice * 2 reps per plot * 2 time points == 24 total samples), BULK metagenomes were produced via direct DNA extraction with the Powersoil kit + Illumina Sequences, and Viromes were produced by first physically separating viral particles from soil and microbes, then extracting DNA from these viral particles and sequencings using Illumina. 6 Viromes failed to produce enough sequences for analysis. All six of these samples were from the July time point and were evenly split between the conventional and organic management practices. The following is a workflow notes for performing metagenome assembly and Genome binning.

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
