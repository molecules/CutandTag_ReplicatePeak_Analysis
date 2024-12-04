# CutandTag_ReplicatePeak_Analysis

![Cut&Tag_ReplicatePeak](/images/Cut&Tag.png)
+ OpenAI. (2024). NewDescription. DALL-E. Retrieved from OpenAI.

# 1) Project Description
CutAndTag_ReplicatePeak_Analysis is a Snakemake pipeline adapted created in order to... This pipeline is designed to take processed Cut-and-Tag sequencing data to analyze of chromatin accessibility and DNA-protein interactions. It includes steps for...

A compact dataset is included within the repository for testing purposes, along with example scripts for analyzing publicly available Cut-and-Tag datasets. This pipeline extends the original protocol, offering a robust framework for both routine analysis and more complex studies.

# 2) Instructions to run on Slurm managed HPC
2A. Clone repository
```
git clone https://github.com/JK-Cobre-Help/CutandTag_ReplicatePeak_Analysis.git
```
2B. Load modules
```
module purge
module load slurm python/3.10 pandas/2.2.3 numpy/1.22.3 matplotlib/3.7.1
```
2C. Modify samples and config file
```
vim samples.csv
vim config.yml
```
2D. Dry Run
```
snakemake -npr
```
2E. Run on HPC with config.yml options
```
sbatch --wrap="snakemake -j 999 --use-envmodules --latency-wait 60 --cluster-config config/cluster_config.yml --cluster 'sbatch -A {cluster.account} -p {cluster.partition} --cpus-per-task {cluster.cpus-per-task}  -t {cluster.time} --mem {cluster.mem} --output {cluster.output}'"
```

# 3) Explanation of samples.csv
Note. Make sure to check sample.csv before each run

The samples.csv file in the config folder has paths to the test fastq files. You must replace those paths with those for your own fastq files. The first column of each row is the sample name. This name will be used for all output files. Columns 2 and 3 are the paths to the paired fastq files. Column 4 is the sample type (either "treatment" or "control"). Column 5 is the name of the corresponding Control sample for each treated sample (use "NA" if the sample is a control).

| sample             | fastq1                        | fastq2                        | sampleType | Control   |
|--------------------|-------------------------------|-------------------------------|------------|-----------|
| K27ac_50_trimmed   | K27ac_50_trimmed_R1.fastq.gz  | K27ac_50_trimmed_R2.fastq.gz  | control    | NA        |
| K27me3_50_trimmed  | K27me3_50_trimmed_R1.fastq.gz | K27me3_50_trimmed_R1.fastq.gz | control    | NA        |


Sample naming recommendation for correct plot output
- "Histone" + "_" + "Replicate" + "Any other identifier"
- Examples:
    + K27ac_50
    + K27me3_5
    + K27ac_50_trimmed
    + H3K27me3_rep1
    + H3K4me3_rep2_set1
    + H3K27ac_rep3_control
    + H3K27ac_rep3_treatment

# 4) Explanation of config.yml
Note. Make sure to check config.yml for the appropriate genome alignment

The config.yml is used to identify the file path of the bowtie2 genome index, specify effective genome size and genome for macs2. There is also information about specific modules and version numbers to maintain dependencies in the snakemake workflow. Running the mm10 genome does not require any modifications to the config.yml. When using the hg38 genome the following need to be modified with the information provided in the config.yml but commented out.

Run hg38 samples in snakemake pipeline
- config.yml 
    + change bowtie2 genome index file path
    + change bamCoverage effective genome size
    + change macs2 genome size

# 5) Citations
Zheng, Y., Ahmad, K., & Henikoff, S. (2019). CUT&Tag for efficient epigenomic profiling of small samples and single cells. Protocols.io, dx.doi.org/10.17504/protocols.io.bjk2kkye
