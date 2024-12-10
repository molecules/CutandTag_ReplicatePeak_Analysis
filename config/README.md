# CutAndTag_ReplicatePeak_Analysis

![ReplicatePeaks](/images/replicatePeaks.png)  
*OpenAI. (2024). Scientific data visualization: Replicate peak analysis in bioinformatics [AI-generated image]. DALL-E. Retrieved from ChatGPT interface.*

## 1) Project Description

**CutAndTag_ReplicatePeak_Analysis** is a Snakemake pipeline designed for downstream peak analysis on processed Cut-and-Tag sequencing data. Instead of starting from raw FASTQ reads, it begins with already aligned and filtered BAM files. The pipeline focuses on:

- Identifying reproducible peaks
- Generating consensus peak sets
- Visualizing overlaps and signal distributions across multiple samples or experimental conditions

**Note**: If you need to process raw FASTQ files (e.g., quality control, trimming, alignment), consider using the [CutandTag_Analysis_Snakemake](https://github.com/JK-Cobre-Help/CutandTag_Analysis_Snakemake) pipeline first. It provides the cleaned and aligned data suitable as input for **CutAndTag_ReplicatePeak_Analysis**.

## Key Features

- **Peak Calling with MACS2**:  
  Calls peaks for each sample using MACS2, a widely-used tool for identifying enriched regions in ChIP- or Cut-and-Tag sequencing data.

- **Merged and Consensus Peak Sets**:  
  Uses sample groupings defined in `samples.csv` (the "Set" column) to merge peaks and generates a consensus peak set by applying a reproducibility threshold.

- **Consensus Peak Conversion**:  
  Converts consensus peak sets into BAM and BigWig formats for efficient genome browser visualization and downstream analysis.

- **Euler Plots of Overlaps**:  
  Creates Euler diagrams to illustrate how consensus peaks emerge from overlapping individual sample peaks within a set.

- **Midpoint and Overlap Analysis**:  
  Identifies peak midpoints and quantifies overlaps for a detailed exploration of peak distribution across samples.

- **Heatmaps for Signal Distribution**:  
  Uses consensus peak midpoints to generate heatmaps that visualize signal intensity and distribution across multiple conditions or sample sets.

## 2) Intended Use Case

This pipeline is ideal for researchers who have already processed their Cut-and-Tag data (alignment, filtering, etc.) and want to:

- Identify reproducible peaks across replicates or conditions.
- Visualize peak overlaps and generate integrated summaries.
- Compare signal intensity profiles around consensus peak midpoints.

By combining the preliminary processing pipeline ([CutandTag_Analysis_Snakemake](https://github.com/JK-Cobre-Help/CutandTag_Analysis_Snakemake)) with this downstream analysis, you create a robust end-to-end workflow.

## 3) Dependencies and Configuration

All parameters (e.g., genome size, MACS2 q-values, minimum overlaps for consensus peaks, tool paths) are specified in `config/config.yml`.

### Explanation of `config.yml`

The `config.yml` file controls genome settings, tool versions, and other workflow parameters.

- **Genome and Effective Genome Size**:  
  For human (hg38): `genome: "hs"` and `effective_genome_size: 2913022398`  
  For mouse (mm10): `genome: "mm"` and `effective_genome_size: 2730871774`

By default, the pipeline is set up for hg38 (`hs`). To run mm10 samples, uncomment and adjust the settings in `config.yml` accordingly.

### Changing Genomes

To switch between mm10 and hg38:

- **Effective Genome Size**: Update `effective_genome_size` (e.g., `2730871774` for mm10 or `2913022398` for hg38).
- **Chrom Sizes File**: Update `chrom_sizes` to point to the correct file (e.g., `resources/mm10.chrom.sizes` or `resources/hg38.chrom.sizes`).
- **MACS2 Genome Parameter**: Set `genome: "mm"` for mouse or `genome: "hs"` for human.

### Tool Versions and Modules

`config.yml` also specifies tool versions (e.g., `deeptools`, `macs2`, `samtools`, `bedtools`, `R`) to ensure reproducibility and consistent results.

## 4) Tools & Modules

The pipeline uses:

- **MACS2** for peak calling  
- **bedtools** and **samtools** for working with peak and alignment formats  
- **deeptools** for coverage, matrix computation, and heatmap generation  
- **R** with Bioconductor for merging peaks, generating consensus sets, and creating Euler diagrams

## 5) Example Data

A small, pre-processed dataset is included to quickly test the pipeline and validate your environment. It demonstrates the workflow’s steps from peak calling through to final visualization.

### Relationship to Previous Protocols

This pipeline builds on [CutandTag_Analysis_Snakemake](https://github.com/JK-Cobre-Help/CutandTag_Analysis_Snakemake), which processes raw FASTQs. **CutAndTag_ReplicatePeak_Analysis** focuses on downstream analysis, starting from aligned BAM files to produce detailed consensus peak sets and visual summaries.

## 6) Explanation of `samples.csv`

`samples.csv` specifies the samples to analyze, their BAM file locations, and how they are grouped into sets. The file has three columns: `sample`, `bam`, and `set`.

**Example `samples.csv`:**
```csv
sample,bam,set
Treatment_Rep1,resources/test1.bam,Set1
Treatment_Rep2,resources/test1A.bam,Set1
Treatment_Rep3,resources/test1B.bam,Set1
Control_Rep1,resources/input1.bam,Set2
Control_Rep2,resources/input1A.bam,Set2
Control_Rep3,resources/input2B.bam,Set2
```

**sample:** Unique sample name (used in output filenames)  
**bam:** Path to the aligned BAM file  
**set:** Sample grouping for consensus peak analysis

Choose descriptive sample names for clarity in outputs and plots.

## 7) Instructions to Run on a Slurm-Managed HPC

**Step A: Clone the repository**
```bash
git clone https://github.com/JK-Cobre-Help/CutandTag_ReplicatePeak_Analysis.git

