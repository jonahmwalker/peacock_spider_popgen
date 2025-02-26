

This snakemake pipeline will produce VCFs in specified genomic regions from raw read files in unaligned cram format.
This was developed to accommodate the very large genomes of peacock spider (5-10Gb) and therefore minimises storage usage.


Steps in this snakemake pipeline:
A1: Download cram-format raw reads from iRODs and unweave forward and reverse reads on the fly (required for trimming) with samtools. Then run QC on raw reads with fastQC.
A2: Alternatively, unweave forward and reverse reads from cram-format raw reads already downloaded with samtools. Then run QC on raw reads with fastQC.
B: Trim reads with fastp (relaxed filter: min length 50bp, remove adapter content). Then run QC on trimmed reads with fastQC.
C: Map filtered reads to the reference genome to produce aligned cram files with BWA-mem2 index. Then index the aligned cram files and get stats with samtools.
D: Remove duplicates and unmapped reads from the aligned cram files with samtools. Then index the filtered cram files and get stats with samtools.
E: Call variants in specified genomic windows with bcftools mpileup | call (relaxed filter: min MAPQ 30, min baseQ 20).


# Setup

You should git clone this repository.
```
git clone https://github.com/jonahmwalker/peacock_spider_popgen/new/main/02_reads_to_vcf
```

It is set up to work on the Sanger farm HPC.
You may need modify the following input files for your purposes:
1. `0_scripts_snakemake/maratus_popgen_reads_to_vcf.snakefile` specifies the location of the reference genome file, the location of the metadata file, and the location of the genomic intervals file. It also prefixes the output files in a manner consistent with the reference genome (in my case, something like 'qqMarSpei1.1'). It also specifies the resources required for each step when submitting jobs to LSF, and how to modify these on restarts if a job fails.
2. `0_scripts_snakemake/maratus_popgen_reads_to_vcf.metadata.tsv` specifies the metadata of the sequenced samples to be included in the pipeline. In this tsv file: 'sequence_id' is a unique identifier of a particular set of reads (I use sequencing-run#_lane#_plex-tag#); 'specimen_id' is the unique name of a particular specimen from which the reads derive; 'library_id' is the unique identifer of a particular illumina library prepared from a particular specimen (I use specimen-id_library#'; and path is the path to the reads file (must be in quotation marks).
3. `0_scripts_snakemake/maratus_popgen_reads_to_vcf.intervals_qqMarSpei1.1.tsv` specifies the genomic regions to use when calling variants in step E. In this tsv file: 'interval' is a unique code for each interval region; and 'coordinates' are the genomic coordiantes of that interval (must be in quotation marks).
4. `0_scripts_snakemake/maratus_popgen_reads_to_vcf.cluster.yaml` specifies the behaviour of the cluster to which the job is submitted.


# Running the pipeline

First, ensure that the `tmp`, `log`, and `0_scripts_snakemake` directories are present.
```
cd 02_reads_to_vcf && ls
```

Next, on the Sanger farm HPC, log into iRODS to access data.
```
iinit
```

Then activate conda environment for snakemake v7.8.0.
The yaml of this conda environment can be found in `0_conda_env`
```
module load conda
conda activate snakemake_7.8.0
```

Do a dry run of the pipeline to check all is well.
```
snakemake -n -s 0_scripts_snakemake
```

Then run the submission command!
It allows max 999 jobs at a time, reruns any unfinished steps from previous runs, and restarts failed jobs up to 2 times (see the snakemake for resource allocation increases on restarts).
```
snakemake -s 0_scripts_snakemake/maratus_popgen_reads_to_vcf.snakefile \
  --rerun-incomplete -j 999 --restart-times 2 \
  --cluster "bsub -G {cluster.account} -R 'select[mem>{resources.mem_mb}] rusage[mem={resources.mem_mb}] span[hosts=1]' -M {resources.mem_mb} -n {resources.cpu} -q {resources.queue} -o {cluster.output} -e {cluster.error}" \
  --cluster-config 0_scripts_snakemake/maratus_popgen_reads_to_vcf.cluster.yaml --cluster-cancel bkill
```

