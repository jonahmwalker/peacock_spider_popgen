# Jonah Walker

# Specify path to all modules used
samtools_badger = "/software/badger/module-builds/samtools/1.21/bin/samtools"
samtools = "/software/treeoflife/shpc/0.1.26/wrapper/quay.io/biocontainers/samtools/1.20--h50ea8bc_0/bin/samtools"
fastp = "/software/treeoflife/shpc/0.1.26/wrapper/quay.io/biocontainers/fastp/0.23.4--hadf994f_3/bin/fastp"
bwa = "/software/treeoflife/shpc/0.1.26/wrapper/quay.io/biocontainers/bwa/0.7.17--h7132678_9/bin/bwa"
bwa_mem2 = "/software/treeoflife/shpc/0.1.26/wrapper/quay.io/biocontainers/bwa-mem2/2.2.1--hd03093a_2/bin/bwa-mem2"
bcftools = "/software/treeoflife/shpc/0.1.26/wrapper/quay.io/biocontainers/bcftools/1.20--h8b25389_0/bin/bcftools"


# Specify reference file
# NOTE — Reference species ID specified in file names ("qqMarSpei1.1")
reference = "/lustre/scratch122/tol/teams/meier/analysis/jonah/maratus_popgen/01_reference/qqMarSpei1.1/GCA_963932465.1.qqMarSpei1.1.renamed.fasta.gz"
# Specify metadata file
import pandas as pd
metadata_table = pd.read_csv('0_scripts_snakemake/maratus_popgen_reads_to_vcf.metadata.tsv', sep='\t', index_col='sequence_id')
sequence_id = metadata_table.index
specimen_id = metadata_table["specimen_id"]
library_id = metadata_table["library_id"]
path = metadata_table["path"]
# Specify intervals file
import pandas as pd
intervals_table = pd.read_csv('0_scripts_snakemake/maratus_popgen_reads_to_vcf.intervals_qqMarSpei1.1.tsv', sep='\t', index_col='interval')
interval = intervals_table.index
coordinates = intervals_table["coordinates"]


# NOTE — Tell snakemake what the final file output will be so that it knows the directionality of the pipeline. Hash out all but one of the options below.
# NOTE — YOU MUST ALSO CHANGE THE OCCURENCES OF 'TEMP()' IN THE OUTPUT LINES OF RULES SO THAT THE CORRECT FILES ARE DELETED!
rule all:
# To finsh with separated raw fastq files:
    # input: expand("A_fastq_raw/{sequence_id}.raw.R1.fq.gz", zip, sequence_id=sequence_id)
# To finish with the trimmed fastq files:
    # input: expand("B_fastq_trimmed/{sequence_id}.trimmed.R1.fq.gz", zip, sequence_id=sequence_id)
# To finish with the aligned cram files:
    # input: expand("C_cram_aligned_qqMarSpei1.1/{specimen_id}.{sequence_id}.qqMarSpei1.1.cram.flagstats", zip, specimen_id=specimen_id, sequence_id=sequence_id)
# To finish with the deduplicated and processed (indexed and stats-ed) aligned cram files:
    # input: expand("D_cram_aligned_qqMarSpei1.1_rm/{specimen_id}.{sequence_id}.qqMarSpei1.1.rm.cram.flagstats", zip, specimen_id=specimen_id, sequence_id=sequence_id)
# To finish with the vcf:
    # input: expand("E_vcf_qqMarSpei1.1/maratus.qqMarSpei1.1.{interval}.vcf.gz", interval=interval)


# NOTE — In all the steps below, you need to specify the cpus, memory on first attempt, and queue to submit for LSF.


# Step A1: Download and convert interleaved cram files from iRODS
rule A1_cram_to_fastq_iRODS:
    params:
        cram = lambda wildcards: path[wildcards.sequence_id]
    output:
        f=temp("A_fastq_raw/{sequence_id}.raw.R1.fq.gz"),
        r=temp("A_fastq_raw/{sequence_id}.raw.R2.fq.gz")
    resources:
        mem_mb=lambda wildcards, attempt: 1000 * 2 * attempt,
        cpu=5,
        queue="normal",
        jobname="A1_cram_to_fastq_iRODS"
    shell:
        """
        {samtools_badger} fastq irods:{params.cram} -N -0 /dev/null -1 {output.f} -2 {output.r};
        fastqc {output.f} {output.r}
        """
# Step A2: Convert pre-downloaded interleaved cram files
# rule A2_cram_to_fastq_farm:
#     input:
#         cram = lambda wildcards: path[wildcards.sequence_id]
#     output:
#         f=temp("A_fastq_raw/{sequence_id}.raw.R1.fq.gz"),
#         r=temp("A_fastq_raw/{sequence_id}.raw.R2.fq.gz")
#     resources:
#         mem_mb=lambda wildcards, attempt: 1000 * 2 * attempt,
#         cpu=5,
#         queue="normal",
#         jobname="A2_cram_to_fastq_farm"
#     shell:
#         """
#         {samtools_badger} fastq {input.cram} -N -0 /dev/null -1 {output.f} -2 {output.r}
#         fastqc {output.f} {output.r}
#         """


# Step B: Filter reads
rule B_filter_fastq:
    input:
        f="A_fastq_raw/{sequence_id}.raw.R1.fq.gz",
        r="A_fastq_raw/{sequence_id}.raw.R2.fq.gz"
    output:
        f=temp("B_fastq_trimmed/{sequence_id}.trimmed.R1.fq.gz"),
        r=temp("B_fastq_trimmed/{sequence_id}.trimmed.R2.fq.gz"),
        html="B_fastq_trimmed/{sequence_id}.fastp.html",
        json="B_fastq_trimmed/{sequence_id}.fastp.json"
    resources:
        mem_mb=lambda wildcards, attempt: 3000 * 2 * attempt,
        cpu=10,
        queue="normal",
        jobname="B_filter_fastq"
    threads: 8
    shell: 
        """
        {fastp} --in1 {input.f} --in2 {input.r} --out1 {output.f} --out2 {output.r} -l 50 -g --detect_adapter_for_pe -w {threads} -h {output.html} -j {output.json};
        fastqc {output.f} {output.r}
        """


# Step C: Map the filtered reads to the reference genome, index, and get stats
rule C_map_fastq:
    params:
        sequence_id = "{sequence_id}",
        specimen_id = "{specimen_id}",
        library_id = lambda wildcards: library_id[wildcards.sequence_id]
    input:
        fq1="B_fastq_trimmed/{sequence_id}.trimmed.R1.fq.gz",
        fq2="B_fastq_trimmed/{sequence_id}.trimmed.R2.fq.gz"
    output:
        cram=temp("C_cram_aligned_qqMarSpei1.1/{specimen_id}.{sequence_id}.qqMarSpei1.1.cram"),
        index=temp("C_cram_aligned_qqMarSpei1.1/{specimen_id}.{sequence_id}.qqMarSpei1.1.cram.crai"),
        stats="C_cram_aligned_qqMarSpei1.1/{specimen_id}.{sequence_id}.qqMarSpei1.1.cram.stats",
        flagstats="C_cram_aligned_qqMarSpei1.1/{specimen_id}.{sequence_id}.qqMarSpei1.1.cram.flagstats"
    resources:
        mem_mb=lambda wildcards, attempt: 45000 * 2 * attempt,
        cpu=20,
        queue="long",
        jobname="C_map_fastq"
    threads: 20
    shell:
        """
        {bwa_mem2} mem {reference} {input.fq1} {input.fq2} -R "@RG\tSM:{params.specimen_id}\tLB:{params.library_id}\tPU:{params.sequence_id}\tPL:ILLUMINA\tID:{params.specimen_id}.{params.library_id}.{params.sequence_id}" -t {threads} | {samtools} sort -O bam -u -T tmp/C_{params.specimen_id}.{params.sequence_id}.sort - | {samtools} view -T {reference} -C -o {output.cram} -;
        {samtools} index {output.cram} > {output.index};
        {samtools} stats {output.cram} > "{output.stats}";
        {samtools} flagstats {output.cram} > "{output.flagstats}"
        """


# Step D: Remove duplicates and unmapped reads from aligned cram files, index, and get stats
rule D_process_cram:
    params:
        specimen_id = "{specimen_id}",
        sequence_id = "{sequence_id}"
    input:
        cram="C_cram_aligned_qqMarSpei1.1/{specimen_id}.{sequence_id}.qqMarSpei1.1.cram",
        index="C_cram_aligned_qqMarSpei1.1/{specimen_id}.{sequence_id}.qqMarSpei1.1.cram.crai"
    output:
        cram="D_cram_aligned_qqMarSpei1.1_rm/{specimen_id}.{sequence_id}.qqMarSpei1.1.rm.cram",
        index="D_cram_aligned_qqMarSpei1.1_rm/{specimen_id}.{sequence_id}.qqMarSpei1.1.rm.cram.crai",
        markdup_stats="D_cram_aligned_qqMarSpei1.1_rm/{specimen_id}.{sequence_id}.qqMarSpei1.1.rm.cram.samtools_markdup_stats",
        stats="D_cram_aligned_qqMarSpei1.1_rm/{specimen_id}.{sequence_id}.qqMarSpei1.1.rm.cram.stats",
        flagstat="D_cram_aligned_qqMarSpei1.1_rm/{specimen_id}.{sequence_id}.qqMarSpei1.1.rm.cram.flagstats"
    resources:
        mem_mb=lambda wildcards, attempt: 17500 * 2 * attempt,
        cpu=20,
        queue="normal",
        jobname="D_process_cram"
    threads: 20
    shell:
        """
        {samtools} collate -@ {threads} -T tmp/D_{params.specimen_id}.{params.sequence_id}.collate -O -u {input.cram} | {samtools} fixmate -@ {threads} -m -u - - | {samtools} sort -@ {threads} -T tmp/D_{params.specimen_id}.{params.sequence_id}.sort -u - | {samtools} markdup -@ {threads} -T tmp/D_{params.specimen_id}.{params.sequence_id}.markdup -r -f {output.markdup_stats} --json -u - - | {samtools} view -@ {threads} -F 4 - -o {output.cram};
        {samtools} index {output.cram} -o {output.index};
        {samtools} stats {output.cram} > {output.stats};
        {samtools} flagstat {output.cram} > {output.flagstat}
        """


# Step E: Call variants from the bam files with bcftools.
rule E_call_variants:
    params:
        coordinates=lambda wildcards: coordinates[wildcards.interval]
    input:
        cram=expand("D_cram_aligned_qqMarSpei1.1_rm/{specimen_id}.{sequence_id}.qqMarSpei1.1.rm.cram", zip, specimen_id=specimen_id, sequence_id=sequence_id),
        index=expand("D_cram_aligned_qqMarSpei1.1_rm/{specimen_id}.{sequence_id}.qqMarSpei1.1.rm.cram.crai", zip, specimen_id=specimen_id, sequence_id=sequence_id)
    output:
        "E_vcf_qqMarSpei1.1/maratus.qqMarSpei1.1.{interval}.vcf.gz"
    resources:
        mem_mb=lambda wildcards, attempt: 10000 * 2 * attempt,
        cpu=1,
        queue="long",
        jobname="E_call_variants"
    shell:
        "{bcftools} mpileup -q 30 -Q 20 -a AD,DP,SP -Ou -f {reference} -r {params.coordinates} {input.cram} | {bcftools} call -f GQ,GP -m -Oz -o {output}"
