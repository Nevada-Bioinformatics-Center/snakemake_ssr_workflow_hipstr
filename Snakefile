import glob

configfile: "config.yaml"

wrappers_version="v2.6.0"

inputdirectory=config["input_data"]
#inputbam=config["input_data_bam"]
ref=config["ref"]
minreads=config["minreads"].split(",")
regions=config["regions"]
hipstr_exec=config["hipstr_exec"]
print(inputdirectory)
SAMPLES, =glob_wildcards(inputdirectory+"/{sample}_R1_001.fastq.gz", followlinks=True)

print(SAMPLES)

##### target rules #####
rule all:
    input: 
       "qc/multiqc_report_bwa.html",
       expand("mapped/{sample}.sorted.bam.bai", sample=SAMPLES),
       expand("hipstr/hipstr_results_minreads{minread}.stats", minread=minreads),
       expand("hipstr/hipstr_results_minreads{minread}.vcf.gz", minread=minreads)


rule bwa_index:
    input:
        ref
    output:
        idx=multiext("genome", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        "logs/bwa_index/genome.log"
    params:
        prefix="genome",
        algorithm="bwtsw"
    resources: time_min=520, mem_mb=20000, cpus=1
    wrapper:
        f"{wrappers_version}/bio/bwa/index"

rule trimmomatic_se:
    input:
        inputdirectory+"/{sample}_R1_001.fastq.gz"
    output:
        "trimmed/{sample}.1.fastq.gz"
    log:
        "logs/trimmomatic/{sample}.log"
    params:
        #trimmer = [f"ILLUMINACLIP:{config['ref']['adapter']}:2:30:10 SLIDINGWINDOW:4:15 LEADING:30 MINLEN:36"],
        trimmer = ["SLIDINGWINDOW:4:15 LEADING:30 MINLEN:36"],
        extra="",
        compression_level="-9"
    threads: 8
    resources: time_min=480, mem_mb=20000, cpus=8
    wrapper:
        f"{wrappers_version}/bio/trimmomatic/se"



rule fastqc_posttrim_r1:
    input:
        "trimmed/{sample}.1.fastq.gz"
    output:
        html="qc/fastqc_posttrim/{sample}_r1.html",
        zip="qc/fastqc_posttrim/{sample}_r1_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: ""
    log:
        "logs/fastqc_posttrim/{sample}_r1.log"
    resources: time_min=520, mem_mb=20000, cpus=1
    threads: 1
    wrapper:
        f"{wrappers_version}/bio/fastqc"

rule fastqc_pretrim_r1:
    input:
        inputdirectory+"/{sample}_R1_001.fastq.gz"
    output:
        html="qc/fastqc_pretrim/{sample}_r1.html",
        zip="qc/fastqc_pretrim/{sample}_r1_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: ""
    log:
        "logs/fastqc_pretrim/{sample}_r1.log"
    resources: time_min=320, mem_mb=8000, cpus=1
    threads: 1
    wrapper:
        f"{wrappers_version}/bio/fastqc"

rule multiqc_pre:
    input:
        expand("qc/fastqc_pretrim/{sample}_r1_fastqc.zip", sample=SAMPLES),
    output:
        "qc/multiqc_report_pretrim.html"
    log:
        "logs/multiqc_pre.log"
    resources: time_min=320, mem_mb=8000, cpus=1
    wrapper:
        f"{wrappers_version}/bio/multiqc"


rule bwa_mem:
    input:
        reads="trimmed/{sample}.1.fastq.gz",
        idx=multiext("genome", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    output:
        "mapped/{sample}.bam"
    log:
        "logs/bwa_mem/{sample}.log"
    params:
        index="genome",
        extra=r"-R '@RG\tID:{sample}\tSM:{sample}\tLB:{sample}'",
        sort="none",             # Can be 'none', 'samtools' or 'picard'.
        sort_order="coordinate",  # Can be 'queryname' or 'coordinate'.
        sort_extra=""            # Extra args for samtools/picard.
    threads: 16
    resources: time_min=1320, mem_mb=20000, cpus=16
    wrapper:
        f"{wrappers_version}/bio/bwa/mem"

rule samtools_sort:
    input:
        "mapped/{sample}.bam"
    output:
        "mapped/{sample}.sorted.bam"
    threads: 8
    resources: time_min=220, mem_mb=10000, cpus=10
    log:
        "logs/samtools_sort/sort_{sample}.log"
    conda:
        "envs/samtools.yaml"
    shell:
        "samtools sort -@ {threads} {input} > {output} 2> {log}"


rule samtools_flagstat:
    input:
        "mapped/{sample}.sorted.bam"
    output:
        "mapped/{sample}.sorted.bam.flagstat"
    threads: 1
    resources: time_min=320, mem_mb=8000, cpus=1
    wrapper:
        f"{wrappers_version}/bio/samtools/flagstat"

rule multiqc_bwa:
    input:
        expand("mapped/{sample}.sorted.bam.flagstat", sample=SAMPLES),
        expand("mapped/{sample}.sorted.bam", sample=SAMPLES),
        expand("logs/trimmomatic/{sample}.log", sample=SAMPLES),
        #expand("qc/fastqc_posttrim/{sample}_r1_fastqc.zip", sample=SAMPLES),
    output:
        "qc/multiqc_report_bwa.html"
    log:
        "logs/multiqc.log"
    threads: 1
    resources: time_min=320, mem_mb=8000, cpus=1
    wrapper:
        f"{wrappers_version}/bio/multiqc"


rule samtools_index:
    input:
        "mapped/{sample}.sorted.bam",
    output:
        "mapped/{sample}.sorted.bam.bai",
    log:
        "logs/samtools_index/{sample}.log",
    params:
        extra="",  # optional params string
    threads: 4  # This value - 1 will be sent to -@
    resources: time_min=320, mem_mb=8000, cpus=4
    wrapper:
        f"{wrappers_version}/bio/samtools/index"


rule multiqc_posttrim:
    input:
        expand("qc/fastqc_posttrim/{sample}_r1_fastqc.zip", sample=SAMPLES),
        expand("mapped/{sample}.sorted.bam", sample=SAMPLES)
    output:
        "qc/multiqc_posttrim_report_all.html"
    log:
        "logs/multiqc_posttrim_all.log"
    threads: 1
    resources: time_min=520, mem_mb=40000, cpus=1
    wrapper:
        f"{wrappers_version}/bio/multiqc"

rule samtools_faidx:
    input:
        ref, 
    output:
        ref + ".fai" 
    threads: 8
    resources: time_min=220, mem_mb=10000, cpus=10
    log:
        "logs/samtools_faidx/faidx.log"
    conda:
        "envs/samtools.yaml"
    shell:
        "samtools faidx {input} 2> {log}"

rule run_hipstr:
    input:
        bams=expand("mapped/{sample}.sorted.bam", sample=SAMPLES),
        bai=expand("mapped/{sample}.sorted.bam.bai", sample=SAMPLES),
        ref=ref, 
        regions=regions,
        fai=ref + ".fai"
    output:
        vcf="hipstr/hipstr_results_minreads{minread}.vcf.gz",
        viz="hipstr/hipstr_results_minreads{minread}.viz.gz",
        stdlog="hipstr/hipstr_minreads{minread}.log",
    params:
        listbams=",".join(expand("mapped/{sample}.sorted.bam", sample=SAMPLES)),
        hipstr=hipstr_exec,
    threads: 2
    resources: time_min=4320, mem_mb=40000, cpus=2   ## 3 days runtime
    log:
        "logs/hipstr/hipstr_minreads{minread}.log",
    shell:
        "{params.hipstr} --bams {params.listbams} --fasta {ref} --regions {regions} --str-vcf {output.vcf} --viz-out {output.viz} --log {output.stdlog} --min-reads {wildcards.minread} --use-unpaired 2> {log}"

rule bcftstats:
    input:
        vcf="hipstr/hipstr_results_minreads{minread}.vcf.gz",
        fa=ref
    output:
        "hipstr/hipstr_results_minreads{minread}.stats"
    log: 
        "logs/filter/bcfstats_minreads{minread}.log"
    resources: time_min=220, mem_mb=8000, cpus=1
    threads: 1
    conda:
        "envs/bcftools.yaml"
    shell:
        "bcftools stats -F {input.fa} -s - {input.vcf} > {output} 2> {log}"

