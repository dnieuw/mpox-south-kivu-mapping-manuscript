import sys, string, shutil, glob
from pandas import read_table

sample_data = read_table(config["sample_data"], dtype={"UniqueID":str,"FASTQ_path":str,"Reference":str,"primers":str,"Sequence_name":str,"gzipped":str,"coverage":str,"min_length":str})
SAMPLES = sample_data["UniqueID"].tolist()
REFERENCE = sample_data["Reference"].tolist()[0]
sample_data = sample_data.set_index('UniqueID')
sample_data = sample_data.to_dict("index")

rule all:
    input:
        consensus="alignment/all_consensus_aligned.fasta",
        stats_raw="readstats/raw.tsv",
        stats_qc="readstats/QC.tsv",
        stats_trimmed="readstats/trimmed.tsv",
        stats_mapped="readstats/mapped.tsv"

rule merge_barcodes:
    input:
        lambda wildcards: sample_data[wildcards.sample_id]["FASTQ_path"]
    output:
        "raw/{sample_id}.fastq"
    threads: 1
    shell:
        """
        zcat {input}/*.fastq.gz > {output}
        """

### PREPROCESSING ###
rule QC_Nanopore_reads:
    input:
        "raw/{sample_id}.fastq"
    output:
        fastq="QC/{sample_id}.fastq",
        report="QC/{sample_id}.html"
    threads: 2
    params:
        min_length=lambda wildcards: sample_data[wildcards.sample_id]["min_length"]
    shell:
        """
        fastp -i {input} -o {output.fastq} -j /dev/null -h {output.report} \
        --disable_trim_poly_g \
        --disable_adapter_trimming \
        --qualified_quality_phred 10 \
        --unqualified_percent_limit 50 \
        --length_required {params.min_length} \
        -w {threads} > /dev/null 2>&1
        """

rule cut_primers:
    input:
        fastq="QC/{sample_id}.fastq"
    output:
        cut="trimmed/{sample_id}_cutadapt.fastq"
    threads: 1
    shell:
        """
        cutadapt -u 30 -u -30 -o {output.cut} {input.fastq} -j {threads} --quiet
        """

rule map_to_primer_reference:
	input:
		fastq="trimmed/{sample_id}_cutadapt.fastq"
	output:
		mapped="trimmed/{sample_id}_mapped.bam"
	params:
		reference=lambda wildcards: sample_data[wildcards.sample_id]["Primer_reference"]
	threads: 4
	shell:
		"""
		minimap2 -Y -t {threads} -x map-ont -a {params.reference} {input.fastq} 2> /dev/null | samtools view -bF 4 - | samtools sort -@ {threads} - > {output.mapped}
		"""

rule trim_primers:
    input:
        mapped="trimmed/{sample_id}_mapped.bam"
    output:
        clipped="trimmed/{sample_id}_clipped.bam",
        trimmed="trimmed/{sample_id}_trimmed.fastq"
    params:
        reference=lambda wildcards: sample_data[wildcards.sample_id]["Primer_reference"],
        primers=lambda wildcards: sample_data[wildcards.sample_id]["primers"],
        min_length=lambda wildcards: sample_data[wildcards.sample_id]["min_length"]
    threads: 1
    shell:
        """
        samtools index {input.mapped}

        scripts/ampliclip.py \
        --infile {input.mapped} \
        --outfile {output.clipped}_ \
        --outfastq {output.trimmed} \
        --primerfile {params.primers} \
        --referencefile {params.reference}\
        -fwd LEFT -rev RIGHT \
        --padding 10 --mismatch 2 --minlength {params.min_length}

        samtools sort {output.clipped}_ > {output.clipped}
        rm {output.clipped}_
        """

rule map_to_reference:
    input:
        fastq="trimmed/{sample_id}_trimmed.fastq",
        reference=lambda wildcards: sample_data[wildcards.sample_id]["Reference"]
    output:
        "mapped/{sample_id}_mapped.bam"
    threads: 4
    shell:
        """
        minimap2 -Y -t {threads} -x map-ont -a {input.reference} {input.fastq} 2> /dev/null | samtools view -bF 4 - | samtools sort -@ {threads} - > {output}
        """

### CONSENSUS CALLING ###
rule create_depth_file:
    input:
        bam="mapped/{sample_id}_mapped.bam",
        reference=lambda wildcards: sample_data[wildcards.sample_id]["Reference"]
    output:
        "depth/{sample_id}.tsv"
    shell:
        """
        samtools mpileup -a -A -Q 0 -d 8000 -f {input.reference} {input.bam} | cut -f 2,3,4 --output-delimiter=',' > {output}
        """

rule create_vcf:
    input:
        bam="mapped/{sample_id}_mapped.bam",
        reference=lambda wildcards: sample_data[wildcards.sample_id]["Reference"]
    output:
        "vcf/{sample_id}.vcf"
    threads: 4
    params:
        coverage=lambda wildcards: sample_data[wildcards.sample_id]["coverage"]
    shell:
        """
        samtools index -@ {threads} {input.bam} > /dev/null 2>&1
        scripts/bam2vcf.py -c {threads} -d {params.coverage} -af 0.1 -r {input.reference} -b {input.bam} -o {output} 
        """

rule filter_vcf:
    input:
        "vcf/{sample_id}.vcf"
    output:
        "vcf/{sample_id}_filtered.vcf"
    threads: 1
    shell:
        """
        scripts/filtervcf.py -i {input} -o {output}
        """

rule create_consensus:
    input:
        vcf="vcf/{sample_id}_filtered.vcf",
        depth="depth/{sample_id}.tsv",
        reference=lambda wildcards: sample_data[wildcards.sample_id]["Reference"]
    output:
        "consensus/{sample_id}_consensus.fasta"
    threads: 1
    params:
        coverage=lambda wildcards: sample_data[wildcards.sample_id]["coverage"],
        sequence_name=lambda wildcards: sample_data[wildcards.sample_id]["Sequence_name"]
    shell:
        """
        scripts/vcf2consensus.py -v {input.vcf} -d {input.depth} -r {input.reference} -o {output} -dp {params.coverage} -n "{params.sequence_name}_consensus"
        """

rule aggregate_consensus:
    input:
        expand("consensus/{sample_id}_consensus.fasta",sample_id=SAMPLES)
    output:
        "alignment/all_consensus.fasta"
    threads: 1
    shell:
        """
        cat {input} > {output}
        """

rule align_consensus:
    input:
        "alignment/all_consensus.fasta"
    output:
        "alignment/all_consensus_aligned.fasta"
    threads: 8
    params:
        reference=REFERENCE
    shell:
        """
        minimap2 -t {threads} -a -x asm20 --sam-hit-only --secondary=no --score-N=0 {params.reference} {input} -o alignment/all_consensus_aligned.sam
        gofasta sam toMultiAlign -s alignment/all_consensus_aligned.sam -o {output}
        """

### STATS GENERATION ###
rule generate_readstats_raw:
    input:
        expand("raw/{sample_id}.fastq", sample_id=SAMPLES)
    output:
        "readstats/raw.tsv"
    threads: 1
    shell:
        """
        seqkit stats -T {input} > {output}
        """

rule generate_readstats_QC:
    input:
        expand("QC/{sample_id}.fastq", sample_id=SAMPLES)
    output:
        "readstats/QC.tsv"
    threads: 1
    shell:
        """
        seqkit stats -T {input} > {output}
        """

rule generate_readstats_trimmed:
    input:
        expand("trimmed/{sample_id}_trimmed.fastq", sample_id=SAMPLES)
    output:
        "readstats/trimmed.tsv"
    threads: 1
    shell:
        """
        seqkit stats -T {input} > {output}
        """

rule generate_readstats_mapped:
    input:
        expand("mapped/{sample_id}_mapped.bam", sample_id=SAMPLES)
    output:
        "readstats/mapped.tsv"
    threads: 1
    shell:
        """
        for file in {input}; do
            samtools fastq $file | seqkit stats -T --stdin-label $file | tail -1
        done > {output}
        """