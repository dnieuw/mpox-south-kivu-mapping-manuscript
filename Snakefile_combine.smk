import sys, string, shutil, glob
from pandas import read_table

sample_data = read_table(config["sample_data"], dtype={"LongID":str,"ShortID":str,"Long_path":str,"Short_path":str,"Sequence_name":str,"Reference":str,"coverage":str}, comment="#")
SAMPLES = sample_data["Sequence_name"].tolist()
REFERENCE = sample_data["Reference"].tolist()[0]
sample_data = sample_data.set_index('Sequence_name')
sample_data = sample_data.to_dict("index")

rule all:
    input:
        consensus="alignment/all_consensus_aligned.fasta",
        stats_merged="readstats/merged.tsv"

rule merge_fastq:
    input:
        fastq1=lambda wildcards: sample_data[wildcards.sample_id]["Long_path"],
        fastq2=lambda wildcards: sample_data[wildcards.sample_id]["Short_path"]
    output:
        "merged/{sample_id}_merged.fastq"
    threads: 1
    shell:
        """
        cat {input.fastq1} {input.fastq2} > {output}
        """

rule map_to_reference:
    input:
        fastq="merged/{sample_id}_merged.fastq",
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
        sequence_name=lambda wildcards: wildcards.sample_id
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
rule generate_readstats_trimmed:
    input:
        expand("merged/{sample_id}_merged.fastq", sample_id=SAMPLES)
    output:
        "readstats/merged.tsv"
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