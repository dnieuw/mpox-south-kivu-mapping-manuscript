# Repository associated with the Mpox South-Kivu mapping manuscript

## Initial data analysis and generation of genomes

Raw sequence data of the two different amplicon-based approaches were preprocessed separately and merged afterwards to generate the final whole genome consensus sequences. Preprocessing was done by fastp with parameters “--disable_trim_polyg  --disable_adapter_trimming --qualified_quality_phred 10 --unqualified_percent_limit 50 --length_required 1000“ for the long or “—length_required 150” for the short  amplicon panel. 

Primer trimming was performed using [ampliclip](https://github.com/dnieuw/Ampliclip). The long amplicon reads were mapped to the reference genome NC_063383.1, based on which the primers for these amplicons were designed, with the parameters “-Y -x map-ont”. Using Ampliclip and the available long amplicon primers, with parameters “--padding 10 --mismatch 2 --minlength 1000”, the mapped reads were softclipped returned as trimmed reads.

For the short amplicon set initial primer trimming was performed by removing 30 nucleotides from both ends of the amplicon using cutadapt with parameters “-u 30 -u -30” to be certain all primers were trimmed. Additional primer trimming was performed using by mapping the reads to a truncated version of ON585033.1 (last 6491 nt removed) using minimap with parameters “-Y -x map-ont”. The reference was truncated to handle the absence of amplicons in this region for this amplicon set. For the short amplicon protocol, the primer sequences are not publicly available; we assumed that these were the 30 nt flanking regions of the amplicons (start and end position of the amplicons are available), which was used as an input for Ampliclip. We used parameters “--padding 10 --mismatch 2 --minlength 150” for the short amplicon set.

The resulting clipped fastq files of the short and long amplicon set were merged for samples for which both datasets were available. Clipped reads were mapped to the NC_003310.1 reference genome using minimap2 with the same parameters as above. 

Coverage depth was calculated using samtools mpileup with parameters “-a -A -Q 0 -d 8000”. 

Variant calling was performed by generating a VCF file using our tool bam2vcf and parameters “-af 0.1” and filtering the VCF file to retain only major variants. Major variants were incorporated in the genome using a custom vcf2consensus script and regions with coverage below the coverage threshold of 30x were masked with an “N”.

After merging a multiple sequence alignment was made using gofasta following the recommendations on their [github page](https://github.com/virus-evolution/gofasta?tab=readme-ov-file#sam-to-fasta-format-conversion)

Additional manual curation was done by scanning through the multiple sequence alignment and verifying mutations that were only present in either the long or the short amplicon set data. The position of the mutations was verified to not match a primer location or be located in a complex homopolymeric or di-, tri-, or multimer repeat region.

## Analysis steps

### Installation

Install the required softwares using conda and following command:

''conda env create -f environment.yml''

### Analyze long amplicon data

Create a directory for the output of the workflow (e.g. results_long).

Create a config file (e.g. sample_data_long.tsv) containing the following columns:

| UniqueID   | FASTQ_path                   | Reference                         | Primers                                            | Primer_reference                     | Sequence_name     | Gzipped | Coverage | Min_length |
|------------|------------------------------|-----------------------------------|----------------------------------------------------|--------------------------------------|-------------------|---------|----------|------------|
| pool1_NB01 | /data/fastq_pass/barcode01   | /data/references/NC_003310.fasta  | /data/primer_sequences/long_amplicon_primers.fasta | /data/references/NC_063383.1.fasta   | 3_SouthKivu_Apr24 | TRUE    | 30       | 1000       |
| pool1_NB02 | /data/fastq_pass/barcode02   | /data/references/NC_003310.fasta  | /data/primer_sequences/long_amplicon_primers.fasta | /data/references/NC_063383.1.fasta   | 5_SouthKivu_Apr24 | TRUE    | 30       | 1000       |
| pool1_NB03 | /data/fastq_pass/barcode03   | /data/references/NC_003310.fasta  | /data/primer_sequences/long_amplicon_primers.fasta | /data/references/NC_063383.1.fasta   | 6_SouthKivu_May   | TRUE    | 30       | 1000       |

Run the following command:

''snakemake --snakefile Snakefile_long.smk --directory results_long --config sample_data_long.tsv --cores 16''

### Analyze short amplicon data

Create a directory for the output of the workflow (e.g. results_short).

Create a config file (e.g. sample_data_short.tsv) containing the following columns:

| UniqueID     | FASTQ_path                   | Reference                         | Primers                                              | Primer_reference                             | Sequence_name     | Gzipped | Coverage | Min_length |
|--------------|------------------------------|-----------------------------------|------------------------------------------------------|----------------------------------------------|-------------------|---------|----------|------------|
| 3_pool1_NB01 | /data/fastq_pass/barcode01   | /data/references/NC_003310.fasta  | /data/references/short_amplicon_primers.fasta.fasta  | /data/references/ON585033.1_truncated.fasta  | 3_SouthKivu_Apr24 | TRUE    | 30       | 150        |
| 5_pool1_NB02 | /data/fastq_pass/barcode02   | /data/references/NC_003310.fasta  | /data/references/short_amplicon_primers.fasta.fasta  | /data/references/ON585033.1_truncated.fasta  | 5_SouthKivu_Apr24 | TRUE    | 30       | 150        |
| 6_pool1_NB03 | /data/fastq_pass/barcode03   | /data/references/NC_003310.fasta  | /data/references/short_amplicon_primers.fasta.fasta  | /data/references/ON585033.1_truncated.fasta  | 6_SouthKivu_May24 | TRUE    | 30       | 150        |

''snakemake --snakefile Snakefile_short.smk --directory results_short --config sample_data_short.tsv --cores 16''

### Analyze combined amplicon data

Create a directory for the output of the workflow (e.g. results_combined).

Create a config file (e.g. sample_data_combined.tsv) containing the following columns:

| LongID     | ShortID      | Long_path                                             | Short_path                                              | Sequence_name       | Reference                        | Coverage |
|------------|--------------|-------------------------------------------------------|---------------------------------------------------------|---------------------|----------------------------------|----------|
| pool1_NB01 | 3_pool1_NB01 | /data/results_long/trimmed/pool1_NB06_trimmed.fastq   | /data/results_short/trimmed/14_pool1_NB06_trimmed.fastq | 14_SouthKivu_May24  | /data/references/NC_003310.fasta | 30       |
| pool1_NB02 | 5_pool1_NB02 | /data/results_long/trimmed/pool1_NB07_trimmed.fastq   | /data/results_short/trimmed/15_pool1_NB07_trimmed.fastq | 15_SouthKivu_Oct23  | /data/references/NC_003310.fasta | 30       |
| pool1_NB03 | 6_pool1_NB03 | /data/results_long/trimmed/pool1_NB08_trimmed.fastq   | /data/results_short/trimmed/16_pool1_NB08_trimmed.fastq | 16_SouthKivu_Mar24  | /data/references/NC_003310.fasta | 30       |

''snakemake --snakefile Snakefile_combine.smk --directory results_combined --config sample_data_combined.tsv --cores 16''