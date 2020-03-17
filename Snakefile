## run all workflows for the project

# config file containing samples and parameters
configfile: "config.yml"

# import workflows
include: "rules/create_alignment_refs.smk"
include: "rules/align_reads.smk"
include: "rules/extract_dge.smk"

# ALL RULE -----------------------------------------------------------------------------------------

# run whole workflow to align reads of all samples and extract dge for specified number of cells
# (samples and cells are read from config file)
rule all:
  input:
    align = expand("results/alignment/{sample}_align_report.html", sample = config["samples"]),
    dge = expand("results/dge/{sample}_dge_report.html", sample = config["samples"])

# create whole-genome and tap-seq alignment references
rule alignment_references:
  input:
    expand("data/alignment_references/{align_ref}/{ref_file}",
      align_ref = config["create_vector_ref"]["vector_fasta"],
      ref_file = ["cropseq_ref.dict", "cropseq_ref.refFlat", "genomeDir"])

# run workflow until read alignment. this allows to determine the number of cells before extracting
# dge data and finishing the workflow
rule align:
  input:
    expand("results/alignment/{sample}_align_report.html", sample = config["samples"])
    
# run workflow until digital gene expression (DGE) extraction. this concludes the data processing
# part of the project
rule dge:
  input:
    align = expand("results/alignment/{sample}_align_report.html", sample = config["samples"]),
    dge = expand("results/dge/{sample}_dge_report.html", sample = config["samples"])
