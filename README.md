# TAP-seq workflow

This repository contains a [snakemake](https://snakemake.readthedocs.io/en/stable/index.html)
workflow that handles all data processing from fastq files to transcript counts and perturbation
status matrices from TAP-seq experiments.

The workflow can be downloaded by simply cloning the repository into a location of choice:
```
git clone https://github.com/argschwind/TAPseq_workflow.git
```

The workflow can be executed through snakemake and conda. This only requires that conda and
snakemake are installed. All other required dependencies will be installed through conda.

The data processing workflow consists of following steps:

## 1. Create TAP-seq alignment references
TAP-seq uses custom alignment references with added CROP-seq vector transcripts to the identify
perturbation of each cell. The rules in workflow `rules/create_alignment_refs.smk` allow creation of
such aligment references. The workflow supports references containing either the whole transcriptome
or only TAP-seq target gene loci. The type of the alignment reference is provided by its name. The 
alignment reference name is made according to following pattern: **`species_type_suffix`**. A
reference named `hg38_tapseq_ref` therefore specifies a TAP-seq, i.e. only containing target gene
loci based on the human genome hg38. `ref` serves as a suffix to distinguish different references of
the same type, e.g. references for different perturbation or target gene pools.

Alignment references for all samples are specifed in the first part of the `config.yml` file.
Parameters for reference creation are also found in the config file. CROP-seq vector and target
gene lists can be provided via files in `meta_data/cropseq_vectors` and
`meta_data/target_gene_panels`. Which lists are used for each specified alignment reference is
defined in the config file under: `step 1: create alignment reference`. See example files and
reference to learn more about this.

This example workflow specifies two human alignment references with the same CROP-seq vectors, but
one containing the whole transcriptome and one only example TAP-seq target gene loci. Alignment
references for e.g. mouse could be created by changing/adding urls for mm10 in
`download_genome_annot` in `config.yml` (for whole transcriptome) and changing/adding the mm10
BSgenome object in `create_tapseq_ref` (for target genes only). Using the species tag in alignment
reference names enables creating references for multiple species/genomes within one project.

All alignment references used by the workflow can be created using the following command. The
`--sjdbOverhang` STAR parameter might have to be adjusted in the `create_genomedir` section in the
config file.

```
# create all aligment references defined in the config file (--jobs = number of threads to use in
# parallel, please adjust; -n = dryrun, remove it to execute)
snakemake --use-conda --jobs 2 alignment_references -n
```

This uses parallel computing, but the number of available threads of course depends on your system.
By default STAR uses up to 5 threads, but this can be adjusted in the `create_genomedir` and
`star align` sections of the config file. Providing 10 cores would therefore mean that 2 processes
can be run entirely in parallel with each STAR process using 5 threads each.

## 2. Align reads
Reads can be aligned to created references using the snakemake rules in `rules/align_reads.smk`.
This is based on the [Drop-seq tools](http://mccarrolllab.org/dropseq/) workflow. Input paired end
fastq files for every sample are specified under `samples` in the config file. This assumes that
files for each sample are located in one directory and follow the naming scheme
`prefix_sample_1_sequence.txt.gz` and `prefix_sample_2_sequence.txt.gz` for read 1 and read 2 files.
This naming scheme will likely have to be adapted by changing the `get_fastq_files` function in
`rules/align_reads.smk`.

The cell number and alignment reference for each sample is specified by `cell_numbers` and
`align_ref` in the config file. Parameters for each workflow step can also be changed via the config
file.

Reads for all samples can be aligned by running following command.  This also creates an alignment
report for every sample in `results/alignment`.

```
# align reads for all samples
snakemake --use-conda --jobs 2 align_reads -n
```

## 3. Extract digital gene expression (DGE)
The last step of data processing consists of extracting transcript counts and the perturbation
status for each cell in every sample. The main output of this step are `dge.txt` and
`perturbation_status.txt` files for every sample containing the transcript counts and detected
CROP-seq vector perturbations per cell. This also applies chimeric read filtering as proposed by
[Dixit et al., 2016](https://www.biorxiv.org/content/10.1101/093237v1.full).

DGE data can be extracted for all samples using following command, which also creates DGE reports in
`results/dge`.

```
snakemake --use-conda --jobs 2 extract_dge -n
```

## Executing the whole workflow
The whole workflow can be exectuted for all samples at once using the snakemake "all" rule by simply
running:

```
snakemake --use-conda
```
