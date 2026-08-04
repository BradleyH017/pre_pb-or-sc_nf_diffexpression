
# Description

This pipeline performs differential gene expression and gene set enrichment in droplet single cell RNA-seq. As input, this pipeline takes an H5AD single cell experiment, or a pre-pseudobulked h5ad (more convenient if testing many different contrasts while correcting for the same set of covariates). An example of how to construct the pseudobulk h5ad can be found in `bin/make_pseudobulk_adata.py`. 

NOTE: All constrasts / covariates must be present in the .obs of the H5AD, and if testing multiple contrasts against a core set of covariates, these can be supplied via the `target_variable_list_file` `.txt` file - one per line, '#' lines are ignored.

For differential gene expression, we support the following packages:
1. [MAST](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0844-5)
2. [edgeR](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2796818/)
3. [DESeq2](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8)

For gene set enrichment, we support the following packages:
1. [fGSEA](https://www.biorxiv.org/content/10.1101/060012v3)

For multiple correction, we support:
1. [Benjamini-Hochberg](https://www.jstor.org/stable/2346101?seq=1#metadata_info_tab_contents)
2. [IHW](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4930141/)

For a full description of possible parameters, see [README-params.md](docs/README-params.md).

For more information on running the pipeline, see:
* [Common problems](docs/README-common_problems.md)
* [Output directory structure](docs/README-outdir_structure.md)


# Quickstart

Quickstart for deploying this pipeline locally and on a high performance compute cluster.


## 1. Set up the environment

See [environment README](env/README.md) to set up environment. Once the environment is set up, activate the `conda` environment:

```bash
source activate sc_diff_expr
```

Alternatively, if using singularity or docker, one can pull the image from [henryjt/sc_nf_diffexpression:1.0.0](https://hub.docker.com/layers/196450988/henryjt/sc_nf_diffexpression/1.0.0/images/sha256-da59d053c402d3ba2f610488a91e5dead9a2821ac0bb565723ca5c9bef4f1d5e?context=repo).


## 2. Prepare the input files

Generate and/or edit input files for the pipeline.

The pipeline takes as input:
1. **--file_anndata**:  H5AD file containing sequencing data. Required.
2. **-params-file**:  YAML file containing analysis parameters. Required.

Examples of these files can be found in `demo/`.


## 3. Run pipeline

**NOTE**: All input file paths should be full paths.

To run:
```bash
nextflow run \
  "/path/to/repo/dir/main.nf" \
  -profile "local" \ # could also be configuration for cluster
  --file_anndata "/path/to/data/data.h5ad" \
  -params-file "/path/to/config/params.yml"
```

Examples:
* [Local with docker](demo/run_differential_expression_demo-docker.sh)
* [Local with singularity](demo/run_differential_expression_demo-singularity.sh)
* [Cluster with conda](demo/run_differential_expression_demo-conda.sh)


Authors: Henry Taylor and Leland Taylor
