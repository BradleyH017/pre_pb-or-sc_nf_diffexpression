### Bradley set up
Get code 
```
git clone https://github.com/andersonlab/sc_nf_diffexpression
```


1. Getting singularity\\
```
module load ISG/singularity/3.11.4
export SINGULARITY_CACHEDIR=$PWD/.singularity_cache
mkdir -p "$SINGULARITY_CACHEDIR"
singularity pull docker://henryjt/sc_nf_diffexpression:1.0.0

# update the anndata
cp sc_nf_diffexpression_1.0.0.sif sc_nf_diffexpression_1.0.0-update.sif
mkdir -p ${PIPELINE_DIR}/packages/python_libs
module load ISG/singularity/3.11.4 
singularity exec -B /lustre:/lustre ./sc_nf_diffexpression_1.0.0-update.sif pip install --target=./packages/python_libs "anndata==0.9.2"


rm -rf /lustre/scratch125/humgen/projects_v2/ibdresponse/analysis/bradley_analysis/sc_nf_diffexpression/packages/python_libs/*

PIPELINE_DIR="/lustre/scratch125/humgen/projects_v2/ibdresponse/analysis/bradley_analysis/sc_nf_diffexpression"
module load ISG/singularity/3.11.4
singularity exec -B /lustre ${PIPELINE_DIR}/sc_nf_diffexpression_1.0.0-update.sif \
    pip install --target=${PIPELINE_DIR}/packages/python_libs \
    "numpy==1.24.4" "pandas==2.0.3" "anndata==0.9.2" "scanpy==1.9.8"


```

2. Then update the `params.yml` file

3. Then run (in a tmux)
```bash
bash rundiff.sh
```


