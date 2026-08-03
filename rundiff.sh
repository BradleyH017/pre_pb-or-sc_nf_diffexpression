#! bin/bash


# Remove old logs but not the most previous run
# rm -r *html.*;
# rm .nextflow.log.*;
# rm flowchart.png.*;
# rm trace.txt.*;

# SPECIFY ANNDATA
ANNDATA_PATH="/lustre/scratch125/humgen/projects_v2/ibdresponse/analysis/bradley_analysis/archetypal_analysis_from_mrvi/results/objects/adata_plus_archetypes.h5ad"

# Add Singularity to path
# PATH=$PATH:/software/singularity/v3.10.0/bin
module load cellgen/singularity
module load nextflow-23.10.0

# Nextflow settings
export NXF_OPTS="-Xms25G -Xmx25G"
# Uncomment this if get strange bus errors
# export NXF_OPTS="${NXF_OPTS} -Dleveldb.mmap=false" # No resume functionality
export NXF_HOME=$(pwd)
export NXF_WORK="${NXF_HOME}/work"
export NXF_TEMP="${NXF_HOME}/nextflow_temp"
export NXF_CONDA_CACHEDIR="${NXF_HOME}/nextflow_conda"
export NXF_SINGULARITY_CACHEDIR="${NXF_HOME}/cache_singularity"

# Farm specific settings
#export QT_QPA_PLATFORM='offscreen'
export LSB_DEFAULT_JOBGROUP="/${USER}/nf"
# export LSB_DEFAULTGROUP="sc-eqtl-ibd"
# export LSB_DEFAULTGROUP="cnv_15x"
export LSB_DEFAULTGROUP="team152"
export PYTHONHASHSEED=0

mkdir $NXF_TEMP
# mkdir results
# cp ${INPUT} results/
# cp ${EQTL_REPO}/conf/extra_confs/* results/

# To stop TclError
#export QT_QPA_PLATFORM='offscreen'

# Run nextflow
nextflow run "./main.nf" -profile "lsf" --file_anndata "$ANNDATA_PATH" -params-file "params.yml" -resume

# bash rundiff.sh
