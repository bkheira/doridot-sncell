#!/bin/bash

################################ Slurm options #################################

### Job name
#SBATCH --job-name=sncell-pipeline

### Output
#SBATCH --output=sncell-%j.out  # both STDOUT and STDERR
##SBATCH -o slurm.%N.%j.out  # STDOUT file with the Node name and the Job ID
##SBATCH -e slurm.%N.%j.err  # STDERR file with the Node name and the Job ID

### Limit run time "days-hours:minutes:seconds"
#SBATCH --time=24:00:00

### Requirements
#SBATCH --partition=fast
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=5GB

################################################################################

echo '########################################'
echo 'Date:' $(date --iso-8601=seconds)
echo 'User:' $USER
echo 'Host:' $HOSTNAME
echo 'Job Name:' $SLURM_JOB_NAME
echo 'Job Id:' $SLURM_JOB_ID
echo 'Directory:' $(pwd)
echo '########################################'
echo 'Pipeline SNCell version: v0.1.0'
echo '-------------------------'
echo 'Main module versions:'


start0=`date +%s`

# modules loading
module purge
module load conda snakemake/6.5.0
module load R
conda --version
python --version
echo 'snakemake' && snakemake --version

echo '-------------------------'
echo 'PATH:'
echo $PATH
echo '########################################'

snakemake --cluster "sbatch -c {threads} --mem={resources.mem_gb}G" --use-conda --jobs=10 --retries 3
