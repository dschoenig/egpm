#!/bin/bash
#SBATCH --account=def-cricrime 
#SBATCH --cpus-per-task=4
#SBATCH --mem=6G
#SBATCH --time=1:30:00
#SBATCH --array=1-100
#SBATCH --mail-user=schonig.daniel@courrier.uqam.ca
#SBATCH --mail-type=ALL

module load StdEnv/2020 gcc/9.3.0 gdal/3.4.3 geos/3.10.2 python/3.10 udunits/2.2.28 r/4.2.1

Rscript 2_egp_100_sam1_som50.R $SLURM_CPUS_PER_TASK $SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_COUNT
