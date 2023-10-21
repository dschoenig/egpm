#!/bin/bash
#SBATCH --account=def-cricrime 
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH --array=1-10
#SBATCH --mail-user=schonig.daniel@courrier.uqam.ca
#SBATCH --mail-type=ALL

module load StdEnv/2020 gcc/9.3.0 gdal/3.5.1 geos/3.10.2 python/3.10 udunits/2.2.28 r/4.3.1

Rscript 1_landscapes_noeffect.R imbalance_low noeff_imbalance_low $SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_COUNT