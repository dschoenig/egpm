#!/bin/bash
#SBATCH --account=def-cricrime 
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=7:00:00
#SBATCH --array=1-100
#SBATCH --mail-user=schonig.daniel@courrier.uqam.ca
#SBATCH --mail-type=ALL

module load StdEnv/2020 gcc/9.3.0 gdal/3.5.1 geos/3.10.2 python/3.10 udunits/2.2.28 r/4.2.2

Rscript 1_landscapes.R imbalance_high 1000 1000 0.3 $SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_COUNT
Rscript 1_landscapes.R test 1 100 0.3 1 1 4
