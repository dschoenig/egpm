#!/bin/bash
#SBATCH --account=def-cricrime 
#SBATCH --cpus-per-task=4
#SBATCH --mem=6G
#SBATCH --time=0:45:00
#SBATCH --mail-user=schonig.daniel@courrier.uqam.ca
#SBATCH --mail-type=ALL

module load StdEnv/2023 gcc/12.3 gdal/3.4.3 geos/3.12.0 python/3.11.5 udunits/2.2.28 r/4.2.1

Rscript 2_egp_1000_sam0.01_som50.R $SLURM_CPUS_PER_TASK 1 100
