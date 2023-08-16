#!/bin/bash
#SBATCH --account=def-cricrime 
#SBATCH --cpus-per-task=32
#SBATCH --mem=24G
#SBATCH --time=6:00:00
#SBATCH --mail-user=schonig.daniel@courrier.uqam.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=brm_binary_high

module load StdEnv/2020 gcc/9.3.0 gdal/3.5.1 geos/3.10.2 python/3.10 udunits/2.2.28 r/4.2.2

Rscript 5_compare_brm.R binary_high
