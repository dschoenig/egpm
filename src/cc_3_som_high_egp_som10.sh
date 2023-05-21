#!/bin/bash
#SBATCH --account=def-cricrime 
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=9:00:00
#SBATCH --mail-user=schonig.daniel@courrier.uqam.ca
#SBATCH --mail-type=ALL

module load StdEnv/2020 gcc/9.3.0 gdal/3.5.1 geos/3.10.2 python/3.10 udunits/2.2.28 r/4.2.2

Rscript 3_embedding_som.R imbalance_high egp_som10 embedding
