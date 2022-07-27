#!/bin/bash
#SBATCH --account=def-cricrime 
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=4:00:00
#SBATCH --mail-user=schonig.daniel@courrier.uqam.ca
#SBATCH --mail-type=ALL

module load StdEnv/2020 gcc/9.3.0 gdal/3.4.3 geos/3.10.2 python/3.10 udunits/2.2.28 r/4.2.1

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
Rscript test_200_sam0.01_som10_e1e4_k50.R 4
