sbatch cc_5_compare_binary_high.sh
sbatch cc_5_compare_binary_low.sh
sbatch cc_5_compare_normal_high.sh
sbatch cc_5_compare_normal_low.sh
sbatch cc_5_compare_tweedie_high.sh
sbatch cc_5_compare_tweedie_low.sh


scp schoed@narval.computecanada.ca:scratch/egpm/results/comparison/mod.normal_low.rds ../results/comparison/
scp schoed@narval.computecanada.ca:scratch/egpm/results/comparison/mod.normal_high.rds ../results/comparison/
scp schoed@narval.computecanada.ca:scratch/egpm/results/comparison/mod.binary_low.rds ../results/comparison/
scp schoed@narval.computecanada.ca:scratch/egpm/results/comparison/mod.binary_high.rds ../results/comparison/
scp schoed@narval.computecanada.ca:scratch/egpm/results/comparison/mod.tweedie_low.rds ../results/comparison/
scp schoed@narval.computecanada.ca:scratch/egpm/results/comparison/mod.tweedie_high.rds ../results/comparison/
