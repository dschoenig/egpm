#!/bin/bash
# Rerun block 2 (see rerun.txt): Tweedie/binary summaries and models, Gaussian marginals.
# Start only after all block 1 jobs have finished.
set -euo pipefail
cd "$(dirname "$0")"

scripts=(
  cc_3_summary_tweedie_low.sh
  cc_3_summary_tweedie_high.sh
  cc_3_summary_binary_low.sh
  cc_3_summary_binary_high.sh
  cc_3_summary_noeff_tweedie_low.sh
  cc_3_summary_noeff_tweedie_high.sh
  cc_3_summary_noeff_binary_low.sh
  cc_3_summary_noeff_binary_high.sh
  cc_2_models_tweedie_low_egp_som25.sh
  cc_2_models_tweedie_low_egp_som10.sh
  cc_2_models_tweedie_low_egp_som50.sh
  cc_2_models_tweedie_low_egp_som25_unequal.sh
  cc_2_models_tweedie_high_egp_som25.sh
  cc_2_models_tweedie_high_egp_som10.sh
  cc_2_models_tweedie_high_egp_som50.sh
  cc_2_models_tweedie_high_egp_som25_unequal.sh
  cc_2_models_binary_low_egp_som25.sh
  cc_2_models_binary_low_egp_som10.sh
  cc_2_models_binary_low_egp_som50.sh
  cc_2_models_binary_low_egp_som25_unequal.sh
  cc_2_models_binary_high_egp_som25.sh
  cc_2_models_binary_high_egp_som10.sh
  cc_2_models_binary_high_egp_som50.sh
  cc_2_models_binary_high_egp_som25_unequal.sh
  cc_2_models_noeff_tweedie_low_egp_som25.sh
  cc_2_models_noeff_tweedie_high_egp_som25.sh
  cc_2_models_noeff_binary_low_egp_som25.sh
  cc_2_models_noeff_binary_high_egp_som25.sh
  cc_3_marginal_low_egp_som25.sh
  cc_3_marginal_low_egp_som10.sh
  cc_3_marginal_low_egp_som50.sh
  cc_3_marginal_low_egp_som25_unequal.sh
  cc_3_marginal_high_egp_som25.sh
  cc_3_marginal_high_egp_som10.sh
  cc_3_marginal_high_egp_som50.sh
  cc_3_marginal_high_egp_som25_unequal.sh
  cc_3_marginal_noeff_low_egp_som25.sh
  cc_3_marginal_noeff_high_egp_som25.sh
)

for s in "${scripts[@]}"; do
  [ -f "$s" ] || { echo "Missing script: $s" >&2; exit 1; }
done

# 3_marginal_egp.R resumes from an existing .inter.rds file, which would mix
# old results into the rerun. Remove them for the Gaussian marginals rerun here.
for f in ../results/imbalance_{low,high}/egp_som{25,10,50,25_unequal}.inter.rds \
         ../results/noeff_imbalance_{low,high}/egp_som25.inter.rds; do
  if [ -f "$f" ]; then
    echo "Removing $f"
    rm -f "$f"
  fi
done

for s in "${scripts[@]}"; do
  echo "Submitting $s"
  sbatch "$s"
done
