#!/bin/bash
# Rerun block 3 (see rerun.txt): Tweedie/binary marginals.
# Start only after all block 2 jobs have finished.
set -euo pipefail
cd "$(dirname "$0")"

scripts=(
  cc_3_marginal_tweedie_low_egp_som25.sh
  cc_3_marginal_tweedie_low_egp_som10.sh
  cc_3_marginal_tweedie_low_egp_som50.sh
  cc_3_marginal_tweedie_low_egp_som25_unequal.sh
  cc_3_marginal_tweedie_high_egp_som25.sh
  cc_3_marginal_tweedie_high_egp_som10.sh
  cc_3_marginal_tweedie_high_egp_som50.sh
  cc_3_marginal_tweedie_high_egp_som25_unequal.sh
  cc_3_marginal_binary_low_egp_som25.sh
  cc_3_marginal_binary_low_egp_som10.sh
  cc_3_marginal_binary_low_egp_som50.sh
  cc_3_marginal_binary_low_egp_som25_unequal.sh
  cc_3_marginal_binary_high_egp_som25.sh
  cc_3_marginal_binary_high_egp_som10.sh
  cc_3_marginal_binary_high_egp_som50.sh
  cc_3_marginal_binary_high_egp_som25_unequal.sh
  cc_3_marginal_noeff_tweedie_low_egp_som25.sh
  cc_3_marginal_noeff_tweedie_high_egp_som25.sh
  cc_3_marginal_noeff_binary_low_egp_som25.sh
  cc_3_marginal_noeff_binary_high_egp_som25.sh
)

for s in "${scripts[@]}"; do
  [ -f "$s" ] || { echo "Missing script: $s" >&2; exit 1; }
done

# 3_marginal_egp.R resumes from an existing .inter.rds file, which would mix
# old results into the rerun. Remove them for the marginals rerun here.
for f in ../results/{tweedie,binary}_imbalance_{low,high}/egp_som{25,10,50,25_unequal}.inter.rds \
         ../results/noeff_{tweedie,binary}_imbalance_{low,high}/egp_som25.inter.rds; do
  if [ -f "$f" ]; then
    echo "Removing $f"
    rm -f "$f"
  fi
done

for s in "${scripts[@]}"; do
  echo "Submitting $s"
  sbatch "$s"
done
