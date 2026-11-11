#!/bin/bash
# Rerun block 5 of 6 (508 jobs): binary high models, marginals for the models of block 4.
set -euo pipefail
cd "$(dirname "$0")"
source rerun_lib.sh

scripts=(
  cc_2_models_binary_high_egp_som25.sh
  cc_2_models_binary_high_egp_som10.sh
  cc_2_models_binary_high_egp_som50.sh
  cc_2_models_binary_high_egp_som25_unequal.sh
  cc_3_marginal_binary_low_egp_som25.sh
  cc_3_marginal_binary_low_egp_som10.sh
  cc_3_marginal_binary_low_egp_som50.sh
  cc_3_marginal_binary_low_egp_som25_unequal.sh
  cc_3_marginal_noeff_tweedie_low_egp_som25.sh
  cc_3_marginal_noeff_tweedie_high_egp_som25.sh
  cc_3_marginal_noeff_binary_low_egp_som25.sh
  cc_3_marginal_noeff_binary_high_egp_som25.sh
)
prev_scripts=(
  cc_2_models_binary_low_egp_som25.sh
  cc_2_models_binary_low_egp_som10.sh
  cc_2_models_binary_low_egp_som50.sh
  cc_2_models_binary_low_egp_som25_unequal.sh
  cc_2_models_noeff_tweedie_low_egp_som25.sh
  cc_2_models_noeff_tweedie_high_egp_som25.sh
  cc_2_models_noeff_binary_low_egp_som25.sh
  cc_2_models_noeff_binary_high_egp_som25.sh
  cc_3_summary_noeff_tweedie_low.sh
  cc_3_summary_noeff_tweedie_high.sh
  cc_3_summary_noeff_binary_low.sh
  cc_3_summary_noeff_binary_high.sh
  cc_3_marginal_noeff_low_egp_som25.sh
  cc_3_marginal_noeff_high_egp_som25.sh
  cc_3_marginal_tweedie_high_egp_som25.sh
  cc_3_marginal_tweedie_high_egp_som10.sh
  cc_3_marginal_tweedie_high_egp_som50.sh
  cc_3_marginal_tweedie_high_egp_som25_unequal.sh
)
inter_files=(
  ../results/binary_imbalance_low/egp_som25.inter.rds
  ../results/binary_imbalance_low/egp_som10.inter.rds
  ../results/binary_imbalance_low/egp_som50.inter.rds
  ../results/binary_imbalance_low/egp_som25_unequal.inter.rds
  ../results/noeff_tweedie_imbalance_low/egp_som25.inter.rds
  ../results/noeff_tweedie_imbalance_high/egp_som25.inter.rds
  ../results/noeff_binary_imbalance_low/egp_som25.inter.rds
  ../results/noeff_binary_imbalance_high/egp_som25.inter.rds
)

submit_block
