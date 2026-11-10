#!/bin/bash
# Rerun block 2 of 6 (815 jobs): remaining Gaussian models, Tweedie low models,
# landscape summaries, marginals for the Gaussian models of block 1.
set -euo pipefail
cd "$(dirname "$0")"
source rerun_lib.sh

scripts=(
  cc_2_models_high_egp_som25_unequal.sh
  cc_2_models_noeff_low_egp_som25.sh
  cc_2_models_noeff_high_egp_som25.sh
  cc_2_models_tweedie_low_egp_som25.sh
  cc_2_models_tweedie_low_egp_som10.sh
  cc_2_models_tweedie_low_egp_som50.sh
  cc_2_models_tweedie_low_egp_som25_unequal.sh
  cc_3_summary_tweedie_low.sh
  cc_3_summary_tweedie_high.sh
  cc_3_summary_binary_low.sh
  cc_3_summary_binary_high.sh
  cc_3_summary_noeff_tweedie_low.sh
  cc_3_summary_noeff_tweedie_high.sh
  cc_3_summary_noeff_binary_low.sh
  cc_3_summary_noeff_binary_high.sh
  cc_3_marginal_low_egp_som25.sh
  cc_3_marginal_low_egp_som10.sh
  cc_3_marginal_low_egp_som50.sh
  cc_3_marginal_low_egp_som25_unequal.sh
  cc_3_marginal_high_egp_som25.sh
  cc_3_marginal_high_egp_som10.sh
  cc_3_marginal_high_egp_som50.sh
)
prev_scripts=(
  cc_1_landscapes_low_tweedie.sh
  cc_1_landscapes_high_tweedie.sh
  cc_1_landscapes_low_binary.sh
  cc_1_landscapes_high_binary.sh
  cc_1_landscapes_noeff_low_tweedie.sh
  cc_1_landscapes_noeff_high_tweedie.sh
  cc_1_landscapes_noeff_low_binary.sh
  cc_1_landscapes_noeff_high_binary.sh
  cc_2_models_low_egp_som25.sh
  cc_2_models_low_egp_som10.sh
  cc_2_models_low_egp_som50.sh
  cc_2_models_low_egp_som25_unequal.sh
  cc_2_models_high_egp_som25.sh
  cc_2_models_high_egp_som10.sh
  cc_2_models_high_egp_som50.sh
)
inter_files=(
  ../results/imbalance_low/egp_som25.inter.rds
  ../results/imbalance_low/egp_som10.inter.rds
  ../results/imbalance_low/egp_som50.inter.rds
  ../results/imbalance_low/egp_som25_unequal.inter.rds
  ../results/imbalance_high/egp_som25.inter.rds
  ../results/imbalance_high/egp_som10.inter.rds
  ../results/imbalance_high/egp_som50.inter.rds
)

submit_block
