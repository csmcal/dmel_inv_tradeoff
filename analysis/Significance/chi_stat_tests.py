
# main_help = ""
"""
# For calculating a statistic (chi squared) of the differences in inversion frequencies between experimental cohorts,
#   and comparing that to an expected distribution of read counts generated from a model
#   using a maximum-likelihood estimated paternal population (not paternal sample) allele count
# CONSIDER 
# 1) alternate statistics like abs(obs-exp)/exp, or percentage point difference,
#    though these are likely monotonic with frequency distance and so give similar p-values
# 2) use the exact possible allele counts as the potential chi-value counts,
#    which would give a little jitter/error difference, but would calculate with much less memory/processor overhead
#
### A common testing call is:

  python chi2_pval_test.py --count_file ../Alignment/inv_freqs/toy_freq_data.csv --library_data ../Alignment/read_metadata.csv --od testing_temp/ --out_prefix test --calc_chi2 --verbose

### For testing multiple frequency change sets:

  python chi2_pval_test.py --stat_target two_changes_combined_cohorts --comp_cohorts embryo_comb_off_two_way --bootstrap 1000 --count_file ../Alignment/inv_freqs/toy_freq_data.csv --library_data ../Alignment/read_metadata.csv --od testing_temp/ --out_prefix test --verbose

### Main runs on the server:

python chi2_pval_test.py --stat_target two_changes_combined_cohorts --comp_cohorts embryo_all_comb_two_way --bootstrap 1000000 --count_file ../inv_freqs/freq_data_combined.csv --library_data ../read_metadata.csv --od comb_2way/ --out_prefix 2023-07-11

python chi2_pval_test.py --stat_target one_change_combined_cohorts --comp_cohorts all_including_combinations --bootstrap 1000000 --count_file ../inv_freqs/freq_data_combined.csv --library_data ../read_metadata.csv --od comb_1way/ --out_prefix 2023-07-11

python chi2_pval_test.py --stat_target one_change_combined_cohorts --comp_cohorts ecl_sex_comb --bootstrap 1000000 --count_file ../inv_freqs/freq_data_combined.csv --library_data ../read_metadata.csv --od comb_1way/ --out_prefix 2023-07-23

python chi2_pval_test.py --stat_target one_change_combined_cohorts --comp_cohorts ecl_sex_comb --count_file ../inv_freqs/freq_data_combined.csv --library_data ../read_metadata.csv --od comb_1way/ --out_prefix 2023-07-23 --bootstrap 1000000
"""

# Import packages
import numpy as np
import argparse
import os
import scipy as sp
from scipy.stats import binom
from scipy.stats import hypergeom
from functools import partial



# A function running the basic chi statistic calculation on an observation and expectation
def chi(observed,expected):
    return ((observed-expected)**2)/expected

# calculates the chi statistic for a specific observation, with vectorized opertaions
def calc_chi_statistic(read_counts,read_totals,expected_counts,allele_totals=None,count_mode="read",
        combination_index=None,combine=False):
    # Check the inputs
    if allele_totals is None and count_mode == "allele":
        raise ValueError("can only calculate an allele-count chi statistic when allele_totals are provided")
    if combination_index is None and combine:
        raise ValueError("can only combine counts when the combination_index at which to combine into two groups is provided")
    if combine and count_mode == "read":
        raise ValueError("only expected to combine the 1...n other counts if considering allele_counts")
    if count_mode != "read" and count_mode != "allele":
        raise ValueError("chi statistic calculation count_mode must only be one of \"read\" or \"allele\"")

    # Generate the array of counts to compare, either reads or alleles,
    #  and combine the allele counts on either side of the combination index if specified
    if count_mode == "read":
        pos_counts = np.array(read_counts)
        neg_counts = np.subtract(read_totals,pos_counts)
        exp_pos_counts = np.array(expected_counts)
        exp_neg_counts = np.subtract(read_totals,exp_pos_counts)
    elif count_mode == "allele":
        expected_freqs = np.divide(expected_counts,read_totals)
        pos_proportion = np.divide(read_counts,read_totals)
        pos_counts = np.multiply(pos_proportion,allele_totals)
        neg_counts = np.subtract(allele_totals,pos_counts)
        exp_pos_counts = np.multiply(expected_freqs,allele_totals)
        exp_neg_counts = np.subtract(allele_totals,exp_pos_counts)
        if combine:
            last_dim = len(pos_counts.shape)-1
            pos_counts = np.stack((np.sum(pos_counts[...,:combination_index],axis=last_dim),
                np.sum(pos_counts[...,combination_index:],axis=last_dim)),axis=last_dim)
            neg_counts = np.stack((np.sum(neg_counts[...,:combination_index],axis=last_dim),
                np.sum(neg_counts[...,combination_index:],axis=last_dim)),axis=last_dim)
            last_dim_exp = len(exp_pos_counts.shape)-1
            exp_pos_counts = np.stack((np.sum(exp_pos_counts[...,:combination_index],axis=last_dim_exp),
                np.sum(exp_pos_counts[...,combination_index:],axis=last_dim_exp)),axis=last_dim_exp)
            exp_neg_counts = np.stack((np.sum(exp_neg_counts[...,:combination_index],axis=last_dim_exp),
                np.sum(exp_neg_counts[...,combination_index:],axis=last_dim_exp)),axis=last_dim_exp)
    else:
        raise ValueError("Unexpected count_mode value: "+str(count_mode))

    # Calculate the sum of the chi statistic for all relevant comparisons
    last_dim = len(pos_counts.shape)-1
    chi_sum = np.add(np.sum(chi(pos_counts,exp_pos_counts),axis=last_dim),np.sum(chi(neg_counts,exp_neg_counts),axis=last_dim))

    return chi_sum



