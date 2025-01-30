
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

# Define and collect input arguments
def parse_args():
    parser = argparse.ArgumentParser()

    # Create a new group to store required arguments
    required_args = parser.add_argument_group('required named arguments')

    # Add required arguments to the parser
    required_args.add_argument('--count_file', type=str, required=True,
                        help='The .csv file containing the estimated inverted and standard read counts, etc.')
                        # ../Alignment/inv_freqs/freq_data.csv
    required_args.add_argument('--library_data', type=str, required=True,
                        help='The .csv file with the library metadata')
                        # ../Alignment/read_metadata.csv
    required_args.add_argument('--od', type=str, required=True,
                        help='The output directory to write reshaped data, confidence interval data, and bootstrap replicates to')
                        # ./

    # Add optional arguments to the parser
    parser.add_argument('--varN', type=int, default=316,
                        help='The binomial sampling size expected to account for sample-independent library-prep-dependent variance')
    parser.add_argument('--alpha', type=float, default=0.05,
                        help='The single-tailed global alpha level for multiple-test correction')
    parser.add_argument('--count', type=str, default='allele', choices=('allele','read'),
                        help='A mode variable determining the form of count used in the statistic calculations')
    parser.add_argument('--stat_target', type=str, default='one_change_combined_cohorts', 
                        choices=('one_change_combined_cohorts',
                            'one_change_lone_cohorts',\
                            'two_changes_combined_cohorts',\
                            'one_and_two_changes_combined_cohorts'),
                        help='A mode variable determining the form of frequncy change(s) that will be used')
    parser.add_argument('--stat_tail', type=str, default='two', choices=('one',
                            'two',\
                            'one_and_two'),
                        help='A mode variable determining whether to test for more extreme stats '+\
                        'in the direction of the observation, or all directions')
    parser.add_argument('--cust_mt', type=str, 
                        default=("Benjamini-Yekutieli"), nargs='+',
                        # default=("Bonferroni","Dunn-Sidak","Benjamini-Hochberg","Benjamini-Yekutieli"), nargs='+',
                        help='When using custom multiple testing mode, the script will only attempt MT based on the named tests')
    parser.add_argument('--comp_cohorts', type=str, default='embryo_combined_offspring',
                        choices=('embryo_combined_offspring','embryo_combined_sexes','embryo_combined_eclosion','embryo_all',\
                            'embryo_all_combinations','all_individual','all_including_combinations','ecl_sex_comb','all_exper_comb',\
                            'embryo_comb_off_two_way','embryo_all_comb_two_way','custom'),
                        help='A mode variable determining the set of cohort comparisons to make within the data')
    parser.add_argument('--cust_comps', type=str, default=None, nargs='+',
                        help='When using custom comparison mode, the comparisons to make: expects comparisons in the form group1-group2')
    parser.add_argument('--cust_lines', type=str, default=None, nargs='+',
                        help='When using custom comparison mode, the experimental lines to use in the comparison')
    parser.add_argument('--cust_invs', type=str, default=None, nargs='+',
                        help='When using custom comparison mode, the experimental inversions to test within an experiment')
    parser.add_argument('--source_pat_from_obs', default=False, action='store_true',
                        help='--source_pat_from_obs causes the script to use the observed paternal allele count, instead of '+\
                            'estimating a maximum likelihood paternal allele count '+\
                            'from the neutral experimental model including the offspring cohorts to be tested. '+\
                            'Estimating the neutral model paternal frequencies from tested cohorts is more conservative '+\
                            'and avoids having improper chi-statistic calculations due to selective effects in the observed data not being tested at the time. '+\
                            'Calculation of selection coefficients ignores this and always uses the count estimates from the '+\
                            'observed paternal sequencing.')

    parser.add_argument('--out_prefix', type=str, default='',
                        help='The output file prefix to use')
    parser.add_argument('--write_sep_no_mt', default=False, action='store_true',
                        help='A flag telling the script to stop after pval calculation and write each pval to a unique text file')
    parser.add_argument('--mt_from_sep', default=False, action='store_true',
                        help='A flag telling the script to read separate files from the output directory for each p_val expected '+\
                            'instead of calculating them, then run multiple testing analysis')
    parser.add_argument('--run_sep_nohup', default=False, action='store_true',
                        help='A flag telling the script to spawn separate nohup python instances to run each comparison/line/inv combination '+\
                            'and write them to csv separately, instead of calculating each in this instanace of the script/using '+\
                            'the python multiprocessing package. The script will not run analyses for comparisons that already have '+\
                            'output files of the expected form')

    # Bootstrap, test separately that both cohort changes are codirectional in each direction: the probability that the chi statistic is  and that both have chi stats at least as 
    parser.add_argument('--bootstrap', type=int, default=None,
                        help='Set the script to approximately calculate p-values by random sampling the given number of saples from '+\
                            'the estimated neutral distribution, rather than considering the probabilities of all observations across that distribution')
    parser.add_argument('--bs_chunk', type=int, default=10000,
                        help='The number of bootstrap samples to process in a vectorized manner at once; '+\
                            'if larger than the total bootstrap sample size, calculates all samples at once')
    parser.add_argument('--bs_precalc_dist', default=False, action='store_true',
                        help='--bs_precalc_dist causes the script to calculate and store the probabilities '+\
                            'for the entire set of potential counts for use in sampling')
    # Selection coefficient ABC acceptance tolerance in the difference in observed vs expected frequency
    #    sel. coef. ABC uses the --bootstrap value for the bootstrap size
    parser.add_argument('--s_tol', type=float, default=0.01,
                        help='The maximum difference between any sampled frequency and the observed frequency '+\
                        'in order to accept the sample during selection coefficient ABC sampling')
    parser.add_argument('--s_quantile', type=float, default=0.05,
                        help='The proportion of the distribution to cut off in presenting the high and low quantiles of accepted selection coefficient estimates; '+\
                        'the default is to use a 95% interval, specified by the default --s_quantile of 0.05, for 0.025 and 0.975 quantile values')

    parser.add_argument('--calc_s', default=False, action='store_true',
                        help='--calc_s causes the script to calculate selection coefficients for the comparisons '+\
                        'using a simple Wright Fisher model and more realistic ABC estimations. '+\
                        'If no specific calculation is flagged, the script defaults to attempting to calculate all.')
    parser.add_argument('--calc_prod', default=False, action='store_true',
                        help='--calc_prod causes the script to calculate p-values based on the product '+\
                        'of frequency changes in multiple-change comparisons; requires multiple comparison. '+\
                        'If no specific calculation is flagged, the script defaults to attempting to calculate all.')
    parser.add_argument('--calc_chi2', default=False, action='store_true',
                        help='--calc_chi2 causes the script to calculate p-values for the comparisons '+\
                        'based on the chi-squared statistic of the inverted and non-inverted reads. '+\
                        'If no specific calculation is flagged, the script defaults to attempting to calculate all.')

    parser.add_argument('--verbose', default=False, action='store_true',
                        help='--verbose causes the script to print out descriptions of the current computation')
    parser.add_argument('--nomp', default=True, action='store_false',
                        help='--nomp causes the script to forgo using the multiprocessing python package to compute p-values asynchronously')
    parser.add_argument('--num_proc', type=int, default=0,
                        help='The number of processors to use in multiprocessing, overrides prop_proc when set > 0')
    parser.add_argument('--prop_proc', type=float, default=0.34,
                        help='The proportion of available processors to use in multiprocessing (scaled down), default is 0.34')

    parser.add_argument('--mt_comps_comb', default=False, action='store_true',
                        help='--mt_comps_comb causes the script to perform multiple test corrections across categories of comparison instead of separating them')

    args = parser.parse_args()
    return args



# Define the (internal shorthand) names for expermental cohorts and their groupings under consideration
cohort_group_dict = {
    "pat"   : [0],
    "emb"   : [1],
    "eaf"   : [2],
    "laf"   : [3],
    "eam"   : [4],
    "lam"   : [5],
    "foff"  : [2,3],
    "moff"  : [4,5],
    "eoff"  : [2,4],
    "loff"  : [3,5],
    "off"   : [2,3,4,5]
}

# A dictionary to convert shorthand cohort/group names to output names
group_full_name_dict = {
    "pat"   : "Paternal",
    "emb"   : "Embryo",
    "eaf"   : "Early-Eclosing Adult Female Offspring",
    "laf"   : "Late-Eclosing Adult Female Offspring",
    "eam"   : "Early-Eclosing Adult Male Offspring",
    "lam"   : "Late-Eclosing Adult Male Offspring",
    "foff"  : "Adult Female Offspring",
    "moff"  : "Adult Male Offspring",
    "eoff"  : "Early-Eclosing Adult Offspring",
    "loff"  : "Late-Eclosing Adult Offspring",
    "off"   : "Adult Offspring"
    # "off"   : "All Adult Offspring"
}

# A dictionary to convert shorthand statistic names to output file statistic headers
stat_full_name_dict = {
    ("chi2_pval","one_change_combined_cohorts")      : "Chi_squared_on_combined_",
    ("chi2_pval","one_change_lone_cohorts")          : "Chi_squared_on_separate_",
    ("chi2_pval","two_changes_combined_cohorts")     : "Chi_squared_on_two_cohorts_testing_first_decrease_second_increase_with_combined_",
    ("prod_pval","two_changes_combined_cohorts")     : "Frequency_product_on_two_cohorts_testing_first_decrease_second_increase_with_combined_",
    ("chi2_pval","chi2_combined_two_way_ben_cost")   : "Chi_squared_on_two_cohorts_testing_first_increase_second_decrease_with_combined_",
    ("sel_coeff","one_change_combined_cohorts")      : "Selection_coefficients_from_combined_",
    ("sel_coeff","one_change_lone_cohorts")          : "Selection_coefficients_from_separate_",
    ("sel_coeff","two_changes_combined_cohorts")     : "Selection_coefficients_from_two_frequency_changes_with_combined_"
}

# A list of stat modes with multiple comparisons needed, for reference in the script
mult_comp_stats = ['two_changes_combined_cohorts','chi2_combined_two_way_ben_cost']

# Define the set of paternal census sizes, as a dictionary on inbred maternal line
pat_census_dict = {
    "192"  : 322,
    "251"  : 528,
    "254"  : 666,
    "418"  : 474
}

# Define the set of accepted maternal line names, and their internal representations
line_name_dict = {
    "192"     : "192",
    "251"     : "251",
    "254"     : "254",
    "418"     : "418",
    "ZI192N"  : "192",
    "ZI251N"  : "251",
    "ZI254N"  : "254",
    "ZI418N"  : "418",
    "192N"    : "192",
    "251N"    : "251",
    "254N"    : "254",
    "418N"    : "418"
}

# Define the set of accepted inversion names, and their internal representations
inv_name_dict = {
    "In(2L)t"   : "t",
    "In(2R)NS"  : "NS",
    "In(3L)Ok"  : "Ok",
    "In(3R)K"   : "K",
    "2Lt"       : "t",
    "2RNS"      : "NS",
    "3LOk"      : "Ok",
    "3RK"       : "K",
    "t"         : "t",
    "NS"        : "NS",
    "Ok"        : "Ok",
    "K"         : "K",
    "2L"        : "t",
    "2R"        : "NS",
    "3L"        : "Ok",
    "3R"        : "K"
}

# Define the chromosome on which each inversion resides
inv_chrom_lookup = {
    "t" : "2L",
    "NS" : "2R",
    "Ok" : "3L",
    "K" : "3R"
}


# Define the comparison pairs keyed to the input mode, 
comp_mode_dict = {
    "embryo_combined_offspring"  : [("pat","emb"),("emb","off")],
    "embryo_combined_sexes"      : [("pat","emb"),("emb","foff"),("emb","moff")],
    "embryo_combined_eclosion"   : [("pat","emb"),("emb","eoff"),("emb","loff")],
    "all_exper_comb"             : [("pat","emb"),("emb","off"),("eoff","loff"),("foff","moff")],
    "embryo_all"                 : [("pat","emb"),("emb","eaf"),("emb","laf"),("emb","eam"),("emb","lam")],
    "embryo_all_combinations"    : [("pat","emb"),("emb","eaf"),("emb","laf"),("emb","eam"),("emb","lam"),
                                    ("emb","off"),("emb","foff"),("emb","moff"),("emb","eoff"),("emb","loff")],
    "no_pat_off_all_indiv"       : [("pat","emb"),("emb","eaf"),("emb","laf"),("emb","eam"),("emb","lam"),
                                    ("eaf","laf"),("eaf","eam"),("eaf","lam"),
                                    ("laf","eam"),("laf","lam")],
    "all_individual"             : [("pat","emb"),("emb","eaf"),("emb","laf"),("emb","eam"),("emb","lam"),
                                    ("pat","eaf"),("pat","laf"),("pat","eam"),("pat","lam"),
                                    ("eaf","laf"),("eaf","eam"),("eaf","lam"),
                                    ("laf","eam"),("laf","lam"),
                                    ("eam","lam")],
    "all_including_combinations" : [("pat","emb"),("emb","eaf"),("emb","laf"),("emb","eam"),("emb","lam"),
                                    ("pat","eaf"),("pat","laf"),("pat","eam"),("pat","lam"),
                                    ("eaf","laf"),("eaf","eam"),("eaf","lam"),
                                    ("laf","eam"),("laf","lam"),
                                    ("eam","lam"),
                                    ("emb","off"),
                                    ("emb","foff"),("emb","moff"),
                                    ("emb","eoff"),("emb","loff"),
                                    ("pat","off"),
                                    ("pat","foff"),("pat","moff"),
                                    ("pat","eoff"),("pat","loff"),
                                    ("eoff","loff"),("foff","moff")],
    "ecl_sex_comb"               : [("eoff","loff"),("foff","moff")],
    "embryo_comb_off_two_way"    : [(("pat","emb"),("emb","off"))],
    "embryo_all_comb_two_way"    : [(("pat","emb"),("emb","off")),
                                    (("pat","emb"),("emb","eaf")),(("pat","emb"),("emb","laf")),
                                    (("pat","emb"),("emb","eam")),(("pat","emb"),("emb","lam")),
                                    (("pat","emb"),("emb","moff")),(("pat","emb"),("emb","foff")),
                                    (("pat","emb"),("emb","eoff")),(("pat","emb"),("emb","loff")),
                                    (("emb","eaf"),("emb","laf")),(("emb","eam"),("emb","lam")),
                                    (("emb","eaf"),("emb","eam")),(("emb","laf"),("emb","lam")),
                                    (("emb","moff"),("emb","foff")),(("emb","eoff"),("emb","loff"))],
}

# For counting chromosome proportions to compare to expected proportions
class HashableTuple(object): 
 
 def __init__(self, val):
    self.val = val 
 
 def __hash__(self):
    return hash(str(self.val)) 
 
 def __repr__(self):
    # define this method to get clean output
    return str(self.val) 
 
 def __eq__(self, other):
    return str(self.val) == str(other.val)


# For printing a sample of generated gametes
def __printGameteSample(individual,N):
    for i in np.arange(N):
        print(individual.genGamete())
    return

# A helper function to check if the supplied custom comparison is acceptable
def check_cust_comp(comp):
    comparison = tuple([c.strip() for c in comp.split('-')])
    if len(comparison) != 2:
        raise ValueError("Script is only configured to generate chi square statistics that compare two groupings "+\
            "of experimental cohorts at a time, "+string(len(comparison))+\
            " found in custom comparison \""+comp+"\"\n")
    if comparison[0] not in cohort_group_dict or comparison[1] not in cohort_group_dict:
        raise ValueError("Custom comparisons only configured to compare certain named cohorts or groupings:\n"+\
            string(cohort_group_dict.keys())+"\nInput custom comparison \""+comp+"\" does not match\n")
    return comparison


# Extracts data about the paired-end read pools to be used in analysis
# Includes file names to analyze
def parse_library_metadata_to_dict(read_metadata_file_path):
    read_metadata = {}
    with open(read_metadata_file_path) as read_meta_file:
        lines = read_meta_file.readlines()
        for line in lines[1:]:
            data = line.strip().split(',')
            replicate = str(data[0])
            age_stage = str(data[1])
            chrom = str(data[2])
            inv = str(data[3])
            amp_name = str(data[4])
            amp_start = int(data[5])
            amp_end = int(data[6])
            library_ID = str(data[7])
            DNA_ID = str(data[8])
            sequencing_name = str(data[9])
            UDI_ID = int(data[10])
            i5_barcode = str(data[11])
            i7_barcode_rc = str(data[12])
            sample_size = int(data[13])
            DNA_concentration = float(data[14])
            paired_read_file_1 = str(data[15])
            paired_read_file_2 = str(data[16])
            # alt_read_file_1 = str(data[17])
            # alt_read_file_2 = str(data[18])

            inv_name = "In("+chrom+")"+inv
            std_name = "In("+chrom+")Std"

            read_metadata[(replicate,chrom,age_stage)] = [replicate,age_stage,chrom,inv,inv_name,std_name,
                amp_name,library_ID,sequencing_name,UDI_ID,i5_barcode,i7_barcode_rc,
                sample_size,DNA_concentration,paired_read_file_1,paired_read_file_2,amp_start,amp_end]
    return read_metadata


# Extracts the called read count and frequency data
#   removes half the called reads, exclusively from the standard set (known std mothers)
#   returns a dictionary with library name as key and the metadata, adjusted proportion and counts as value
def parse_count_data_to_dict(count_file_path):
    count_data = {}
    with open(count_file_path) as count_file:
        lines = count_file.readlines()
        for line in lines[1:]:
            data = line.strip().split(',')
            sequencing_name = str(data[0]).strip()
            line = str(data[1]).strip()
            cohort = str(data[2]).strip()
            chrom = str(data[3]).strip()
            inv = str(data[4]).strip()
            sample_ind_N = int(data[5])
            allele_N = 2*sample_ind_N

            inv_prop = float(data[6])
            std_prop = float(data[7])
            num_inv = int(data[8])
            num_std = int(data[9])
            num_recomb = int(data[10])
            num_mismap = int(data[11])
            num_multimap = int(data[12])

            if cohort != 'Parental ♂':
                total_called = num_inv + num_std
                adjusted_std_count = num_std - total_called//2
                adjusted_total = num_inv + adjusted_std_count
                adjusted_inv_proportion = num_inv/adjusted_total
                adjusted_allele_N = sample_ind_N

                if adjusted_std_count < 1:
                    print("Experimental discrepancy: expected >= half of called reads to be standard orientation:")
                    print(sequencing_name+" has "+str(adjusted_std_count)+" called standard reads out of "+str(adjusted_total)+" after adjustment")

                count_data[(line,chrom,cohort)] = [sequencing_name,line,cohort,chrom,inv,adjusted_allele_N,adjusted_inv_proportion,
                    num_inv,adjusted_std_count]
            else:
                count_data[(line,chrom,cohort)] = [sequencing_name,line,cohort,chrom,inv,allele_N,inv_prop,
                    num_inv,num_std]
    return count_data


# # A helper function to order D. mel chromosome names, 
# #   assumes only X, 2L, 2R, 3L, 3R will be passed
# def chrom_order(chrom):
#     chrom_val = 0
#     if chrom == '2L':
#         chrom_val == 1
#     if chrom == '2R':
#         chrom_val == 2
#     elif chrom == '3L':
#         chrom_val == 3
#     elif chrom == '3R':
#         chrom_val == 4
#     return chrom_val

# # A helper function to order the combination of line and chromosome
# def experimental_order(line,chrom):
#     line_val = int(line)
#     chrom_val = chrom_order(chrom)
#     sort_val = 10*line_val+chrom_val
#     return sort_val

# A helper list with the cohort strings in defined order
cohort_order = ['Parental ♂','Embryo','Early ♀','Late ♀','Early ♂','Late ♂']

"""
# To combine separate sequencing libraries into experimental sets
Returns experimental data lines as:
0    line
1    chrom
2    inv
3    ["Parental ♂", pat_ind_N, pat_freq, pat_read_N]
4    ["Embryo", emb_ind_N, emb_freq, emb_read_N]
5    ["Early ♀", eaf_ind_N, eaf_freq, eaf_read_N]
6    ["Late ♀", laf_ind_N, laf_freq, laf_read_N]
7    ["Early ♂", eam_ind_N, eam_freq, eam_read_N]
8    ["Late ♂", lam_ind_N, lam_freq, lam_read_N]
"""
def reshape_read_count_data(count_data,library_metadata):
    experiment_data = []
    lines = []
    for lib in library_metadata:
        if lib[0] not in lines:
            lines += [lib[0]]
    lines.sort()
    chroms = []
    inv_dict = {}
    for lib in library_metadata:
        if lib[1] not in chroms:
            chroms += [lib[1]]
            inv_dict[lib[1]] = library_metadata[lib][3]
    chroms.sort()
    for line in lines:
        for chrom in chroms:
            experiment = [line,chrom,inv_dict[chrom]]
            contains_cohorts = False
            for cohort in cohort_order:
                allele_count = None
                inv_freq = None
                read_count = None
                if (line,chrom,cohort) in count_data:
                    contains_cohorts = True
                    # Expecting count_data formatted as:
                    #   [sequencing_name,line,cohort,chrom,inv,allele_N,inv_prop,num_inv,num_std]
                    allele_count = count_data[(line,chrom,cohort)][5]
                    inv_freq = count_data[(line,chrom,cohort)][6]
                    inv_read_count = count_data[(line,chrom,cohort)][7]
                    std_read_count = count_data[(line,chrom,cohort)][8]
                    read_total = inv_read_count + std_read_count
                    experiment += [[cohort,allele_count,inv_freq,inv_read_count,std_read_count,read_total]]
                else:
                    print("Missing cohort: "+str((line,chrom,cohort)))
                    experiment += [[cohort,None,None,None,None,None]]
            if contains_cohorts: 
                experiment_data += [experiment]
    return experiment_data

# Takes two *frequencies* and estimates a selection coefficient in an infinite size Wright-Fisher model:
#  Solves p_1 * (1+s)/(p_1*(1+s)+(1-p_1)*1) = p_2
#     (p_1 = p, p_2 = q)
#     p+ps = q(p+ps+1-p)
#     ps-qps = q-p
#     s = (q-p)/(p-qp)
#  Vectorized with numpy
def WF_estimate_sel_coeff(p_1,p_2):
    # print("Frequencies used in W-F est. of sel. coeff:")
    # print(p_1)
    # print(p_2)
    return np.divide(np.subtract(p_2,p_1),np.subtract(p_1,np.multiply(p_1,p_2)))


# For generating the p-value as the accepted bootstrap samples out of the total bootstrap samples,
#    for comparisons with only one cohort comparison involved, though a cohort may be allele-combined
def calc_sing_WF_estimate_sel_coeff(read_counts,read_totals,allele_totals,combination_index,
        mode="off-off",combine=True):
    # print("calculating WF estimate for a single transition")
    read_counts = np.array(read_counts)
    read_totals = np.array(read_totals)
    # print(read_counts)
    # print(read_totals)
    allele_totals = np.array(allele_totals)
    if len(read_counts) > 2:
        if combine:
            last_dim = len(read_counts.shape)-1
            allele_counts_1 = np.multiply(np.divide(read_counts[...,:combination_index],
                read_totals[...,:combination_index]),allele_totals[...,:combination_index])
            freq_1 = np.divide(np.sum(allele_counts_1,axis=last_dim),
                np.sum(allele_totals[...,:combination_index],axis=last_dim))
            allele_counts_2 = np.multiply(np.divide(read_counts[...,combination_index:],
                read_totals[...,combination_index:]),allele_totals[...,combination_index:])
            freq_2 = np.divide(np.sum(allele_counts_2,axis=last_dim),
                np.sum(allele_totals[...,combination_index:],axis=last_dim))
        else:
            raise ValueError("can only calculate a single selection coefficient if "+\
                "given allele totals and combining based on alleles")
    else:
        assert len(read_counts) == 2
        freq_1 = np.divide(read_counts[...,0],read_totals[...,0])
        freq_2 = np.divide(read_counts[...,1],read_totals[...,1])

    # #ALREADY ADJUSTED FREQUENCY
    # if mode == "pat-off":
    #     freq_2 = freq_2*2

    print("Frequencies to test with: "+str((freq_1,freq_2)))

    return WF_estimate_sel_coeff(freq_1,freq_2)

# For generating the p-value as the accepted bootstrap samples out of the total bootstrap samples,
#    for comparisons with only one cohort comparison involved
def calc_mult_WF_estimate_sel_coeff(read_counts,read_totals,allele_totals,index_translations,modes,
        combine=True):

    if verbose: print("Performing direct estimate of selection coefficients from the WF equation "+\
        "for change in frequency: Calculating for both comparisons, using 2 * the offspring frequency "+\
        "for paternal-offspring comparisons")

    read_counts_1  = [read_counts[i] for i in index_translations[0]+index_translations[1]]
    read_totals_1 = [read_totals[i] for i in index_translations[0]+index_translations[1]]
    allele_totals_1 = [allele_totals[i] for i in index_translations[0]+index_translations[1]]
    comb_ind_1 = len(index_translations[0])


    read_counts_2  = [read_counts[i] for i in index_translations[2]+index_translations[3]]
    read_totals_2 = [read_totals[i] for i in index_translations[2]+index_translations[3]]
    allele_totals_2 = [allele_totals[i] for i in index_translations[2]+index_translations[3]]
    comb_ind_2 = len(index_translations[2])

    assert(all([m == "off" for m in modes[1:]])), "Only configured to handle cohort groups containing <= 1 paternal"
    if modes[0] == "pat":
        if index_translations[0][0] == 0:
            assert(len(index_translations[0])==1),"Given a grouping including a paternal but also >1 cohort "
            mode_1 = "pat-off"
        else:
            print("paternal cohort given but not as first comparison - UNEXPECTED")
            mode_1 = "off-off"
        if index_translations[2][0] == 0:
            assert(len(index_translations[2])==1),"Given a grouping including a paternal but also >1 cohort "
            mode_2 = "pat-off"
        else:
            mode_2 = "off-off"
    else:
        assert(modes[0] == "off")
        mode_1 = "off-off"
        mode_2 = "off-off"

    s_1 = calc_sing_WF_estimate_sel_coeff(read_counts_1,read_totals_1,
        allele_totals_1,comb_ind_1,mode_1,combine)
    s_2 = calc_sing_WF_estimate_sel_coeff(read_counts_2,read_totals_2,
            allele_totals_2,comb_ind_2,mode_2,combine)

    sel_coefficients = [s_1,s_2]

    return sel_coefficients

# DO: estimate s, a confidence interval 
#   try to get 1000 acceptances for a resampling
#   DFE is unclear, scale of effects unclear, use uniform -1 to 1 (rerun with a zero )
#       just use s 0-1, direction already inferred (from the associated test)
#       coestimate both s values simultaneously


# Helper function to generate the running log-probability of the 
def __gen_var_probs(allele_total,potential_pat_allele_count,pat_allele_total,var_est_N,mode="off"):
    # print(str((allele_total,potential_pat_allele_count,pat_allele_total,var_est_N,mode)))
    pot_allele_counts = np.arange(allele_total+1)
    pat_freq_est = potential_pat_allele_count/pat_allele_total
    if mode == "pat":
        log_prob_allele_counts = hypergeom.logpmf(pot_allele_counts,pat_allele_total,
            potential_pat_allele_count,allele_total)
    elif mode == "off":
        log_prob_allele_counts = binom.logpmf(pot_allele_counts,allele_total,pat_freq_est)
    else:
        raise ValueError("count emission model mode must be only one of \"off\" or \"pat\"")

    pot_allele_freqs = np.divide(pot_allele_counts,allele_total)
    pot_var_counts = np.arange(var_est_N+1)
    log_prob_var_counts = binom.logpmf(pot_var_counts[np.newaxis,:],var_est_N,pot_allele_freqs[:,np.newaxis])

    running_log_sum = sp.special.logsumexp(np.add(log_prob_allele_counts[:,np.newaxis],log_prob_var_counts),axis=0)

    return(running_log_sum)


def logProb_obs_count_from_pat_est(read_count,read_total,allele_total,
        potential_pat_allele_count,pat_allele_total,var_est_N,mode="off"):
    # pot_allele_counts = np.arange(allele_total+1)
    # pat_freq_est = potential_pat_allele_count/pat_allele_total
    # if mode == "pat":
    #     log_prob_allele_counts = hypergeom.logpmf(pot_allele_counts,pat_allele_total,
    #         potential_pat_allele_count,allele_total)
    # elif mode == "off":
    #     log_prob_allele_counts = binom.logpmf(pot_allele_counts,allele_total,pat_freq_est)
    # else:
    #     raise ValueError("count emission model mode must be only one of \"off\" or \"pat\"")

    # pot_allele_freqs = np.divide(pot_allele_counts,allele_total)
    # pot_var_counts = np.arange(var_est_N+1)
    # log_prob_var_counts = binom.logpmf(pot_var_counts[np.newaxis,:],var_est_N,pot_allele_freqs[:,np.newaxis])

    # running_log_sum = sp.special.logsumexp(np.add(log_prob_allele_counts[:,np.newaxis],log_prob_var_counts),axis=0)

    running_log_sum = __gen_var_probs(allele_total,potential_pat_allele_count,pat_allele_total,var_est_N,mode)

    pot_var_freqs = np.divide(np.arange(var_est_N+1),var_est_N)

    log_prob_obs_read_count = binom.logpmf(read_count,read_total,pot_var_freqs)

    total_log_prob = sp.special.logsumexp(np.add(running_log_sum,log_prob_obs_read_count))
    return total_log_prob


def logProb_obs_count_combination_from_pat_est(read_counts,read_totals,allele_totals,
        potential_pat_allele_count,pat_allele_total,var_est_N,mode="off-off"):
    log_probs = np.zeros(len(read_counts))
    # Handle paternal cohorts, check for correct mode specification
    if mode == "off-off":
        log_probs[0] = logProb_obs_count_from_pat_est(read_counts[0],read_totals[0],
            allele_totals[0],potential_pat_allele_count,pat_allele_total,var_est_N,mode="off")
    elif mode == "pat-off":
        log_probs[0] = logProb_obs_count_from_pat_est(read_counts[0],read_totals[0],
            allele_totals[0],potential_pat_allele_count,pat_allele_total,var_est_N,mode="pat")
    else:
        raise ValueError("ML model mode must be only one of \"off-off\" or \"pat-off\"")
    # Calculate the log probabilities of the rest of the count observations
    for i in np.arange(1,len(read_counts)):
        log_probs[i] = logProb_obs_count_from_pat_est(read_counts[i],read_totals[i],allele_totals[i],
            potential_pat_allele_count,pat_allele_total,var_est_N,mode="off")
    total_log_prob = np.sum(log_probs)
    # print("log-prob estimate: "+str(total_log_prob))
    return total_log_prob




# Return the log-probability distribution of observing each count (the index in the distribution)
#   from the given ML paternal census allele count
def logProb_count_dist_from_pat_est(read_count,read_total,allele_total,
        potential_pat_allele_count,pat_allele_total,var_est_N,mode="off"):
    running_log_sum = __gen_var_probs(allele_total,potential_pat_allele_count,
        pat_allele_total,var_est_N,mode)
    pot_var_freqs = np.divide(np.arange(var_est_N+1),var_est_N)
    pot_read_counts = np.arange(read_total+1)
    log_prob_read_counts = binom.logpmf(pot_read_counts[np.newaxis,:],read_total,pot_var_freqs[:,np.newaxis])
    running_log_sum = sp.special.logsumexp(np.add(running_log_sum[:,np.newaxis],log_prob_read_counts),axis=0)
    return running_log_sum

# Return a list of log-probability distributions for each count category,
#  of the probability of observing that count (the index in the distribution) from the given ML paternal census allele count
def logProb_count_dists_from_pat_est(read_counts,read_totals,allele_totals,
        potential_pat_allele_count,pat_allele_total,var_est_N,mode="off-off"):
    ML_est_count_probs = []
    # Handle paternal cohorts, check for correct mode specification
    if mode == "off-off":
        ML_est_count_probs += [logProb_count_dist_from_pat_est(read_counts[0],read_totals[0],
            allele_totals[0],potential_pat_allele_count,pat_allele_total,var_est_N,mode="off")]
    elif mode == "pat-off":
        ML_est_count_probs += [logProb_count_dist_from_pat_est(read_counts[0],read_totals[0],
            allele_totals[0],potential_pat_allele_count,pat_allele_total,var_est_N,mode="pat")]
    else:
        raise ValueError("ML model mode must be only one of \"off-off\" or \"pat-off\"")
    # Calculate the log probabilities of the rest of the count observation distributions
    for i in np.arange(1,len(read_counts)):
         ML_est_count_probs += [logProb_count_dist_from_pat_est(read_counts[i],read_totals[i],allele_totals[i],
            potential_pat_allele_count,pat_allele_total,var_est_N,mode="off")]
    return ML_est_count_probs

# Calculates the expectation of a 1d array of values with a matching array of log probabilities
#   :vectorized
def expectation(log_probs,vals):
    E = np.sum(np.multiply(np.exp(log_probs),vals))
    return E

# Generate only the maximum likelihood paternal allele count
def gen_ML_pat_count(read_counts,read_totals,allele_totals,
        paternal_pop_allele_total,var_est_N,mode="off-off"):
    # Generate the likelihood of each possible paternal allele count,
    #   and estimate the maximum likelihood 
    if verbose: print("Estimating maximum likelihood paternal population allele count")
    # pat_est_log_probs = np.zeros(paternal_pop_allele_total-1)
    ML_pat_allele_count = None
    max_log_prob = -float('inf')
    for potential_pat_allele_count in np.arange(1,paternal_pop_allele_total):
        curr_log_prob = logProb_obs_count_combination_from_pat_est(read_counts,read_totals,allele_totals,
            potential_pat_allele_count,paternal_pop_allele_total,var_est_N,mode)
        # pat_est_log_probs[potential_pat_allele_count-1] = curr_log_prob
        if curr_log_prob > max_log_prob:
            max_log_prob = curr_log_prob
            ML_pat_allele_count = potential_pat_allele_count
    if verbose:
        print("ML_pat_allele_count:\t\t\t\t"+str(ML_pat_allele_count))
        print("Paternal allele census size:\t\t\t"+str(paternal_pop_allele_total))
        print("Maximum likelihood paternal allele frequency:\t"+\
            str(ML_pat_allele_count/paternal_pop_allele_total))
    return ML_pat_allele_count

def gen_expected_counts_separately(read_counts,read_totals,allele_totals,
        ML_pat_allele_count,paternal_pop_allele_total,var_est_N,mode="off-off"):
    ML_est_expected_counts = []
    # Handle paternal cohorts, check for correct mode specification
    if mode == "off-off":
        ML_est_expected_counts += [expectation(logProb_count_dist_from_pat_est(read_counts[0],read_totals[0],
            allele_totals[0],ML_pat_allele_count,paternal_pop_allele_total,var_est_N,mode="off"),np.arange(read_totals[0]+1))]
    elif mode == "pat-off":
        ML_est_expected_counts += [expectation(logProb_count_dist_from_pat_est(read_counts[0],read_totals[0],
            allele_totals[0],ML_pat_allele_count,paternal_pop_allele_total,var_est_N,mode="pat"),np.arange(read_totals[0]+1))]
    else:
        raise ValueError("ML model mode must be only one of \"off-off\" or \"pat-off\"")
    # Calculate the log probabilities of the rest of the count observation distributions
    for i in np.arange(1,len(read_counts)):
         ML_est_expected_counts += [expectation(logProb_count_dist_from_pat_est(read_counts[i],read_totals[i],allele_totals[i],
            ML_pat_allele_count,paternal_pop_allele_total,var_est_N,mode="off"),np.arange(read_totals[i]+1))]
    return ML_est_expected_counts


# Estimate the maximum likelihood paternal frequency
#   the first observed count/total is always assumed to be the paternal in "pat-off" mode
def gen_count_prob_dists(read_counts,read_totals,allele_totals,
        ML_pat_allele_count,paternal_pop_allele_total,var_est_N,mode="off-off"):

    # check for correct mode specification
    if mode != "off-off" and mode != "pat-off":
        raise ValueError("ML model mode must be only one of \"off-off\" or \"pat-off\"")

    ML_pat_allele_freq = ML_pat_allele_count/paternal_pop_allele_total

    if verbose:
        print("Generating the log-likelihood distribution of possible "+\
            "observed read counts")
    ML_est_count_probs = logProb_count_dists_from_pat_est(read_counts,read_totals,allele_totals,
        ML_pat_allele_count,paternal_pop_allele_total,var_est_N,mode)

    if verbose:
        print("Generating the expected read counts from the LLH distribution")
    ML_est_expected_counts = [expectation(ML_est_count_probs[i],np.arange(read_totals[i]+1)) for i in np.arange(len(read_totals))]

    return(ML_est_count_probs,ML_est_expected_counts,ML_pat_allele_freq)


# Estimate the maximum likelihood paternal frequency
#   the first observed count/total is always assumed to be the paternal in "pat-off" mode
def gen_ML_pat_count_and_count_prob_dists(read_counts,read_totals,allele_totals,
        paternal_pop_allele_total,var_est_N,mode="off-off"):

    # check for correct mode specification
    if mode != "off-off" and mode != "pat-off":
        raise ValueError("ML model mode must be only one of \"off-off\" or \"pat-off\"")

    ML_pat_allele_count = gen_ML_pat_count(read_counts,read_totals,allele_totals,
        paternal_pop_allele_total,var_est_N,mode)

    ML_pat_allele_freq = ML_pat_allele_count/paternal_pop_allele_total

    if verbose:
        print("Generating the log-likelihood distribution of possible "+\
            "observed read counts")
    ML_est_count_probs = logProb_count_dists_from_pat_est(read_counts,read_totals,allele_totals,
        ML_pat_allele_count,paternal_pop_allele_total,var_est_N,mode)

    if verbose:
        print("Generating the expected read counts from the LLH distribution")
    ML_est_expected_counts = [expectation(ML_est_count_probs[i],np.arange(read_totals[i]+1)) for i in np.arange(len(read_totals))]

    # potential_paternal_allele_counts = np.arange(1,paternal_pop_allele_total)
    # potential_paternal_allele_freqs = np.divide(potential_paternal_allele_counts,paternal_pop_allele_total)

    # # Store the log-probability matrices for each observation set
    # #   (matrices indexed in possible paternal allele counts by possible observed read counts) 
    # log_prob_matrices = []
    # # Calculate the log probabilities of the first set
    # pot_allele_counts_1 = np.arange(allele_totals[0]+1)
    # # running_log_sum = None
    # if mode == "pat-off":
    #     log_prob_allele_counts_1 = hypergeom.logpmf(pot_allele_counts_1[np.newaxis,:],
    #         paternal_pop_allele_total,potential_paternal_allele_counts[:,np.newaxis],allele_totals[0])
    # else:
    #     log_prob_allele_counts_1 = binom.logpmf(pot_allele_counts_1[np.newaxis,:],
    #         allele_totals[0],potential_paternal_allele_freqs[:,np.newaxis])
    # # running_log_sum = log_prob_allele_counts_1


    # pot_allele_freqs_1 = np.divide(pot_allele_counts_1,allele_totals[0])
    # pot_var_counts = np.arange(var_est_N+1)
    # log_prob_var_counts_1 = binom.logpmf(pot_var_counts[np.newaxis,:],var_est_N,pot_allele_freqs_1[:,np.newaxis])

    # running_log_sum_1 = sp.special.logsumexp(np.add(log_prob_allele_counts_1[:,:,np.newaxis],log_prob_var_counts_1[np.newaxis,:,:]),axis=1)

    # pot_var_freqs = np.divide(pot_var_counts,var_est_N)
    # pot_read_counts_1 = np.arange(read_totals[0]+1)
    # log_prob_read_counts_1 = binom.logpmf(pot_read_counts_1[np.newaxis,:],read_totals[0],pot_var_freqs[:,np.newaxis])

    # running_log_sum_1 = sp.special.logsumexp(np.add(running_log_sum_1[:,:,np.newaxis],log_prob_read_counts_1[np.newaxis,:,:]),axis=1)
    # log_prob_matrices += [running_log_sum_1]

    # # calculate the log probabilities of the remaining sets
    # for i in np.arange(1,len(read_counts)):
    #     pot_allele_counts = np.arange(allele_totals[i]+1)
    #     log_prob_allele_counts = binom.logpmf(pot_allele_counts[np.newaxis,:],
    #         allele_totals[i],potential_paternal_allele_freqs[:,np.newaxis])

    #     pot_allele_freqs = np.divide(pot_allele_counts,allele_totals[i])
    #     # pot_var_counts = np.arange(1,var_est_N)
    #     log_prob_var_counts = binom.logpmf(pot_var_counts[np.newaxis,:],var_est_N,pot_allele_freqs[:,np.newaxis])

    #     running_log_sum = sp.special.logsumexp(np.add(log_prob_allele_counts[:,:,np.newaxis],log_prob_var_counts[np.newaxis,:,:]),axis=1)

    #     # pot_var_freqs = np.divide(pot_var_counts,var_est_N)
    #     pot_read_counts = np.arange(read_totals[i]+1)
    #     log_prob_read_counts = binom.logpmf(pot_read_counts[np.newaxis,:],read_totals[i],pot_var_freqs[:,np.newaxis])

    #     running_log_sum = sp.special.logsumexp(np.add(running_log_sum[:,:,np.newaxis],log_prob_read_counts[np.newaxis,:,:]),axis=1)

    #     log_prob_matrices += [running_log_sum]
    

    # # Identify the maximum likelihood paternal population allele count
    # ML_pat_allele_count = None
    # i_ML = None
    # min_nLL = float('inf')
    # assert (log_prob_matrices[0].shape[0] == paternal_pop_allele_total-1)
    # # print(len(log_prob_matrices))
    # # for arr in log_prob_matrices:
    # #     print(arr.shape)
    # # print(read_counts)
    # # print(read_totals)
    # for i in np.arange(log_prob_matrices[0].shape[0]):
    #     curr_neg_log_likelihood = 0
    #     for l in np.arange(len(read_counts)):
    #         curr_neg_log_likelihood += -log_prob_matrices[l][i,read_counts[l]]
    #     if curr_neg_log_likelihood < min_nLL:
    #         ML_pat_allele_count = i+1
    #         i_ML = i
    #         min_nLL = curr_neg_log_likelihood
    # ML_pat_allele_freq = ML_pat_allele_count/paternal_pop_allele_total


    # # Use the ML paternal population allele count to select the probability distribution of count outcomes
    # ML_est_obs_probs = []
    # for i in np.arange(len(read_counts)):
    #     ML_est_obs_probs += [log_prob_matrices[i][i_ML,]]

    return(ML_est_count_probs,ML_est_expected_counts,ML_pat_allele_count,ML_pat_allele_freq)




# A function running the basic chi statistic calculation on an observation and expectation
def chi(observed,expected):
    return ((observed-expected)**2)/expected

# calculates the chi statistic for a specific observation
def calc_chi_statistic_unvectorized(read_counts,read_totals,expected_counts,allele_totals=None,count_mode="read",
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

    # Generate the list of counts to compare, either reads or alleles,
    #  and combining the second on allele counts if specified
    expected_freqs = [expected_counts[i]/read_totals[i] for i in np.arange(len(read_counts))]
    pos_counts = []
    neg_counts = []
    exp_pos_counts = []
    exp_neg_counts = []
    if count_mode == "read":
        for i in np.arange(len(read_counts)):
            pos_counts += [read_counts[i]]
            neg_counts += [read_totals[i]-read_counts[i]]
            expected_pos = expected_freqs[i]*read_totals[i]
            exp_pos_counts += [expected_pos]
            exp_neg_counts += [read_totals[i]-expected_pos]
    elif count_mode == "allele":
        pos_allele_counts = []
        neg_allele_counts = []
        exp_pos_allele_counts = []
        exp_neg_allele_counts = []
        for i in np.arange(len(read_counts)):
            pos_proportion = read_counts[i]/read_totals[i]
            pos_allele_counts += [pos_proportion*allele_totals[i]]
            neg_allele_counts += [(1-pos_proportion)*allele_totals[i]]
            expected_pos = expected_freqs[i]*allele_totals[i]
            exp_pos_allele_counts += [expected_pos]
            exp_neg_allele_counts += [allele_totals[i]-expected_pos]
        if combine:
            # pos_count_tail = 0
            # neg_count_tail = 0
            # exp_pos_tail = 0
            # exp_neg_tail = 0
            # for i in np.arange(1,len(read_counts)):
            #     pos_count_tail += pos_allele_counts[i]
            #     neg_count_tail += neg_allele_counts[i]
            #     exp_pos_tail += exp_pos_allele_counts[i]
            #     exp_neg_tail += exp_neg_allele_counts[i]
            # pos_counts = [pos_allele_counts[0],pos_count_tail]
            # neg_counts = [neg_allele_counts[0],neg_count_tail]
            # exp_pos_counts = [exp_pos_allele_counts[0],exp_pos_tail]
            # exp_neg_counts = [exp_neg_allele_counts[0],exp_neg_tail]
            pos_counts = [sum(pos_allele_counts[:combination_index]),sum(pos_allele_counts[combination_index:])]
            neg_counts = [sum(neg_allele_counts[:combination_index]),sum(neg_allele_counts[combination_index:])]
            exp_pos_counts = [sum(exp_pos_allele_counts[:combination_index]),sum(exp_pos_allele_counts[combination_index:])]
            exp_neg_counts = [sum(exp_neg_allele_counts[:combination_index]),sum(exp_neg_allele_counts[combination_index:])]
        else:
            pos_counts = pos_allele_counts
            neg_counts = neg_allele_counts
            exp_pos_counts = exp_pos_allele_counts
            exp_neg_counts = exp_neg_allele_counts
    else:
        raise ValueError("Unexpected count_mode value: "+str(count_mode))

    # Calculate the sum of the chi statistic for all relevant comparisons
    chi_sum = 0
    for i in np.arange(len(pos_counts)):
        chi_sum += chi(pos_counts[i],exp_pos_counts[i])
        chi_sum += chi(neg_counts[i],exp_neg_counts[i])

    return chi_sum

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


# Generate the chi statistic for the entire landscape of read count combinations
def calc_chi_stat_dist(read_totals,expected_counts,allele_totals=None,count_mode="read",
        combination_index=None,combine=False):
    dimensions = [t+1 for t in read_totals]
    chi_stats = np.zeros(dimensions) # ndtype=
    for hypothetical_read_counts in np.ndindex(*dimensions):
        chi_statistic = calc_chi_statistic(hypothetical_read_counts,read_totals,expected_counts,
            allele_totals,count_mode,combination_index,combine)
        chi_stats[hypothetical_read_counts] = chi_statistic
    return chi_stats

# Helper function to get the log probability of a specific index of read counts from the ML model
def __get_logProb_counts(index,ML_est_obs_logProbs):
    logProb = 0 
    for dimension,dim_i in np.ndenumerate(index):
        logProb += ML_est_obs_logProbs[dimension[0]][dim_i]
    return logProb

# Calculate the sum of the probabilities of counts under the ML model with a statistic greater than the observed count
# Assumes that the probability data is given as a list of ndarrays,
#   where each array is as long as the total number of observed reads,
#   and the log probability in each position is the probability of that count under the ML model
#   So, the dimension of the stats array should equal the length of each array in ML_est_obs_logProbs, in the same order
def calc_log_pval_from_stat(ML_est_obs_logProbs,stats,observed_stat):
    sig_logProb_vals = []
    for index, val in np.ndenumerate(stats):
        if val >= observed_stat:
            logProb = __get_logProb_counts(index,ML_est_obs_logProbs)
            sig_logProb_vals += [logProb]
    total_log_pval = sp.special.logsumexp(sig_logProb_vals)
    return total_log_pval

# # A helper function to calculate the cumulant sum of exponentiated log-probabilities
# #   with the sum of input indexes from 0 to i stored in output index i
# def __cumulant_partial_sum_exp(log_probs):
#     cumulants = np.zeros(len(log_probs))
#     for i in np.arange(len(log_probs)):
#         running_sum_exp = 
#     return cumulants



def calc_chi_stat_pval(read_counts,read_totals,expected_counts,ML_est_obs_logProbs,
        allele_totals=None,count_mode="read",combination_index=None,combine=False):
    # Check the inputs (may be redundant, also done when calculating the chi statistic for the observed counts)
    if allele_totals is None and count_mode == "allele":
        raise ValueError("can only calculate an allele-count chi statistic when allele_totals are provided")
    if combination_index is None and combine:
        raise ValueError("can only combine counts when the combination_index at which to combine into two groups is provided")
    if combine and count_mode == "read":
        raise ValueError("only expected to combine the 1...n other counts if considering allele_counts")
    if count_mode != "read" and count_mode != "allele":
        raise ValueError("chi statistic calculation count_mode must only be one of \"read\" or \"allele\"")

    # expected_freqs = [expected_counts[i]/read_totals[i] for i in np.arange(len(read_counts))]

    if verbose: print("Calculating the chi statistics and comparing to the observed for p-value generation")

    observed_stat = calc_chi_statistic(read_counts,read_totals,expected_counts,allele_totals,
        count_mode,combination_index,combine)
    logProb_obs_counts = __get_logProb_counts(read_counts,ML_est_obs_logProbs)
    if verbose:
        print("Observed Counts:\t"+str(read_counts)+"\nExpected Counts:\t"+str(expected_counts)+\
            "\nTotal Counts:\t\t"+str(read_totals)+"\nLog Prob of Counts:\t"+str(logProb_obs_counts)+\
            "\nChi Square:\t\t"+str(observed_stat))

    scaled_diff_exp_obs = [(expected_counts[i]-read_counts[i])/(read_totals[i]+1) for i in np.arange(len(read_counts))]
    approx_proportion_area_within_obs_stat = np.pi*np.sum(np.power(scaled_diff_exp_obs,2)) # pi*r^2 within a unit cube
    calc_exterior = approx_proportion_area_within_obs_stat > 0.5
    if verbose:
        print("Approximate proportion of the test area closer than observed stat:\n"+\
            str(approx_proportion_area_within_obs_stat))
        if calc_exterior:
            print("Calculating by summing the exterior space - counts with statistic >= observed")
        else:
            print("Calculating by summing the interior space - counts with statistic < observed")

    # Subtract this (approximately maximum) term from each log before summing the exponents, for numerical stability
    logsumexp_factor = logProb_obs_counts
    if calc_exterior:
        logsumexp_factor = __get_logProb_counts(np.array(expected_counts).round().astype('int32'),ML_est_obs_logProbs) # If summing the inside

    if combine and len(read_counts)>2:
        if verbose:
            print("(Cohort-combined statistics) Calculating the p-value and the chi-statistic values "+\
                "independently for each read count combination\n")
        # Simultaneously calculate the p-value and the chi-statistic values
        running_sum_probDiffs = 0
        dimensions = [t+1 for t in read_totals]
        hypothetical_read_counts = np.ndindex(*dimensions[:-1]) # Hold out the last dimension to allow breaking the loop when no longer significant
        if calc_exterior:
            for hyp_counts in hypothetical_read_counts:
                held_count = 0
                max_held_count = read_totals[-1]
                curr_counts = list(hyp_counts)+[held_count]
                potential_chi_val = calc_chi_statistic(curr_counts,read_totals,
                    expected_counts,allele_totals,count_mode,combination_index,combine)
                while held_count < max_held_count and potential_chi_val >= observed_stat:
                    logProb = __get_logProb_counts(curr_counts,ML_est_obs_logProbs)
                    running_sum_probDiffs += np.exp(logProb-logsumexp_factor)
                    held_count += 1
                    curr_counts = list(hyp_counts)+[held_count]
                    potential_chi_val = calc_chi_statistic(curr_counts,read_totals,
                        expected_counts,allele_totals,count_mode,combination_index,combine)
                if held_count == max_held_count:
                    if potential_chi_val >= observed_stat:
                        logProb = __get_logProb_counts(curr_counts,ML_est_obs_logProbs)
                        running_sum_probDiffs += np.exp(logProb-logsumexp_factor)
                else:
                    assert held_count < max_held_count
                    min_held_count = held_count
                    held_count = max_held_count
                    curr_counts = list(hyp_counts)+[held_count]
                    potential_chi_val = calc_chi_statistic(curr_counts,read_totals,
                        expected_counts,allele_totals,count_mode,combination_index,combine)
                    while held_count > min_held_count and potential_chi_val >= observed_stat:
                        logProb = __get_logProb_counts(curr_counts,ML_est_obs_logProbs)
                        running_sum_probDiffs += np.exp(logProb-logsumexp_factor)
                        held_count -= 1
                        curr_counts = list(hyp_counts)+[held_count]
                        potential_chi_val = calc_chi_statistic(curr_counts,read_totals,
                            expected_counts,allele_totals,count_mode,combination_index,combine)
                    if held_count == min_held_count:
                        if potential_chi_val >= observed_stat:
                            logProb = __get_logProb_counts(curr_counts,ML_est_obs_logProbs)
                            running_sum_probDiffs += np.exp(logProb-logsumexp_factor)
        else:
            for hyp_counts in hypothetical_read_counts:
                held_count_upper = int(expected_counts[-1])+1 # Assumes the expectation is never at the extreme value
                max_held_count = read_totals[-1]
                curr_counts = list(hyp_counts)+[held_count_upper]
                potential_chi_val = calc_chi_statistic(curr_counts,read_totals,
                    expected_counts,allele_totals,count_mode,combination_index,combine)
                while held_count_upper < max_held_count and potential_chi_val < observed_stat:
                    logProb = __get_logProb_counts(curr_counts,ML_est_obs_logProbs)
                    running_sum_probDiffs += np.exp(logProb-logsumexp_factor)
                    held_count_upper += 1
                    curr_counts = list(hyp_counts)+[held_count_upper]
                    potential_chi_val = calc_chi_statistic(curr_counts,read_totals,
                        expected_counts,allele_totals,count_mode,combination_index,combine)
                if held_count_upper == max_held_count:
                    if potential_chi_val < observed_stat:
                        logProb = __get_logProb_counts(curr_counts,ML_est_obs_logProbs)
                        running_sum_probDiffs += np.exp(logProb-logsumexp_factor)
                held_count_lower = int(expected_counts[-1]) # Assumes the expectation is never at the extreme value
                min_held_count = 0
                curr_counts = list(hyp_counts)+[held_count_lower]
                potential_chi_val = calc_chi_statistic(curr_counts,read_totals,
                    expected_counts,allele_totals,count_mode,combination_index,combine)
                while held_count_lower > min_held_count and potential_chi_val < observed_stat:
                    logProb = __get_logProb_counts(curr_counts,ML_est_obs_logProbs)
                    running_sum_probDiffs += np.exp(logProb-logsumexp_factor)
                    held_count_lower -= 1
                    curr_counts = list(hyp_counts)+[held_count_lower]
                    potential_chi_val = calc_chi_statistic(curr_counts,read_totals,
                        expected_counts,allele_totals,count_mode,combination_index,combine)
                if held_count_lower == min_held_count:
                    if potential_chi_val < observed_stat:
                        logProb = __get_logProb_counts(curr_counts,ML_est_obs_logProbs)
                        running_sum_probDiffs += np.exp(logProb-logsumexp_factor)
    else:
        if verbose:
            print("(Non-combined statistics) Pre-calculating partial chi squared values "+\
                "for each count category, to memoize part of the chi-statistic calculation "+\
                "when finding the p-value and chi-stat for each read count combination\n")
        # Pre-calculate the partial chi squared values for each count category
        chi_stat_partials = []
        for i in np.arange(len(read_counts)):
            potential_pos_read_counts = np.arange(read_totals[i]+1)
            potential_neg_read_counts = np.arange(read_totals[i],-1,-1)
            if count_mode == "read":
                expected_pos = expected_counts[i]
                expected_neg = read_totals[i] - expected_counts[i]
                potential_pos = potential_pos_read_counts
                potential_neg = potential_neg_read_counts
            elif count_mode == "allele":
                potential_pos_allele_counts = potential_pos_read_counts/read_totals[i]*allele_totals[i]
                potential_neg_allele_counts = potential_neg_read_counts/read_totals[i]*allele_totals[i]
                expected_pos = expected_counts[i]/read_totals[i]*allele_totals[i]
                expected_neg = allele_totals[i] - expected_pos
                potential_pos = potential_pos_allele_counts
                potential_neg = potential_neg_allele_counts
            else:
                raise ValueError("chi statistic calculation count_mode must only be one of \"read\" or \"allele\"")
            chi_stat_pos_partial = chi(potential_pos,expected_pos)
            chi_stat_neg_partial = chi(potential_neg,expected_neg)
            chi_stat_partial = chi_stat_pos_partial + chi_stat_neg_partial
            chi_stat_partials += [chi_stat_partial]
        # Calculate the p-value from these partial chi-statistic values
        running_sum_probDiffs = 0
        dimensions = [t+1 for t in read_totals]
        hypothetical_read_counts = np.ndindex(*dimensions[:-1]) # Hold out the last dimension to allow breaking the loop when no longer significant
        if calc_exterior:
            for hyp_counts in hypothetical_read_counts:
                held_count = 0
                max_held_count = read_totals[-1]
                curr_counts = list(hyp_counts)+[held_count]
                potential_chi_val = np.sum([chi_stat_partials[i][curr_counts[i]] for i in np.arange(len(read_counts))])
                while held_count < max_held_count and potential_chi_val >= observed_stat:
                    logProb = __get_logProb_counts(curr_counts,ML_est_obs_logProbs)
                    running_sum_probDiffs += np.exp(logProb-logsumexp_factor)
                    held_count += 1
                    curr_counts = list(hyp_counts)+[held_count]
                    potential_chi_val = np.sum([chi_stat_partials[i][curr_counts[i]] for i in np.arange(len(read_counts))])
                if held_count == max_held_count:
                    if potential_chi_val >= observed_stat:
                        logProb = __get_logProb_counts(curr_counts,ML_est_obs_logProbs)
                        running_sum_probDiffs += np.exp(logProb-logsumexp_factor)
                else:
                    assert held_count < max_held_count
                    min_held_count = held_count
                    held_count = max_held_count
                    curr_counts = list(hyp_counts)+[held_count]
                    potential_chi_val = np.sum([chi_stat_partials[i][curr_counts[i]] for i in np.arange(len(read_counts))])
                    while held_count > min_held_count and potential_chi_val >= observed_stat:
                        logProb = __get_logProb_counts(curr_counts,ML_est_obs_logProbs)
                        running_sum_probDiffs += np.exp(logProb-logsumexp_factor)
                        held_count -= 1
                        curr_counts = list(hyp_counts)+[held_count]
                        potential_chi_val = np.sum([chi_stat_partials[i][curr_counts[i]] for i in np.arange(len(read_counts))])
                    if held_count == min_held_count:
                        if potential_chi_val >= observed_stat:
                            logProb = __get_logProb_counts(curr_counts,ML_est_obs_logProbs)
                            running_sum_probDiffs += np.exp(logProb-logsumexp_factor)
        else:
            for hyp_counts in hypothetical_read_counts:
                held_count_upper = int(expected_counts[-1])+1 # Assumes the expectation is never at the extreme value
                max_held_count = read_totals[-1]
                curr_counts = list(hyp_counts)+[held_count_upper]
                potential_chi_val = np.sum([chi_stat_partials[i][curr_counts[i]] for i in np.arange(len(read_counts))])
                while held_count_upper < max_held_count and potential_chi_val < observed_stat:
                    logProb = __get_logProb_counts(curr_counts,ML_est_obs_logProbs)
                    running_sum_probDiffs += np.exp(logProb-logsumexp_factor)
                    held_count_upper += 1
                    curr_counts = list(hyp_counts)+[held_count_upper]
                    potential_chi_val = np.sum([chi_stat_partials[i][curr_counts[i]] for i in np.arange(len(read_counts))])
                if held_count_upper == max_held_count:
                    if potential_chi_val < observed_stat:
                        logProb = __get_logProb_counts(curr_counts,ML_est_obs_logProbs)
                        running_sum_probDiffs += np.exp(logProb-logsumexp_factor)
                held_count_lower = int(expected_counts[-1]) # Assumes the expectation is never at the extreme value
                min_held_count = 0
                curr_counts = list(hyp_counts)+[held_count_lower]
                potential_chi_val = np.sum([chi_stat_partials[i][curr_counts[i]] for i in np.arange(len(read_counts))])
                while held_count_lower > min_held_count and potential_chi_val < observed_stat:
                    logProb = __get_logProb_counts(curr_counts,ML_est_obs_logProbs)
                    running_sum_probDiffs += np.exp(logProb-logsumexp_factor)
                    held_count_lower -= 1
                    curr_counts = list(hyp_counts)+[held_count_lower]
                    potential_chi_val = np.sum([chi_stat_partials[i][curr_counts[i]] for i in np.arange(len(read_counts))])
                if held_count_lower == min_held_count:
                    if potential_chi_val < observed_stat:
                        logProb = __get_logProb_counts(curr_counts,ML_est_obs_logProbs)
                        running_sum_probDiffs += np.exp(logProb-logsumexp_factor)

    total_log_pval = logsumexp_factor + np.log(running_sum_probDiffs)
    if verbose:
        if calc_exterior:
            print("Exterior log-p-value:\t\t\t"+str(total_log_pval))
        else:
            print("Interior log-p-value:\t\t\t"+str(total_log_pval))
    if not calc_exterior:
        total_log_pval = np.log(1-np.exp(total_log_pval))
    if verbose: print("Final log-p-value:\t\t\t"+str(total_log_pval))

    return (observed_stat,total_log_pval)



# Helper function to generate the running log-probability of a sampled experimental simulation
#   up to the sampling that accounts for the sample-size independent variation introduced in library prep
def __gen_var_samples(allele_total,ML_pat_allele_count,pat_allele_total,var_est_N,sample_size,mode="off"):
    pat_freq_est = ML_pat_allele_count/pat_allele_total
    if mode == "pat":
        sampled_allele_counts = hypergeom.rvs(pat_allele_total,ML_pat_allele_count,
            allele_total,size=sample_size)
    elif mode == "off":
        sampled_allele_counts = binom.rvs(allele_total,pat_freq_est,size=sample_size)
    else:
        raise ValueError("count emission model mode must be only one of \"off\" or \"pat\"")

    sampled_allele_freqs = np.divide(sampled_allele_counts,allele_total)
    sampled_var_counts = binom.rvs(var_est_N,sampled_allele_freqs)

    return(sampled_var_counts)


# A function to sample the (modeled) neutral distribution of read counts 
def sample_cohort_with_pat_freq(allele_total,ML_pat_allele_count,pat_allele_total,
        var_est_N,read_total,sample_size,mode="off"):

    # if verbose:
    #     print("Running neutral sampling of a cohort from a given paternal freq")
    #     print((allele_total,ML_pat_allele_count,pat_allele_total,
    #         var_est_N,read_total,sample_size,mode))

    sampled_var_counts = __gen_var_samples(allele_total,ML_pat_allele_count,pat_allele_total,
        var_est_N,sample_size,mode)

    sampled_var_freqs = np.divide(sampled_var_counts,var_est_N)

    sampled_read_counts = binom.rvs(read_total,sampled_var_freqs)

    return sampled_read_counts


# Helper function to generate the running log-probability of a sampled experimentsl simulation
#   up to the sampling that accounts for the sample-size independent variation introduced in library prep
#   while incorporating SELECTION modeled as a 
def __gen_var_samples_with_selection(allele_total,ML_pat_allele_count,pat_allele_total,var_est_N,
        sample_size,sel_coefficients,mode="off"):
    if mode == "pat":
        sampled_allele_counts = hypergeom.rvs(pat_allele_total,ML_pat_allele_count,
            allele_total,size=sample_size)
    elif mode == "off":
        # pat_freq_est = ML_pat_allele_count/pat_allele_total
        # sel_mult = 1
        # if isinstance(sel_coefficients,(int,float)):
        #     sel_mult = 1+sel_coefficients
        # else:
        #     for s in sel_coefficients:
        #         sel_mult *= (1+s)
        # # prob_post_sel = pat_freq_est*(1+sel_coeff)/(pat_freq_est*(1+sel_coeff)+(1-pat_freq_est))
        # prob_post_sel = pat_freq_est*sel_mult/(pat_freq_est*sel_mult+(1-pat_freq_est))
        # sampled_allele_counts = binom.rvs(allele_total,prob_post_sel,size=sample_size)

        pat_freq_est = ML_pat_allele_count/pat_allele_total
        sel_mults = 1+sel_coefficients
        freqs_post_sel = pat_freq_est*sel_mults/(pat_freq_est*sel_mults+(1-pat_freq_est))
        sampled_allele_counts = binom.rvs(allele_total,freqs_post_sel)
    else:
        raise ValueError("count emission model mode must be only one of \"off\" or \"pat\"")

    sampled_allele_freqs = np.divide(sampled_allele_counts,allele_total)
    sampled_var_counts = binom.rvs(var_est_N,sampled_allele_freqs)

    return(sampled_var_counts)


# A function to sample the (modeled) neutral distribution of read counts 
def sample_cohort_with_pat_freq_and_selection(allele_total,ML_pat_allele_count,pat_allele_total,
        var_est_N,read_total,sample_size,sel_coefficients,mode="off"):
    sampled_var_counts = __gen_var_samples_with_selection(allele_total,ML_pat_allele_count,pat_allele_total,
        var_est_N,sample_size,sel_coefficients,mode)

    sampled_var_freqs = np.divide(sampled_var_counts,var_est_N)
    sampled_read_counts = binom.rvs(read_total,sampled_var_freqs)

    return sampled_read_counts

# Helper function to run bootstrap sampling for single-comparison stats
def __bootstrap_sample(read_counts,read_totals,ML_pat_allele_count,pat_allele_total,
        var_est_N,modes,expected_counts,allele_totals,
        tail,expected_direction,
        count_mode,combination_index,combine,
        observed_stat,observed_direction,
        bs_sample_size):
    # np.random.choice treats an integer parameter n as np.arange(n)
    sampled_counts = np.array([sample_cohort_with_pat_freq(allele_totals[c],ML_pat_allele_count,pat_allele_total,
        var_est_N,read_totals[c],bs_sample_size,modes[c]) for c in np.arange(len(read_counts))]).transpose()

    sampled_stats = calc_chi_statistic(sampled_counts,read_totals,expected_counts,allele_totals,
        count_mode,combination_index,combine)

    if tail == "two":
        accepted_samples = sampled_stats >= observed_stat
    if tail == "one":
        sampled_dirs = get_direction(sampled_counts,read_totals,allele_totals,count_mode,
            combination_index,combine)
        if observed_direction == expected_direction:
            accepted_samples = np.logical_and(sampled_dirs == observed_direction,sampled_stats >= observed_stat)
        else:
            accepted_samples = np.logical_or(np.logical_and(sampled_dirs == observed_direction,sampled_stats <= observed_stat),
                sampled_dirs != observed_direction)

    number_GEQ = np.sum(accepted_samples)

    return number_GEQ
    

# Helper function to run bootstrap sampling for single-comparison stats
def __bootstrap_sample_precalc_dist(read_counts,read_totals,ML_est_read_count_logProbs,
        expected_counts,allele_totals,
        tail,expected_direction,
        count_mode,combination_index,combine,
        observed_stat,observed_direction,
        bs_sample_size):
    # np.random.choice treats an integer parameter n as np.arange(n)
    sampled_counts = np.array([np.random.choice(read_totals[c]+1, p=np.exp(ML_est_read_count_logProbs[c]), 
        size=bs_sample_size) for c in np.arange(len(read_counts))]).transpose()

    sampled_stats = calc_chi_statistic(sampled_counts,read_totals,expected_counts,allele_totals,
        count_mode,combination_index,combine)

    # number_GEQ = np.sum(sampled_stats >= observed_stat)

    if tail == "two":
        accepted_samples = sampled_stats >= observed_stat
    if tail == "one":
        sampled_dirs = get_direction(sampled_counts,read_totals,allele_totals,count_mode,
            combination_index,combine)
        if observed_direction == expected_direction:
            accepted_samples = np.logical_and(sampled_dirs == observed_direction,sampled_stats >= observed_stat)
        else:
            accepted_samples = np.logical_or(np.logical_and(sampled_dirs == observed_direction,sampled_stats <= observed_stat),
                sampled_dirs != observed_direction)

    number_GEQ = np.sum(accepted_samples)

    return number_GEQ

# For generating the p-value as the accepted bootstrap samples out of the total bootstrap samples,
#    for comparisons with only one cohort comparison involved
def calc_chi_stat_pval_bootstrap(read_counts,read_totals,expected_counts,ML_est_read_count_logProbs,ML_pat_allele_count,pat_allele_total,
        var_est_N,modes,allele_totals=None,count_mode="read",combination_index=None,combine=False,tail="two",expected_direction=None):
    # Check the inputs (may be redundant, also done when calculating the chi statistic for the observed counts)
    if allele_totals is None and count_mode == "allele":
        raise ValueError("can only calculate an allele-count chi statistic when allele_totals are provided")
    if combination_index is None and combine:
        raise ValueError("can only combine counts when the combination_index at which to combine into two groups is provided")
    if combine and count_mode == "read":
        raise ValueError("only expected to combine the 1...n other counts if considering allele_counts")
    if count_mode != "read" and count_mode != "allele":
        raise ValueError("chi statistic calculation count_mode must only be one of \"read\" or \"allele\"")
    if tail == "one" and expected_direction is None:
        raise ValueError("one tailed statistic tests require a direction to polarize the test")

    # expected_freqs = [expected_counts[i]/read_totals[i] for i in np.arange(len(read_counts))]

    if verbose: print("Calculating the chi statistics for random samplings from the estimated neutral distribution "+\
        "and comparing to the observed for p-value generation")

    observed_stat = calc_chi_statistic(read_counts,read_totals,expected_counts,allele_totals,
        count_mode,combination_index,combine)

    if tail=="one":
        observed_direction = get_direction(read_counts,read_totals,allele_totals,count_mode,combination_index,combine)
    else:
        observed_direction = None


    if verbose:
        if bs_precalc_dist:
            logProb_obs_counts = __get_logProb_counts(read_counts,ML_est_read_count_logProbs)
        else:
            logProb_obs_counts = "N/A - bootstrap sampling w/o keeping log probs of the complete distribution"
        print("Observed Counts:\t"+str(read_counts)+"\nExpected Counts:\t"+str(expected_counts)+\
            "\nTotal Counts:\t\t"+str(read_totals)+"\nLog Prob of Counts:\t"+str(logProb_obs_counts)+\
            "\nChi Square:\t\t"+str(observed_stat))

    # poss_counts = []
    # ML_est_read_count_Probs = []
    # for c in np.arange(len(read_counts)):
    #     poss_counts += [np.arange(read_totals[c])]
    #     ML_est_read_count_Probs += [np.exp(ML_est_read_count_logProbs[c])]
    #     print("Log prob conversion check for cohort "+str(c)+": "+str(sum(ML_est_read_count_Probs[c])))

    if verbose:
        if bs_precalc_dist:
            for c in np.arange(len(read_counts)):
                ML_est_read_count_probs = np.exp(ML_est_read_count_logProbs[c])
                print("Log prob conversion check for cohort "+str(c)+": "+str(sum(ML_est_read_count_probs))+\
                    " total probability in "+str(len(ML_est_read_count_probs))+" values")

    number_GEQ = 0
    num_std_bs_chunks = N_bootstrap//bs_chunk
    remainder_bs_N = N_bootstrap%bs_chunk
    # Bootstrap sample in smaller chunks to lower the total memory requirement,
    #    but save execution time by vectorization
    for i in np.arange(num_std_bs_chunks):
        if bs_precalc_dist:
            number_GEQ += __bootstrap_sample_precalc_dist(read_counts,read_totals,ML_est_read_count_logProbs,
                expected_counts,allele_totals,tail,expected_direction,count_mode,combination_index,combine,
                observed_stat,observed_direction,bs_chunk)
        else:
            number_GEQ += __bootstrap_sample(read_counts,read_totals,ML_pat_allele_count,pat_allele_total,
                var_est_N,modes,expected_counts,allele_totals,tail,expected_direction,count_mode,combination_index,combine,
                observed_stat,observed_direction,bs_chunk)
    # Sample the remainder
    if bs_precalc_dist:
        number_GEQ += __bootstrap_sample_precalc_dist(read_counts,read_totals,ML_est_read_count_logProbs,
            expected_counts,allele_totals,tail,expected_direction,count_mode,combination_index,combine,
            observed_stat,observed_direction,remainder_bs_N)
    else:
        number_GEQ += __bootstrap_sample(read_counts,read_totals,ML_pat_allele_count,pat_allele_total,
            var_est_N,modes,expected_counts,allele_totals,tail,expected_direction,count_mode,combination_index,combine,
            observed_stat,observed_direction,remainder_bs_N)

    pval = number_GEQ/N_bootstrap

    # log_pval = np.log(pval)

    return (observed_stat,pval,number_GEQ,N_bootstrap)


# A helper function for extracting the direction of frequency change, vectorized
#   only defined for combined statistics or where only 2 cohorts are considered
#  RETURNS 1 if the second cohort set has frequency > the first (a strict increase), 0 otherwise
def get_direction(read_counts,read_totals,allele_totals=None,count_mode="read",
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


    if not combine and len(allele_totals) > 2:
        raise ValueError("can only calculate a direction change for combined cohorts or two cohorts")

    # Generate the array of counts to compare, either reads or alleles,
    #  and combine the allele counts on either side of the combination index if specified
    allele_totals = np.array(allele_totals)
    if count_mode == "read":
        pos_counts = np.array(read_counts)
        # neg_counts = np.subtract(read_totals,pos_counts)
        pos_freqs = np.divide(pos_counts,read_totals)
        # neg_freqs = np.divide(neg_counts,read_totals)
    elif count_mode == "allele":
        # expected_freqs = np.divide(expected_counts,read_totals)
        pos_proportion = np.divide(read_counts,read_totals)
        pos_counts = np.multiply(pos_proportion,allele_totals)
        # neg_counts = np.subtract(allele_totals,pos_counts)
        if combine:
            last_dim = len(pos_counts.shape)-1
            pos_counts = np.stack((np.sum(pos_counts[...,:combination_index],axis=last_dim),
                np.sum(pos_counts[...,combination_index:],axis=last_dim)),axis=last_dim)
            # neg_counts = np.stack((np.sum(neg_counts[...,:combination_index],axis=last_dim),
            #     np.sum(neg_counts[...,combination_index:],axis=last_dim)),axis=last_dim)
            last_allele_dim = len(allele_totals.shape)-1
            total_allele_counts = np.stack((np.sum(allele_totals[...,:combination_index],axis=last_allele_dim),
                np.sum(allele_totals[...,combination_index:],axis=last_allele_dim)),axis=last_allele_dim)

            pos_freqs = np.divide(pos_counts,total_allele_counts)
            # neg_freqs = np.divide(neg_counts,total_allele_counts)
        else:
            pos_freqs = np.divide(pos_counts,allele_totals)
            # neg_freqs = np.divide(neg_counts,allele_totals)           
    else:
        raise ValueError("Unexpected count_mode value: "+str(count_mode))

    # Calculate the difference
    last_dim = len(pos_freqs.shape)-1
    # Asserted that there are only 2 frequencies at this point, but can again
    assert(pos_freqs.shape[-1]==2),"Should only ever have 2 uncombined or 2 combined frequencies per entry at this point"
    freq_difference = np.subtract(pos_freqs[...,1],pos_freqs[...,0])
    strict_increase = freq_difference > 0

    if verbose:
        strict_decrease = freq_difference < 0
        num_neither = np.sum(np.logical_and(np.logical_not(strict_increase),np.logical_not(strict_decrease)))
        assert(num_neither == np.sum(freq_difference == 0))
        if num_neither > 0:
            print("Encountered "+str(num_neither)+" direction checks with no frequency change at all/in either direction")
            print(pos_freqs)
            print(freq_difference)

    return strict_increase


# Helper function to run bootstrap sampling for mult-comparison stats
def __bootstrap_sample_mult_comp(read_counts,read_totals,allele_totals,ML_est_read_count_logProbs,ML_pat_allele_count,
        pat_allele_total,var_est_N,index_translations,modes,read_totals_1,expected_counts_1,allele_totals_1,comb_ind_1,
        read_totals_2,expected_counts_2,allele_totals_2,comb_ind_2,count_mode,combine,tail,expected_directions,
        observed_stat_1,observed_direction_1,observed_stat_2,observed_direction_2,bs_sample_size):

    if bs_precalc_dist:
        # np.random.choice treats an integer parameter n as np.arange(n)
        sampled_counts = [np.random.choice(read_totals[c]+1, p=np.exp(ML_est_read_count_logProbs[c]), 
            size=bs_sample_size) for c in np.arange(len(read_counts))]
    else:
        sampled_counts = [sample_cohort_with_pat_freq(allele_totals[c],ML_pat_allele_count,pat_allele_total,
            var_est_N,read_totals[c],bs_sample_size,modes[c]) for c in np.arange(len(read_counts))]

    sampled_counts_1 = np.array([sampled_counts[c] for c in index_translations[0]+index_translations[1]]).transpose()
    sampled_stats_1 = calc_chi_statistic(sampled_counts_1,read_totals_1,expected_counts_1,allele_totals_1,
        count_mode,comb_ind_1,combine)
    sampled_dirs_1 = get_direction(sampled_counts_1,read_totals_1,allele_totals_1,count_mode,
        comb_ind_1,combine)

    sampled_counts_2 = np.array([sampled_counts[c] for c in index_translations[2]+index_translations[3]]).transpose()
    sampled_stats_2 = calc_chi_statistic(sampled_counts_2,read_totals_2,expected_counts_2,allele_totals_2,
        count_mode,comb_ind_2,combine)
    sampled_dirs_2 = get_direction(sampled_counts_2,read_totals_2,allele_totals_2,count_mode,
        comb_ind_2,combine)

    if tail == "two":
        if observed_direction_1 == int(not bool(expected_directions[0])):
            accept_comp_1_rev = np.logical_and(sampled_dirs_1 == observed_direction_1,sampled_stats_1 >= observed_stat_1)
        else:
            accept_comp_1_rev = np.logical_or(np.logical_and(sampled_dirs_1 == observed_direction_1,sampled_stats_1 <= observed_stat_1),
                sampled_dirs_1 != observed_direction_1)

        if observed_direction_2 == int(not bool(expected_directions[1])):
            accept_comp_2_rev = np.logical_and(sampled_dirs_2 == observed_direction_2,sampled_stats_2 >= observed_stat_2)
        else:
            accept_comp_2_rev = np.logical_or(np.logical_and(sampled_dirs_2 == observed_direction_2,sampled_stats_2 <= observed_stat_2),
                sampled_dirs_2 != observed_direction_2)

    if observed_direction_1 == expected_directions[0]:
        accept_comp_1 = np.logical_and(sampled_dirs_1 == observed_direction_1,sampled_stats_1 >= observed_stat_1)
    else:
        accept_comp_1 = np.logical_or(np.logical_and(sampled_dirs_1 == observed_direction_1,sampled_stats_1 <= observed_stat_1),
            sampled_dirs_1 != observed_direction_1)

    if observed_direction_2 == expected_directions[1]:
        accept_comp_2 = np.logical_and(sampled_dirs_2 == observed_direction_2,sampled_stats_2 >= observed_stat_2)
    else:
        accept_comp_2 = np.logical_or(np.logical_and(sampled_dirs_2 == observed_direction_2,sampled_stats_2 <= observed_stat_2),
            sampled_dirs_2 != observed_direction_2)

    # print(tail)
    # print(accept_comp_1)
    # print(accept_comp_2)

    # number_GEQ = np.sum(np.logical_and(accept_comp_1,accept_comp_2))
    if tail == "two":
        number_GEQ = np.sum(np.logical_or(np.logical_or(accept_comp_1,accept_comp_2),np.logical_or(accept_comp_1_rev,accept_comp_2_rev)))
    elif tail == "one":
        number_GEQ = np.sum(np.logical_or(accept_comp_1,accept_comp_2))

    return number_GEQ

# For generating the p-value as the accepted bootstrap samples out of the total bootstrap samples,
#    for comparisons with only one cohort comparison involved
def calc_mult_comb_allele_cohort_chi_stat_pval_bootstrap(read_counts,read_totals,expected_counts,ML_est_read_count_logProbs,
        ML_pat_allele_count,pat_allele_total,var_est_N,allele_totals,index_translations,modes,tail,expected_directions,
        count_mode="allele",combine=True):

    if verbose: print("Performing bootstrap sampling to estimate the p-value: Calculating the chi statistics for both comparisons "+\
        "in each sample from the estimated neutral distribution and comparing to the observed chi stats for acceptance")

    read_counts_1  = [read_counts[i] for i in index_translations[0]+index_translations[1]]
    read_totals_1 = [read_totals[i] for i in index_translations[0]+index_translations[1]]
    expected_counts_1 = [expected_counts[i] for i in index_translations[0]+index_translations[1]]
    allele_totals_1 = [allele_totals[i] for i in index_translations[0]+index_translations[1]]
    comb_ind_1 = len(index_translations[0])

    observed_stat_1 = calc_chi_statistic(read_counts_1,read_totals_1,expected_counts_1,allele_totals_1,
        count_mode,comb_ind_1,combine)
    observed_direction_1 = get_direction(read_counts_1,read_totals_1,allele_totals_1,count_mode,
        comb_ind_1,combine)


    read_counts_2  = [read_counts[i] for i in index_translations[2]+index_translations[3]]
    read_totals_2 = [read_totals[i] for i in index_translations[2]+index_translations[3]]
    expected_counts_2 = [expected_counts[i] for i in index_translations[2]+index_translations[3]]
    allele_totals_2 = [allele_totals[i] for i in index_translations[2]+index_translations[3]]
    comb_ind_2 = len(index_translations[2])

    observed_stat_2 = calc_chi_statistic(read_counts_2,read_totals_2,expected_counts_2,allele_totals_2,
        count_mode,comb_ind_2,combine)
    observed_direction_2 = get_direction(read_counts_2,read_totals_2,allele_totals_2,count_mode,
        comb_ind_2,combine)


    if verbose:
        if bs_precalc_dist:
            logProb_obs_counts = __get_logProb_counts(read_counts,ML_est_read_count_logProbs)
        else:
            logProb_obs_counts = "N/A - bootstrap sampling w/o keeping log probs of the complete distribution"
        print("Observed Counts:\t"+str(read_counts)+"\nExpected Counts:\t"+str(expected_counts)+\
            "\nPossible Totals:\t\t"+str(read_totals)+"\nLog Prob of Counts:\t"+str(logProb_obs_counts)+\
            "\nChi Square 1:\t\t"+str(observed_stat_1)+"\nIncrease 1:\t\t"+str(bool(observed_direction_1))+\
            "\nChi Square 2:\t\t"+str(observed_stat_2)+"\nIncrease 2:\t\t"+str(bool(observed_direction_2)))
        if bs_precalc_dist:
            for c in np.arange(len(read_counts)):
                ML_est_read_count_probs = np.exp(ML_est_read_count_logProbs[c])
                print("Log prob conversion check for cohort "+str(c)+": "+str(sum(ML_est_read_count_probs))+\
                    " total probability in "+str(len(ML_est_read_count_probs))+" values")

    number_GEQ = 0
    num_std_bs_chunks = N_bootstrap//bs_chunk
    remainder_bs_N = N_bootstrap%bs_chunk
    # Bootstrap sample in smaller chunks to lower the total memory requirement,
    #    but save execution time by vectorization
    for i in np.arange(num_std_bs_chunks):
        number_GEQ += __bootstrap_sample_mult_comp(read_counts,read_totals,allele_totals,
            ML_est_read_count_logProbs,ML_pat_allele_count,pat_allele_total,var_est_N,
            index_translations,modes,read_totals_1,expected_counts_1,allele_totals_1,comb_ind_1,
            read_totals_2,expected_counts_2,allele_totals_2,comb_ind_2,count_mode,combine,
            tail,expected_directions,observed_stat_1,observed_direction_1,observed_stat_2,observed_direction_2,
            bs_chunk)
    # Sample the remainder
    if remainder_bs_N > 0:
        number_GEQ += __bootstrap_sample_mult_comp(read_counts,read_totals,allele_totals,
            ML_est_read_count_logProbs,ML_pat_allele_count,pat_allele_total,var_est_N,
            index_translations,modes,read_totals_1,expected_counts_1,allele_totals_1,comb_ind_1,
            read_totals_2,expected_counts_2,allele_totals_2,comb_ind_2,count_mode,combine,
            tail,expected_directions,observed_stat_1,observed_direction_1,observed_stat_2,observed_direction_2,
            remainder_bs_N)

    pval = number_GEQ/N_bootstrap

    observed_stats = [observed_stat_1,observed_stat_2]

    return (observed_stats,pval,number_GEQ,N_bootstrap)

# To return a boolean array, where an element is true iff
#   the corresponding index in the input array is within the given tolerance of the expectation
def __test_sample_acceptance(data,expectation,census_totals,tolerance,mode="absolute"):

    if verbose: print("Testing the sampled read acceptance:\n"+\
        "Tolerance: "+str(tolerance)+"\n"+\
        "Observed read counts: "+str(expectation)+"\n"+\
        "Census totals: "+str(census_totals)+"\n"+\
        "Accepted lower bound: "+str(np.add(expectation,np.multiply(census_totals,-tolerance)))+"\n"+\
        "Accepted upper bound: "+str(np.add(expectation,np.multiply(census_totals,tolerance)))+"\n"+\
        "Samples: "+str(data)+"\n")


    assert (tolerance > 0), "sampling tolerance given to __test_sample_acceptance must be positive"
    if mode == "proportion":
        return np.all(np.logical_and(data > np.add(expectation,np.multiply(expectation,1-tolerance)),
            data < np.add(expectation,np.multiply(expectation,1+tolerance))),axis=1)
    elif mode == "absolute":
        return np.all(np.logical_and(data > np.add(expectation,np.multiply(census_totals,-tolerance)),
            data < np.add(expectation,np.multiply(census_totals,tolerance))),axis=1)
    else:
        raise ValueError("__test_sample_acceptance only takes modes of \"proportion\" and \"absolute\" when calculating "+\
            "the acceptance of given data within a given (proportion or absolute value) tolerance of the given expectation")
    return


def __ABC_sel_coeff_sample_single_comp(read_counts,read_totals,pat_allele_count,pat_allele_total,
        var_est_N,modes,allele_totals,tail,expected_direction,count_mode,combination_index,combine,
        observed_direction,s_tolerance,WF_est_sel_coeff,bs_chunk):
    
    # sel_coefficients = np.arange()
    rng = np.random.default_rng()
    selection_magnitude = abs(WF_est_sel_coeff)*10
    selection_coefficients = rng.uniform(low=-1.0, high=max(1.0,selection_magnitude), size=bs_chunk)

    # np.random.choice treats an integer parameter n as np.arange(n)
    sampled_sel_counts = np.array([sample_cohort_with_pat_freq_and_selection(allele_totals[c],pat_allele_count,pat_allele_total,
        var_est_N,read_totals[c],bs_chunk,selection_coefficients,modes[c]) for c in np.arange(len(read_counts))]).transpose()

    expected_sel_counts = np.array(read_counts).transpose()
    census_totals = np.array(read_totals).transpose()
    # print((sampled_sel_counts,expected_sel_counts,s_tolerance))
    accepted_samples = __test_sample_acceptance(sampled_sel_counts,expected_sel_counts,census_totals,s_tolerance)
    # print(selection_coefficients)
    # print(accepted_samples)

    accepted_sel_coeff = selection_coefficients[accepted_samples]

    return accepted_sel_coeff



# For generating the selection coefficient as the average selection coefficient in accepted simulations
#    generated from a [0,1) uniform distribution on selection magnitudes,
#    with paternal frequency taken directly from paternal sample frequency.
#    For comparisons with multiple cohorts involved
# NOTE - some of the opposing direction calls will have such values because of the 
def calc_single_comp_sel_coeff_ABC(read_counts,read_totals,pat_allele_count,
        pat_allele_total,var_est_N,modes,s_tolerance,quantile,WF_dir_estimate_s,
        allele_totals=None,count_mode="read",combination_index=None,combine=False,tail="two",expected_direction=None):
    if verbose: print("Performing ABC estimate of selection coefficients from simulations "+\
        "sampling from a uniform [0,1) prior on selection strength "+\
        "for paternal-offspring comparisons")

    # if verbose:
    #     if bs_precalc_dist:
    #         logProb_obs_counts = __get_logProb_counts(read_counts,ML_est_read_count_logProbs)
    #     else:
    #         logProb_obs_counts = "N/A - bootstrap sampling w/o keeping log probs of the complete distribution"
    #     print("Observed Counts:\t"+str(read_counts)+"\nExpected Counts:\t"+str(expected_counts)+\
    #         "\nTotal Counts:\t\t"+str(read_totals)+"\nLog Prob of Counts:\t"+str(logProb_obs_counts)+\
    #         "\nChi Square:\t\t"+str(observed_stat))

    # if verbose:
    #     if bs_precalc_dist:
    #         for c in np.arange(len(read_counts)):
    #             ML_est_read_count_probs = np.exp(ML_est_read_count_logProbs[c])
    #             print("Log prob conversion check for cohort "+str(c)+": "+str(sum(ML_est_read_count_probs))+\
    #                 " total probability in "+str(len(ML_est_read_count_probs))+" values")


    if tail=="one":
        observed_direction = get_direction(read_counts,read_totals,allele_totals,count_mode,combination_index,combine)
    else:
        observed_direction = None

    # Bootstrap sample in smaller chunks to lower the total memory requirement,
    #    but save execution time by vectorization
    num_std_bs_chunks = N_bootstrap//bs_chunk
    remainder_bs_N = N_bootstrap%bs_chunk
    accepted_s = np.array([])
    for i in np.arange(num_std_bs_chunks):
        a = __ABC_sel_coeff_sample_single_comp(read_counts,read_totals,pat_allele_count,pat_allele_total,
            var_est_N,modes,allele_totals,tail,expected_direction,count_mode,combination_index,combine,
            observed_direction,s_tolerance,WF_dir_estimate_s,bs_chunk)
        accepted_s = np.append(accepted_s,a)
    # Sample the remainder
    if remainder_bs_N > 0:
        a = __ABC_sel_coeff_sample_single_comp(read_counts,read_totals,pat_allele_count,pat_allele_total,
            var_est_N,modes,allele_totals,tail,expected_direction,count_mode,combination_index,combine,
            observed_direction,s_tolerance,WF_dir_estimate_s,remainder_bs_N)
        accepted_s = np.append(accepted_s,a)

    if verbose:
        print("Calculating upper and lower "+str(quantile)+" quantiles from "+str(len(accepted_s))+\
            " accepted selection coefficients")
        print(len(accepted_s))
        print(accepted_s)


    if len(a) == 0:
        print("WARNING: no accepted selection coefficient simulations\n"+\
            "Tolerance: "+str(s_tolerance)+"\n"+\
            "Read counts: "+str(read_counts)+"\n"+\
            "Census totals: "+str(read_totals)+"\n"+\
            "Accepted lower bound: "+str(np.add(read_counts,np.multiply(read_totals,-s_tolerance)))+"\n"+\
            "Accepted upper bound: "+str(np.add(read_counts,np.multiply(read_totals,s_tolerance)))+"\n"+\
            "Combination index: "+str(combination_index)+"\n"+\
            "Paternal/Offspring modes: "+str(modes)+"\n")
        
        sel_coefficient = [None,None,None]
    else:
        s = np.mean(accepted_s)
        lower_quant = np.quantile(accepted_s,quantile/2)
        upper_quant = np.quantile(accepted_s,1-quantile/2)

        sel_coefficient = [s, upper_quant, lower_quant]

    return sel_coefficient


def __ABC_sel_coeff_sample_mult_comp(read_counts,read_totals,allele_totals,
            pat_allele_count,pat_allele_total,var_est_N,index_translations,
            modes,expected_directions,#is_sequential_sel,
            s_tolerance,WF_est_sel_coeffs,bs_chunk):

    # if is_sequential_sel:
    #     assert(index_translations[1] == index_translations[2]), "Sequential selection requires "+\
    #         "that the second comparison starts with the second cohort from the first comparison"
    is_sequential_sel = index_translations[1] == index_translations[2]
    
    selection_coefficients = []
    rng = np.random.default_rng()
    assert(expected_directions is not None)
    selection_magnitude_1 = abs(WF_est_sel_coeffs[0])*10
    selection_magnitude_2 = abs(WF_est_sel_coeffs[1])*10
    if expected_directions is None:
        selection_coefficients += [rng.uniform(low=-1.0, high=max(1.0,selection_magnitude_1), size=bs_chunk)]
        selection_coefficients += [rng.uniform(low=-1.0, high=max(1.0,selection_magnitude_2), size=bs_chunk)]
    else:
        dir_1 = expected_directions[0]
        if bool(dir_1):
            selection_coefficients += [rng.uniform(low=0, high=max(1.0,selection_magnitude_1), size=bs_chunk)]
        else:
            selection_coefficients += [rng.uniform(low=-1.0, high=0.0, size=bs_chunk)]

        dir_2 = expected_directions[0]
        if bool(dir_1):
            selection_coefficients += [rng.uniform(low=0, high=max(1.0,selection_magnitude_2), size=bs_chunk)]
        else:
            selection_coefficients += [rng.uniform(low=-1.0, high=0.0, size=bs_chunk)]

        # for direction in expected_directions:
        #     if bool(direction):
        #         selection_coefficients += [rng.uniform(size=bs_chunk)]
        #     else:
        #         selection_coefficients += [rng.uniform(low=-1.0, high=0.0, size=bs_chunk)]

    # print(index_translations)
    # print(len(read_counts))
    # print(selection_coefficients)


    sampled_sel_counts = []
    expected_sel_counts = []
    census_totals = []
    for c in index_translations[1]:
        expected_sel_counts += [read_counts[c]]
        census_totals += [read_totals[c]]
        sampled_sel_counts += [sample_cohort_with_pat_freq_and_selection(allele_totals[c],pat_allele_count,pat_allele_total,
            var_est_N,read_totals[c],bs_chunk,selection_coefficients[0],modes[c])]
    for c in index_translations[3]:
        expected_sel_counts += [read_counts[c]]
        census_totals += [read_totals[c]]
        if is_sequential_sel:
            sampled_sel_counts += [sample_cohort_with_pat_freq_and_selection(allele_totals[c],pat_allele_count,pat_allele_total,
                var_est_N,read_totals[c],bs_chunk,np.multiply(selection_coefficients[0],selection_coefficients[1]),modes[c])]
        else:
            sampled_sel_counts += [sample_cohort_with_pat_freq_and_selection(allele_totals[c],pat_allele_count,pat_allele_total,
                var_est_N,read_totals[c],bs_chunk,selection_coefficients[1],modes[c])]

    sampled_sel_counts = np.array(sampled_sel_counts).transpose()
    expected_sel_counts = np.array(expected_sel_counts).transpose()
    # expected_sel_counts = np.array(read_counts).transpose()
    accepted_samples = __test_sample_acceptance(sampled_sel_counts,expected_sel_counts,census_totals,s_tolerance)
    accepted_sel_coeff = [s[accepted_samples] for s in selection_coefficients]
    # number_accepted = len(accepted_sel_coeff[0])

    return accepted_sel_coeff



# For generating the selection coefficient as the average selection coefficient in accepted simulations
#    generated from a [0,1) uniform distribution on selection magnitudes,
#    with paternal frequency taken directly from paternal sample frequency.
#    For comparisons with multiple cohorts involved
# NOTE - some of the opposing direction calls will have such values because of the 
def calc_mult_comp_sel_coeff_ABC(read_counts,read_totals,pat_allele_count,
        pat_allele_total,var_est_N,allele_totals,index_translations,modes,expected_directions,
        s_tolerance,quantile,WF_est_sel_coeffs):
    if verbose: print("Performing ABC estimate of selection coefficients from simulations "+\
        "sampling from a uniform [0,1) prior on selection strength "+\
        "for paternal-offspring comparisons")

    # read_counts_1  = [read_counts[i] for i in index_translations[0]+index_translations[1]]
    # read_totals_1 = [read_totals[i] for i in index_translations[0]+index_translations[1]]
    # allele_totals_1 = [allele_totals[i] for i in index_translations[0]+index_translations[1]]
    # comb_ind_1 = len(index_translations[0])


    # read_counts_2  = [read_counts[i] for i in index_translations[2]+index_translations[3]]
    # read_totals_2 = [read_totals[i] for i in index_translations[2]+index_translations[3]]
    # allele_totals_2 = [allele_totals[i] for i in index_translations[2]+index_translations[3]]
    # comb_ind_2 = len(index_translations[2])

    assert(all([m == "off" for m in modes[1:]])), "Only configured to handle cohort groups containing <= 1 paternal"
    if modes[0] == "pat":
        if index_translations[0][0] == 0:
            assert(len(index_translations[0])==1),"A group with a paternal cohort must contain no other cohorts"
            mode_1 = "pat-off"
        else:
            print("paternal cohort given but not as first comparison - UNEXPECTED") # This should never trigger
            mode_1 = "off-off"
        if index_translations[2][0] == 0:
            assert(len(index_translations[2])==1),"A group with a paternal cohort must contain no other cohorts"
            mode_2 = "pat-off"
        else:
            mode_2 = "off-off"
    else:
        assert(modes[0] == "off")
        mode_1 = "off-off"
        mode_2 = "off-off"

    num_std_bs_chunks = N_bootstrap//bs_chunk
    remainder_bs_N = N_bootstrap%bs_chunk
    # Bootstrap sample in smaller chunks to lower the total memory requirement,
    #    but save execution time by vectorization
    accepted_s1 = np.array([])
    accepted_s2 = np.array([])
    if verbose: print("Calculating accepted selection coefficient draws")
    for i in np.arange(num_std_bs_chunks):
        [a1, a2] = __ABC_sel_coeff_sample_mult_comp(read_counts,read_totals,allele_totals,
            pat_allele_count,pat_allele_total,var_est_N,index_translations,
            modes,expected_directions,s_tolerance,WF_est_sel_coeffs,bs_chunk)
        print([a1, a2])
        accepted_s1 = np.append(accepted_s1,a1)
        accepted_s2 = np.append(accepted_s2,a2)
    # Sample the remainder
    if remainder_bs_N > 0:
        [a1, a2] = __ABC_sel_coeff_sample_mult_comp(read_counts,read_totals,allele_totals,
            pat_allele_count,pat_allele_total,var_est_N,index_translations,
            modes,expected_directions,s_tolerance,WF_est_sel_coeffs,remainder_bs_N)
        print([a1, a2])
        accepted_s1 = np.append(accepted_s1,a1)
        accepted_s2 = np.append(accepted_s2,a2)

    # (accepted_s1, accepted_s2) = zip(*accepted_s_pairs)
    # s_1, upper_quant_1, lower_quant_1 = __get_sel_coeff_and_quants(accepted_s1,quantile)
    # s_2, upper_quant_2, lower_quant_2 = __get_sel_coeff_and_quants(accepted_s2,quantile)

    assert(len(a1) == len(a2)), "With two selection coefficients, there should be equally many accepted values/simulations"
    if len(a1) == 0 and len(a2) == 0:
        print("WARNING: no accepted selection coefficient simulations\n"+\
            "Tolerance: "+str(s_tolerance)+"\n"+\
            "Read counts: "+str(read_counts)+"\n"+\
            "Census totals: "+str(read_totals)+"\n"+\
            "Read count indexes: "+str(index_translations)+"\n"+\
            "Paternal/Offspring modes: "+str(modes)+"\n")
        
        sel_coefficients = [None,None,None,None,None,None]
    else:
        s_1 = np.mean(accepted_s1)
        lower_quant_1 = np.quantile(accepted_s1,quantile/2)
        upper_quant_1 = np.quantile(accepted_s1,1-quantile/2)
        s_2 = np.mean(accepted_s2)
        lower_quant_2 = np.quantile(accepted_s2,quantile/2)
        upper_quant_2 = np.quantile(accepted_s2,1-quantile/2)
        sel_coefficients = [s_1, upper_quant_1, lower_quant_1, s_2, upper_quant_2, lower_quant_2]

    return sel_coefficients


# calculates the product of the two frequency changes in a comparison,
#  as a statistic for assessing the significance of the frequency change, with vectorized opertaions
#  ASSUMES that either 2 read counts are given, or counts are combined into two allele count categories
def calc_freq_change(read_counts,read_totals,expected_counts,allele_totals=None,count_mode="read",
        combination_index=None,combine=False,freq_change_metric="absolute"):
    # Check the inputs
    if allele_totals is None and count_mode == "allele":
        raise ValueError("can only calculate an allele-count frequency change product statistic when allele_totals are provided")
    if combination_index is None and combine:
        raise ValueError("can only combine counts when the combination_index at which to combine into two groups is provided")
    if combine and count_mode != "allele":
        raise ValueError("only expected to combine the 1...n other counts if considering allele_counts")
    if count_mode != "read" and count_mode != "allele":
        raise ValueError("frequency change product statistic calculation count_mode must only be one of \"read\" or \"allele\"")

    # Generate the array of counts to compare, either reads or alleles,
    #  and combine the allele counts on either side of the combination index if specified

    # expected_freqs = np.divide(expected_counts,read_totals)
    pos_proportion = np.divide(read_counts,read_totals)
    if count_mode == "allele" and combine:
        pos_counts = np.multiply(pos_proportion,allele_totals)
        # exp_pos_counts = np.multiply(expected_freqs,allele_totals)
        allele_tots = np.array(allele_totals)

        last_dim = len(pos_counts.shape)-1
        comb_pos_counts = np.stack((np.sum(pos_counts[...,:combination_index],axis=last_dim),
            np.sum(pos_counts[...,combination_index:],axis=last_dim)),axis=last_dim)
        # last_dim_exp = len(exp_pos_counts.shape)-1
        # comb_exp_pos_counts = np.stack((np.sum(exp_pos_counts[...,:combination_index],axis=last_dim_exp),
        #     np.sum(exp_pos_counts[...,combination_index:],axis=last_dim_exp)),axis=last_dim_exp)
        last_dim_allele = len(allele_tots.shape)-1
        comb_allele_totals = np.stack((np.sum(allele_tots[...,:combination_index],axis=last_dim_allele),
            np.sum(allele_tots[...,combination_index:],axis=last_dim_allele)),axis=last_dim_allele)

        # expected_freqs = np.divide(comb_exp_pos_counts,comb_allele_totals)
        pos_proportion = np.divide(comb_pos_counts,comb_allele_totals)
    else:
        raise ValueError("Unexpected count_mode value: "+str(count_mode))

    # Calculate the frequency change statistic (at the relevant vectorized level)
    last_dim = len(pos_proportion.shape)-1
    if freq_change_metric == "absolute":
        # change = np.subtract(pos_proportion[...,1],pos_proportion[...,0],axis=last_dim)
        change = np.subtract(pos_proportion[...,1],pos_proportion[...,0])
    if freq_change_metric == "proportion":
        # change = np.divide(pos_proportion[...,1],pos_proportion[...,0],axis=last_dim)
        change = np.divide(pos_proportion[...,1],pos_proportion[...,0])

    return change


# Helper function to run bootstrap sampling for mult-comparison frequency change product stats
def __bootstrap_sample_mult_comp_freq_change_prod(read_counts,read_totals,allele_totals,ML_est_read_count_logProbs,ML_pat_allele_count,
        pat_allele_total,var_est_N,index_translations,modes,read_totals_1,expected_counts_1,allele_totals_1,comb_ind_1,
        read_totals_2,expected_counts_2,allele_totals_2,comb_ind_2,count_mode,combine,tail,expected_directions,observed_freq_change_product,
        observed_direction_1,observed_direction_2,bs_sample_size):

    if bs_precalc_dist:
        # np.random.choice treats an integer parameter n as np.arange(n)
        sampled_counts = [np.random.choice(read_totals[c]+1, p=np.exp(ML_est_read_count_logProbs[c]), 
            size=bs_sample_size) for c in np.arange(len(read_counts))]
    else:
        sampled_counts = [sample_cohort_with_pat_freq(allele_totals[c],ML_pat_allele_count,pat_allele_total,
            var_est_N,read_totals[c],bs_sample_size,modes[c]) for c in np.arange(len(read_counts))]

    sampled_counts_1 = np.array([sampled_counts[c] for c in index_translations[0]+index_translations[1]]).transpose()
    sampled_freq_changes_1 = calc_freq_change(sampled_counts_1,read_totals_1,expected_counts_1,allele_totals_1,
        count_mode,comb_ind_1,combine)
    sampled_dirs_1 = get_direction(sampled_counts_1,read_totals_1,allele_totals_1,count_mode,
        comb_ind_1,combine)

    sampled_counts_2 = np.array([sampled_counts[c] for c in index_translations[2]+index_translations[3]]).transpose()
    sampled_freq_changes_2 = calc_freq_change(sampled_counts_2,read_totals_2,expected_counts_2,allele_totals_2,
        count_mode,comb_ind_2,combine)
    sampled_dirs_2 = get_direction(sampled_counts_2,read_totals_2,allele_totals_2,count_mode,
        comb_ind_2,combine)

    sampled_freq_change_products = np.multiply(sampled_freq_changes_1,sampled_freq_changes_2)


    if tail == "two":
        number_GEQ = np.sum(sampled_freq_change_products<=observed_freq_change_product)
    else:
        assert(tail == "one")
        number_GEQ = np.sum(np.logical_and(sampled_dirs_1==expected_directions[0],sampled_dirs_2==expected_directions[1],sampled_freq_change_products<=observed_freq_change_product))

    return number_GEQ

# For generating the p-value as the accepted bootstrap samples out of the total bootstrap samples,
#    for comparisons with only one cohort comparison involved
def calc_mult_comb_allele_cohort_freq_change_product_stat_pval_bootstrap(read_counts,read_totals,
        expected_counts,ML_est_read_count_logProbs,ML_pat_allele_count,pat_allele_total,var_est_N,
        allele_totals,index_translations,modes,tail,expected_directions,count_mode="allele",combine=True):

    if verbose: print("Performing bootstrap sampling to estimate the p-value: Calculating the chi statistics for both comparisons "+\
        "in each sample from the estimated neutral distribution and comparing to the observed chi stats for acceptance")

    read_counts_1  = [read_counts[i] for i in index_translations[0]+index_translations[1]]
    read_totals_1 = [read_totals[i] for i in index_translations[0]+index_translations[1]]
    expected_counts_1 = [expected_counts[i] for i in index_translations[0]+index_translations[1]]
    allele_totals_1 = [allele_totals[i] for i in index_translations[0]+index_translations[1]]
    comb_ind_1 = len(index_translations[0])

    observed_freq_change_1 = calc_freq_change(read_counts_1,read_totals_1,expected_counts_1,allele_totals_1,
        count_mode,comb_ind_1,combine)
    observed_direction_1 = get_direction(read_counts_1,read_totals_1,allele_totals_1,count_mode,
        comb_ind_1,combine)


    read_counts_2  = [read_counts[i] for i in index_translations[2]+index_translations[3]]
    read_totals_2 = [read_totals[i] for i in index_translations[2]+index_translations[3]]
    expected_counts_2 = [expected_counts[i] for i in index_translations[2]+index_translations[3]]
    allele_totals_2 = [allele_totals[i] for i in index_translations[2]+index_translations[3]]
    comb_ind_2 = len(index_translations[2])

    observed_freq_change_2 = calc_freq_change(read_counts_2,read_totals_2,expected_counts_2,allele_totals_2,
        count_mode,comb_ind_2,combine)
    observed_direction_2 = get_direction(read_counts_2,read_totals_2,allele_totals_2,count_mode,
        comb_ind_2,combine)

    observed_freq_change_product = np.multiply(observed_freq_change_1,observed_freq_change_2)

    if verbose:
        if bs_precalc_dist:
            logProb_obs_counts = __get_logProb_counts(read_counts,ML_est_read_count_logProbs)
        else:
            logProb_obs_counts = "N/A - bootstrap sampling w/o keeping log probs of the complete distribution"
        print("Observed Counts:\t"+str(read_counts)+"\nExpected Counts:\t"+str(expected_counts)+\
            "\nPossible Totals:\t\t"+str(read_totals)+"\nLog Prob of Counts:\t"+str(logProb_obs_counts)+\
            "\nFreq Change 1:\t\t"+str(observed_freq_change_1)+"\nIncrease 1:\t\t"+str(bool(observed_direction_1))+\
            "\nFreq Change 2:\t\t"+str(observed_freq_change_2)+"\nIncrease 2:\t\t"+str(bool(observed_direction_2))+\
            "\nFreq Change Product:\t"+str(observed_freq_change_product))
        if bs_precalc_dist:
            for c in np.arange(len(read_counts)):
                ML_est_read_count_probs = np.exp(ML_est_read_count_logProbs[c])
                print("Log prob conversion check for cohort "+str(c)+": "+str(sum(ML_est_read_count_probs))+\
                    " total probability in "+str(len(ML_est_read_count_probs))+" values")

    number_GEQ = 0
    num_std_bs_chunks = N_bootstrap//bs_chunk
    remainder_bs_N = N_bootstrap%bs_chunk
    # Bootstrap sample in smaller chunks to lower the total memory requirement,
    #    but save execution time by vectorization
    for i in np.arange(num_std_bs_chunks):
        number_GEQ += __bootstrap_sample_mult_comp_freq_change_prod(read_counts,read_totals,allele_totals,
            ML_est_read_count_logProbs,ML_pat_allele_count,pat_allele_total,var_est_N,
            index_translations,modes,read_totals_1,expected_counts_1,allele_totals_1,comb_ind_1,
            read_totals_2,expected_counts_2,allele_totals_2,comb_ind_2,count_mode,combine,tail,
            expected_directions,observed_freq_change_product,observed_direction_1,observed_direction_2,
            bs_chunk)
    # Sample the remainder
    if remainder_bs_N > 0:
        number_GEQ += __bootstrap_sample_mult_comp_freq_change_prod(read_counts,read_totals,allele_totals,
            ML_est_read_count_logProbs,ML_pat_allele_count,pat_allele_total,var_est_N,
            index_translations,modes,read_totals_1,expected_counts_1,allele_totals_1,comb_ind_1,
            read_totals_2,expected_counts_2,allele_totals_2,comb_ind_2,count_mode,combine,tail,
            expected_directions,observed_freq_change_product,observed_direction_1,observed_direction_2,
            remainder_bs_N)

    pval = number_GEQ/N_bootstrap

    return (observed_freq_change_product,pval,number_GEQ,N_bootstrap)


# Helper function to reset the global variables when import in multiprocessing fails
def __reset_all_globals(N_bootstrap_temp,bs_chunk_temp,bs_precalc_dist_temp,verbose_temp,s_tolerance_temp,s_quantile_temp,
        calc_chi2_pval_temp,calc_prod_pval_temp,calc_sel_coeff_temp,count_mode_temp,target_mode_temp,tail_mode_temp):
    global N_bootstrap
    N_bootstrap = N_bootstrap_temp
    global bs_chunk
    bs_chunk = bs_chunk_temp
    global bs_precalc_dist
    bs_precalc_dist = bs_precalc_dist_temp

    global verbose
    verbose = verbose_temp

    global s_tolerance
    s_tolerance = s_tolerance_temp
    global s_quantile
    s_quantile = s_quantile_temp

    global calc_chi2_pval
    calc_chi2_pval = calc_chi2_pval_temp
    global calc_prod_pval
    calc_prod_pval = calc_prod_pval_temp
    global calc_sel_coeff
    calc_sel_coeff = calc_sel_coeff_temp

    global count_mode
    count_mode = count_mode_temp
    global target_mode
    target_mode = target_mode_temp
    global tail_mode
    tail_mode = tail_mode_temp
    return

def __reset_globals(N_bootstrap_temp,bs_chunk_temp,bs_precalc_dist_temp,verbose_temp,s_tolerance_temp,s_quantile_temp):
    global N_bootstrap
    N_bootstrap = N_bootstrap_temp
    global bs_chunk
    bs_chunk = bs_chunk_temp
    global bs_precalc_dist
    bs_precalc_dist = bs_precalc_dist_temp

    global verbose
    verbose = verbose_temp

    global s_tolerance
    s_tolerance = s_tolerance_temp
    global s_quantile
    s_quantile = s_quantile_temp
    return

# Manages a p-value calculation for a single comparison, relies on globally defined variables like cohort_group_dict
def comparison_pval(comparison,experiment,paternal_pop_allele_total,var_est_N,N_bootstrap,bs_chunk,bs_precalc_dist,verbose,
        s_tolerance,s_quantile,calc_chi2_pval,calc_prod_pval,calc_sel_coeff,calc_pat_from_off,
        count_mode="allele",target_mode="one_change_combined_cohorts",tail_mode="two",expected_direction=None):
    
    # Multiprocessing will deadlock if not using the "spawn" start method,
    #   and using "spawn" blanks the global variables, so they must be reset
    # __reset_all_globals(N_bootstrap,bs_chunk,bs_precalc_dist,verbose,s_tolerance,s_quantile,
    #     calc_chi2_pval,calc_prod_pval,calc_sel_coeff,count_mode,target_mode,tail_mode)
    __reset_globals(N_bootstrap,bs_chunk,bs_precalc_dist,verbose,s_tolerance,s_quantile)

    # Extract the relevant experimental data
    # experimental_cohort_data data formatted as [line,chrom,inv,pat_data,emb_data,eaf_data,laf_data,eam_data,lam_data]
    # Each *_data within experimental_cohort_data formatted as [cohort,allele_count,inv_freq,inv_read_count,std_read_count,read_total]
    line = experiment[0]
    chrom = experiment[1]
    inv = experiment[2]
    experimental_cohort_data = experiment[3:]

    if verbose:
        out_msg = "\nComparing "+str(comparison)+" for In("+str(chrom)+")"+str(inv)+\
            " differences in replicate Zi"+str(line)+"N\n\np-value calculation will asses "
        if target_mode=="one_change_combined_cohorts":
            out_msg += "the chi-square statistic on a 2x2 table with each group of counts combined to one, as a sum of "+\
                str(count_mode)+" counts assessing when the chi-square is greater than the observed"
        elif target_mode=="one_change_lone_cohorts":
            out_msg += "the chi-square statistic on a nx2 table with every cohort in either compared group considered separately by "+\
                str(count_mode)+" counts assessing when the chi-square is greater than the observed"
        elif target_mode=="two_changes_combined_cohorts":
            out_msg += "two separate chi-square statistics on two 2x2 tables with each group of counts combined to one, as a sum of "+\
                str(count_mode)+" counts assessing when the chi-squares are both greater than the observed, "+\
                "with the first comparison involving a frequency decrease and the second an increase"
        elif target_mode=="chi2_combined_two_way_ben_cost":
            out_msg += "two separate chi-square statistics on two 2x2 tables with each group of counts combined to one, as a sum of "+\
                str(count_mode)+" counts assessing when the chi-squares are both greater than the observed, "+\
                "with the first comparison involving a frequency decrease and the second an increase"
        if tail_mode=="two":
            out_msg += " as a two-tailed test "
        elif tail_mode=="one":
            out_msg += " as a one-tailed test "
        if expected_direction is not None:
            out_msg += " with expected direction "+str(expected_direction)+" where 1 is an increase, 0 a decrease"
        print(out_msg)

    if target_mode in mult_comp_stats:
        # Shouldn't need to check comparisons again for correct tuple use
        ((g11_name,g12_name),(g21_name,g22_name)) = comparison
        (g11,g12,g21,g22) = (cohort_group_dict[g11_name],cohort_group_dict[g12_name],
            cohort_group_dict[g21_name],cohort_group_dict[g22_name])
        cohorts_present = []
        index_translations = []
        for cohorts in (g11,g12,g21,g22):
            group_index_translations = []
            for cohort in cohorts:
                if len(cohorts_present) == 0:
                    cohorts_present += [cohort]
                    group_index_translations += [0]
                else:
                    j = 0
                    while j < len(cohorts_present) and cohorts_present[j] != cohort:
                        j += 1
                    if j == len(cohorts_present):
                        cohorts_present += [cohort]
                        group_index_translations += [j]
                    else:
                        assert(cohorts_present[j] == cohort),"Should always break the loop on one of the two conditions"
                        group_index_translations += [j]
            index_translations += [group_index_translations]

                 
        mode = None
        modes = None
        assert("pat" not in (g12_name,g21_name,g22_name)),"Script is only configured to interpret "+\
            "multiple comparison stats where the first given cohort is paternal (\"pat\") or none are."
        if g11_name == "pat":
            mode = "pat-off"
            modes = ['pat']+['off']*(len(cohorts_present)-1)
        else:
            mode = "off-off"
            modes = ['off']*len(cohorts_present)


        experimental_cohort_data = experiment[3:]
        # Each cohort_data formatted as [cohort,allele_count,inv_freq,inv_read_count,std_read_count,read_total]

        read_counts = []
        read_totals = []
        allele_totals = []
        expected_count_freqs = []
        for i in cohorts_present:
            read_counts += [experimental_cohort_data[i][3]]
            read_totals += [experimental_cohort_data[i][5]]
            allele_totals += [experimental_cohort_data[i][1]]
            expected_count_freqs += [experimental_cohort_data[i][3]/experimental_cohort_data[i][5]]

    else:
        # Generate the count data and grouping necessary for chi/pval calculation
        (group_1_name,group_2_name) = comparison
        (group_1,group_2) = (cohort_group_dict[group_1_name],cohort_group_dict[group_2_name])
        mode = None
        if group_1_name == "pat":
            mode = "pat-off"
            modes = ['pat']+['off']*(len(group_1)+len(group_2)-1)
        else:
            mode = "off-off"
            modes = ['off']*(len(group_1)+len(group_2))
        grouping_index = len(group_1)
        experimental_cohort_data = experiment[3:]
        # Each cohort_data formatted as [cohort,allele_count,inv_freq,inv_read_count,std_read_count,read_total]

        read_counts = []
        read_totals = []
        allele_totals = []
        expected_count_freqs = []
        for i in list(group_1)+list(group_2):
            read_counts += [experimental_cohort_data[i][3]]
            read_totals += [experimental_cohort_data[i][5]]
            allele_totals += [experimental_cohort_data[i][1]]
            expected_count_freqs += [experimental_cohort_data[i][3]/experimental_cohort_data[i][5]]

    
    # calc_chi2_pval = which_calc%2 == 1
    # calc_sel_coeff = (which_calc//2)%2 == 1
    # calc_prod_pval = (which_calc//4)%2 == 1
    if calc_sel_coeff or not calc_pat_from_off:
        assert (experimental_cohort_data[0][0] == 'Parental ♂'), "Need to give a paternal inv frequency to estimate selection coefficients"
        pat_pop_allele_count_from_obs = int(experimental_cohort_data[0][2]*paternal_pop_allele_total)

    if calc_prod_pval:
        assert (N_bootstrap is not None and target_mode in mult_comp_stats), "Frequency change product comparisons are only implemented for bootstrapped two-way multiple comparisons"
    if verbose:
        print("Parameters:\nread_counts\t"+str(read_counts)+"\nread_totals\t"+\
            str(read_totals)+"\nallele_totals\t"+str(allele_totals)+\
            "\npaternal_pop_allele_total\t"+str(paternal_pop_allele_total)+\
            "\nvar_est_N\t\t\t"+str(var_est_N)+\
            "\nmode\t\t\t\t"+str(mode)+\
            "\ncalc chi2-based p-values\t"+str(calc_chi2_pval)+\
            "\ncalc freq∆-prod-based p-values\t"+str(calc_prod_pval)+\
            "\ncalc sel coef\t\t\t"+str(calc_sel_coeff))


    if target_mode in mult_comp_stats:
        allele_counts = []
        for i in np.arange(len(read_counts)):
            allele_counts += [read_counts[i]/read_totals[i]*allele_totals[i]]


        combined_allele_freqs_1 = [sum([allele_counts[i] for i in index_translations[0]])/\
                sum([allele_totals[i] for i in index_translations[0]]),
            sum([allele_counts[i] for i in index_translations[1]])/\
                sum([allele_totals[i] for i in index_translations[1]])]
        combined_allele_freqs_2 = [sum([allele_counts[i] for i in index_translations[2]])/\
                sum([allele_totals[i] for i in index_translations[2]]),
            sum([allele_counts[i] for i in index_translations[3]])/\
                sum([allele_totals[i] for i in index_translations[3]])]

        combined_freq_diff_1 = combined_allele_freqs_1[1]-combined_allele_freqs_1[0]
        combined_freq_diff_2 = combined_allele_freqs_2[1]-combined_allele_freqs_2[0]

        combined_allele_freqs = [combined_allele_freqs_1,combined_allele_freqs_2] 
        combined_freq_diff = [combined_freq_diff_1,combined_freq_diff_2]
    else:
        # Generate the frequency difference between the two groups based on summed allele counts
        allele_counts = []
        for i in np.arange(len(read_counts)):
            allele_counts += [read_counts[i]/read_totals[i]*allele_totals[i]]
        combined_allele_freqs = [sum(allele_counts[:grouping_index])/sum(allele_totals[:grouping_index]),
            sum(allele_counts[grouping_index:])/sum(allele_totals[grouping_index:])]
        combined_freq_diff = combined_allele_freqs[1]-combined_allele_freqs[0]

    if verbose: print("Combined allele frequencies:\t\t"+str(combined_allele_freqs)+\
        "\nCombined allele frequency difference:\t"+str(combined_freq_diff))


    # Store the output calculations in a dictionary
    calculated_data = {}


    # Check whether to calculate a p-value, or skip and instead calculate another statistic
    if not calc_chi2_pval:
        obs_stat = None
        pval = None
        num_GEQ = None
        num_total = None
    else:
        # Use the data to estimate a maximum likelihood paternal allele count/frequency,
        #    and the likelihoods of observed counts in the associated 'neutral' model
        if N_bootstrap is not None and not bs_precalc_dist:
            if calc_pat_from_off:
                ML_pat_allele_count = gen_ML_pat_count(read_counts,read_totals,allele_totals,
                    paternal_pop_allele_total,var_est_N,mode)
            else:
                ML_pat_allele_count = pat_pop_allele_count_from_obs

            ML_est_expected_counts = gen_expected_counts_separately(read_counts,read_totals,allele_totals,
                ML_pat_allele_count,paternal_pop_allele_total,var_est_N,mode)
            ML_est_obs_logProbs = None
            ML_pat_allele_freq = ML_pat_allele_count/paternal_pop_allele_total
        else:
            if calc_pat_from_off:
                # Calculate the log-probability of each possible count in its own distribution,
                #   after estimating the maximum likelihood paternal allele count
                ML_model = gen_ML_pat_count_and_count_prob_dists(read_counts,read_totals,
                    allele_totals,paternal_pop_allele_total,var_est_N,mode)
                (ML_est_obs_logProbs,ML_est_expected_counts,ML_pat_allele_count,ML_pat_allele_freq) = ML_model
            else:
                ML_pat_allele_count = pat_pop_allele_count_from_obs
                ML_model = gen_count_prob_dists(read_counts,read_totals,allele_totals,
                    pat_pop_allele_count_from_obs,paternal_pop_allele_total,var_est_N,mode)
                (ML_est_obs_logProbs,ML_est_expected_counts,ML_pat_allele_freq) = ML_model

        # ML_est_expected_freqs = [ML_est_expected_counts[i]/read_totals[i] for i in np.arange(len(read_totals))]
        
        if N_bootstrap is None:
            num_GEQ = None
            num_total = None

            if tail_mode != "two":
                raise ValueError("Script is only configured to run non-Bootstrap single change chi-square comparison tests as two-tailed tests")

            # Calculate the log-probability that any count combination under the model has a statistic >= to the observed
            if target_mode == "one_change_combined_cohorts":
                # Calculate the log-probability that any count combination under the model has a statistic greater than the observed
                (obs_stat,log_pval) = calc_chi_stat_pval(read_counts,read_totals,ML_est_expected_counts,ML_est_obs_logProbs,allele_totals,
                    count_mode,grouping_index,combine=True)
            elif target_mode == "one_change_lone_cohorts":
                # Same, without combining the cohorts
                (obs_stat,log_pval) = calc_chi_stat_pval(read_counts,read_totals,ML_est_expected_counts,ML_est_obs_logProbs,allele_totals,
                    count_mode,combination_index=None,combine=False)
            elif target_mode in mult_comp_stats:
                raise ValueError("Script is only configured to run non-bootstrap tests for single change chi-square comparison tests")

            # Exponentiate the log-probability
            pval = np.exp(log_pval)

        else:
            # Bootstrap sample to approximate the log-probability that any count combination under the model has a statistic >= to the observed
            if target_mode == "one_change_combined_cohorts":
                assert(count_mode=="allele")
                # Calculate the log-probability that any count combination under the model has a statistic greater than the observed
                (obs_stat,pval,num_GEQ,num_total) = calc_chi_stat_pval_bootstrap(read_counts,read_totals,ML_est_expected_counts,ML_est_obs_logProbs,
                    ML_pat_allele_count,paternal_pop_allele_total,var_est_N,modes,allele_totals,count_mode,grouping_index,combine=True,
                    tail=tail_mode,expected_direction=expected_direction)
            elif target_mode == "one_change_lone_cohorts":
                # Same, without combining the cohorts
                (obs_stat,pval,num_GEQ,num_total) = calc_chi_stat_pval_bootstrap(read_counts,read_totals,ML_est_expected_counts,ML_est_obs_logProbs,
                    ML_pat_allele_count,paternal_pop_allele_total,var_est_N,modes,allele_totals,count_mode,combination_index=None,combine=False,
                    tail=tail_mode,expected_direction=expected_direction)
            elif target_mode=="two_changes_combined_cohorts":
                assert(count_mode=="allele")
                (obs_stat,pval,num_GEQ,num_total) = calc_mult_comb_allele_cohort_chi_stat_pval_bootstrap(read_counts,read_totals,
                    ML_est_expected_counts,ML_est_obs_logProbs,ML_pat_allele_count,paternal_pop_allele_total,var_est_N,
                    allele_totals,index_translations,modes,
                    tail=tail_mode,expected_directions=expected_direction,count_mode=count_mode,combine=True)
            # elif target_mode=="chi2_combined_two_way_ben_cost":
            #     assert(count_mode=="allele")
            #     (obs_stat,pval,num_GEQ,num_total) = calc_mult_comb_allele_cohort_chi_stat_pval_bootstrap(read_counts,read_totals,
            #         ML_est_expected_counts,ML_est_obs_logProbs,ML_pat_allele_count,paternal_pop_allele_total,var_est_N,
            #         allele_totals,index_translations,modes,tail,[1,0],count_mode,combine=True)

        if verbose: print("Raw chi-2 based p-value:\t\t\t\t"+str(pval))
        # calc_key_chi2 = HashableTuple((comparison,chrom,inv,line,target_mode,count_mode,"chi2_pval",tail_mode,expected_direction))
        calc_key_chi2 = (comparison,chrom,inv,line,target_mode,count_mode,"chi2_pval",tail_mode,expected_direction)
        calculated_data[calc_key_chi2] = [pval,obs_stat,combined_allele_freqs,combined_freq_diff,ML_pat_allele_count,
            paternal_pop_allele_total,ML_est_expected_counts,ML_est_obs_logProbs,num_GEQ,num_total]



    # Check whether to calculate a p-value, or skip and instead calculate another statistic
    if not calc_prod_pval:
        obs_stat_prod = None
        pval_prod = None
        num_GEQ_prod = None
        num_total_prod = None
    else:
        # ASSUMES asserted N_bootstrap is not None and target_mode in mult_comp_stats already above
        assert(N_bootstrap is not None and target_mode in mult_comp_stats)


        # if already did precalc for chi2, assume no need to redo
        if not calc_chi2_pval:
            # Use the data to estimate a maximum likelihood paternal allele count/frequency,
            #    and the likelihoods of observed counts in the associated 'neutral' model
            if not bs_precalc_dist:
                if calc_pat_from_off:
                    ML_pat_allele_count = gen_ML_pat_count(read_counts,read_totals,allele_totals,
                        paternal_pop_allele_total,var_est_N,mode)
                else:
                    ML_pat_allele_count = pat_pop_allele_count_from_obs

                ML_est_expected_counts = gen_expected_counts_separately(read_counts,read_totals,allele_totals,
                    ML_pat_allele_count,paternal_pop_allele_total,var_est_N,mode)
                ML_est_obs_logProbs = None
                ML_pat_allele_freq = ML_pat_allele_count/paternal_pop_allele_total
            else:
                if calc_pat_from_off:
                    # Calculate the log-probability of each possible count in its own distribution,
                    #   after estimating the maximum likelihood paternal allele count
                    ML_model = gen_ML_pat_count_and_count_prob_dists(read_counts,read_totals,
                        allele_totals,paternal_pop_allele_total,var_est_N,mode)
                    (ML_est_obs_logProbs,ML_est_expected_counts,ML_pat_allele_count,ML_pat_allele_freq) = ML_model
                else:
                    ML_pat_allele_count = pat_pop_allele_count_from_obs
                    ML_model = gen_count_prob_dists(read_counts,read_totals,allele_totals,
                        pat_pop_allele_count_from_obs,paternal_pop_allele_total,var_est_N,mode)
                    (ML_est_obs_logProbs,ML_est_expected_counts,ML_pat_allele_freq) = ML_model

        
        if target_mode=="two_changes_combined_cohorts":
            assert(count_mode=="allele")
            (obs_stat_prod,pval_prod,num_GEQ_prod,num_total_prod) = \
                calc_mult_comb_allele_cohort_freq_change_product_stat_pval_bootstrap(read_counts,read_totals,
                    ML_est_expected_counts,ML_est_obs_logProbs,ML_pat_allele_count,paternal_pop_allele_total,
                    var_est_N,allele_totals,index_translations,modes,tail_mode,expected_direction,count_mode,combine=True)
        # elif target_mode=="chi2_combined_two_way_ben_cost":
        #     assert(count_mode=="allele")
        #     (obs_stat_prod,pval_prod,num_GEQ_prod,num_total_prod) = \
        #         calc_mult_comb_allele_cohort_freq_change_product_stat_pval_bootstrap(read_counts,read_totals,
        #             ML_est_expected_counts,ML_est_obs_logProbs,ML_pat_allele_count,paternal_pop_allele_total,
        #             var_est_N,allele_totals,index_translations,modes,tail_mode,[1,0],count_mode,combine=True)

        if verbose: print("Raw frequency change product p-value:\t\t\t\t"+str(pval_prod))
        # calc_key_prod = HashableTuple((comparison,chrom,inv,line,target_mode,count_mode,"prod_pval",tail_mode,expected_direction))
        calc_key_prod = (comparison,chrom,inv,line,target_mode,count_mode,"prod_pval",tail_mode,expected_direction)
        calculated_data[calc_key_prod] = [pval_prod,obs_stat_prod,combined_allele_freqs,combined_freq_diff,ML_pat_allele_count,
            paternal_pop_allele_total,ML_est_expected_counts,ML_est_obs_logProbs,num_GEQ_prod,num_total_prod]


    # Check whether to calculate a selection coefficient, or skip and instead calculate another statistic
    if not calc_sel_coeff:
        WF_dir_estimate_s = None
        ABC_estimate_s_and_quant = None
        # confidence_interval = None
    else:
        # Calculate selection coefficient estimates from frequency alone using a simple Wright-Fisher model
        if target_mode in mult_comp_stats:
            WF_dir_estimate_s = calc_mult_WF_estimate_sel_coeff(read_counts,read_totals,allele_totals,
                index_translations,modes,combine=True)
        else:
            WF_dir_estimate_s = calc_sing_WF_estimate_sel_coeff(read_counts,read_totals,allele_totals,
                grouping_index,mode,combine=True)
        # elif target_mode == "one_change_lone_cohorts":
        #     WF_dir_estimate_s = calc_sing_WF_estimate_sel_coeff(read_counts,read_totals,allele_totals,
        #         grouping_index,mode,combine=True)
        if verbose: print("Simple W-F selection coefficient estimate:\t"+str(WF_dir_estimate_s))


        # Calculate selection coefficients via an ABC simulation with uniformly sampled selection coefficients
        if target_mode in mult_comp_stats:
            if target_mode=="two_changes_combined_cohorts":
                assert(count_mode=="allele")
                ABC_estimate_s_and_quant = calc_mult_comp_sel_coeff_ABC(read_counts,read_totals,pat_pop_allele_count_from_obs,
                    paternal_pop_allele_total,var_est_N,allele_totals,index_translations,modes,expected_direction,
                    s_tolerance,s_quantile,WF_dir_estimate_s)
            # elif target_mode=="chi2_combined_two_way_ben_cost":
            #     assert(count_mode=="allele")
            #     ABC_estimate_s_and_quant = calc_mult_comp_sel_coeff_ABC(read_counts,read_totals,pat_pop_allele_count_from_obs,
            #         paternal_pop_allele_total,var_est_N,allele_totals,index_translations,modes,[1,0],
            #         s_tolerance,quantile,WF_dir_estimate_s)
        elif N_bootstrap is not None:
            if target_mode == "one_change_combined_cohorts":
                # Calculate the log-probability that any count combination under the model has a statistic greater than the observed
                ABC_estimate_s_and_quant = calc_single_comp_sel_coeff_ABC(read_counts,read_totals,pat_pop_allele_count_from_obs,
                    paternal_pop_allele_total,var_est_N,modes,s_tolerance,s_quantile,WF_dir_estimate_s,
                    allele_totals,count_mode,grouping_index,combine=True,
                    tail=tail_mode,expected_direction=expected_direction)
            elif target_mode == "one_change_lone_cohorts":
                # Same, without combining the cohorts
                ABC_estimate_s_and_quant = calc_single_comp_sel_coeff_ABC(read_counts,read_totals,pat_pop_allele_count_from_obs,
                    paternal_pop_allele_total,var_est_N,modes,s_tolerance,s_quantile,WF_dir_estimate_s,
                    allele_totals,count_mode,combination_index=None,combine=False,
                    tail=tail_mode,expected_direction=expected_direction)
        else:
            ABC_estimate_s_and_quant = [None,None,None]

        if verbose: print("Approximate Bayesian Computation (ABC) selection coefficient estimate:\n"+str(ABC_estimate_s_and_quant))
        # calc_key_sel = HashableTuple((comparison,chrom,inv,line,target_mode,count_mode,"sel_coeff",s_tolerance,expected_direction,quantile))
        calc_key_sel = (comparison,chrom,inv,line,target_mode,count_mode,"sel_coeff",s_tolerance,expected_direction,s_quantile)
        calculated_data[calc_key_sel] = [WF_dir_estimate_s,ABC_estimate_s_and_quant,combined_allele_freqs,
            combined_freq_diff,paternal_pop_allele_total]

    # comparison_data = [comparison,chrom,inv,line,pval,obs_stat,combined_allele_freqs,combined_freq_diff,
    #     ML_pat_allele_count,paternal_pop_allele_total,ML_est_expected_counts,ML_est_obs_logProbs,num_GEQ,num_total]
    # comparison_data = [comparison,chrom,inv,line,pval,obs_stat,combined_allele_freqs,combined_freq_diff,
    #     ML_pat_allele_count,paternal_pop_allele_total,ML_est_expected_counts,ML_est_obs_logProbs,num_GEQ,num_total,
    #     WF_dir_estimate_s,ABC_estimate_s_and_quant,pval_prod,obs_stat_prod,num_GEQ_prod,num_total_prod]
    return calculated_data



# For each experimental comparison, perform maximum likelihood comparison
def calc_all_pval(experiment_data,pat_census_dict,var_est_N,cust_lines,cust_invs,comparison_tuples,calc_pat_from_off):
        # comparison_tuples=[(,True,False,True,"allele","one_change_combined_cohorts","two",None)]):
        # calc_chi2_pval=True,calc_prod_pval=False,calc_sel_coeff=True,
        # count_mode="allele",target_mode="one_change_combined_cohorts",tail_mode="one"):

    comparison_data_sets = {}

    if multiprocessing:
        # Prepare the multiprocessing pool
        from multiprocessing import set_start_method
        set_start_method("spawn")
        import multiprocessing as mp

        if num_proc > 0:
            pool_size = num_proc
        else:
            pool_size = int(mp.cpu_count()*prop_proc)
        pool = mp.Pool(pool_size)

        # Collect and write the p-values apply_async()
        async_objects = []
        print('Performing p-value calculations in asynchronous parallel:\nMultiprocessing pool size:\t\t'+\
            str(pool_size)+' pocessor cores')
        # Experimental data formatted as [line,chrom,inv,pat_data,emb_data,eaf_data,laf_data,eam_data,lam_data]
        # Each *_data formatted as [cohort,allele_count,inv_freq,inv_read_count,std_read_count,read_total]
        for experiment in experiment_data:
            if (cust_lines is None or str(experiment[0]) in cust_lines) and (cust_invs is None or str(experiment[2]) in cust_invs):
                if verbose:
                    print("\n\nBeginning analysis of In("+str(experiment[1])+")"+str(experiment[2])+\
                        " differences in replicate Zi"+str(experiment[0])+"N")
                for (comparison,calc_chi2_pval,calc_prod_pval,calc_sel_coeff,
                        count_mode,target_mode,tail_mode,expected_direction) in comparison_tuples:
                    paternal_pop_allele_total = 2*pat_census_dict[experiment[0]]
                    async_objects += [pool.apply_async(comparison_pval,args=(comparison,experiment,
                        paternal_pop_allele_total,var_est_N,N_bootstrap,bs_chunk,bs_precalc_dist,verbose,
                        s_tolerance,s_quantile,calc_chi2_pval,calc_prod_pval,calc_sel_coeff,calc_pat_from_off,
                        count_mode,target_mode,tail_mode,expected_direction))]

        # Wait for the results with .get() -- could use .wait instead, and a for loop, so as to avoid a None list output
        # comparison_data_sets = [result.get() for result in async_objects]
        comparison_data_sets = {k: v for result in async_objects for k, v in result.get().items()}

        # Close the pool and return the comparison outputs
        pool.close()
        pool.join()
    else:
        # Experimental data formatted as [line,chrom,inv,pat_data,emb_data,eaf_data,laf_data,eam_data,lam_data]
        # Each *_data formatted as [cohort,allele_count,inv_freq,inv_read_count,std_read_count,read_total]
        for experiment in experiment_data:
            if (cust_lines is None or str(experiment[0]) in cust_lines) and (cust_invs is None or str(experiment[2]) in cust_invs):
                if verbose:
                    print("\n\nBeginning analysis of In("+str(experiment[1])+")"+str(experiment[2])+\
                        " differences in replicate Zi"+str(experiment[0])+"N")
                # for comparison in comparisons:
                for (comparison,calc_chi2_pval,calc_prod_pval,calc_sel_coeff,
                        count_mode,target_mode,tail_mode,expected_direction) in comparison_tuples:
                    paternal_pop_allele_total = 2*pat_census_dict[experiment[0]]
                    comparison_pval_data = comparison_pval(comparison,experiment,paternal_pop_allele_total,
                        var_est_N,N_bootstrap,bs_chunk,bs_precalc_dist,verbose,s_tolerance,s_quantile,
                        calc_chi2_pval,calc_prod_pval,calc_sel_coeff,calc_pat_from_off,
                        count_mode,target_mode,tail_mode,expected_direction)
                    # comparison_data_sets += [comparison_pval_data]
                    comparison_data_sets.update(comparison_pval_data)

    return comparison_data_sets







# # For each experimental comparison, perform maximum likelihood comparison
# def calc_all_sel_coeff(experiment_data,pat_census_dict,var_est_N,
#         comparisons,cust_lines,cust_invs,
#         count_mode="allele",target_mode="one_change_combined_cohorts",tail_mode="one"):

#     comparison_data_sets = []

#     if multiprocessing:
#         # Prepare the multiprocessing pool
#         from multiprocessing import set_start_method
#         set_start_method("spawn")
#         import multiprocessing as mp

#         if num_proc > 0:
#             pool_size = num_proc
#         else:
#             pool_size = int(mp.cpu_count()*prop_proc)
#         pool = mp.Pool(pool_size)

#         # Collect and write the selection coefficient calcs from apply_async()
#         async_objects = []
#         print('Performing selection coefficient calculations in asynchronous parallel:')
#         # Experimental data formatted as [line,chrom,inv,pat_data,emb_data,eaf_data,laf_data,eam_data,lam_data]
#         # Each *_data formatted as [cohort,allele_count,inv_freq,inv_read_count,std_read_count,read_total]
#         for experiment in experiment_data:
#             if (cust_lines is None or str(experiment[0]) in cust_lines) and (cust_invs is None or str(experiment[2]) in cust_invs):
#                 if verbose:
#                     print("\n\nBeginning analysis of In("+str(experiment[1])+")"+str(experiment[2])+\
#                         " differences in replicate Zi"+str(experiment[0])+"N")
#                 for comparison in comparisons:
#                     paternal_pop_allele_total = 2*pat_census_dict[experiment[0]]
#                     async_objects += [pool.apply_async(comparison_sel_coeff,args=(comparison,experiment,
#                         paternal_pop_allele_total,var_est_N,N_bootstrap,bs_chunk,bs_precalc_dist,verbose,
#                         count_mode,target_mode,tail_mode))]

#         # Wait for the results with .get() -- could use .wait instead, and a for loop, so as to avoid a None list output
#         comparison_sel_coeff_sets = [result.get() for result in async_objects]

#         # Close the pool and return the comparison outputs
#         pool.close()
#         pool.join()
#     else:
#         # Experimental data formatted as [line,chrom,inv,pat_data,emb_data,eaf_data,laf_data,eam_data,lam_data]
#         # Each *_data formatted as [cohort,allele_count,inv_freq,inv_read_count,std_read_count,read_total]
#         for experiment in experiment_data:
#             if (cust_lines is None or str(experiment[0]) in cust_lines) and (cust_invs is None or str(experiment[2]) in cust_invs):
#                 if verbose:
#                     print("\n\nBeginning analysis of In("+str(experiment[1])+")"+str(experiment[2])+\
#                         " selection coefficients in replicate Zi"+str(experiment[0])+"N")
#                 for comparison in comparisons:
#                     paternal_pop_allele_total = 2*pat_census_dict[experiment[0]]
#                     comparison_s_data = comparison_sel_coeff(comparison,experiment,paternal_pop_allele_total,
#                         var_est_N,N_bootstrap,bs_chunk,bs_precalc_dist,verbose,count_mode,target_mode,tail_mode)
#                     comparison_sel_coeff_sets += [comparison_s_data]

#     return comparison_sel_coeff_sets










# Helper functions to factor out reused code
def __get_comparison_str(comparison):
        if type(comparison[0]) is tuple:
            comparison_str = '-'.join(comparison[0])+':'+'-'.join(comparison[1])
        elif type(comparison[0]) is not str:
            raise ValueError("Unexpected form of comparison, should only be (str,str) or "+\
                "((str,str),(str,str)); given comparison \""+string(comparison))+"\""
        else:
            comparison_str = '-'.join(comparison)
        return comparison_str

def get_sep_pval_csv_filename(out_file_prefix,comparison,line,inversion,calc_name,target_mode,count_mode,tail_mode,expected_direction):
    print((out_file_prefix,comparison,line,inversion,calc_name,target_mode,count_mode,tail_mode,expected_direction))
    comparison_str = __get_comparison_str(comparison)
    out_file_name = out_file_prefix+'_'+comparison_str+'_'+line+'_'+inversion+'_'+calc_name+\
            '_'+target_mode+'_'+count_mode+'_'+tail_mode
    if tail_mode == "one":
        out_file_name+='_'+str(expected_direction)+'.csv'
    else:
        out_file_name += '.csv'
    return out_file_name

def __gen_stat_name(target_mode,tail_mode,count_mode,calc_name):
    stat_full_name = ''
    if (calc_name,target_mode) not in stat_full_name_dict:
        raise ValueError("script only configured for managing certain statistical tests:\n"+\
            string(stat_full_name_dict.keys())+"\nInput statistical mode \""+str(target_mode)+"\" does not match\n")
    else:
        stat_full_name += stat_full_name_dict[(calc_name,target_mode)]+count_mode+"_counts"
    return stat_full_name

def __header_1(stat_full_name,isMultComp=False):
    if stat_full_name[:17] == 'Frequency_product':
        return "Comparison,Cohort_1-1,Cohort_1-2,Cohort_2-1,Cohort_2-2,Line,Chrom,Inv,Pval,"+stat_full_name+','
    elif isMultComp:
        return "Comparison,Cohort_1-1,Cohort_1-2,Cohort_2-1,Cohort_2-2,Line,Chrom,Inv,Pval,"+stat_full_name+'_1,'+\
                stat_full_name+'_2,'
    else:
        return "Comparison,Cohort_1,Cohort_2,Line,Chrom,Inv,Pval,"+stat_full_name+','

def __header_2(stat_full_name,isMultComp=False):
    header = ''
    if N_bootstrap is not None:
        header = "Bootstrap_Successes,Bootstrap_Replicates,"
    if isMultComp:
        return header+"Allele_Freq_1-1,Allele_Freq_1-2,Allele_Freq_2-1,Allele_Freq_2-2,"+\
            "Allele_Freq_Diff_1,Allele_Freq_Diff_2,ML_Paternal_Allele_Count,Paternal_Census_Size\n"
    else:
        return header+"Allele_Freq_1,Allele_Freq_2,Allele_Freq_Diff,ML_Paternal_Allele_Count,Paternal_Census_Size\n"

def __out_line_1(comp_key,comp_data,isMultComp=False):
    if comp_key[6] == 'prod_pval':
        return '-'.join(comp_key[0][0])+':'+'-'.join(comp_key[0][1])+','+group_full_name_dict[comp_key[0][0][0]]+','+\
            group_full_name_dict[comp_key[0][0][1]]+','+group_full_name_dict[comp_key[0][1][0]]+','+group_full_name_dict[comp_key[0][1][1]]+','+\
            str(comp_key[3])+','+str(comp_key[1])+','+str(comp_key[2])+','+str(comp_data[0])+','+str(comp_data[1])+','
    elif isMultComp:
        return '-'.join(comp_key[0][0])+':'+'-'.join(comp_key[0][1])+','+group_full_name_dict[comp_key[0][0][0]]+','+\
            group_full_name_dict[comp_key[0][0][1]]+','+group_full_name_dict[comp_key[0][1][0]]+','+group_full_name_dict[comp_key[0][1][1]]+','+\
            str(comp_key[3])+','+str(comp_key[1])+','+str(comp_key[2])+','+str(comp_data[0])+','+str(comp_data[1][0])+','+str(comp_data[1][1])+','
    else:
        return '-'.join(comp_key[0])+','+group_full_name_dict[comp_key[0][0]]+','+group_full_name_dict[comp_key[0][1]]+','+\
            str(comp_key[3])+','+str(comp_key[1])+','+str(comp_key[2])+','+str(comp_data[0])+','+str(comp_data[1])+','

def __out_line_2(comp_data,isMultComp=False):
    line = ''
    if N_bootstrap is not None:
        assert(N_bootstrap==comp_data[9]),"Expecting the bootstrap sample size "+\
            "to be equal to the --bootstrap parameter"
        line += str(comp_data[8])+','+str(comp_data[9])+','
    if isMultComp:
        return line+str(comp_data[2][0][0])+','+str(comp_data[2][0][1])+','+str(comp_data[2][1][0])+','+str(comp_data[2][1][1])+','+\
                    str(comp_data[3][0])+','+str(comp_data[3][1])+','+str(comp_data[4])+','+str(comp_data[5])+'\n'
    else:
        return line+str(comp_data[2][0])+','+str(comp_data[2][1])+','+str(comp_data[3])+','+str(comp_data[4])+','+str(comp_data[5])+'\n'


def __header_sel_coeff(stat_full_name,isMultComp=False):
    header = ''
    if isMultComp:
        header += "Comparison,Cohort_1-1,Cohort_1-2,Cohort_2-1,Cohort_2-2,Line,Chrom,Inv,WF_est_"+stat_full_name+'_1,WF_est_'+\
                stat_full_name+'_2,ABC_est_1,ABC_est_1_upper_quantile,ABC_est_1_lower_quantile,'+\
                'ABC_est_2,ABC_est_2_upper_quantile,ABC_est_2_lower_quantile,quantile,s_tolerance,'
    else:
        header += "Comparison,Cohort_1,Cohort_2,Line,Chrom,Inv,WF_est_"+stat_full_name+\
                ',ABC_est,ABC_est_upper_quantile,ABC_est_lower_quantile,quantile,s_tolerance,'
    # if N_bootstrap is not None:
    #     header += "Bootstrap_Successes,Bootstrap_Replicates,"
    if isMultComp:
        return header+"Allele_Freq_1-1,Allele_Freq_1-2,Allele_Freq_2-1,Allele_Freq_2-2,"+\
            "Allele_Freq_Diff_1,Allele_Freq_Diff_2,Paternal_Census_Size\n"
    else:
        return header+"Allele_Freq_1,Allele_Freq_2,Allele_Freq_Diff,Paternal_Census_Size\n"


def __out_line_sel_coeff(comp_key,comp_data,isMultComp=False):
    line = ''
    if isMultComp:
        line += '-'.join(comp_key[0][0])+':'+'-'.join(comp_key[0][1])+','+group_full_name_dict[comp_key[0][0][0]]+','+\
            group_full_name_dict[comp_key[0][0][1]]+','+group_full_name_dict[comp_key[0][1][0]]+','+group_full_name_dict[comp_key[0][1][1]]+','+\
            str(comp_key[3])+','+str(comp_key[1])+','+str(comp_key[2])+','+str(comp_data[0][0])+','+str(comp_data[0][1])+','+\
            str(comp_data[1][0])+','+str(comp_data[1][1])+','+str(comp_data[1][2])+','+\
            str(comp_data[1][3])+','+str(comp_data[1][4])+','+str(comp_data[1][5])+','+\
            str(comp_key[9])+','+str(comp_key[7])+','
    else:
        line += '-'.join(comp_key[0])+','+group_full_name_dict[comp_key[0][0]]+','+group_full_name_dict[comp_key[0][1]]+','+\
            str(comp_key[3])+','+str(comp_key[1])+','+str(comp_key[2])+','+str(comp_data[0])+','+\
            str(comp_data[1][0])+','+str(comp_data[1][1])+','+str(comp_data[1][2])+','+\
            str(comp_key[9])+','+str(comp_key[7])+','
    # if N_bootstrap is not None:
    #     assert(N_bootstrap==comp_data[13]),"Expecting the bootstrap sample size "+\
    #         "to be equal to the --bootstrap parameter"
    #     line += str(comp_data[8])+','+str(comp_data[9])+','
    if isMultComp:
        return line+str(comp_data[2][0][0])+','+str(comp_data[2][0][1])+','+str(comp_data[2][1][0])+','+str(comp_data[2][1][1])+','+\
                    str(comp_data[3][0])+','+str(comp_data[3][1])+','+str(comp_data[4])+'\n'
    else:
        return line+str(comp_data[2][0])+','+str(comp_data[2][1])+','+str(comp_data[3])+','+str(comp_data[4])+'\n'

# A standardized file output of nohup parallel execution to read from when combining the results
def write_sel_coeff_results(sel_comps,sel_comp_data,out_file_prefix):
    out_file_name = out_file_prefix+'-'+'sel_coeff'+'.csv'
    if verbose:
        print("Writing selection coefficient calculations to "+str(out_file_name))

    comp = sel_comps[0]
    calc_name = comp[6]
    assert(calc_name == 'sel_coeff'), "writing selection coefficients should only be called for selection coefficient calculation results"
    target_mode = comp[4]
    count_mode = comp[5]
    tail_mode = comp[7]

    stat_full_name = __gen_stat_name(target_mode,tail_mode,count_mode,'sel_coeff')
    isMultComp = target_mode in mult_comp_stats
    header = __header_sel_coeff(stat_full_name,isMultComp)

    with open(out_file_name, 'w') as out_file:
        out_file.write(header)

        for i in np.arange(len(sel_comps)):
            comp = sel_comps[i]
            calc_name = comp[6]
            assert(calc_name == 'sel_coeff'), "writing selection coefficients should only be called for selection coefficient calculation results"
            target_mode = comp[4]
            count_mode = comp[5]
            tail_mode = comp[7]
            """
            for 'sel_coeff':
                (comparison,chrom,inv,line,target_mode,count_mode,"sel_coeff",s_tolerance,expected_direction,quantile) 
            """
            comp_data = sel_comp_data[i]
            """
            for 'sel_coeff':
                [WF_dir_estimate_s,ABC_estimate_s_and_quant,combined_allele_freqs,
                combined_freq_diff,paternal_pop_allele_total]
            """

            if verbose:
                print("Writing selection calculation output for "+str(comp[0]))
                print(comp)
                print(comp_data)
            line = __out_line_sel_coeff(comp,comp_data,isMultComp)
            out_file.write(line)
    return

# A standardized file output of nohup parallel execution to read from when combining the results
def write_sep_pval_results(comparison_data,out_file_prefix):
        # target_mode,tail_mode,count_mode):
    for comp in comparison_data:
        calc_name = comp[6]
        target_mode = comp[4]
        count_mode = comp[5]
        tail_mode = comp[7]
        expected_direction = comp[8]
        """
        comp arranged by calc_name:
        for 'chi2_pval':
            (comparison,chrom,inv,line,target_mode,count_mode,"chi2_pval",tail_mode,expected_direction)
        for 'prod_pval':
            (comparison,chrom,inv,line,target_mode,count_mode,"prod_pval",tail_mode,expected_direction)
        for 'sel_coeff':
            (comparison,chrom,inv,line,target_mode,count_mode,"sel_coeff",s_tolerance,expected_direction,quantile) 
        """
        comp_data = comparison_data[comp]
        """
        comp_data arranged by calc_name:
        for 'chi2_pval':
            [pval,obs_stat,combined_allele_freqs,combined_freq_diff,ML_pat_allele_count,
            paternal_pop_allele_total,ML_est_expected_counts,ML_est_obs_logProbs,num_GEQ,num_total]
        for 'prod_pval':
            [pval_prod,obs_stat_prod,combined_allele_freqs,combined_freq_diff,ML_pat_allele_count,
            paternal_pop_allele_total,ML_est_expected_counts,ML_est_obs_logProbs,num_GEQ_prod,num_total_prod]
        for 'sel_coeff':
            [WF_dir_estimate_s,ABC_estimate_s_and_quant,combined_allele_freqs,
            combined_freq_diff,paternal_pop_allele_total]
        """
        out_file_name = get_sep_pval_csv_filename(out_file_prefix,comp[0],comp[3],comp[2],
                                    calc_name,target_mode,count_mode,tail_mode,expected_direction)

        if verbose:
            print("Writing output for "+str(comp[0])+" to "+str(out_file_name))
            print(comp)
            print(comp_data)

        stat_full_name = __gen_stat_name(target_mode,tail_mode,count_mode,calc_name)
        isMultComp = target_mode in mult_comp_stats

        if calc_name == "chi2_pval" or calc_name == "prod_pval":

            header = __header_1(stat_full_name,isMultComp)+__header_2(stat_full_name,isMultComp)

            with open(out_file_name, 'w') as out_file:
                out_file.write(header)

                # Expecting comparison data as [comparison,chrom,inv,line,p_val,obs_stat,combined_allele_freqs,
                #  combined_freq_diff,ML_est_pat_count,paternal_pop_allele_total,ML_est_expected_counts,ML_est_obs_logProbs]
                line = __out_line_1(comp,comp_data,isMultComp)+__out_line_2(comp_data,isMultComp)

                out_file.write(line)

        elif calc_name == "sel_coeff":
            header = __header_sel_coeff(stat_full_name,isMultComp)
            with open(out_file_name, 'w') as out_file:
                out_file.write(header)
                line = __out_line_sel_coeff(comp,comp_data,isMultComp)
                out_file.write(line)





        # if target_mode in mult_comp_stats:
        #     header = "Comparison,Cohort_1-1,Cohort_1-2,Cohort_2-1,Cohort_2-2,Line,Chrom,Inv,Pval,"+stat_full_name+'_1,'+\
        #         stat_full_name+'_2,'+"Allele_Freq_1-1,Allele_Freq_1-2,Allele_Freq_2-1,Allele_Freq_2-2,"+\
        #         "Allele_Freq_Diff_1,Allele_Freq_Diff_2,ML_Paternal_Allele_Count,Paternal_Census_Size\n"
            
        #     with open(out_file_name, 'w') as out_file:
        #         out_file.write(header)

        #         # Expecting comparison data as [comparison,chrom,inv,line,p_val,obs_stat,combined_allele_freqs,
        #         #  combined_freq_diff,ML_est_pat_count,paternal_pop_allele_total,ML_est_expected_counts,ML_est_obs_logProbs]
        #         line = '('+'-'.join(comparison[0])+','+'-'.join(comparison[1])+')'+','+group_full_name_dict[comp[0][0][0]]+','+\
        #             group_full_name_dict[comp[0][0][1]]+','+group_full_name_dict[comp[0][1][0]]+','+group_full_name_dict[comp[0][1][1]]+','+\
        #             str(comp[3])+','+str(comp[1])+','+str(comp[2])+','+str(comp[4])+','+str(comp[5][0])+','+str(comp[5][1])+','+\
        #             str(comp[6][0][0])+','+str(comp[6][0][1])+','+str(comp[6][1][0])+','+str(comp[6][1][1])+','+\
        #             str(comp[7][0])+','+str(comp[7][1])+','+str(comp[8])+','+str(comp[9])+'\n'

        #         out_file.write(line)
        # else:
        #     header = "Comparison,Cohort_1,Cohort_2,Line,Chrom,Inv,Pval,"+stat_full_name+','+\
        #         "Allele_Freq_1,Allele_Freq_2,Allele_Freq_Diff,ML_Paternal_Allele_Count,Paternal_Census_Size\n"
            
        #     with open(out_file_name, 'w') as out_file:
        #         out_file.write(header)

        #         # Expecting comparison data as [comparison,chrom,inv,line,p_val,obs_stat,combined_allele_freqs,
        #         #  combined_freq_diff,ML_est_pat_count,paternal_pop_allele_total,ML_est_expected_counts,ML_est_obs_logProbs]
        #         line = '-'.join(comp[0])+','+group_full_name_dict[comp[0][0]]+','+group_full_name_dict[comp[0][1]]+','+\
        #             str(comp[3])+','+str(comp[1])+','+str(comp[2])+','+str(comp[4])+','+str(comp[5])+','+\
        #             str(comp[6][0])+','+str(comp[6][1])+','+str(comp[7])+','+str(comp[8])+','+str(comp[9])+'\n'

        #         out_file.write(line)
    return


# For generating a log file unique to each nohup process
def get_nohup_log_filename(out_file_prefix,comparison_str,line,inv,calc_chi2_pval,calc_prod_pval,calc_sel_coeff,
    target_mode,count_mode,tail_mode):
    log_file_name = out_file_prefix+'_nohup_'+comparison_str+'_'+line+'_'+inv+'_'+\
            'c'+str(int(calc_chi2_pval))+'p'+str(int(calc_prod_pval))+'s'+str(int(calc_sel_coeff))+'_'+\
            target_mode+'_'+count_mode+'_'+tail_mode+'.log'
    return log_file_name

# For running each single comparison as a new nohup script and generating an associated output csv,
#   has no memory or cpu management
def run_comps_separately_with_nohup(called_read_count_file,library_metadata_file_path,output_prefix,output_directory,
    global_alpha,var_est_N,comparison_tuples,cust_lines,cust_invs,calc_pat_from_off,experiment_data,out_file_prefix):
    # global_alpha,var_est_N,count_mode,target_mode,tail_mode,comparisons,cust_lines,cust_invs,experiment_data,out_file_prefix):

    from subprocess import Popen

    call_args = []
    call_args += ["nohup","python","chi2_pval_test.py","--library_data",library_metadata_file_path,"--count_file",
        called_read_count_file,"--od",output_directory]
    if output_prefix != '':
        call_args += ["--out_prefix",output_prefix]
    call_args += ["--varN",str(var_est_N),
        "--alpha",str(global_alpha),"--s_tol",str(s_tolerance),"--s_quantile",str(s_quantile),
        "--bs_chunk",str(bs_chunk),"--write_sep_no_mt","--nomp"]
    
    if verbose:
        call_args += ["--verbose"]

    if N_bootstrap is not None:
        call_args += ["--bootstrap",str(N_bootstrap)]

    if bs_precalc_dist:
        call_args += ["--bs_precalc_dist"]

    if not calc_pat_from_off:
        call_args += ["--source_pat_from_obs"]


    print('\nSpawning nohup python processes to calculate and write p-value data for individual comparisons:')
    # Experimental data formatted as [line,chrom,inv,pat_data,emb_data,eaf_data,laf_data,eam_data,lam_data]
    # Each *_data within is formatted as [cohort,allele_count,inv_freq,inv_read_count,std_read_count,read_total]
    from os.path import exists
    for experiment in experiment_data:
        line = str(experiment[0])
        inv = str(experiment[2])

        if (cust_lines is None or line in cust_lines) and (cust_invs is None or inv in cust_invs):
            if verbose:
                print("\n\nBeginning separate proccesses for analysis of In("+str(experiment[1])+")"+str(inv)+\
                    " differences in replicate Zi"+str(line)+"N")
            for (comparison,calc_chi2_pval,calc_prod_pval,calc_sel_coeff,
                    count_mode,target_mode,tail_mode,expected_direction) in comparison_tuples:
                # Get the comparison string to work with
                comparison_str = __get_comparison_str(comparison)

                # Check for the presence of the expected comparison output file:
                log_file_name = get_nohup_log_filename(out_file_prefix,comparison_str,line,
                                                    inv,calc_chi2_pval,calc_prod_pval,calc_sel_coeff,
                                                    target_mode,count_mode,tail_mode)
                # log_file_name = get_sep_pval_csv_filename(out_file_prefix,comparison,line,inv,'nohup')[:-3]+'log'
                if exists(log_file_name):
                    print("Warning - nohup_p-value.log runtime file already present for "+str(comparison)+' '+\
                        str(line)+' '+str(inv)+', SKIPPING')
                else:
                    # Add further run-specific comparisons to the arguments
                    curr_comp_args = list(call_args)
                    curr_comp_args += ["--count",count_mode,"--stat_target",target_mode,
                                        "--stat_tail",tail_mode,
                                        "--comp_cohorts","custom","--cust_comps",comparison_str,
                                        "--cust_lines",line,"--cust_invs",inv
                                        ]

                    if calc_chi2_pval:
                        curr_comp_args += ["--calc_chi2"]
                    if calc_prod_pval:
                        curr_comp_args += ["--calc_prod"]
                    if calc_sel_coeff:
                        curr_comp_args += ["--calc_s"]

                    # Add the nohup output location to the arguments
                    nohup_out_file_name = log_file_name
                    nohup_shell_tail = [">>",nohup_out_file_name,"2>&1","&"]

                    # Call the single-comparison python run,
                    #   but not if filling out the 
                    if verbose: print(' '.join(curr_comp_args+nohup_shell_tail))
                    Popen(curr_comp_args,
                        stdout=open(nohup_out_file_name, 'a'),
                        stderr=open(nohup_out_file_name, 'a'),
                        shell=False)
    return

# Helper function to read an output line and interpret it into comparison data
def __read_sep_pval_line(data_line,comparison,
        target_mode,count_mode,calc_mode,tail_mode,expected_direction,
        isMultComp):
    """
    calc_mode == "prod_pval", 2way, bootstrap:
        Comparison,Cohort_1-1,Cohort_1-2,Cohort_2-1,Cohort_2-2,Line,Chrom,Inv,Pval,"+stat_full_name+",
        Bootstrap_Successes,Bootstrap_Replicates,
        Allele_Freq_1-1,Allele_Freq_1-2,Allele_Freq_2-1,Allele_Freq_2-2,
        Allele_Freq_Diff_1,Allele_Freq_Diff_2,ML_Paternal_Allele_Count,Paternal_Census_Size

    calc_mode == "chi2_pval", 2way, bootstrap:
        [0]Comparison,Cohort_1-1,Cohort_1-2,Cohort_2-1,Cohort_2-2,Line,Chrom,Inv,Pval,"+stat_full_name_1+","+stat_full_name_2+",
        [11]Bootstrap_Successes,Bootstrap_Replicates,
        Allele_Freq_1-1,Allele_Freq_1-2,Allele_Freq_2-1,Allele_Freq_2-2,
        Allele_Freq_Diff_1,Allele_Freq_Diff_2,ML_Paternal_Allele_Count,Paternal_Census_Size
        

    calc_mode == "chi2_pval", 1way:
        [0]Comparison,Cohort_1,Cohort_2,Line,Chrom,Inv,Pval,"+stat_full_name+','
        Bootstrap_Successes,Bootstrap_Replicates,
        Allele_Freq_1,Allele_Freq_2,Allele_Freq_Diff,ML_Paternal_Allele_Count,Paternal_Census_Size

    calc_mode == "sel_coeff", 2way:
        Comparison,Cohort_1-1,Cohort_1-2,Cohort_2-1,Cohort_2-2,Line,Chrom,Inv,
        WF_est_"+stat_full_name+"_1,WF_est_"+stat_full_name+"_2,
        ABC_est_1,ABC_est_1_upper_quantile,ABC_est_1_lower_quantile,
        ABC_est_2,ABC_est_2_upper_quantile,ABC_est_2_lower_quantile,
        quantile,s_tolerance,Allele_Freq_1-1,Allele_Freq_1-2,Allele_Freq_2-1,Allele_Freq_2-2,
        Allele_Freq_Diff_1,Allele_Freq_Diff_2,Paternal_Census_Size

    calc_mode == "sel_coeff", 1way:
        Comparison,Cohort_1,Cohort_2,Line,Chrom,Inv,
        WF_est_Selection_coefficients_from_combined_allele_counts,ABC_est,ABC_est_upper_quantile,ABC_est_lower_quantile,
        quantile,s_tolerance,Allele_Freq_1,Allele_Freq_2,Allele_Freq_Diff,Paternal_Census_Size

    """
    data = data_line.strip().split(',')
    if calc_mode == "sel_coeff":
        raise ValueError
        if isMultComp:
            assert(all([data[i+1]==group_full_name_dict[(comparison[0]+comparison[1])[i]] for i in np.arange(4)]))
            chrom = str(data[6])
            inv = str(data[7])
            line = str(data[5])
            WF_est = [float(data[8]),float(data[9])]
            ABC_estimate_s_and_quant = [float(i) for i in np.arange(10,16)]
            quantile = float(data[16])
            s_tolerance = float(data[17])
            combined_allele_freqs = [[float(data[18]),float(data[19])],[float(data[20]),float(data[21])]]
            combined_freq_diff = [float(data[22]),float(data[23])]
            paternal_pop_allele_total = int(data[24])

        else:
            assert(all([data[i+1]==group_full_name_dict[comparison[i]] for i in np.arange(2)]))
            chrom = str(data[4])
            inv = str(data[5])
            line = str(data[3])
            WF_est = float(data[6])
            if N_bootstrap is not None:
                ABC_estimate_s_and_quant = [float(i) for i in np.arange(7,10)]
            else:
                ABC_estimate_s_and_quant = [None,None,None]
            quantile = float(data[10])
            s_tolerance = float(data[11])
            combined_allele_freqs = [float(data[12]),float(data[13])]
            combined_freq_diff = float(data[14])
            paternal_pop_allele_total = int(data[15])


        # Expecting: (comparison,chrom,inv,line,target_mode,count_mode,"sel_coeff",s_tolerance,expected_direction,s_quantile)
        comp_key = (comparison,chrom,inv,line,target_mode,count_mode,calc_mode,s_tolerance,expected_direction,s_quantile)
        # [WF_dir_estimate_s,ABC_estimate_s_and_quant,combined_allele_freqs,combined_freq_diff,paternal_pop_allele_total]
        comp_data = [WF_est,ABC_estimate_s_and_quant,combined_allele_freqs,combined_freq_diff,paternal_pop_allele_total]

    else:
        if isMultComp:
            assert(all([data[i+1]==group_full_name_dict[(comparison[0]+comparison[1])[i]] for i in np.arange(4)]))
            chrom = str(data[6])
            inv = str(data[7])
            line = str(data[5])
            p_val = float(data[8])
            if calc_mode == "chi2_pval":
                obs_stat = [float(data[9]),float(data[10])]
                rem_data = data[11:]
            elif calc_mode == "prod_pval":
                obs_stat = float(data[9])
                rem_data = data[10:]
            if N_bootstrap is not None:
                num_GEQ = int(rem_data[0])
                num_total = int(rem_data[1])
                assert(N_bootstrap == num_total),"Script requires specifying a "+\
                    "--bootstrap sample size equal to that used for data read in"
                rem_data = rem_data[2:]
            else:
                num_GEQ = None
                num_total = None
            combined_allele_freqs = [[float(rem_data[0]),float(rem_data[1])],[float(rem_data[2]),float(rem_data[3])]]
            combined_freq_diff = [float(rem_data[4]),float(rem_data[5])]
            ML_est_pat_count = float(rem_data[6])
            paternal_pop_allele_total = int(rem_data[7])
            ML_est_expected_counts = None
            ML_est_obs_logProbs = None
        else:
            assert(all([data[i+1]==group_full_name_dict[comparison[i]] for i in np.arange(2)]))
            chrom = str(data[4])
            inv = str(data[5])
            line = str(data[3])
            p_val = float(data[6])
            obs_stat = float(data[7])
            if N_bootstrap is not None:
                num_GEQ = int(data[8])
                num_total = int(data[9])
                assert(N_bootstrap == num_total),"Script requires specifying a "+\
                    "--bootstrap sample size equal to that used for data read in"
                rem_data = data[10:]
            else:
                num_GEQ = None
                num_total = None
            combined_allele_freqs = [float(rem_data[0]),float(rem_data[1])]
            combined_freq_diff = float(rem_data[2])
            ML_est_pat_count = float(rem_data[3])
            paternal_pop_allele_total = int(rem_data[4])
            ML_est_expected_counts = None
            ML_est_obs_logProbs = None

        # comparison_pval_data = [comparison,chrom,inv,line,p_val,obs_stat,combined_allele_freqs,
        #             combined_freq_diff,ML_est_pat_count,paternal_pop_allele_total,
        #             ML_est_expected_counts,ML_est_obs_logProbs,num_GEQ,num_total]

        # Expecting: (comparison,chrom,inv,line,target_mode,count_mode,calc_mode,tail_mode,expected_direction)
        comp_key = (comparison,chrom,inv,line,target_mode,count_mode,calc_mode,tail_mode,expected_direction)
        # [pval,obs_stat,combined_allele_freqs,combined_freq_diff,ML_pat_allele_count,
        #       paternal_pop_allele_total,ML_est_expected_counts,ML_est_obs_logProbs,num_GEQ,num_total]
        comp_data = [p_val,obs_stat,combined_allele_freqs,combined_freq_diff,ML_est_pat_count,paternal_pop_allele_total,
            ML_est_expected_counts,ML_est_obs_logProbs,num_GEQ,num_total]

    comparison_pval_data = {comp_key:comp_data}

    return comparison_pval_data

# Helper factoring out code to read a full single-entry file
def __read_sep_pval_file(out_file_name,comparison,target_mode,count_mode,
            calc_mode,tail_mode,expected_direction,isMultComp):
    comparison_pval_data = []
    with open(out_file_name) as pval_file:
        lines = pval_file.readlines()
        assert len(lines) == 2, "collecting separate results expects only one result line per file"
        isMultComp = target_mode in mult_comp_stats
        comparison_pval_data = __read_sep_pval_line(lines[1],comparison,target_mode,count_mode,
            calc_mode,tail_mode,expected_direction,isMultComp)
    return comparison_pval_data



def __file_check_and_read_sep_pval_file(out_file_prefix,comparison,line,inv,calc_mode,target_mode,
            count_mode,tail_mode,expected_direction,isMultComp):
    from os.path import exists
    if verbose:
        print("\n\nBeginning reading pval data of In("+inv_chrom_lookup[inv]+")"+str(inv)+\
            " differences in replicate Zi"+str(line)+"N")
    # Check for the presence of the expected comparison output file,
    out_file_name = get_sep_pval_csv_filename(out_file_prefix,comparison,line,inv,calc_mode,
                                target_mode,count_mode,tail_mode,expected_direction)
    if not exists(out_file_name): # p-value.csv?
        print("Warning - .csv data file NOT present for "+str(comparison)+' '+\
            str(line)+' '+str(inv)+' '+str(calc_mode)+', SKIPPING file:\n'+\
            out_file_name+"\n")
    else:
        # out_file_name = get_sep_pval_csv_filename(out_file_prefix,comparison,experiment[0],experiment[2],calc_mode)
        # out_file_name = get_sep_pval_csv_filename(out_file_prefix,comparison,line,inv,calc_mode)
        if verbose: print(out_file_name)
        if calc_mode == "chi2_pval" or calc_mode == "prod_pval":
            return __read_sep_pval_file(out_file_name,comparison,target_mode,count_mode,
                    calc_mode,tail_mode,expected_direction,isMultComp)
    return

# For each experimental comparison, read in the maximum likelihood pval from an expected csv file
def collect_pval_results(experiment_data,out_file_prefix,
        comparisons,cust_lines,cust_invs,comparison_tuples):

    comparison_data_sets = {}

    print('Collecting p-values from pre-calculated individual comparison files:')
    # Experimental data formatted as [line,chrom,inv,pat_data,emb_data,eaf_data,laf_data,eam_data,lam_data]
    # Each *_data formatted as [cohort,allele_count,inv_freq,inv_read_count,std_read_count,read_total]
    for experiment in experiment_data:
        if (cust_lines is None or str(experiment[0]) in cust_lines) and (cust_invs is None or str(experiment[2]) in cust_invs):
            line = experiment[0]
            chrom = experiment[1]
            inv = experiment[2]
            for (comparison,calc_chi2_pval,calc_prod_pval,calc_sel_coeff,
                    count_mode,target_mode,tail_mode,expected_direction) in comparison_tuples:
                isMultComp = target_mode in mult_comp_stats
                if calc_chi2_pval:
                    calc_mode = "chi2_pval"
                    comparison_pval_data = __file_check_and_read_sep_pval_file(out_file_prefix,comparison,line,inv,calc_mode,
                        target_mode,count_mode,tail_mode,expected_direction,isMultComp)
                    if comparison_pval_data is not None:
                        comparison_data_sets.update(comparison_pval_data)
                if calc_prod_pval:
                    calc_mode = "prod_pval"
                    comparison_pval_data = __file_check_and_read_sep_pval_file(out_file_prefix,comparison,line,inv,calc_mode,
                        target_mode,count_mode,tail_mode,expected_direction,isMultComp)
                    if comparison_pval_data is not None:
                        comparison_data_sets.update(comparison_pval_data)
    return comparison_data_sets

# Helper function for running Fisher's combined p-value calculation on a list of p-values
def fishers_combined(pvals):
    chi_square_score = np.sum(-2*np.log(pvals))
    degree_f = 2*len(pvals)
    comb_pval = sp.stats.distributions.chi2.sf(chi_square_score, degree_f)
    return comb_pval



# For calculating a single p-value for the entire set of replicate lines,
#   using Fisher's method:
#      comparing to X^2 distribution of 2k degrees of freedom
#      -2*sum(ln(p_i)) for k independent tests of p-values p_i
#   Consider: the tests across lines are likely to be independent, and are the intended target of combination,
#   but p-values for paternal-to-embryo are likely to be positively correlated with those of embryo-to-any adult offspring,
#   and correlation between inversions is plausible but unclear
#   *This combination step should come before multiple test correction*
def combine_pvals_by_line_and_dir(comps,comp_data,comparisons,cust_lines,cust_invs,calc_mode="chi2_pval"):
    stat_index = {
        "chi2":4,
        "prod":16
    }
    comb_pval_data = []
    inversions = list(dict.fromkeys(list(inv_name_dict.values())))
    if cust_invs is not None:
        inversions = cust_invs
    target_mode = comps[0][4]
    tail_mode = comps[0][7]
    if target_mode in mult_comp_stats:
        directions = ((0,1),(1,0))
    else:
        directions = (0,1)

    for comparison in comparisons:
        for inv in inversions:
            if tail_mode == "one":
                for direction in directions:
                    pvals = []
                    chrom = ''
                    for i in np.arange(len(comps)):
                        if ((comps[i][0] == comparison) and (comps[i][2] == inv) and \
                                (cust_lines is None or comps[i][3] in cust_lines) and (comps[i][8] == direction)):
                            chrom = comps[i][1]
                            pvals += [comp_data[i][0]]
                    if len(pvals) > 0:
                        comb_pval = fishers_combined(pvals)
                        comb_pval_data += [(comparison,chrom,inv,direction,comb_pval)]
            else:
                pvals = []
                chrom = ''
                for i in np.arange(len(comps)):
                    if ((comps[i][0] == comparison) and (comps[i][2] == inv) and \
                            (cust_lines is None or comps[i][3] in cust_lines)):
                        chrom = comps[i][1]
                        pvals += [comp_data[i][0]]
                if len(pvals) > 0:
                    comb_pval = fishers_combined(pvals)
                    comb_pval_data += [(comparison,chrom,inv,None,comb_pval)]
    return comb_pval_data

    # comp: (comparison,chrom,inv,line,target_mode,count_mode,calc_mode,tail_mode,expected_direction)
    # comp_data: [pval,obs_stat,combined_allele_freqs,combined_freq_diff,ML_pat_allele_count,
    #       paternal_pop_allele_total,ML_est_expected_counts,ML_est_obs_logProbs,num_GEQ,num_total]


# def combine_pvals_by_line(comps,comp_data,comparisons,cust_lines,cust_invs,calc_mode="chi2_pval"):
#     stat_index = {
#         "chi2":4,
#         "prod":16
#     }
#     comb_pval_data = []
#     inversions = list(dict.fromkeys(list(inv_name_dict.values())))
#     if cust_invs is not None:
#         inversions = cust_invs
#     # print(inversions)
#     for comp in comparisons:
#         for inv in inversions:
#             # indexes = []
#             pvals = []
#             chrom = ''
#             for i in np.arange(len(comps)):
#                 if ((comps[i][0] == comp) and (comps[i][2] == inv) and (cust_lines is None or comps[i][3] in cust_lines)):
#                     # indexes += [i]
#                     chrom = comps[i][1]
#                     # pvals += [comp_data[i][stat_index[calc_mode]]]
#                     pvals += [comp_data[i][0]]
#             # num_tests = len(indexes)
#             num_tests = len(pvals)
#             if num_tests > 0:
#                 chi_square_score = np.sum(-2*np.log(pvals))
#                 degree_f = 2*num_tests
#                 comb_pval = sp.stats.distributions.chi2.sf(chi_square_score, degree_f)
#                 comb_pval_data += [(comp,chrom,inv,comb_pval)]
#     return comb_pval_data

# Further helper functions for the Fisher's combined data
def __comb_header(stat_full_name,isMultComp=False):
    if isMultComp:
        return "Comparison,Cohort_1-1,Cohort_1-2,Cohort_2-1,Cohort_2-2,Chrom,Inv,Fishers_Pval_on_"+\
            stat_full_name+'\n'
    else:
        return "Comparison,Cohort_1,Cohort_2,Chrom,Inv,Fishers_Pval_on_"+stat_full_name+'\n'

def __comb_out_line(comp_data,isMultComp=False):
    if isMultComp:
        return '-'.join(comp_data[0][0])+':'+'-'.join(comp_data[0][1])+','+group_full_name_dict[comp_data[0][0][0]]+','+\
            group_full_name_dict[comp_data[0][0][1]]+','+group_full_name_dict[comp_data[0][1][0]]+','+group_full_name_dict[comp_data[0][1][1]]+','+\
            str(comp_data[1])+','+str(comp_data[2])+','+str(comp_data[3])+'\n'
    else:
        return '-'.join(comp_data[0])+','+group_full_name_dict[comp_data[0][0]]+','+group_full_name_dict[comp_data[0][1]]+','+\
            str(comp_data[1])+','+str(comp_data[2])+','+str(comp_data[3])+'\n'

# For writing line-combined significance data
def write_combined_pval_data(comb_pval_data,out_file_name,target_mode,tail_mode,count_mode):

    if verbose: print("Writing combined p-value output to "+str(out_file_name))
    
    stat_full_name = __gen_stat_name(target_mode,tail_mode,count_mode)

    isMultComp = target_mode in mult_comp_stats
    header = __comb_header(stat_full_name,isMultComp)
    if tail_mode == "one":
        header = header[:-1]+',Direction\n'
    if N_bootstrap is not None:
        header = header[:-1]+',Bootstrap_Sample_Size\n'

    with open(out_file_name, 'w') as out_file:
        out_file.write(header)
        
        # Expecting combined p-value comparison data as [comparison,chrom,inv,combined_p_val]
        for comp_data in comb_pval_data:
            line = __comb_out_line(comp_data,isMultComp)
            if tail_mode == "one":
                line = line[:-1]+','+'-'.join(direction)+'\n'
            if N_bootstrap is not None:
                line = line[:-1]+','+str(N_bootstrap)+'\n'
            out_file.write(line)
    return



# Run a Bonferroni test correction for the Family-Wise Error Rate 
def bonferroni(alpha,pvals):
    if verbose: print("Running Bonferroni Correction")
    num_tests = len(pvals)
    adjustment_factor = 1/num_tests
    adjusted_pvals = np.divide(pvals,adjustment_factor)
    null_rejected = adjusted_pvals < alpha
    if verbose:
        print(pvals)
        print(adjusted_pvals)
        print(adjustment_factor)
        print(null_rejected)
    return (adjusted_pvals,null_rejected)

# Run a Dunn-Sidak test correction for the Family-Wise Error Rate 
def dunn_sidak(alpha,pvals):
    if verbose: print("Running Dunn-Šidák Correction")
    num_tests = len(pvals)
    sidak_alpha = 1-(1-alpha)**(1/num_tests)
    null_rejected = np.less(pvals,sidak_alpha)
    adjustment_factor = sidak_alpha/alpha
    adjusted_pvals = np.divide(pvals,adjustment_factor)
    if verbose:
        print(pvals)
        print(adjusted_pvals)
        print(adjustment_factor)
        print(null_rejected)
    return (adjusted_pvals,null_rejected)

# Run a  Benjamini-Hochberg test correction for the False Discovery Rate 
def benjamini_hochberg(alpha,pvals):
    if verbose: print("Running Benjamini-Hochberg Correction")
    num_tests = len(pvals)
    null_rejected = np.zeros(num_tests,dtype=bool)
    ascending_order = np.argsort(pvals)
    k = 1
    curr_index = ascending_order[k-1]
    curr_pval = pvals[curr_index]
    while curr_pval <= k/num_tests*alpha and k < num_tests:
        null_rejected[curr_index] = True
        k += 1
        curr_index = ascending_order[k-1]
        curr_pval = pvals[curr_index]
    if k == num_tests:
        if curr_pval <= k/num_tests*alpha:
            null_rejected[curr_index] = True
            k += 1
    if k == 1: k+= 1 # Not clear what to do here, adjusted p-values my not really be appropriate for a step-up
    ben_hoch_alpha = (k-1)/num_tests*alpha
    adjustment_factor = ben_hoch_alpha/alpha
    adjusted_pvals = np.divide(pvals,adjustment_factor)
    if verbose:
        print(pvals)
        print(adjusted_pvals)
        print(adjustment_factor)
        print(null_rejected)
    return (adjusted_pvals,null_rejected)

# Calculates the i'th harmonic number, or the sum 1 + 1/2 + 1/3 + ... 1/i
def harmonic_num(i):
    return sum([1/j for j in np.arange(1,i+1)])

# Run a Benjamini-Yekutieli test correction for the False Discovery Rate 
def benjamini_yekutieli(alpha,pvals):
    if verbose: print("Running Benjamini-Yekutieli Correction")
    num_tests = len(pvals)
    null_rejected = np.zeros(num_tests,dtype=bool)
    ascending_order = np.argsort(pvals)
    k = 1
    curr_index = ascending_order[k-1]
    curr_pval = pvals[curr_index]
    while curr_pval <= (k/(num_tests*harmonic_num(num_tests)))*alpha and k < num_tests:
        null_rejected[curr_index] = True
        k += 1
        curr_index = ascending_order[k-1]
        curr_pval = pvals[curr_index]
    if k == num_tests:
        if curr_pval <= (k/(num_tests*harmonic_num(num_tests)))*alpha:
            null_rejected[curr_index] = True
            k += 1
    if k == 1: k+= 1 # Not clear what to do here, adjusted p-values my not really be appropriate for a step-up
    ben_yekut_alpha = ((k-1)/(num_tests*harmonic_num(num_tests)))*alpha
    adjustment_factor = ben_yekut_alpha/alpha
    adjusted_pvals = np.divide(pvals,adjustment_factor)
    if verbose:
        print(pvals)
        print(adjusted_pvals)
        print(adjustment_factor)
        print(null_rejected)
    return (adjusted_pvals,null_rejected)

# Run a number of multiple testing correction algorithms at once
def run_mult_mult_test_corrections(global_alpha,pvals,test_list):
    test_name_dict = {
        "Benjamini-Yekutieli"  : partial(benjamini_yekutieli,global_alpha),
        "Benjamini-Hochberg"   : partial(benjamini_hochberg,global_alpha),
        "Dunn-Sidak"           : partial(dunn_sidak,global_alpha),
        "Dunn-Šidák"           : partial(dunn_sidak,global_alpha),
        "Šidák"                : partial(dunn_sidak,global_alpha),
        "Sidak"                : partial(dunn_sidak,global_alpha),
        "Bonferroni"           : partial(bonferroni,global_alpha)
    }
    if verbose:
        print("\n\nBeginning calculation of significance under multiple testing corrections "+str(test_list))
    multiple_test_corrections = []
    print(test_list)
    for test_name in test_list:
        print(test_name)
        (adjusted_pvals,null_rejected) = test_name_dict[test_name](pvals)
        multiple_test_corrections += [(test_name,adjusted_pvals,null_rejected)]
    return multiple_test_corrections


def write_comparisons_with_mult_test_data(comps,comp_data,mult_test_corrections,
        target_mode,tail_mode,count_mode,calc_mode,out_file_name):

    if verbose: print("Writing line-separate output with multiple test corrections to "+str(out_file_name))
    
    # stat_full_name = ''
    # if target_mode == "one_change_combined_cohorts":
    #      stat_full_name += "Chi_Squared_on_Combined_"+count_mode+"_Counts"
    # elif target_mode == "one_change_lone_cohorts":
    #      stat_full_name += "Chi_Squared_on_Separate_"+count_mode+"_Counts"

    stat_full_name = __gen_stat_name(target_mode,tail_mode,count_mode,calc_mode)

    isMultComp = target_mode in mult_comp_stats
    header = __header_1(stat_full_name,isMultComp)

    if tail_mode == "one":
        header +='Direction,'

    if mult_test_corrections is not None:
        for (test_name,adjusted_pvals,null_rejected) in mult_test_corrections:
            header += test_name+"_Sig,"+test_name+"_Adjusted_Pval,"
    
    header += __header_2(stat_full_name,isMultComp)


    with open(out_file_name, 'w') as out_file:
        out_file.write(header)

        # Expecting comparison data as [comparison,chrom,inv,line,p_val,obs_stat,combined_allele_freqs,
        #  combined_freq_diff,ML_est_pat_count,paternal_pop_allele_total,
        #  ML_est_expected_counts,ML_est_obs_logProbs,num_GEQ]
        for i in np.arange(len(comp_data)):
            comp = comps[i]
            comp_dat = comp_data[i]

            line = __out_line_1(comp,comp_dat,isMultComp)

            if tail_mode == "one":
                direction = comp[8]
                if isMultComp:
                    assert len(direction) == 2
                    line += str(direction[0])+'-'+str(direction[1])+','
                else:
                    line += str(direction)+','
            
            if mult_test_corrections is not None:
                for (test_name,adjusted_pvals,null_rejected) in mult_test_corrections:
                    line += str(null_rejected[i])+','+str(adjusted_pvals[i])+','

            line += __out_line_2(comp_dat,isMultComp)

            out_file.write(line)
    return



# For writing line-combined comparison data with multiple-test-corrected p-values
def write_combined_comparisons_with_mult_test_data(comparison_data,mult_test_corrections,
        target_mode,tail_mode,count_mode,calc_name,out_file_name):

    if verbose: print("Writing combined p-value output with multiple test corrections to "+str(out_file_name))
    
    stat_full_name = __gen_stat_name(target_mode,tail_mode,count_mode,calc_name)

    isMultComp = target_mode in mult_comp_stats

    header = __comb_header(stat_full_name,isMultComp)[:-1]+","
    if tail_mode == "one":
        header +='Direction,'
    if mult_test_corrections is not None:
        for (test_name,adjusted_pvals,null_rejected) in mult_test_corrections:
            header += test_name+"_Sig,"+test_name+"_Adjusted_Pval,"
    if N_bootstrap is not None:
        header += "Bootstrap_Sample_Size,"
    header = header[:-1]+"\n"

    with open(out_file_name, 'w') as out_file:
        out_file.write(header)
        # Expecting combined p-value comparison data as [comparison,chrom,inv,combined_p_val]
        for i in np.arange(len(comparison_data)):
            (comparison,chrom,inv,direction,comb_pval) = comparison_data[i]
            line = __comb_out_line((comparison,chrom,inv,comb_pval),isMultComp)[:-1]+","
            if tail_mode == "one":
                if isMultComp:
                    assert len(direction) == 2
                    line += str(direction[0])+'-'+str(direction[1])+','
                else:
                    line += str(direction)+','
            if mult_test_corrections is not None:
                for (test_name,adjusted_pvals,null_rejected) in mult_test_corrections:
                    line += str(null_rejected[i])+','+str(adjusted_pvals[i])+','
            if N_bootstrap is not None:
                line += str(N_bootstrap)+','
            line = line[:-1]+'\n'
            out_file.write(line)
    return


def run_mt_corr_fishers_comb_write_output(global_alpha,sep_pvals,mult_tests,
        comps,comp_data,target_mode,tail_mode,count_mode,calc_mode,exp_dir,
        comp_mode,comparisons,cust_lines,cust_invs,out_file_prefix):
    print(tail_mode)
    print(comp_mode)
    sep_mult_test_corrections = run_mult_mult_test_corrections(global_alpha,sep_pvals,mult_tests)

    # Write the first multiple test corrected output to file
    print(comp_mode)
    out_file_name = out_file_prefix+'_'+comp_mode+'_sep-multi-test-corrections.csv'
    write_comparisons_with_mult_test_data(comps,comp_data,sep_mult_test_corrections,
        target_mode,tail_mode,count_mode,calc_mode,out_file_name)

    # Combine the p-values by line using Fisher's combined p-value method
    line_comb_comp_data = combine_pvals_by_line_and_dir(comps,comp_data,comparisons,cust_lines,cust_invs)

    # Calculate the multiple testing corrections/significances with combined lines
    comb_pvals = [comp[4] for comp in line_comb_comp_data]
    comb_mult_test_corrections = run_mult_mult_test_corrections(global_alpha,comb_pvals,mult_tests)

    # Write the line-combined p-value data and multiple-test coreections to file
    comb_out_file_name = out_file_prefix+'_'+comp_mode+'_line-combined-pval-with-mt.csv'
    write_combined_comparisons_with_mult_test_data(line_comb_comp_data,comb_mult_test_corrections,
        target_mode,tail_mode,count_mode,calc_mode,comb_out_file_name)
    return



# # Helper to cut computation off early if mis-configured
# def __test_comparison_and_stat_compatibility(sample_comp,target_mode,tail_mode):
#     if target_mode not in mult_comp_stats:
#         if type(sample_comp[0]) is tuple:
#             raise ValueError("Script is only configured to handle comparison tuples like \""+str(sample_comp)+\
#                 "\" when in a simultaneous chi-square statistical mode on multiple comparisons. Instead "+\
#                 "the target_mode is \""+str(target_mode)+"\"")
#     else:
#         if type(sample_comp) is not tuple:
#             raise ValueError("When given a simultaneous chi-square statistical mode for running multiple "+\
#                 "comparisons, comparisons are expected to be provided in the form "+\
#                 "\"(cohort1.1-cohort1.2,cohort2.1-cohort2.2)\"; was given \""+str(sample_comp)+"\"")

# Helper to cut computation off early if mis-configured
def __test_comparison_and_stat_compatibility(sample_comp,count_mode,target_mode,tail_mode):
    if target_mode not in mult_comp_stats:
        if type(sample_comp[0]) is tuple:
            print("Script is only configured to handle comparison tuples like \""+str(sample_comp)+\
                "\" when in a simultaneous chi-square statistical mode on multiple comparisons. Instead "+\
                "the target_mode is \""+str(target_mode)+"\"")
            return False 
    else:
        if count_mode != 'allele':
            print("Script is only configured to handle multiple comparison targets like \""+str(sample_comp)+\
                "\" in \""+str(target_mode)+"\" when in allele count mode. Instead given \""+str(count_mode)+"\"")
            return False
        if type(sample_comp) is not tuple:
            print("When given a simultaneous chi-square statistical mode for running multiple "+\
                "comparisons, comparisons are expected to be provided in the form "+\
                "\"(cohort1.1-cohort1.2,cohort2.1-cohort2.2)\"; was given \""+str(sample_comp)+\
                "\" for comparison target "+str(target_mode))
            return False
    return True

# Helper to organize tests into compatible compute jobs, potentially save compute time, remove incompatible ones
def __compose_comparison_tuples(comparisons,calc_chi2_pval,calc_prod_pval,calc_sel_coeff,count_mode,target_modes,tail_modes):
    comparison_tuples = []
    for comparison in comparisons:
        # draft_tuple = [comparison]
        for target_mode in target_modes:
            if target_mode in mult_comp_stats:
                directions = ((0,1),(1,0))
            else:
                directions = (0,1)
            for tail_mode in tail_modes:
                # if tail_mode == "one" and target_mode not in mult_comp_stats:
                #     print("WARNING: removing a comparison job requesting one-tailed tests for single-change comparisons; "+\
                #         "script isn't configured to test such")
                # else:
                if __test_comparison_and_stat_compatibility(comparison,count_mode,target_mode,tail_mode):
                    # if tail_mode == "two" and :
                    #     if target_mode in mult_comp_stats:
                    #         comparison_tuples += [(comparison,calc_chi2_pval,calc_prod_pval,
                    #             calc_sel_coeff,count_mode,target_mode,tail_mode,None)]
                    #     else:
                    #         comparison_tuples += [(comparison,calc_chi2_pval,False,
                    #             calc_sel_coeff,count_mode,target_mode,tail_mode,None)]
                    # else:
                    if target_mode in mult_comp_stats:
                        for direction in directions:
                            comparison_tuples += [(comparison,calc_chi2_pval,calc_prod_pval,
                                calc_sel_coeff,count_mode,target_mode,tail_mode,direction)]
                    else:
                        if tail_mode == "two":
                            comparison_tuples += [(comparison,calc_chi2_pval,False,
                                calc_sel_coeff,count_mode,target_mode,tail_mode,None)]
                        else:
                            for direction in directions:
                                comparison_tuples += [(comparison,calc_chi2_pval,False,
                                    calc_sel_coeff,count_mode,target_mode,tail_mode,direction)]
    return comparison_tuples




        # if tail_mode != "two":
        #     print("Script is only configured to run single change chi-square comparison tests as two-tailed tests")
        #     return False

        #     if tail_mode != "two":
        #         raise ValueError("Script is only configured to run single change chi-square comparison tests as two-tailed tests")

        #     elif target_mode in mult_comp_stats:
        #         raise ValueError("Script is only configured to run non-bootstrap tests for single change chi-square comparison tests")

        #     if calc_sel_coeff:
        #         assert (experimental_cohort_data[0][0] == 'Parental ♂'), "Need to give a paternal inv frequency to estimate selection coefficients"
        #         pat_pop_allele_count_from_obs = experimental_cohort_data[0][2]*paternal_pop_allele_total

        #     if calc_prod_pval:
        #         assert (N_bootstrap is not None and target_mode in mult_comp_stats), "Frequency change product comparisons are only implemented for bootstrapped two-way multiple comparisons"




# Main function, for running the pipeline when the script is called
def main():
    # Parse arguments
    args = parse_args()

    # Assign input arguments to variables
    called_read_count_file = args.count_file
    library_metadata_file_path = args.library_data
    output_prefix = args.out_prefix
    output_directory = args.od

    global global_alpha
    global var_est_N
    global_alpha = args.alpha
    var_est_N = args.varN

    # global count_mode
    # global target_mode
    # global tail_mode
    count_mode = args.count
    target_mode = args.stat_target
    target_modes = []
    if target_mode == 'one_and_two_changes_combined_cohorts':
        target_modes = ['two_changes_combined_cohorts','one_change_combined_cohorts']
    else:
        target_modes = [target_mode]
    tail_mode = args.stat_tail
    tail_modes = []
    if tail_mode == 'one_and_two':
        tail_modes = ['two','one']
    else:
        tail_modes = [tail_mode]

    calc_pat_from_off = not args.source_pat_from_obs

    mult_tests = args.cust_mt
    if isinstance(mult_tests,str):
        mult_tests = [mult_tests]

    comp_mode = args.comp_cohorts
    # print(comp_mode)
    cust_comps = args.cust_comps
    # print(cust_comps)
    cust_lines = args.cust_lines
    if cust_lines is not None:
        cust_lines = [line_name_dict[line] for line in cust_lines]
    # print(cust_lines)
    cust_invs = args.cust_invs
    if cust_invs is not None:
        cust_invs = [inv_name_dict[inv] for inv in cust_invs]
    # print(cust_invs)

    comparisons = []
    if comp_mode == "custom":
        for comp in cust_comps:
            mult_comp_parts = comp.split(':')
            if target_mode in mult_comp_stats:
                if len(mult_comp_parts) < 2:
                    raise ValueError("When given a simultaneous chi-square statistical mode for running multiple "+\
                        "comparisons, comparisons are expected to be provided in the form "+\
                        "\"cohort1.1-cohort1.2:cohort2.1-cohort2.2\"; script was given cohort \""+str(comp)+"\"")
                # Currently only configured to take two, so takes te first two split entries
                comparison = (check_cust_comp(mult_comp_parts[0]),check_cust_comp(mult_comp_parts[1]))
            else:
                if len(mult_comp_parts) >= 2:
                    raise ValueError("Script is only configured to handle comparison tuples like "+\
                        "\"cohort1.1-cohort1.2:cohort2.1-cohort2.2\" or the given comparison \""+str(comp)+\
                        "\" when in a simultaneous chi-square statistical mode taking multiple comparisons. "+\
                        "Instead the target_mode is \""+str(target_mode)+"\"")
                comparison = check_cust_comp(comp)
            # comparison = tuple(comp.split('-'))
            # if len(comparison) != 2:
            #     raise ValueError("Script is only configured to generate chi square statistics that compare two groupings "+\
            #         "of experimental cohorts at a time, "+string(len(comparison))+\
            #         " found in custom comparison \""+comp+"\"\n")
            # if comparison[0] not in cohort_group_dict or comparison[1] not in cohort_group_dict:
            #     raise ValueError("Custom comparisons only configured to compare certain named cohorts or groupings:\n"+\
            #         string(cohort_group_dict.keys())+"\nInput custom comparison \""+comp+"\" does not match\n")
            comparisons += [comparison]
    elif comp_mode not in comp_mode_dict.keys():
        raise ValueError("comp_mode given as "+str(comp_mode)+", mode expected to be only one of "+str(comp_mode_dict.keys()))
    else:
        comparisons = comp_mode_dict[comp_mode]
    # __test_comparison_and_stat_compatibility(comparisons[0],target_mode,tail_mode)

    global N_bootstrap
    N_bootstrap = args.bootstrap
    assert(target_mode not in mult_comp_stats or N_bootstrap is not None),"Script is only configured "+\
        "to run non-bootstrap tests for single chi-square comparison tests; instead given target_mode \""+str(target_mode)+\
        "\" and N_bootstrap = "+str(N_bootstrap)
    global bs_chunk
    bs_chunk = args.bs_chunk
    if N_bootstrap is not None and bs_chunk > N_bootstrap:
        bs_chunk = N_bootstrap
    global bs_precalc_dist
    bs_precalc_dist = args.bs_precalc_dist

    global s_tolerance
    s_tolerance = args.s_tol
    global s_quantile
    s_quantile = args.s_quantile

    # no_pval_calc = args.no_pval_calc
    # calc_s = args.calc_s
    # calc_prod_pval = args.calc_prod_pval
    # assert(calc_s or not no_pval_calc), "Nothing to calculate: neither selection coefficient or p-values requested"
    s_spec = args.calc_s
    chi2_spec = args.calc_chi2
    prod_spec = args.calc_prod
    calc_all = not s_spec and not chi2_spec and not prod_spec
    # global calc_chi2_pval
    # global calc_prod_pval
    # global calc_sel_coeff
    calc_chi2_pval = calc_all or chi2_spec
    calc_prod_pval = calc_all or prod_spec
    calc_sel_coeff = calc_all or s_spec
    assert(calc_chi2_pval or calc_prod_pval or calc_sel_coeff), "Nothing to calculate: neither selection coefficient or p-values requested"
    # global which_calc
    # which_calc = int(not no_pval_calc)+2*int(calc_s)+4*int(calc_prod_pval and not no_pval_calc)
    # assert(which_calc > 0)


    comparison_tuples = __compose_comparison_tuples(comparisons,calc_chi2_pval,calc_prod_pval,calc_sel_coeff,
                                                    count_mode,target_modes,tail_modes)
    print("\nPrinting intended comparisons as tuples of (comparison,calc_chi2_pval,calc_prod_pval,calc_sel_coeff,"+\
            "count_mode,target_mode,tail_mode,direction):")
    print(comparison_tuples)
    print()

    global verbose
    verbose = args.verbose

    global multiprocessing
    multiprocessing = args.nomp
    global num_proc
    global prop_proc
    prop_proc = args.prop_proc
    num_proc = args.num_proc

    write_sep_no_mt = args.write_sep_no_mt
    mt_from_sep = args.mt_from_sep
    run_sep_nohup = args.run_sep_nohup

    mt_across_comparisons = args.mt_comps_comb

    assert sum([write_sep_no_mt,mt_from_sep,run_sep_nohup]) < 2, "Undefined behavior when "+\
        "more than one of the execution mode flags \"--write_sep_no_mt\", \"--mt_from_sep\", "+\
        "\"--run_sep_nohup\" are set at the same time"

    # Ensure the output directory exists and the string is sanitized
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    if output_directory[-1] != '/':
        output_directory += '/'


    out_file_prefix = ''

    # # Prepare the output file prefix (located here to avoid code duplication)
    # out_file_prefix = target_mode+'_'+count_mode+'_'+tail_mode+'_tailed'
    # if N_bootstrap is not None:
    #     out_file_prefix += '_bootstrap'
    # if calc_chi2_pval:
    #     out_file_prefix += '_chi2Stat'
    # if calc_prod_pval:
    #     out_file_prefix += '_prodStat'
    # if calc_sel_coeff:
    #     out_file_prefix += '_selEst'

    if output_prefix == '':
        out_file_prefix = output_directory+out_file_prefix
    else:
        if out_file_prefix == '':
            out_file_prefix = output_directory+output_prefix
        else:
            out_file_prefix = output_directory+output_prefix+'_'+out_file_prefix

    library_metadata = parse_library_metadata_to_dict(library_metadata_file_path)
    count_data = parse_count_data_to_dict(called_read_count_file)
    experiment_data = reshape_read_count_data(count_data,library_metadata)

    # If spawning separate processes, do so and exit
    if run_sep_nohup:
        run_comps_separately_with_nohup(called_read_count_file,library_metadata_file_path,output_prefix,output_directory,
            global_alpha,var_est_N,comparison_tuples,cust_lines,cust_invs,calc_pat_from_off,
            # global_alpha,var_est_N,count_mode,target_mode,tail_mode,comparisons,cust_lines,cust_invs,
            experiment_data,out_file_prefix)
        return

    # CALCULATING (or reading in) CHI SQUARED STATISTIC P-VALUES
    if mt_from_sep:
        comparison_data_sets = collect_pval_results(experiment_data,out_file_prefix,
            comparisons,cust_lines,cust_invs,comparison_tuples)
        print(comparison_data_sets)
    else:
        # Calculate all the relevant pvals
        comparison_data_sets = calc_all_pval(experiment_data,pat_census_dict,
            var_est_N,cust_lines,cust_invs,comparison_tuples,calc_pat_from_off)
            #var_est_N,comparisons,cust_lines,cust_invs,count_mode,target_mode,tail_mode)
    if calc_sel_coeff:
        sel_comps = []
        for comp in comparison_data_sets:
            if comp[6] == 'sel_coeff':
                sel_comps += [comp]
        if len(sel_comps) == 0:
            print("Set to calculate selection coefficients with --calc_s, but no comparisons with 'sel_coeff' data in comparison_data_sets")
        else:
            sel_comp_data = []
            for comp in sel_comps:
                sel_comp_data += [comparison_data_sets.pop(comp)]
            write_sel_coeff_results(sel_comps,sel_comp_data,out_file_prefix)
    if write_sep_no_mt:
        # Write the output to file individually
        write_sep_pval_results(comparison_data_sets,out_file_prefix)
        return


    # mult_tests = ("Bonferroni","Dunn-Sidak","Benjamini-Hochberg","Benjamini-Yekutieli")

    # Calculate the multiple testing corrections/significances with separate lines

    # print(comparison_data_sets)
    # calc_key_sel = HashableTuple((comparison,chrom,inv,line,target_mode,count_mode,"sel_coeff",s_tolerance,expected_direction,quantile))
    #     calculated_data[calc_key_sel] = [WF_dir_estimate_s,ABC_estimate_s_and_quant,combined_allele_freqs,
    #         combined_freq_diff,paternal_pop_allele_total]



    if calc_chi2_pval or calc_prod_pval:
        testTypes = set()
        comparisonTypes = set()
        for comp in comparison_data_sets:
                testTypes.add(tuple(comp[4:8]))
                comparisonTypes.add(comp[0])
        for testType in testTypes:
            if mt_across_comparisons:
                this_mt_set_comps = []
                this_mt_set_comp_data = []
                this_mt_set_pvals = []
                for comp in comparison_data_sets:
                    if tuple(comp[4:8]) == testType:
                        this_mt_set_comps += [comp]
                        this_mt_set_comp_data += [comparison_data_sets[comp]]
                        this_mt_set_pvals += [comparison_data_sets[comp][0]]
                if len(this_mt_set_comps) > 0:
                    # Prepare the multiple testing output file prefix
                    mt_prefix = '_'.join(testType)
                    if N_bootstrap is not None:
                        mt_prefix += '_bootstrap'
                    if output_prefix == '':
                        mt_prefix = output_directory+mt_prefix
                    else:
                        mt_prefix = output_directory+output_prefix+'_'+mt_prefix

                    run_mt_corr_fishers_comb_write_output(global_alpha,this_mt_set_pvals,mult_tests,
                        this_mt_set_comps,this_mt_set_comp_data,testType[0],testType[3],testType[1],testType[2],None,
                        comp_mode,comparisons,cust_lines,cust_invs,mt_prefix)
                else:
                    print("CHECK - no tests of class "+str(testType))
            else:
                for compType in comparisonTypes:
                    this_mt_set_comps = []
                    this_mt_set_comp_data = []
                    this_mt_set_pvals = []
                    for comp in comparison_data_sets:
                        if comp[0] == compType and tuple(comp[4:8]) == testType:
                            this_mt_set_comps += [comp]
                            this_mt_set_comp_data += [comparison_data_sets[comp]]
                            this_mt_set_pvals += [comparison_data_sets[comp][0]]
                    if len(this_mt_set_comps) > 0:
                        # Prepare the multiple testing output file prefix
                        mt_prefix = __get_comparison_str(compType)+'_'+'_'.join(testType) ########## Better name?
                        if N_bootstrap is not None:
                            mt_prefix += '_bootstrap'
                        if output_prefix == '':
                            mt_prefix = output_directory+mt_prefix
                        else:
                            mt_prefix = output_directory+output_prefix+'_'+mt_prefix

                        run_mt_corr_fishers_comb_write_output(global_alpha,this_mt_set_pvals,mult_tests,
                            this_mt_set_comps,this_mt_set_comp_data,testType[0],testType[3],testType[1],testType[2],None,
                            __get_comparison_str(compType),[compType],cust_lines,cust_invs,mt_prefix)
                    else:
                        print("CHECK - no tests of class "+str(testType)+" for comparison "+str(compType))

    # comp: (comparison,chrom,inv,line,target_mode,count_mode,calc_mode,tail_mode,expected_direction)
    # comp_data: [pval,obs_stat,combined_allele_freqs,combined_freq_diff,ML_pat_allele_count,
    #       paternal_pop_allele_total,ML_est_expected_counts,ML_est_obs_logProbs,num_GEQ,num_total]




    # if calc_chi2_pval:
    #     # sep_pvals = [comp[4] for comp in comparison_data_sets]
    #     chi2_comps = []
    #     chi2_comp_data = []
    #     sep_chi2_pvals = []
    #     testTypes = set()
    #     comparisonTypes = set
    #     for comp in comparison_data_sets:
    #         if comp[6] == 'chi2_pval':
    #             chi2_comps += [comp]
    #             chi2_comp_data += [comparison_data_sets[comp]]
    #             sep_chi2_pvals += [comparison_data_sets[comp][0]]
    #             testTypes.add((comp[]))
    #     if len(chi2_comps) > 0:
    #         run_mt_corr_fishers_comb_write_output(global_alpha,sep_chi2_pvals,mult_tests,
    #             chi2_comps,chi2_comp_data,target_mode,tail_mode,count_mode,comp[6],None,
    #             comp_mode,comparisons,cust_lines,cust_invs,out_file_prefix)
    #     else:
    #         print("WARNING - calc_chi2_pval is set to true, but no comparisons have 'chi2_pval' calc_mode")


    # if calc_prod_pval:
    #     prod_comps = []
    #     prod_comp_data = []
    #     sep_prod_pvals = []
    #     for comp in comparison_data_sets:
    #         if comp[6] == 'prod_pval':
    #             prod_comps += [comp]
    #             prod_comp_data += [comparison_data_sets[comp]]
    #             sep_prod_pvals += [comparison_data_sets[comp][0]]
    #     if len(prod_comps) > 0:
    #         run_mt_corr_fishers_comb_write_output(global_alpha,sep_prod_pvals,mult_tests,
    #             prod_comps,prod_comp_data,target_mode,tail_mode,count_mode,comp[6],None,
    #             comp_mode,comparisons,cust_lines,cust_invs,out_file_prefix)
    #     else:
    #         print("WARNING - calc_prod_pval is set to true, but no comparisons have 'prod_pval' calc_mode")

    return



# Run when this script is called directly
if __name__ == '__main__':
    main()





# # For calculating a single p-value for the entire set of replicate lines,
# #   using Fisher's method:
# #      comparing to X^2 distribution of 2k degrees of freedom
# #      -2*sum(ln(p_i)) for k independent tests of p-values p_i
# #   Consider: the tests across lines are likely to be independent, and are the intended target of combination,
# #   but p-values for paternal-to-embryo are likely to be positively correlated with those of embryo-to-any adult offspring,
# #   and correlation between inversions is plausible but unclear
# #   IS THERE AN ISSUE WITH COMBINING multiple-hypothesis-testing-corrected p-values in this way?
# def combine_pvals_by_line(comparison_data,comparisons,lines,cust_invs,mult_test_corrections):
#     comb_pval_dict = {}
#     inversions = list(dict.fromkeys(list(inv_name_dict.values())))
#     if cust_invs is not None:
#         inversions = cust_invs
#     # print(inversions)
#     for comp in comparisons:
#         for inv in inversions:
#             indexes = []
#             chrom = ''
#             for i in np.arange(len(comparison_data)):
#                 if ((comparison_data[i][0] == comp) and (comparison_data[i][2] == inv) and (lines is None or comparison_data[i][3] in lines)):
#                     indexes += [i]
#                     chrom = comparison_data[i][1]
#             num_tests = len(indexes)
#             combined_pval_data = []
#             for (test_name,adjusted_pvals,null_rejected) in mult_test_corrections:
#                 corrected_pvals = []
#                 for i in indexes:
#                     corrected_pvals += [adjusted_pvals[i]]
#                 chi_square_score = np.sum(-2*np.log(corrected_pvals))
#                 degree_f = 2*num_tests
#                 comb_pval = sp.stats.distributions.chi2.sf(chi_square_score, degree_f)
#                 combined_pval_data += [(test_name,comb_pval)]

#             comb_pval_dict[(comp,chrom,inv)] = combined_pval_data
#     return comb_pval_dict

# # For writing line-combined significance data
# def write_combined_pval_data(comb_pval_dict,out_file_name,target_mode,tail_mode,count_mode):

#     if verbose: print("Writing combined p-value output to "+str(out_file_name))
    
#     # stat_full_name = ''
#     # if target_mode == "one_change_combined_cohorts":
#     #      stat_full_name += "Chi_Squared_on_Combined_"+count_mode+"_Counts"
#     # elif target_mode == "one_change_lone_cohorts":
#     #      stat_full_name += "Chi_Squared_on_Separate_"+count_mode+"_Counts"

#     header = "Comparison,Cohort_1,Cohort_2,Chrom,Inv,"
#     for (test_name,comb_pval) in list(comb_pval_dict.values())[0]:
#         header += "Fishers_Pval_on_"+test_name+"_Adjusted_Pval,"
#     header = header[:-1]+'\n'


#     with open(out_file_name, 'w') as out_file:
#         out_file.write(header)

#         # Expecting comparison data as [comparison,chrom,inv,line,p_val,obs_stat,combined_allele_freqs,
#         #  combined_freq_diff,ML_est_pat_count,paternal_pop_allele_total,ML_est_expected_counts,ML_est_obs_logProbs]
#         for (comp,chrom,inv) in comb_pval_dict:
#             line = '-'.join(comp)+','+group_full_name_dict[comp[0]]+','+group_full_name_dict[comp[1]]+','+\
#                 str(chrom)+','+str(inv)+','
#             for (test_name,comb_pval) in comb_pval_dict[(comp,chrom,inv)]:
#                 line += str(comb_pval)+','
#             line = line[:-1]+'\n'

#             out_file.write(line)
#     return
