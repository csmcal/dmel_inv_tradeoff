
# For calculating percentage point differences in inversion frequencies between experimental cohorts,
#   and bootstrap resampling confidence intervals on those differences
# CONSIDER - generating percent change data as well

# Import packages
import numpy as np
import argparse
import os
import scipy as sp
from scipy.stats import binom
from functools import partial


# Define and collect input arguments
def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('--repN', type=int, default=100000,
                        help='The number of replicate simulations of the experimental pipeline to perform')
    parser.add_argument('--varN', type=int, default=316,
                        help='The binomial sampling size expected to account for sample-independent library-prep-dependent variance')
    parser.add_argument('--rejT', type=float, default=0.01,
                        help='The threshold difference between sampled and observed read frequency sum, used to reject samples')
    parser.add_argument('--alpha', type=float, default=0.05,
                        help='The single-tailed global alpha level to generate a two-tailed bonferroni-corrected alpha from for each confidence interval calculation')
    parser.add_argument('--count_file', type=str, required=True,
                        help='The .csv file containing the estimated inverted and standard read counts, etc.')
                        # inv_freqs/freq_data.csv
    parser.add_argument('--library_data', type=str, required=True,
                        help='The .csv file with the library metadata')
                        # read_metadata.csv
    parser.add_argument('--od', type=str, required=True,
                        help='The output directory to write reshaped data, confidence interval data, and bootstrap replicates to')
                        # bootstrap/
    parser.add_argument('--out_prefix', type=str, default='',
                        help='The output file prefix to use')
    parser.add_argument('--calc_CI', dest='calc_mode', default=0,
                        action='store_const', const=1,
                        help='A flag to only calculate the resampled confidence intervals, based on no sample-independent variance')
    parser.add_argument('--calc_rej_pval', dest='calc_mode', default=0,
                        action='store_const', const=2,
                        help='A flag to only calculate rejection-sampled p-vals')
    parser.add_argument('--calc_dir_pval', dest='calc_mode', default=0,
                        action='store_const', const=3,
                        help='A flag to only calculate both directional p-vals')
    parser.add_argument('--calc_both_pval', dest='calc_mode', default=0,
                        action='store_const', const=4,
                        help='A flag to only calculate both directional p-vals')

    args = parser.parse_args()
    return args






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

                if adjusted_std_count < 1:
                    print("Experimental discrepancy: expected >= half of called reads to be standard orientation:")
                    print(sequencing_name+" has "+str(adjusted_std_count)+" called standard reads out of "+str(adjusted_total)+" after adjustment")

                count_data[(line,chrom,cohort)] = [sequencing_name,line,cohort,chrom,inv,allele_N,adjusted_inv_proportion,
                    num_inv,adjusted_std_count]
            else:
                count_data[(line,chrom,cohort)] = [sequencing_name,line,cohort,chrom,inv,allele_N,inv_prop,
                    num_inv,num_std]
    return count_data


# A helper function to order D. mel chromosome names, 
#   assumes only X, 2L, 2R, 3L, 3R will be passed
def chrom_order(chrom):
    chrom_val = 0
    if chrom == '2L':
        chrom_val == 1
    if chrom == '2R':
        chrom_val == 2
    elif chrom == '3L':
        chrom_val == 3
    elif chrom == '3R':
        chrom_val == 4
    return chrom_val

# A helper function to order the combination of line and chromosome
def experimental_order(line,chrom):
    line_val = int(line)
    chrom_val = chrom_order(chrom)
    sort_val = 10*line_val+chrom_val
    return sort_val

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
            for cohort in cohort_order:
                sample_size = library_metadata[(line,chrom,cohort)][12]
                inv_freq = None
                read_count = None
                if (line,chrom,cohort) in count_data:
                    inv_freq = count_data[(line,chrom,cohort)][6]
                    read_count = count_data[(line,chrom,cohort)][7]+count_data[(line,chrom,cohort)][8]
                experiment += [[cohort,sample_size,inv_freq,read_count]]
            experiment_data += [experiment]
    return experiment_data



# Run a single bootstrap re-sampling of the paternal pool
#   because the fathers were collected the individuals are here assumed to be the population
def resample_paternal(N_reads,N_rep):

    return

# Run a single bootstrap re-sampling of an offspring pool
def resample_offspring(N_reads,N_indiv,N_rep):

    return 


# Bootstrap resample and carry the frequencies forward for:
#   paternal reads, offspring alleles, offspring reads
#   then return the frequency difference
def resample_paternal_to_offspring(pat_read_freq,pat_read_N,off_allele_N,off_read_N):
    pat_read_freqs = np.array([pat_read_freq,1-pat_read_freq])
    # Resample paternal allele frequencies from reads
    pat_read_resample = np.random.multinomial(pat_read_N,pat_read_freqs)
    pat_allele_freqs = pat_read_resample/sum(pat_read_resample)
    # Resample offspring allele frequencies
    off_allele_resample = np.random.multinomial(off_allele_N,pat_allele_freqs)
    off_allele_freqs = off_allele_resample/sum(off_allele_resample)
    # Resample offspring read frequencies
    off_read_resample = np.random.multinomial(off_read_N,off_allele_freqs)
    off_read_freqs = off_read_resample/sum(off_read_resample)
    
    sampled_off_freq = off_read_freqs[0]
    return (sampled_off_freq)


# Bootstrap resample and carry the frequencies forward for:
#   paternal reads, offspring alleles, offspring reads
#   then return the frequency difference
def resample_offspring_to_paternal(off_read_freq,off_read_N,off_allele_N,pat_read_N):
    pat_read_freqs = np.array([pat_read_freq,1-pat_read_freq])
    # Resample paternal allele frequencies from reads
    pat_read_resample = np.random.multinomial(pat_read_N,pat_read_freqs)
    pat_allele_freqs = pat_read_resample/sum(pat_read_resample)
    # Resample offspring allele frequencies
    off_allele_resample = np.random.multinomial(off_allele_N,pat_allele_freqs)
    off_allele_freqs = off_allele_resample/sum(off_allele_resample)
    # Resample offspring read frequencies
    off_read_resample = np.random.multinomial(off_read_N,off_allele_freqs)
    off_read_freqs = off_read_resample/sum(off_read_resample)
    
    sampled_off_freq = off_read_freqs[0]
    return (sampled_off_freq)


# Bootstrap resample and carry the frequencies forward for:
#   offspring 1 reads, off 1 alleles, off 2 alleles, off 2 reads
#   then return the frequency difference
def resample_offspring_to_offspring(off_1_obs_freq,off_1_allele_N,off_1_read_N,
        off_2_allele_N,off_2_read_N):
    off_1_obs_freqs = np.array([off_1_obs_freq,1-off_1_obs_freq])
    # Resample first offspring read frequencies
    off_1_read_resample = np.random.multinomial(off_1_read_N,off_1_obs_freqs)
    off_1_read_freqs = off_1_read_resample/sum(off_1_read_resample)
    # Resample first offspring allele frequencies
    off_1_allele_resample = np.random.multinomial(off_1_allele_N,off_1_read_freqs)
    off_1_allele_freqs = off_1_allele_resample/sum(off_1_allele_resample)
    # Resample second offspring allele frequencies
    off_2_allele_resample = np.random.multinomial(off_2_allele_N,off_1_allele_freqs)
    off_2_allele_freqs = off_2_allele_resample/sum(off_2_allele_resample)
    # Resample second offspring read frequencies
    off_2_read_resample = np.random.multinomial(off_1_read_N,off_1_obs_freqs)
    off_2_read_freqs = off_2_read_resample/sum(off_2_read_resample)
    
    sampled_off_2_freq = off_2_read_freqs[0]
    return (sampled_off_2_freq)


def calc_pat_off_CIs(pat_data,off_data,N_reps,alpha):
    for value in pat_data+off_data:
        if value is None:
            return (None,None,None,None)

    pat_freq = pat_data[2]
    pat_read_N = pat_data[3]

    off_size = off_data[1]
    off_allele_N = 2*off_size
    off_freq = off_data[2]
    off_read_N = off_data[3]

    perc_point_diff = off_freq - pat_freq
    resampled_diffs = []
    for i in np.arange(N_reps):
        sampled_off_freq = resample_paternal_to_offspring(pat_freq,pat_read_N,off_allele_N,off_read_N)
        sampled_diff = sampled_off_freq - pat_freq
        resampled_diffs += [sampled_diff]

    (low_CI,high_CI) = np.quantile(resampled_diffs,[alpha,1-alpha])
    is_outside = perc_point_diff < low_CI or perc_point_diff > high_CI
    return (perc_point_diff,low_CI,high_CI,is_outside)

def calc_off_off_CIs(off_1_data,off_2_data,N_reps,alpha):
    for value in off_1_data+off_2_data:
        if value is None:
            return (None,None,None,None)

    off_1_size = off_1_data[1]
    off_1_allele_N = 2*off_1_size
    off_1_freq = off_1_data[2]
    off_1_read_N = off_1_data[3]

    off_2_size = off_2_data[1]
    off_2_allele_N = 2*off_2_size
    off_2_freq = off_2_data[2]
    off_2_read_N = off_2_data[3]

    perc_point_diff = off_1_freq - off_2_freq
    resampled_diffs = []
    for i in np.arange(N_reps):
        sampled_off_2_freq = resample_offspring_to_offspring(off_1_freq,
                off_1_allele_N,off_1_read_N,off_2_allele_N,off_2_read_N)
        sampled_diff = sampled_off_2_freq - off_1_freq
        resampled_diffs += [sampled_diff]

    (low_CI,high_CI) = np.quantile(resampled_diffs,[alpha,1-alpha])
    is_outside = perc_point_diff < low_CI or perc_point_diff > high_CI
    return (perc_point_diff,low_CI,high_CI,is_outside)


def calculate_confidence_intervals(experiment_metadata,N_reps,alpha):
    pat_data = experiment_metadata[3]
    emb_data = experiment_metadata[4]
    eaf_data = experiment_metadata[5]
    laf_data = experiment_metadata[6]
    eam_data = experiment_metadata[7]
    lam_data = experiment_metadata[8]

    pat_emb_diff_CI = list(calc_pat_off_CIs(pat_data,emb_data,N_reps,alpha))
    pat_emb_diff_CI += ['Parental ♂','Embryo']
    temp_pat_data = list(pat_data)
    temp_pat_data[1] = emb_data[1]
    rev_p_e_diff_CI = list(calc_pat_off_CIs(emb_data,temp_pat_data,N_reps,alpha))
    if any([val is None for val in rev_p_e_diff_CI]):
        pat_emb_diff_CI += [None,None,None,None,None]
    else:
        rev_p_e_low = -rev_p_e_diff_CI[2]
        rev_p_e_high = -rev_p_e_diff_CI[1]
        rev_sig = rev_p_e_diff_CI[3]
        diff_p_e_low = pat_emb_diff_CI[1] - rev_p_e_low
        diff_p_e_high = pat_emb_diff_CI[2] - rev_p_e_high
        pat_emb_diff_CI += [rev_p_e_low,rev_p_e_high,rev_sig,diff_p_e_low,diff_p_e_high]

    emb_eaf_diff_CI = list(calc_off_off_CIs(emb_data,eaf_data,N_reps,alpha))
    emb_eaf_diff_CI += ['Embryo','Early ♀']
    rev_emb_eaf_diff_CI = list(calc_off_off_CIs(eaf_data,emb_data,N_reps,alpha))
    if any([val is None for val in rev_emb_eaf_diff_CI]):
        emb_eaf_diff_CI += [None,None,None,None,None]
    else:
        rev_m_ef_low = -rev_emb_eaf_diff_CI[2]
        rev_m_ef_high = -rev_emb_eaf_diff_CI[1]
        rev_sig = rev_emb_eaf_diff_CI[3]
        diff_m_ef_low = emb_eaf_diff_CI[1] - rev_m_ef_low
        diff_m_ef_high = emb_eaf_diff_CI[2] - rev_m_ef_high
        emb_eaf_diff_CI += [rev_m_ef_low,rev_m_ef_high,rev_sig,diff_m_ef_low,diff_m_ef_high]

    emb_laf_diff_CI = list(calc_off_off_CIs(emb_data,laf_data,N_reps,alpha))
    emb_laf_diff_CI += ['Embryo','Late ♀']
    rev_emb_laf_diff_CI = list(calc_off_off_CIs(laf_data,emb_data,N_reps,alpha))
    if any([val is None for val in rev_emb_laf_diff_CI]):
        emb_laf_diff_CI += [None,None,None,None,None]
    else:
        rev_m_lf_low = -rev_emb_laf_diff_CI[2]
        rev_m_lf_high = -rev_emb_laf_diff_CI[1]
        rev_sig = rev_emb_laf_diff_CI[3]
        diff_m_lf_low = emb_laf_diff_CI[1] - rev_m_lf_low
        diff_m_lf_high = emb_laf_diff_CI[2] - rev_m_lf_high
        emb_laf_diff_CI += [rev_m_lf_low,rev_m_lf_high,rev_sig,diff_m_lf_low,diff_m_lf_high]

    emb_eam_diff_CI = list(calc_off_off_CIs(emb_data,eam_data,N_reps,alpha))
    emb_eam_diff_CI += ['Embryo','Early ♂']
    rev_emb_eam_diff_CI = list(calc_off_off_CIs(eam_data,emb_data,N_reps,alpha))
    if any([val is None for val in rev_emb_eam_diff_CI]):
        emb_eam_diff_CI += [None,None,None,None,None]
    else:
        rev_m_em_low = -rev_emb_eam_diff_CI[2]
        rev_m_em_high = -rev_emb_eam_diff_CI[1]
        rev_sig = rev_emb_eam_diff_CI[3]
        diff_m_em_low = emb_eam_diff_CI[1] - rev_m_em_low
        diff_m_em_high = emb_eam_diff_CI[2] - rev_m_em_high
        emb_eam_diff_CI += [rev_m_em_low,rev_m_em_high,rev_sig,diff_m_em_low,diff_m_em_high]

    emb_lam_diff_CI = list(calc_off_off_CIs(emb_data,lam_data,N_reps,alpha))
    emb_lam_diff_CI += ['Embryo','Late ♂']
    rev_emb_lam_diff_CI = list(calc_off_off_CIs(lam_data,emb_data,N_reps,alpha))
    if any([val is None for val in rev_emb_lam_diff_CI]):
        emb_lam_diff_CI += [None,None,None,None,None]
    else:
        rev_m_lm_low = -rev_emb_lam_diff_CI[2]
        rev_m_lm_high = -rev_emb_lam_diff_CI[1]
        rev_sig = rev_emb_lam_diff_CI[3]
        diff_m_lm_low = emb_lam_diff_CI[1] - rev_m_lm_low
        diff_m_lm_high = emb_lam_diff_CI[2] - rev_m_lm_high
        emb_lam_diff_CI += [rev_m_lm_low,rev_m_lm_high,rev_sig,diff_m_lm_low,diff_m_lm_high]

    CIs = [pat_emb_diff_CI,emb_eaf_diff_CI,emb_laf_diff_CI,emb_eam_diff_CI,emb_lam_diff_CI]
    return CIs

def print_significant(experiment_data_with_CI):
    for exp_data in experiment_data_with_CI:
        line = exp_data[0]
        chrom = exp_data[1]
        inv = exp_data[2]

def write_data_with_CIs(experiment_data_with_CI,out_file_name):
    header = "Line,Chrom,Inv,N_Pat_Ind,Pat_Inv_Freq,N_Pat_Reads,N_Emb_Ind,Emb_Inv_Freq,N_Emb_Reads,"
    header += "N_Eaf_Ind,Eaf_Inv_Freq,N_Eaf_Reads,N_Laf_Ind,Laf_Inv_Freq,N_Laf_Reads,"
    header += "N_Eam_Ind,Eam_Inv_Freq,N_Eam_Reads,N_Lam_Ind,Lam_Inv_Freq,N_Lam_Reads,"
    header += "Pat_Emb_Diff,Pat_Emb_CI_Low,Pat_Emb_CI_High,"
    header += "Emb_Eaf_Diff,Emb_Eaf_CI_Low,Emb_Eaf_CI_High,Emb_Laf_Diff,Emb_Laf_CI_Low,Emb_Laf_CI_High,"
    header += "Emb_Eam_Diff,Emb_Eam_CI_Low,Emb_Eam_CI_High,Emb_Lam_Diff,Emb_Lam_CI_Low,Emb_Lam_CI_High\n"
    with open(out_file_name, 'w') as out_file:
        out_file.write(header)
        for exp_data in experiment_data_with_CI:
            line = str(exp_data[0])+','+str(exp_data[1])+','+str(exp_data[2])
            for count_data in exp_data[3:9]:
                for value in count_data[1:]:
                    line += ','+str(value)
            for diff_CI_data in exp_data[9:14]:
                for value in diff_CI_data[:3]:
                    line += ','+str(value)
            line += '\n'
            out_file.write(line)
    return

def write_comparisons(experiment_data_with_CI,comparison_file):
    with open(comparison_file,'w') as out_file:
        header = "Line,Chrom,Inv,Cohort_1,Cohort_2,Freq_Diff,Low_CI,High_CI,Outside_CI,"
        header += "Rev_Low_CI,Rev_High_CI,Outside_Rev_CI,Diff_Low_CI,Diff_High_CI\n"
        out_file.write(header)
        for exp_data in experiment_data_with_CI:
            line_prefix = str(exp_data[0])+','+str(exp_data[1])+','+str(exp_data[2])+','
            for comparison_data in exp_data[9:]:
                line = line_prefix+str(comparison_data[4])+','+str(comparison_data[5])+\
                    ','+str(comparison_data[0])+','+str(comparison_data[1])+','+\
                    str(comparison_data[2])+','+str(comparison_data[3])+','+\
                    str(comparison_data[6])+','+str(comparison_data[7])+','+\
                    str(comparison_data[8])+','+str(comparison_data[9])+','+\
                    str(comparison_data[10])+'\n'
                out_file.write(line)
    return


# Rejection sampling - sampled frequencies must be within 0.01 of sum of the frequencies
def rejection_test(rej_threshold,obs_freq_sum,sampled_freq_sum):
    if abs(obs_freq_sum - sampled_freq_sum) > rej_threshold:
        return False
    else:
        return True

def rejection_pval_pat_off(pat_read_freq,pat_read_N,off_allele_N,off_read_N,off_read_freq,var_est_N,rej_N,rej_threshold):
    count_accepted = 0
    count_greater_diff = 0
    obs_freq_sum = pat_read_freq + off_read_freq
    vect_rej_test = np.vectorize(partial(rejection_test,rej_threshold,obs_freq_sum))
    obs_freq_diff = abs(pat_read_freq - off_read_freq)

    while count_accepted < rej_N:
        sampled_freqs = np.random.uniform(size=rej_N)

        # print(sampled_freqs)
        # print(off_allele_N)
        sampled_off_allele_counts = np.random.binomial(off_allele_N,sampled_freqs)
        sampled_off_allele_freqs = np.divide(sampled_off_allele_counts,off_allele_N)
        sampled_off_var_counts = np.random.binomial(var_est_N,sampled_off_allele_freqs)
        sampled_off_var_freqs = np.divide(sampled_off_var_counts,var_est_N)
        sampled_off_read_counts = np.random.binomial(off_read_N,sampled_off_var_freqs)
        sampled_off_read_freqs = np.divide(sampled_off_read_counts,off_read_N)

        # print(sampled_off_read_freqs)

        sampled_pat_var_counts = np.random.binomial(var_est_N,sampled_freqs)
        sampled_pat_var_freqs = np.divide(sampled_pat_var_counts,var_est_N)
        sampled_pat_read_counts = np.random.binomial(pat_read_N,sampled_pat_var_freqs)
        sampled_pat_read_freqs = np.divide(sampled_pat_read_counts,pat_read_N)

        # print(sampled_pat_read_freqs)

        sampled_abs_diffs = np.absolute(np.subtract(sampled_pat_read_freqs,sampled_off_read_freqs))
        # print(sampled_abs_diffs)
        sampled_sums = np.add(sampled_pat_read_freqs,sampled_off_read_freqs)
        # print(sampled_sums)
        accepted_vect_test = vect_rej_test(sampled_sums)
        # print(accepted_vect_test)
        accepted = np.isclose(sampled_sums,obs_freq_sum,rtol=0,atol=rej_threshold)
        # print(accepted)
        if np.any(np.logical_xor(accepted_vect_test,accepted)):
            print("ERROR FLAG - vectorized and numpy acceptance tests disagree")
        count_accepted += np.count_nonzero(accepted)
        # print(count_accepted)
        accepted_diffs = sampled_abs_diffs[accepted]
        count_greater_diff += np.count_nonzero(np.greater_equal(accepted_diffs,obs_freq_diff))

    prob = count_greater_diff/count_accepted
    print('Prob:\t\t\t'+str(prob))
    print('Acceptance count:\t\t'+str(count_accepted))
    return prob

def calc_pat_off_rej_pvals(pat_data,off_data,var_est_N,rej_N,rej_threshold,alpha):
    for value in pat_data+off_data:
        if value is None:
            return (None,None,None)

    pat_read_freq = pat_data[2]
    pat_read_N = pat_data[3]

    off_size = off_data[1]
    off_allele_N = 2*off_size
    off_read_freq = off_data[2]
    off_read_N = off_data[3]

    print("Pat freq: "+str(pat_read_freq))
    print("Off freq: "+str(off_read_freq))
    perc_point_diff = off_read_freq - pat_read_freq
    
    pval = rejection_pval_pat_off(pat_read_freq,pat_read_N,off_allele_N,off_read_N,off_read_freq,var_est_N,rej_N,rej_threshold)
    passes_alpha = alpha > pval

    return (perc_point_diff,pval,passes_alpha)


def rejection_pval_off_off(off1_read_freq,off1_read_N,off1_allele_N,off2_allele_N,off2_read_N,off2_read_freq,var_est_N,rej_N,rej_threshold):
    count_accepted = 0
    count_greater_diff = 0
    obs_freq_sum = off1_read_freq + off2_read_freq
    vect_rej_test = np.vectorize(partial(rejection_test,rej_threshold,obs_freq_sum))
    obs_freq_diff = abs(off1_read_freq - off2_read_freq)

    while count_accepted < rej_N:
        sampled_freqs = np.random.uniform(size=rej_N)

        sampled_off1_allele_counts = np.random.binomial(off1_allele_N,sampled_freqs)
        sampled_off1_allele_freqs = np.divide(sampled_off1_allele_counts,off1_allele_N)
        sampled_off1_var_counts = np.random.binomial(var_est_N,sampled_off1_allele_freqs)
        sampled_off1_var_freqs = np.divide(sampled_off1_var_counts,var_est_N)
        sampled_off1_read_counts = np.random.binomial(off1_read_N,sampled_off1_var_freqs)
        sampled_off1_read_freqs = np.divide(sampled_off1_read_counts,off1_read_N)

        sampled_off2_allele_counts = np.random.binomial(off2_allele_N,sampled_freqs)
        sampled_off2_allele_freqs = np.divide(sampled_off2_allele_counts,off2_allele_N)
        sampled_off2_var_counts = np.random.binomial(var_est_N,sampled_off2_allele_freqs)
        sampled_off2_var_freqs = np.divide(sampled_off2_var_counts,var_est_N)
        sampled_off2_read_counts = np.random.binomial(off2_read_N,sampled_off2_var_freqs)
        sampled_off2_read_freqs = np.divide(sampled_off2_read_counts,off2_read_N)

        sampled_abs_diffs = np.absolute(np.subtract(sampled_off1_read_freqs,sampled_off2_read_freqs))
        sampled_sums = np.add(sampled_off1_read_freqs,sampled_off2_read_freqs)
        accepted_vect_test = vect_rej_test(sampled_sums)
        accepted = np.isclose(sampled_sums,obs_freq_sum,rtol=0,atol=rej_threshold)
        if np.any(np.logical_xor(accepted_vect_test,accepted)):
            print("ERROR FLAG - vectorized and numpy acceptance tests disagree")
        count_accepted += np.count_nonzero(accepted)
        accepted_diffs = sampled_abs_diffs[accepted]
        count_greater_diff += np.count_nonzero(np.greater_equal(accepted_diffs,obs_freq_diff))

    prob = count_greater_diff/count_accepted
    print('Prob:\t\t\t'+str(prob))
    print('Acceptance count:\t\t'+str(count_accepted))
    return prob

def calc_off_off_rej_pvals(off_1_data,off_2_data,var_est_N,rej_N,rej_threshold,alpha):
    for value in off_1_data+off_2_data:
        if value is None:
            return (None,None,None)

    off1_size = off_1_data[1]
    off1_allele_N = 2*off1_size
    off1_read_freq = off_1_data[2]
    off1_read_N = off_1_data[3]

    off2_size = off_2_data[1]
    off2_allele_N = 2*off2_size
    off2_read_freq = off_2_data[2]
    off2_read_N = off_2_data[3]

    print("Off1 freq: "+str(off1_read_freq))
    print("Off2 freq: "+str(off2_read_freq))
    perc_point_diff = off2_read_freq - off1_read_freq
    
    pval = rejection_pval_off_off(off1_read_freq,off1_read_N,off1_allele_N,off2_allele_N,off2_read_N,off2_read_freq,var_est_N,rej_N,rej_threshold)
    passes_alpha = alpha > pval

    return (perc_point_diff,pval,passes_alpha)


def calculate_rej_pvals(experiment_metadata,var_est_N,rej_N,rej_threshold,alpha):
    pat_data = experiment_metadata[3]
    emb_data = experiment_metadata[4]
    eaf_data = experiment_metadata[5]
    laf_data = experiment_metadata[6]
    eam_data = experiment_metadata[7]
    lam_data = experiment_metadata[8]

    print("Calculating for experimental set - "+str(experiment_metadata[:3]))
    # print("Paternal data:")
    # print(pat_data)
    # print("Embryo data:")
    # print(emb_data)

    pat_emb_diff_pval = ['Parental ♂','Embryo'] + list(calc_pat_off_rej_pvals(pat_data,emb_data,var_est_N,rej_N,rej_threshold,alpha))
    if any([val is None for val in pat_emb_diff_pval]):
        print("\'None\' found in pat_emb pval diff. Why?")

    emb_eaf_diff_pval = ['Embryo','Early ♀'] + list(calc_off_off_rej_pvals(emb_data,eaf_data,var_est_N,rej_N,rej_threshold,alpha))
    if any([val is None for val in emb_eaf_diff_pval]):
        print("\'None\' found in emb_eaf pval diff. Why?")

    emb_laf_diff_pval = ['Embryo','Late ♀'] + list(calc_off_off_rej_pvals(emb_data,laf_data,var_est_N,rej_N,rej_threshold,alpha))
    if any([val is None for val in emb_laf_diff_pval]):
        print("\'None\' found in emb_laf pval diff. Why?")

    emb_eam_diff_pval = ['Embryo','Early ♂'] + list(calc_off_off_rej_pvals(emb_data,eam_data,var_est_N,rej_N,rej_threshold,alpha))
    if any([val is None for val in emb_eam_diff_pval]):
        print("\'None\' found in emb_eam pval diff. Why?")

    emb_lam_diff_pval = ['Embryo','Late ♂'] + list(calc_off_off_rej_pvals(emb_data,lam_data,var_est_N,rej_N,rej_threshold,alpha))
    if any([val is None for val in emb_lam_diff_pval]):
        print("\'None\' found in emb_lam pval diff. Why?")

    pvals = [pat_emb_diff_pval,emb_eaf_diff_pval,emb_laf_diff_pval,emb_eam_diff_pval,emb_lam_diff_pval]
    return pvals


def write_rej_pvals(experiment_data_with_dir_pvals,comparison_file):
    with open(comparison_file,'w') as out_file:
        header = "Line,Chrom,Inv,Cohort_1,Cohort_2,Freq_Diff,Pval,Passes_Alpha\n"
        out_file.write(header)
        for exp_data in experiment_data_with_dir_pvals:
            line_prefix = str(exp_data[0])+','+str(exp_data[1])+','+str(exp_data[2])+','
            for comparison_data in exp_data[9:]:
                line = line_prefix+str(comparison_data[0])+','+str(comparison_data[1])+\
                    ','+str(comparison_data[2])+','+str(comparison_data[3])+','+\
                    str(comparison_data[4])+'\n'
                out_file.write(line)
    return







# 2-tailed p-value calculated for the probability
#    of the observed or more extreme OFFSPRING read count,
#    given the observed PATERNAL read count,
#    multiplied by two for a two-tailed approximation
def pval_pat_to_off(pat_read_freq,pat_read_N,off_allele_N,off_read_N,off_read_freq,var_est_N):
    # print(set_freqs)
    # print(read_Ns)
    # print(inv_counts)
    # print(scaling_N)
    pot_pat_read_counts = np.arange(0,pat_read_N+1)
    # print(pot_pat_read_counts)
    # log_prob_read_count = binom.logpmf(pot_pat_read_counts,pat_read_N,set_freqs.reshape((len(set_freqs),1)))
    log_prob_patR_count = binom.logpmf(pot_pat_read_counts,pat_read_N,pat_read_freq)
    # print(log_prob_read_count)

    pot_pat_read_freqs = np.divide(pot_pat_read_counts,pat_read_N)
    pot_var_counts = np.arange(0,var_est_N+1)
    log_prob_var1_count = binom.logpmf(pot_var_counts.reshape((len(pot_var_counts),1)),var_est_N,pot_pat_read_freqs)
    running_log_sum = sp.special.logsumexp(np.add(log_prob_patR_count,log_prob_var1_count),axis=1)

    pot_var_freqs = np.divide(pot_var_counts,var_est_N)
    pot_off_allele_counts = np.arange(0,off_allele_N+1)
    log_prob_offA_count = binom.logpmf(pot_off_allele_counts.reshape((len(pot_off_allele_counts),1)),off_allele_N,pot_var_freqs)
    running_log_sum = sp.special.logsumexp(np.add(running_log_sum,log_prob_offA_count),axis=1)

    pot_off_allele_freqs = np.divide(pot_off_allele_counts,off_allele_N)
    # pot_var_counts = np.arange(1,var_est_N)
    log_prob_var2_count = binom.logpmf(pot_var_counts.reshape((len(pot_var_counts),1)),var_est_N,pot_off_allele_freqs)
    running_log_sum = sp.special.logsumexp(np.add(running_log_sum,log_prob_var2_count),axis=1)

    # pot_var_freqs = np.divide(pot_var_counts,var_est_N)
    # Just consider the frequencies where the difference is larger than observed:
    if pat_read_freq > off_read_freq:
        pot_off_read_counts = np.arange(0,int(off_read_N*off_read_freq)+1)
    else:
        pot_off_read_counts = np.arange(int(off_read_N*off_read_freq),off_read_N+1)
    log_prob_offR_count = binom.logpmf(pot_off_read_counts.reshape((len(pot_off_read_counts),1)),off_read_N,pot_var_freqs)
    # running_log_sum = sp.special.logsumexp(np.add(running_log_sum,log_prob_offR_count),axis=0)

    final_prob = np.exp(sp.special.logsumexp(np.add(running_log_sum,log_prob_offR_count)))

    two_tailed_prob = 2*final_prob

    # print(str(scaling_N)+'\t\t'+str(-total_log_prob))
    print('Prob:\t\t\t'+str(final_prob))
    print('2 Tailed Prob:\t\t'+str(two_tailed_prob))
    return two_tailed_prob

# 2-tailed p-value calculated for the probability
#    of the observed or more extreme PATERNAL read count,
#    given the observed OFFSPRING read count,
#    multiplied by two for a two-tailed approximation
def pval_off_to_pat(pat_read_freq,pat_read_N,off_allele_N,off_read_N,off_read_freq,var_est_N):
    pot_off_read_counts = np.arange(0,off_read_N+1)
    log_prob_offR_count = binom.logpmf(pot_off_read_counts,off_read_N,off_read_freq)
    # print(log_prob_offR_count.shape)

    pot_off_read_freqs = np.divide(pot_off_read_counts,off_read_N)
    pot_var_counts = np.arange(0,var_est_N+1)
    log_prob_var1_count = binom.logpmf(pot_var_counts.reshape((len(pot_var_counts),1)),var_est_N,pot_off_read_freqs)
    # print(log_prob_var1_count.shape)
    running_log_sum = sp.special.logsumexp(np.add(log_prob_offR_count,log_prob_var1_count),axis=1)
    # print(running_log_sum.shape)
    running_log_prob_sum = np.add(log_prob_offR_count,log_prob_var1_count)
    # print(running_log_prob_sum.shape)
    # print()

    pot_var_freqs = np.divide(pot_var_counts,var_est_N)
    pot_off_allele_counts = np.arange(0,off_allele_N+1)
    log_prob_offA_count = binom.logpmf(pot_off_allele_counts.reshape((len(pot_off_allele_counts),1)),off_allele_N,pot_var_freqs)
    # print(log_prob_offA_count.shape)
    running_log_sum = sp.special.logsumexp(np.add(running_log_sum,log_prob_offA_count),axis=1)
    # print(running_log_sum.shape)
    # running_log_prob_sum = np.add(running_log_prob_sum.reshape(1,*running_log_prob_sum.shape),
    #     log_prob_offA_count.reshape(*log_prob_offA_count.shape,1))
    # print(running_log_prob_sum.shape)
    # print()

    pot_off_allele_freqs = np.divide(pot_off_allele_counts,off_allele_N)
    log_prob_var2_count = binom.logpmf(pot_var_counts.reshape((len(pot_var_counts),1)),var_est_N,pot_off_allele_freqs)
    # print(log_prob_var2_count.shape)
    running_log_sum = sp.special.logsumexp(np.add(running_log_sum,log_prob_var2_count),axis=1)
    # print(running_log_sum.shape)
    # running_log_prob_sum = np.add(running_log_prob_sum.reshape(1,*running_log_prob_sum.shape),
    #     log_prob_var2_count.reshape(*log_prob_var2_count.shape,1,1))
    # print(running_log_prob_sum.shape)
    # print()

    # Just consider the frequencies where the difference is larger than observed:
    if off_read_freq > pat_read_freq:
        pot_pat_read_counts = np.arange(0,int(pat_read_N*pat_read_freq)+1)
    else:
        pot_pat_read_counts = np.arange(int(pat_read_N*pat_read_freq),pat_read_N+1)
    log_prob_patR_count = binom.logpmf(pot_pat_read_counts.reshape((len(pot_pat_read_counts),1)),pat_read_N,pot_var_freqs)

    final_prob = np.exp(sp.special.logsumexp(np.add(running_log_sum,log_prob_patR_count)))
    # alt_final_prob = np.exp(sp.special.logsumexp(np.add(running_log_prob_sum.reshape(1,*running_log_prob_sum.shape),
    #     log_prob_patR_count.reshape(*log_prob_patR_count.shape,1,1,1))))
    # print('Alt Prob:\t\t\t'+str(alt_final_prob))

    two_tailed_prob = 2*final_prob
    print('Prob:\t\t\t'+str(final_prob))
    print('2 Tailed Prob:\t\t'+str(two_tailed_prob))

    return two_tailed_prob


# 2-tailed p-value calculated for the probability
#    of the observed or more extreme SECOND OFFSPRING read count,
#    given the observed FIRST OFFSPRING read count,
#    multiplied by two for a two-tailed approximation
def pval_off_to_off(off1_read_freq,off1_read_N,off1_allele_N,off2_allele_N,off2_read_N,off2_read_freq,var_est_N):
    # print(set_freqs)
    # print(read_Ns)
    # print(inv_counts)
    # print(scaling_N)

    pot_off1_read_counts = np.arange(1,off1_read_N)
    log_prob_off1R_count = binom.logpmf(pot_off1_read_counts,off1_read_N,off1_read_freq)

    pot_off1_read_freqs = np.divide(pot_off1_read_counts,off1_read_N)
    pot_var_counts = np.arange(1,var_est_N)
    log_prob_var1_count = binom.logpmf(pot_var_counts.reshape((len(pot_var_counts),1)),var_est_N,pot_off1_read_freqs)
    running_log_sum = sp.special.logsumexp(np.add(log_prob_off1R_count,log_prob_var1_count),axis=1)
    
    pot_var_freqs = np.divide(pot_var_counts,var_est_N)
    pot_off1_allele_counts = np.arange(1,off1_allele_N)
    log_prob_off1A_count = binom.logpmf(pot_off1_allele_counts.reshape((len(pot_off1_allele_counts),1)),off1_allele_N,pot_var_freqs)
    running_log_sum = sp.special.logsumexp(np.add(running_log_sum,log_prob_off1A_count),axis=1)
    
    pot_off1_allele_freqs = np.divide(pot_off1_allele_counts,off1_allele_N)
    pot_off2_allele_counts = np.arange(1,off1_allele_N)
    log_prob_off2A_count = binom.logpmf(pot_off2_allele_counts.reshape((len(pot_off2_allele_counts),1)),off2_allele_N,pot_off1_allele_freqs)
    running_log_sum = sp.special.logsumexp(np.add(running_log_sum,log_prob_off2A_count),axis=1)

    pot_off2_allele_freqs = np.divide(pot_off2_allele_counts,off2_read_N)
    log_prob_var2_count = binom.logpmf(pot_var_counts.reshape((len(pot_var_counts),1)),var_est_N,pot_off2_allele_freqs)
    running_log_sum = sp.special.logsumexp(np.add(running_log_sum,log_prob_var2_count),axis=1)

    # Just consider the frequencies where the difference is larger than observed:
    if off1_read_freq > off2_read_freq:
        pot_off2_read_counts = np.arange(1,int(off2_read_N*off2_read_freq))
    else:
        pot_off2_read_counts = np.arange(int(off2_read_N*off2_read_freq),off2_read_N)
    log_prob_off2R_count = binom.logpmf(pot_off2_read_counts.reshape((len(pot_off2_read_counts),1)),off2_read_N,pot_var_freqs)

    final_prob = np.exp(sp.special.logsumexp(np.add(running_log_sum,log_prob_off2R_count)))

    two_tailed_prob = 2*final_prob

    print('Prob:\t\t\t'+str(final_prob))
    print('2 Tailed Prob:\t\t'+str(two_tailed_prob))
    return two_tailed_prob






def calc_pat_off_2tailed_dir_pvals(pat_data,off_data,var_est_N,alpha):
    for value in pat_data+off_data:
        if value is None:
            return (None,None,None,None,None,None)

    pat_read_freq = pat_data[2]
    pat_read_N = pat_data[3]

    off_size = off_data[1]
    off_allele_N = 2*off_size
    off_read_freq = off_data[2]
    off_read_N = off_data[3]

    perc_point_diff = off_read_freq - pat_read_freq
    
    pval = pval_pat_to_off(pat_read_freq,pat_read_N,off_allele_N,off_read_N,off_read_freq,var_est_N)
    passes_alpha = alpha > pval

    rev_pval = pval_off_to_pat(pat_read_freq,pat_read_N,off_allele_N,off_read_N,off_read_freq,var_est_N)
    rev_passes_alpha = alpha > rev_pval

    pval_ratio = pval/rev_pval

    return (perc_point_diff,pval,passes_alpha,rev_pval,rev_passes_alpha,pval_ratio)

def calc_off_off_2tailed_dir_pvals(off_1_data,off_2_data,var_est_N,alpha):
    for value in off_1_data+off_2_data:
        if value is None:
            return (None,None,None,None,None,None)

    off1_size = off_1_data[1]
    off1_allele_N = 2*off1_size
    off1_read_freq = off_1_data[2]
    off1_read_N = off_1_data[3]

    off2_size = off_2_data[1]
    off2_allele_N = 2*off2_size
    off2_read_freq = off_2_data[2]
    off2_read_N = off_2_data[3]

    perc_point_diff = off2_read_freq - off1_read_freq
    
    pval = pval_off_to_off(off1_read_freq,off1_read_N,off1_allele_N,off2_allele_N,off2_read_N,off2_read_freq,var_est_N)
    passes_alpha = alpha > pval

    rev_pval = pval_off_to_off(off2_read_freq,off2_read_N,off2_allele_N,off1_allele_N,off1_read_N,off1_read_freq,var_est_N)
    rev_passes_alpha = alpha > rev_pval

    pval_ratio = pval/rev_pval

    return (perc_point_diff,pval,passes_alpha,rev_pval,rev_passes_alpha,pval_ratio)



def calculate_2tailed_dir_pvals(experiment_metadata,var_est_N,alpha):
    pat_data = experiment_metadata[3]
    emb_data = experiment_metadata[4]
    eaf_data = experiment_metadata[5]
    laf_data = experiment_metadata[6]
    eam_data = experiment_metadata[7]
    lam_data = experiment_metadata[8]

    print("Calculating for experimental set - "+str(experiment_metadata[:3]))

    pat_emb_diff_pval = ['Parental ♂','Embryo'] + list(calc_pat_off_2tailed_dir_pvals(pat_data,emb_data,var_est_N,alpha))
    if any([val is None for val in pat_emb_diff_pval]):
        print("\'None\' found in pat_emb pval diff. Why?")

    emb_eaf_diff_pval = ['Embryo','Early ♀'] + list(calc_off_off_2tailed_dir_pvals(emb_data,eaf_data,var_est_N,alpha))
    if any([val is None for val in emb_eaf_diff_pval]):
        print("\'None\' found in emb_eaf pval diff. Why?")

    emb_laf_diff_pval = ['Embryo','Late ♀'] + list(calc_off_off_2tailed_dir_pvals(emb_data,laf_data,var_est_N,alpha))
    if any([val is None for val in emb_laf_diff_pval]):
        print("\'None\' found in emb_laf pval diff. Why?")

    emb_eam_diff_pval = ['Embryo','Early ♂'] + list(calc_off_off_2tailed_dir_pvals(emb_data,eam_data,var_est_N,alpha))
    if any([val is None for val in emb_eam_diff_pval]):
        print("\'None\' found in emb_eam pval diff. Why?")

    emb_lam_diff_pval = ['Embryo','Late ♂'] + list(calc_off_off_2tailed_dir_pvals(emb_data,lam_data,var_est_N,alpha))
    if any([val is None for val in emb_lam_diff_pval]):
        print("\'None\' found in emb_lam pval diff. Why?")

    pvals = [pat_emb_diff_pval,emb_eaf_diff_pval,emb_laf_diff_pval,emb_eam_diff_pval,emb_lam_diff_pval]
    return pvals


def write_2tailed_dir_pval_comparisons(experiment_data_with_dir_pvals,comparison_file):
    with open(comparison_file,'w') as out_file:
        header = "Line,Chrom,Inv,Cohort_1,Cohort_2,Freq_Diff,For_Pval,For_Passes_Alpha,"
        header += "Rev_Pval,Rev_Passes_Alpha,Pval_Ratio\n"
        out_file.write(header)
        for exp_data in experiment_data_with_dir_pvals:
            line_prefix = str(exp_data[0])+','+str(exp_data[1])+','+str(exp_data[2])+','
            for comparison_data in exp_data[9:]:
                line = line_prefix+str(comparison_data[0])+','+str(comparison_data[1])+\
                    ','+str(comparison_data[2])+','+str(comparison_data[3])+','+\
                    str(comparison_data[4])+','+str(comparison_data[5])+','+\
                    str(comparison_data[6])+','+str(comparison_data[7])+'\n'
                out_file.write(line)
    return



# Main function, for running the pipeline when the script is called
def main():

    # Parse arguments
    args = parse_args()

    # Assign input arguments to variables
    called_read_count_file = args.count_file
    library_metadata_file_path = args.library_data
    output_prefix = args.out_prefix
    output_directory = args.od

    global_alpha = args.alpha
    N_reps = args.repN
    var_est_N = args.varN
    rej_threshold = args.rejT

    calc_mode = args.calc_mode
    calc_CI = calc_mode == 0 or calc_mode == 1
    calc_rej = calc_mode == 0 or calc_mode == 2 or calc_mode == 4
    calc_dir = calc_mode == 0 or calc_mode == 3 or calc_mode == 4

    # Ensure the output directory exists and the string is sanitized
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    if output_directory[-1] != '/':
        output_directory += '/'



    library_metadata = parse_library_metadata_to_dict(library_metadata_file_path)
    count_data = parse_count_data_to_dict(called_read_count_file)
    experiment_data = reshape_read_count_data(count_data,library_metadata)


    # CALCULATING CONFIDENCE INTERVALS
    if calc_CI:
        num_comparisons = 4*4*5 # lines*invs*compar, Not immediately easy to calculate this algorithmically
        CI_alpha = global_alpha / (2*num_comparisons) # Two-tailed bonferroni correction

        experiment_data_with_CI = []
        for experiment in experiment_data:
            confidence_intervals = calculate_confidence_intervals(experiment,N_reps,CI_alpha)
            experiment_with_CI = experiment + confidence_intervals
            experiment_data_with_CI += [experiment_with_CI]

        data_file = str(N_reps)+'rep_CI_data.csv'
        if output_prefix == '':
            data_file = output_directory+data_file
        else:
            data_file = output_directory+output_prefix+'_'+data_file
        write_data_with_CIs(experiment_data_with_CI,data_file)

        comparison_file = str(N_reps)+'rep_comparisons.csv'
        if output_prefix == '':
            comparison_file = output_directory+comparison_file
        else:
            comparison_file = output_directory+output_prefix+'_'+comparison_file
        write_comparisons(experiment_data_with_CI,comparison_file)



    # CALCULATING REJECTION SAMPLED P-VALS
    if calc_rej:
        num_comparisons = 4*4*5 # lines*invs*compar, Not immediately easy to calculate this algorithmically
        rej_pval_alpha = global_alpha / (2*num_comparisons) # Two-tailed bonferroni correction

        experiment_data_with_rej_pvals = []
        for experiment in experiment_data:
            rej_pvals = calculate_rej_pvals(experiment,var_est_N,N_reps,rej_threshold,rej_pval_alpha)
            experiment_with_rej_pvals = experiment + rej_pvals
            experiment_data_with_rej_pvals += [experiment_with_rej_pvals]

        data_file = str(N_reps)+'rep_rej_pval_data.csv'
        if output_prefix == '':
            data_file = output_directory+data_file
        else:
            data_file = output_directory+output_prefix+'_'+data_file
        write_rej_pvals(experiment_data_with_rej_pvals,data_file)


    # CALCULATING DIRECTIONAL P-VALS
    if calc_dir:
        num_comparisons = 4*4*5 # lines*invs*compar, Not immediately easy to calculate this algorithmically
        dir_pval_alpha = global_alpha / (2*num_comparisons) # Two-tailed bonferroni correction

        experiment_data_with_dir_pvals = []
        for experiment in experiment_data:
            dir_pvals = calculate_2tailed_dir_pvals(experiment,var_est_N,dir_pval_alpha)
            experiment_with_dir_pvals = experiment + dir_pvals
            experiment_data_with_dir_pvals += [experiment_with_dir_pvals]

        comparison_file = '2tail_directional_pval_comparisons.csv'
        if output_prefix == '':
            comparison_file = output_directory+comparison_file
        else:
            comparison_file = output_directory+output_prefix+'_'+comparison_file
        write_2tailed_dir_pval_comparisons(experiment_data_with_dir_pvals,comparison_file)

    return


# Run when this script is called directly
if __name__ == '__main__':
    main()


