
# For extimating a scaling factor to apply to the sample size in resampling,
#   for use in bootstrap resampling confidence intervals on frequency estimates from illumina sequenced amplicons

# Example use:
#   python scaling_factor_estimation.py --count_file inv_freqs/ctrl_freq_data.csv --library_data var_control_meta.csv --od est_factors/ --out_prefix test 

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

    parser.add_argument('--n', type=int, default=1000000,
                        help='The number of replicate simulations of the experimental pipeline to perform')
    parser.add_argument('--count_file', type=str,
                        help='The .csv file containing the estimated inverted and standard read counts, etc.')
                        # inv_freqs/freq_data.csv
    parser.add_argument('--library_data', type=str,
                        help='The .csv file with the library metadata')
                        # var_control_meta.csv
    parser.add_argument('--od', type=str,
                        help='The output directory to write reshaped data, confidence interval data, and bootstrap replicates to')
                        # est_factors/
    parser.add_argument('--out_prefix', type=str, default='',
                        help='The output file prefix to use')
    parser.add_argument('--scan', action='store_true',
                        help='Perform a scan across scalings for graphing')

    args = parser.parse_args()
    return args






# Extracts data about the paired-end read pools to be used in analysis
# Includes file names to analyze
def parse_metadata_to_dict(read_metadata_file_path):
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
def parse_control_count_data_to_dict(count_file_path):
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

            inv_prop = float(data[6])
            std_prop = float(data[7])
            num_inv = int(data[8])
            num_std = int(data[9])
            num_recomb = int(data[10])
            num_mismap = int(data[11])
            num_multimap = int(data[12])

            count_data[(line,chrom,cohort)] = count_data.get((line,chrom,cohort),[]) + [[sample_ind_N,inv_prop,num_inv+num_std,
                num_inv,num_std]]

    return count_data





# # Bootstrap resample a new read count from a scaled allele count,
# #    returns the frequency of the allele in the resample
# def resample_at_scale(set_freq,set_scaling,allele_N,read_N):
#     set_freqs = np.array([set_freq,1-set_freq])
#     # Rescale the allele N
#     rescaled_allele_N = set_scaling*allele_N
#     # Resample the read frequencies
#     read_resample = np.random.multinomial(rescaled_allele_N,set_freqs)
#     read_freqs = read_resample/sum(read_resample)
#     # Return the frequency of the allele in question
#     sampled_freq = read_freqs[0]
#     return (sampled_freq)


# Returns the avg variance of num_reps new draws of read frequencies
#   of the given set of read counts based on the given freqency, allele count rescaling, and allele count
def avg_resample_variance(set_freq,allele_N,read_N_list,num_reps,set_scaling,bess_correction=True):
    # Use numpy to quickly resample the reads
    cast_N = np.array(read_N_list).reshape((len(read_N_list),1))
    rescaled_cast_N = np.multiply(cast_N,set_scaling).astype('int64')
    cast_p = np.repeat(set_freq,num_reps).reshape((1,num_reps))
    sampled_counts = np.random.binomial(rescaled_cast_N,cast_p)
    sampled_freqs = np.divide(sampled_counts,rescaled_cast_N)
    sampled_vars = 0
    if bess_correction:
        sampled_vars = np.var(sampled_freqs,axis=0,ddof=1)
    else:
        sampled_vars = np.var(sampled_freqs,axis=0)
    mean_sampled_var = np.mean(sampled_vars)
    # print("Mean variance of resampled read frequencies:\n"+str(mean_sampled_var))


    # Use the analytic solution for the expectation of variance as num_reps -> INF
    #   = E(1/read_N)*set_freq*(1-set_freq)
    float_rescaled_Ns = rescaled_cast_N.astype('float64')
    mean_inverse_read_N = np.mean(np.power(float_rescaled_Ns,-1))
    analytic_test_var = mean_inverse_read_N*set_freq*(1-set_freq)
    # print("Analytically calculated expectation:\n"+str(analytic_test_var))
    # print("Ratio: "+str(mean_sampled_var/analytic_test_var))

    # # Use the analytic solution for the expectation of variance as num_reps -> INF
    # #   = E(read_N^2)*set_freq^2 + E(read_N)*set_freq*(1-set_freq) + E(read_N)^2*set_freq^2
    # test_read_Ns = np.array(read_N_list)
    # mean_read_N = np.mean(test_read_Ns)
    # mean_squared_read_N = np.mean(np.power(test_read_Ns,2))
    # analytic_test_var = mean_squared_read_N*set_freq**2 + mean_read_N*set_freq*(1-set_freq) + mean_read_N**2*set_freq**2
    # print("Analytically calculated expectation:\n"+str(analytic_test_var))

    return mean_sampled_var

# Generates the average sampled variance for a cohort for each rescaling factor in scaling_factors
def scan_scaling(obs_var,set_freq,allele_N,read_N_list,num_reps,scaling_factors=np.arange(0.025,1,0.025)):
    sampled_avg_vars = []
    for scaling_factor in scaling_factors:
        print("Scaling factor  "+str(scaling_factor))
        sampled_avg_vars += [avg_resample_variance(set_freq,allele_N,read_N_list,num_reps,scaling_factor)]
    var_diffs = []
    for var in sampled_avg_vars:
        var_diffs += [var-obs_var]
    print(sampled_avg_vars)
    print(var_diffs)
    return(sampled_avg_vars,var_diffs)

def scan_scaling_factors_by_extraction(count_data,num_reps,scaling_factors=np.arange(0.025,1,0.025),expected_freq_dict=None):
    scaling_factor_scan_dict = {}
    if expected_freq_dict is not None:
        for cohort_key in expected_freq_dict:
            read_N_list = []
            obs_freqs = []
            sample_ind_N_list = []
            for run_data in count_data[cohort_key]:
                read_N_list += [run_data[2]]
                obs_freqs += [run_data[1]]
                sample_ind_N_list += [run_data[0]]
            ind_N = sample_ind_N_list[0]
            if not all([N == ind_N for N in sample_ind_N_list]):
                print("Not all sample individual numbers where the same in cohort "+\
                    str(cohort_key[0])+' '+str(cohort_key[1])+' '+str(cohort_key[2])+':\n'+\
                    str(sample_ind_N_list))
            allele_N = 2*ind_N
            obs_var = np.var(obs_freqs,ddof=1)
            set_freq = expected_freq_dict[cohort_key]
            print("\n\nScanning scaling factors for "+str(cohort_key[0])+' '+str(cohort_key[1])+' '+\
                str(cohort_key[2])+':\nObserved variance:  '+str(obs_var))
            (sampled_avg_vars,var_diffs) = scan_scaling(obs_var,set_freq,allele_N,read_N_list,num_reps,scaling_factors)
            scaling_factor_scan_dict[cohort_key] = (scaling_factors,sampled_avg_vars,var_diffs,obs_var,allele_N)
    else:
        for cohort_key in count_data:
            read_N_list = []
            obs_freqs = []
            sample_ind_N_list = []
            for run_data in count_data[cohort_key]:
                read_N_list += [run_data[2]]
                obs_freqs += [run_data[1]]
                sample_ind_N_list += [run_data[0]]
            ind_N = sample_ind_N_list[0]
            if not all([N == ind_N for N in sample_ind_N_list]):
                print("Not all sample individual numbers where the same in cohort "+\
                    str(cohort_key[0])+' '+str(cohort_key[1])+' '+str(cohort_key[2])+':\n'+\
                    str(sample_ind_N_list))
            allele_N = 2*ind_N
            obs_var = np.var(obs_freqs,ddof=1)
            set_freq = np.mean(obs_freqs)
            print("\n\nScanning scaling factors for "+str(cohort_key[0])+' '+str(cohort_key[1])+' '+\
                str(cohort_key[2])+':\nObserved variance:  '+str(obs_var))
            (sampled_avg_vars,var_diffs) = scan_scaling(obs_var,set_freq,allele_N,read_N_list,num_reps,scaling_factors)
            scaling_factor_scan_dict[cohort_key] = (scaling_factors,sampled_avg_vars,var_diffs,obs_var,allele_N)
    return scaling_factor_scan_dict


# Uses a simple deterministic local walk to estimate the scaling factor
#   that minimizes the difference between observed and avg sampled variance
def estimate_scaling(obs_var,set_freq,allele_N,read_N_list,num_reps,delta=1e-8):
    lower_est = 1e-4 #ideally we never have an estimate where n < 2 (or even n < 10)
    lower_var = avg_resample_variance(set_freq,allele_N,read_N_list,num_reps,lower_est)
    upper_est = 1
    upper_var = avg_resample_variance(set_freq,allele_N,read_N_list,num_reps,upper_est)
    # print('{0:8.7f}\t\t{1:8.7f}\t\t{2:8.7f}\t\t{3:8.7f}'.format(lower_est,lower_var,upper_var,upper_est))
    while abs(lower_est-upper_est) > delta:

        # slope = (abs(obs_var-curr_var)-abs(obs_var-prev_var))/(curr_est-prev_est)
        slope = (upper_var-lower_var)/(upper_est-lower_est)

        # Guess halfway to the linear intercept between the two, in both directions
        guess_upper_est = upper_est + 0.5*(obs_var-upper_var)/slope
        guess_upper_var = avg_resample_variance(set_freq,allele_N,read_N_list,num_reps,guess_upper_est)

        guess_lower_est = lower_est + 0.5*(obs_var-lower_var)/slope
        guess_lower_var = avg_resample_variance(set_freq,allele_N,read_N_list,num_reps,guess_lower_est)
        
        # Replace one or both of the bounds with the guess
        if guess_upper_var > obs_var and guess_lower_var > obs_var:
            lower_est = guess_upper_est
            lower_var = guess_upper_var
        elif guess_upper_var < obs_var and guess_lower_var < obs_var:
            upper_est = guess_lower_est
            upper_var = guess_lower_var
        elif guess_upper_var < obs_var and guess_lower_var > obs_var:
            lower_est = guess_lower_est
            lower_var = guess_lower_var
            upper_est = guess_upper_est
            upper_var = guess_upper_var
        else:
            print("Mismatched bounds in estimation:")
            print('{0:8.7f}\t\t{1:8.7f}\t\t{2:8.7f}\t\t{3:8.7f}'.format(lower_est,lower_var,upper_var,upper_est))
            break

        # print('{0:8.7f}\t\t{1:8.7f}\t\t{2:8.7f}\t\t{3:8.7f}'.format(lower_est,lower_var,upper_var,upper_est))

    estimated_scaling_factor = np.mean([lower_est,upper_est])
    var_from_est_scaling = avg_resample_variance(set_freq,allele_N,read_N_list,num_reps,estimated_scaling_factor)
    print("Estimated scaling and associated resampled variance as:\n"+\
        str(estimated_scaling_factor)+'\t\t'+str(var_from_est_scaling))

    return (estimated_scaling_factor,var_from_est_scaling)

# Uses a simple deterministic local walk to estimate the scaling factor
#   that minimizes the difference between observed and avg sampled variance
def estimate_scaling_single_bound(obs_var,set_freq,allele_N,read_N_list,num_reps,delta=1e-8):
    lower_est = 1e-4 #ideally we never have an estimate where n < 2 (or even n < 10)
    lower_var = avg_resample_variance(set_freq,allele_N,read_N_list,num_reps,lower_est)
    upper_est = 1
    upper_var = avg_resample_variance(set_freq,allele_N,read_N_list,num_reps,upper_est)
    print('{0:8.7f}\t\t{1:8.7f}\t\t{2:8.7f}\t\t{3:8.7f}'.format(lower_est,lower_var,upper_var,upper_est))
    while abs(lower_est-upper_est) > delta:

        # slope = (abs(obs_var-curr_var)-abs(obs_var-prev_var))/(curr_est-prev_est)
        slope = (upper_var-lower_var)/(upper_est-lower_est)
        print(slope)

        # Guess the linear intercept between the two
        guess_est = upper_est + (obs_var-upper_var)/slope
        print(guess_est)
        guess_var = avg_resample_variance(set_freq,allele_N,read_N_list,num_reps,guess_est)
        
        # Replace one of the bounds with the guess
        if guess_var > obs_var:
            lower_est = guess_est
            lower_var = guess_var
        else:
            upper_est = guess_est
            upper_var = guess_var

        print('{0:8.7f}\t\t{1:8.7f}\t\t{2:8.7f}\t\t{3:8.7f}'.format(lower_est,lower_var,upper_var,upper_est))

    estimated_scaling_factor = np.mean([lower_est,upper_est])
    var_from_est_scaling = avg_resample_variance(set_freq,allele_N,read_N_list,num_reps,estimated_scaling_factor)

    return (estimated_scaling_factor,var_from_est_scaling)


# Uses a simple deterministic local walk to estimate the scaling factor
#   that minimizes the difference between observed and avg sampled variance
def estimate_scaling_semi_newton(obs_var,set_freq,allele_N,read_N_list,num_reps,delta=1e-8):
    prev_est = 0.25
    prev_var = avg_resample_variance(set_freq,allele_N,read_N_list,num_reps,prev_est)
    print('{0:8.7f}\t\t{1:8.7f}'.format(prev_est,prev_var))
    curr_est = 0.75
    curr_var = avg_resample_variance(set_freq,allele_N,read_N_list,num_reps,curr_est)
    print('{0:05}\t\t{1:05}'.format(curr_est,curr_var))
    # while abs(obs_var-curr_var) > delta:
    while abs(curr_est-prev_est) > delta:
        print('{0:05}\t\t{1:05}'.format(curr_est,curr_var))
        if curr_est == prev_est:
            break
        # slope = (abs(obs_var-curr_var)-abs(obs_var-prev_var))/(curr_est-prev_est)
        slope = (curr_var-prev_var)/(curr_est-prev_est)
        prev_est=curr_est
        prev_var=curr_var
        curr_est=curr_est + (obs_var-curr_var)/slope # Linear estimate of where the sampled variance = observed variance
        curr_est=max(1e-4,curr_est)
        curr_est=min(1,curr_est)
        curr_var=avg_resample_variance(set_freq,allele_N,read_N_list,num_reps,curr_est)
    return (curr_est, curr_var)

def estimate_scaling_factors_by_extraction(count_data,num_reps,expected_freq_dict=None):
    scaling_factor_dict = {}
    if expected_freq_dict is not None:
        for cohort_key in expected_freq_dict:
            read_N_list = []
            obs_freqs = []
            sample_ind_N_list = []
            for run_data in count_data[cohort_key]:
                read_N_list += [run_data[2]]
                obs_freqs += [run_data[1]]
                sample_ind_N_list += [run_data[0]]
            ind_N = sample_ind_N_list[0]
            if not all([N == ind_N for N in sample_ind_N_list]):
                print("Not all sample individual numbers where the same in cohort "+\
                    str(cohort_key[0])+' '+str(cohort_key[1])+' '+str(cohort_key[2])+':\n'+\
                    str(sample_ind_N_list))
            allele_N = 2*ind_N
            obs_var = np.var(obs_freqs,ddof=1)
            set_freq = expected_freq_dict[cohort_key]
            print("\n\nEstimating scaling for "+str(cohort_key[0])+' '+str(cohort_key[1])+' '+\
                        str(cohort_key[2])+':\nObserved variance:  '+str(obs_var))
            (estimated_scaling_factor,var_from_est_scaling) = estimate_scaling(obs_var,set_freq,allele_N,read_N_list,num_reps)
            scaling_factor_dict[cohort_key] = [estimated_scaling_factor,var_from_est_scaling,obs_var,allele_N]
    else:
        for cohort_key in count_data:
            read_N_list = []
            obs_freqs = []
            sample_ind_N_list = []
            for run_data in count_data[cohort_key]:
                read_N_list += [run_data[2]]
                obs_freqs += [run_data[1]]
                sample_ind_N_list += [run_data[0]]
            ind_N = sample_ind_N_list[0]
            if not all([N == ind_N for N in sample_ind_N_list]):
                print("Not all sample individual numbers where the same in cohort "+\
                    str(cohort_key[0])+' '+str(cohort_key[1])+' '+str(cohort_key[2])+':\n'+\
                    str(sample_ind_N_list))
            allele_N = 2*ind_N
            obs_var = np.var(obs_freqs,ddof=1)
            set_freq = np.mean(obs_freqs)
            print("\n\nEstimating scaling for "+str(cohort_key[0])+' '+str(cohort_key[1])+' '+\
                        str(cohort_key[2])+':\nObserved variance:  '+str(obs_var))
            (estimated_scaling_factor,var_from_est_scaling) = estimate_scaling(obs_var,set_freq,allele_N,read_N_list,num_reps)
            scaling_factor_dict[cohort_key] = [estimated_scaling_factor,var_from_est_scaling,obs_var,allele_N]
    return scaling_factor_dict


# Uses matplotlib.pyplot to plot a scaling factor scan
def plot_scaling_factor_scan(scaling_factors,sampled_avg_vars,var_diffs,obs_var):
    import matplotlib.pyplot as plt

    # Make a simple plot
    temp_factor=[0,1]
    temp_obs_var=[obs_var,obs_var]
    plt.plot(temp_factor,temp_obs_var,label="Observed Variance")
    plt.plot(scaling_factors, sampled_avg_vars, label="Sampled Variance")
    plt.legend()
    plt.show()

    return

def plot_all_cohort_scaling_factor_scan_diffs_together(scaling_factor_scan_dict):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()

    # Make a simple plot containing all the scans at once
    for cohort_key in scaling_factor_scan_dict:
        scaling_factors = scaling_factor_scan_dict[cohort_key][0]
        var_diffs = scaling_factor_scan_dict[cohort_key][2]
        ax.plot(scaling_factors, var_diffs, label="Cohort "+str(cohort_key[0])+' '+\
                        str(cohort_key[1])+' '+str(cohort_key[2]))

    # Scale the x-axis (scaling factor) to be exponential
    exp = lambda x: 10**(x)
    log = lambda x: np.log(x)
    ax.set_xscale('function', functions=(log, exp))
    ax.set_yscale('function', functions=(exp, log))
    ax.set(xlim=(1e-6,0.9))
    ax.set_xticks(list(np.power(10,np.arange(-4.5,0,.5))))

    # Scale the x-axis (scaling factor) to be exponential
    # exp = lambda x: 10**(x)
    # log = lambda x: np.log(x)
    # ax.set_yscale('function', functions=(exp, log))
    ax.set(ylim=(-2e-4,2e-3))
    # ax.set_yticks(list(np.power(10,np.arange(-4.5,0,.5))))

    # Add a line at y=0, where the difference between observed and sampled is 0
    ax.axhline(dashes=(2,3),color='lightgrey')

    ax.legend()
    plt.show()
    return



def avg_var_dist_to_obs_among_resamples(obs_vars,set_freqs,allele_N_list,read_N_sets,num_reps,scaling_factor):
    var_dist = []
    for i in np.arange(len(obs_vars)):
        sample_var = avg_resample_variance(set_freqs[i],allele_N_list[i],read_N_sets[i],num_reps,scaling_factor)
        var_dist = abs(obs_vars[i]-sample_var)
        var_dist += [var_dist]
    avg_var_dist = np.mean(var_dist)
    return avg_var_dist

def avg_var_among_resamples(set_freqs,allele_N_list,read_N_sets,num_reps,scaling_factor):
    var_list = []
    for i in np.arange(len(set_freqs)):
        sample_var = avg_resample_variance(set_freqs[i],allele_N_list[i],read_N_sets[i],num_reps,scaling_factor)
        var_list += [sample_var]
    avg_var = np.mean(var_list)
    return avg_var


def coestimate_scaling_factors(cohort_list,count_data,num_reps,expected_freq_dict=None,delta=1e-8):
    print("\n\nCoestimating scaling factor for:")
    obs_vars = []
    set_freqs = []
    allele_N_list = []
    read_N_sets = []

    for cohort_key in cohort_list:
        read_N_list = []
        obs_freqs = []
        sample_ind_N_list = []

        for run_data in count_data[cohort_key]:
            read_N_list += [run_data[2]]
            obs_freqs += [run_data[1]]
            sample_ind_N_list += [run_data[0]]
        ind_N = sample_ind_N_list[0]
        if not all([N == ind_N for N in sample_ind_N_list]):
            print("Not all sample individual numbers where the same in cohort "+\
                str(cohort_key[0])+' '+str(cohort_key[1])+' '+str(cohort_key[2])+':\n'+\
                str(sample_ind_N_list))
        allele_N = 2*ind_N
        obs_var = np.var(obs_freqs,ddof=1)
        set_freq = np.mean(obs_freqs)

        print(str(cohort_key[0])+' '+str(cohort_key[1])+' '+str(cohort_key[2]))

        obs_vars += [obs_var]
        set_freqs += [set_freq]
        allele_N_list += [allele_N]
        read_N_sets += [read_N_list]         

    avg_obs_var = np.mean(obs_vars)
    print("Observed average variance:\n"+str(avg_obs_var))

    lower_est = 1e-4 #ideally we never have an estimate where n < 2 (or even n < 10)
    lower_avg_var = avg_var_among_resamples(set_freqs,allele_N_list,read_N_sets,num_reps,lower_est)
    upper_est = 1
    upper_avg_var = avg_var_among_resamples(set_freqs,allele_N_list,read_N_sets,num_reps,upper_est)
    print('{0:8.7f}\t\t{1:8.7f}\t\t{2:8.7f}\t\t{3:8.7f}'.format(lower_est,lower_avg_var,upper_avg_var,upper_est))
    while abs(lower_est-upper_est) > delta:
        # Calculate the slope
        slope = (upper_avg_var-lower_avg_var)/(upper_est-lower_est)

        # Guess halfway to the linear intercept between the two, in both directions
        guess_upper_est = upper_est + 0.5*(avg_obs_var-upper_avg_var)/slope
        guess_upper_var = avg_var_among_resamples(set_freqs,allele_N_list,read_N_sets,num_reps,guess_upper_est)

        guess_lower_est = lower_est + 0.5*(avg_obs_var-lower_avg_var)/slope
        guess_lower_var = avg_var_among_resamples(set_freqs,allele_N_list,read_N_sets,num_reps,guess_lower_est)
        
        # Replace one or both of the bounds with the guess
        if guess_upper_var > avg_obs_var and guess_lower_var > avg_obs_var:
            lower_est = guess_upper_est
            lower_avg_var = guess_upper_var
        elif guess_upper_var < avg_obs_var and guess_lower_var < avg_obs_var:
            upper_est = guess_lower_est
            upper_avg_var = guess_lower_var
        elif guess_upper_var < avg_obs_var and guess_lower_var > avg_obs_var:
            lower_est = guess_lower_est
            lower_avg_var = guess_lower_var
            upper_est = guess_upper_est
            upper_avg_var = guess_upper_var
        else:
            print("Mismatched bounds in estimation:")
            print('{0:8.7f}\t\t{1:8.7f}\t\t{2:8.7f}\t\t{3:8.7f}'.format(lower_est,lower_avg_var,upper_avg_var,upper_est))
            break

        print('{0:8.7f}\t\t{1:8.7f}\t\t{2:8.7f}\t\t{3:8.7f}'.format(lower_est,lower_avg_var,upper_avg_var,upper_est))

    estimated_scaling_factor = np.mean([lower_est,upper_est])
    avg_var_from_est_scaling = avg_var_among_resamples(set_freqs,allele_N_list,read_N_sets,num_reps,estimated_scaling_factor)
    print("Estimated scaling and associated average resampled variance as:\n"+\
        str(estimated_scaling_factor)+'\t\t'+str(avg_var_from_est_scaling))

    return (estimated_scaling_factor,avg_var_from_est_scaling)


# Helper function to convert freq, read, count sets to broadcastable arrays
def convert_to_broadcastable(set_freqs,read_N_sets,inv_count_sets):

    set_freqs_temp = []
    read_N_temp = []
    inv_count_temp = []

    for i in np.arange(len(set_freqs)):
        read_N_temp += read_N_sets[i]
        inv_count_temp += inv_count_sets[i]
        set_freqs_temp += [set_freqs[i]]*len(read_N_sets[i])

    set_freqs_arr = np.array(set_freqs_temp)
    read_N_arr = np.array(read_N_temp)
    inv_count_arr = np.array(inv_count_temp)

    return(set_freqs_arr,read_N_arr,inv_count_arr)



# Returns a function that takes a scaling parameter,
#   and uses it as the sample size in a binomial sampling based on the inital freqs,
#   then a further resampling to the observed inv read count from the number of reads observed
def neg_log_likelihood_func_of_scaling_N(set_freqs,read_Ns,inv_counts,scaling_N):
    # print(set_freqs)
    # print(read_Ns)
    # print(inv_counts)
    # print(scaling_N)
    potential_counts = np.arange(1,scaling_N)
    # print(potential_counts)
    log_prob_scaling_count = binom.logpmf(potential_counts,scaling_N,set_freqs.reshape((len(set_freqs),1)))
    # print(log_prob_scaling_count)
    interm_freqs = np.divide(potential_counts,scaling_N)
    # print(interm_freqs)
    log_prob_inv_count = binom.logpmf(inv_counts.reshape((len(inv_counts),1)),read_Ns.reshape((len(read_Ns),1)),interm_freqs)
    # print(log_prob_inv_count)
    log_prob_both_counts = np.add(log_prob_scaling_count,log_prob_inv_count)
    # print(log_prob_both_counts)
    log_prob_each_draw = sp.special.logsumexp(log_prob_both_counts,axis=1)
    # print(log_prob_each_draw)
    total_log_prob = np.sum(log_prob_each_draw)
    print(str(scaling_N)+'\t\t'+str(-total_log_prob))
    return -total_log_prob

def ML_estimate_scaling_N(set_freqs_arr,read_N_arr,inv_count_arr):
    import scipy as sp
    from functools import partial

    (set_freqs_arr,read_N_arr,inv_count_arr) = convert_to_broadcastable(set_freqs,read_N_sets,inv_count_sets)

    nLL_func = partial(neg_log_likelihood_func_of_scaling_N,set_freqs_arr,read_N_arr,inv_count_arr)
    init_guess = 600
    result = sp.optimize.minimize(nLL_func,init_guess,method='nelder-mead',options={'disp': True})

    return result.x[0]

def ML_estimate_scaling_N_scan(set_freqs_arr,read_N_arr,inv_count_arr,min_par=2,max_par=5000):
    scaling_Ns = np.arange(min_par,max_par+1)
    vect_nLL = np.vectorize(partial(neg_log_likelihood_func_of_scaling_N,set_freqs_arr,read_N_arr,inv_count_arr))
    nLLs = vect_nLL(scaling_Ns)
    return(scaling_Ns,nLLs)




def ML_coestimate_cohort_scale_param(cohort,control_count_data,method='MLE_scan'):
    using_known_freqs = isinstance(cohort,dict)
    if using_known_freqs:
        print("\n\nCoestimating scaling factor using known expected haplotype frequencies:")
    else:
        print("\n\nCoestimating scaling factor using experimentally observed haplotype frequencies:")

    set_freqs = []
    read_N_sets = []
    inv_count_sets = []
    # obs_vars = []

    for extraction in cohort:
        obs_freqs = []
        read_N_list = []
        inv_count_list = []


        for run_data in control_count_data[extraction]:
            obs_freqs += [run_data[1]]
            read_N_list += [run_data[2]]
            inv_count_list += [run_data[3]]

        # obs_var = np.var(obs_freqs,ddof=1)

        set_freq = None
        if using_known_freqs:
            set_freq = cohort[extraction]
        else:
            set_freq = np.mean(obs_freqs)

        print(str(extraction[0])+' '+str(extraction[1])+' '+str(extraction[2]))
        if using_known_freqs:
            print("expected frequency of inversion: "+str(set_freq)+"\n")

        # obs_vars += [obs_var]
        set_freqs += [set_freq]
        read_N_sets += [read_N_list]
        inv_count_sets += [inv_count_list]

    # avg_obs_var = np.mean(obs_vars)
    # print("Observed average variance:\n"+str(avg_obs_var))

    (set_freqs_arr,read_N_arr,inv_count_arr) = convert_to_broadcastable(set_freqs,read_N_sets,inv_count_sets)
    if method == 'MLE_scipy_minimize':
        coestimated_scaling_parameter = ML_estimate_scaling_N(set_freqs_arr,read_N_arr,inv_count_arr)
    elif method == 'MLE_scan':
        (scaling_Ns,nLLs) = ML_estimate_scaling_N_scan(set_freqs_arr,read_N_arr,inv_count_arr)

        do_plotting = False
        if do_plotting:
            # Uses matplotlib.pyplot to plot a scaling factor scan
            import matplotlib.pyplot as plt

            # Make a simple plot
            fig, ax = plt.subplots()

            # ax.plot(scaling_Ns, nLLs, label="-log(likelihood)")
            ax.plot(scaling_Ns, nLLs)
            # ax.legend()

            # Scale the x-axis (scaling factor) to be exponential
            exp = lambda x: 10**(x)
            log = lambda x: np.log(x)
            ax.set_xscale('function', functions=(log, exp))
            ax.set(xlim=(2,5000))

            plt.show()

        coestimated_scaling_parameter = scaling_Ns[np.argmin(nLLs)]

    return coestimated_scaling_parameter


def ML_estimate_separately(cohort,control_count_data):

    return




# Returns the avg variance of num_reps new draws of read frequencies
#   of the given set of read counts based on the given freqency, allele count rescaling, and allele count
def avg_resample_dist_nonN_scale(set_freq,read_N,num_reps,scale_param):

    # Use numpy to quickly resample the reads
    cast_N = np.array(read_N_list).reshape((len(read_N_list),1))
    rescaled_cast_N = np.multiply(cast_N,set_scaling).astype('int64')
    cast_p = np.repeat(set_freq,num_reps).reshape((1,num_reps))
    sampled_counts = np.random.binomial(rescaled_cast_N,cast_p)
    sampled_freqs = np.divide(sampled_counts,rescaled_cast_N)
    sampled_vars = 0
    if bess_correction:
        sampled_vars = np.var(sampled_freqs,axis=0,ddof=1)
    else:
        sampled_vars = np.var(sampled_freqs,axis=0)
    mean_sampled_var = np.mean(sampled_vars)
    # print("Mean variance of resampled read frequencies:\n"+str(mean_sampled_var))


    # Use the analytic solution for the expectation of variance as num_reps -> INF
    #   = E(1/read_N)*set_freq*(1-set_freq)
    float_rescaled_Ns = rescaled_cast_N.astype('float64')
    mean_inverse_read_N = np.mean(np.power(float_rescaled_Ns,-1))
    analytic_test_var = mean_inverse_read_N*set_freq*(1-set_freq)
    # print("Analytically calculated expectation:\n"+str(analytic_test_var))
    # print("Ratio: "+str(mean_sampled_var/analytic_test_var))




# Returns a function that takes a scaling parameter,
#   and uses it as the sample size in a binomial sampling based on the inital freqs,
#   then a further resampling to the observed inv read count from the number of reads observed
def avg_freq_dist(set_freqs,read_Ns,inv_counts,scaling_N):

    # print(set_freqs)
    # print(read_Ns)
    # print(inv_counts)
    # print(scaling_N)
    potential_counts = np.arange(1,scaling_N)
    # print(potential_counts)
    log_prob_scaling_count = binom.logpmf(potential_counts,scaling_N,set_freqs.reshape((len(set_freqs),1)))
    # print(log_prob_scaling_count)
    interm_freqs = np.divide(potential_counts,scaling_N)
    # print(interm_freqs)
    log_prob_inv_count = binom.logpmf(inv_counts.reshape((len(inv_counts),1)),read_Ns.reshape((len(read_Ns),1)),interm_freqs)
    # print(log_prob_inv_count)
    log_prob_both_counts = np.add(log_prob_scaling_count,log_prob_inv_count)
    # print(log_prob_both_counts)
    log_prob_each_draw = sp.special.logsumexp(log_prob_both_counts,axis=1)
    # print(log_prob_each_draw)
    total_log_prob = np.sum(log_prob_each_draw)
    print(str(scaling_N)+'\t\t'+str(-total_log_prob))
    return -total_log_prob

def mean_diff_estimate_scaling_N_scan(set_freqs_arr,read_N_arr,inv_freq_arr,min_par=2,max_par=5000):
    scaling_Ns = np.arange(min_par,max_par+1)
    vect_nLL = np.vectorize(partial(neg_log_likelihood_func_of_scaling_N,set_freqs_arr,read_N_arr,inv_count_arr))
    nLLs = vect_nLL(scaling_Ns)
    return(scaling_Ns,nLLs)


def mean_diff_coestimate_cohort_scale_param(cohort,control_count_data,num_reps,method='scan'):
    using_known_freqs = isinstance(cohort,dict)
    if using_known_freqs:
        print("\n\nCoestimating scaling factor using known expected haplotype frequencies:")
    else:
        print("\n\nCoestimating scaling factor using experimentally observed haplotype frequencies:")

    set_freqs = []
    read_N_sets = []
    inv_freq_sets = []
    # obs_vars = []

    for extraction in cohort:
        obs_freqs = []
        read_N_list = []


        for run_data in control_count_data[extraction]:
            obs_freqs += [run_data[1]]
            read_N_list += [run_data[2]]

        # obs_var = np.var(obs_freqs,ddof=1)

        set_freq = None
        if using_known_freqs:
            set_freq = cohort[extraction]
        else:
            set_freq = np.mean(obs_freqs)

        print(str(extraction[0])+' '+str(extraction[1])+' '+str(extraction[2]))
        if using_known_freqs:
            print("expected frequency of inversion: "+str(set_freq)+"\n")

        # obs_vars += [obs_var]
        set_freqs += [set_freq]
        read_N_sets += [read_N_list]
        inv_freq_sets += [obs_freqs]

    # avg_obs_var = np.mean(obs_vars)
    # print("Observed average variance:\n"+str(avg_obs_var))

    (set_freqs_arr,read_N_arr,inv_freq_arr) = convert_to_broadcastable(set_freqs,read_N_sets,inv_freq_sets)
    if method == 'scan':
        (scaling_Ns,nLLs) = mean_diff_estimate_scaling_N_scan(set_freqs_arr,read_N_arr,inv_freq_arr)

        do_plotting = False
        if do_plotting:
            # Uses matplotlib.pyplot to plot a scaling factor scan
            import matplotlib.pyplot as plt

            # Make a simple plot
            fig, ax = plt.subplots()

            # ax.plot(scaling_Ns, nLLs, label="-log(likelihood)")
            ax.plot(scaling_Ns, nLLs)
            # ax.legend()

            # Scale the x-axis (scaling factor) to be exponential
            exp = lambda x: 10**(x)
            log = lambda x: np.log(x)
            ax.set_xscale('function', functions=(log, exp))
            ax.set(xlim=(2,5000))

            plt.show()

        coestimated_scaling_parameter = scaling_Ns[np.argmin(nLLs)]

    return coestimated_scaling_parameter






def write_all_cohort_scaling_scan(scaling_factor_scan_dict,scaling_factor_file):
    with open(scaling_factor_file,'w') as out_file:
        header = "Cohort,Chrom,Inv,Allele_N,Obs_Freq_Var,Scaling,Sampled_Variance\n"
        out_file.write(header)
        for cohort_key in scaling_factor_scan_dict:
            scan_data = scaling_factor_scan_dict[cohort_key]
            line_prefix = str(cohort_key[0])+','+str(cohort_key[1])+','+str(cohort_key[2])+','+\
                str(scan_data[4])+','+str(scan_data[3])+','
            for i in np.arange(len(scan_data[0])):
                line = line_prefix+str(scan_data[0][i])+','+str(scan_data[1][i])+'\n'
                out_file.write(line)
    return

def write_cohort_scaling_estimates(scaling_factor_dict,scaling_factor_file):
    with open(scaling_factor_file,'w') as out_file:
        header = "Cohort,Chrom,Inv,Allele_N,Obs_Freq_Var,Estimated_Scaling,Sampled_Variance\n"
        out_file.write(header)
        for cohort_key in scaling_factor_dict:
            line = str(cohort_key[0])+','+str(cohort_key[1])+','+str(cohort_key[2])+','+\
                str(scaling_factor_dict[cohort_key][3])+','+str(scaling_factor_dict[cohort_key][2])+','+\
                str(scaling_factor_dict[cohort_key][0])+','+str(scaling_factor_dict[cohort_key][1])+'\n'
            out_file.write(line)
    return

# Main function, for running the pipeline when the script is called
def main():

    # Parse arguments
    args = parse_args()

    # Assign input arguments to variables
    N_reps = args.n
    called_read_count_file = args.count_file
    library_metadata_file_path = args.library_data
    output_prefix = args.out_prefix
    output_directory = args.od

    will_scan = args.scan
    read_rescale_var_est = False

    # Ensure the output directory exists and the string is sanitized
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    if output_directory[-1] != '/':
        output_directory += '/'

    # Generate input data
    control_count_data = parse_control_count_data_to_dict(called_read_count_file)
    expected_freq_dict = {
        ('Cntrl 2','2L','Adult'):0.5,
        ('Cntrl 2','2R','Adult'):0.25,
        ('Cntrl 2','3R','Adult'):0.5,
        ('Cntrl 1','2L','Adult'):0.975,
        ('Cntrl 6','2L','Adult'):0.1,
        ('Cntrl 3','2L','Adult'):0.375
    }
    expected_2L_freq_dict = {
        ('Cntrl 2','2L','Adult'):0.5,
        ('Cntrl 1','2L','Adult'):0.975,
        ('Cntrl 6','2L','Adult'):0.1,
        ('Cntrl 3','2L','Adult'):0.375
    }


    if will_scan:
        # Generate variance and difference from observed variance values for each cohort based on observed means
        scaling_factors=np.power(10,np.arange(-4.5,0,.1))
        scaling_factor_scan_dict = scan_scaling_factors_by_extraction(control_count_data,N_reps,scaling_factors)

        # Write the scan data to file
        data_file = str(N_reps)+'rep_scaling_scan.csv'
        if output_prefix == '':
            data_file = output_directory+data_file
        else:
            data_file = output_directory+output_prefix+'_'+data_file
        write_all_cohort_scaling_scan(scaling_factor_scan_dict,data_file)

        # Plot the scan data
        plot_all_cohort_scaling_factor_scan_diffs_together(scaling_factor_scan_dict)


    if read_rescale_var_est:
        # Generate estimates for each cohort
        scaling_factor_dict = estimate_scaling_factors_by_extraction(control_count_data,N_reps)
        print(scaling_factor_dict)

        # Write the estimated optimal scaling factor data to file
        data_file = str(N_reps)+'rep_scaling_est.csv'
        if output_prefix == '':
            data_file = output_directory+data_file
        else:
            data_file = output_directory+output_prefix+'_'+data_file
        write_cohort_scaling_estimates(scaling_factor_dict,data_file)

        # Coestimate a specific set of cohorts
        cohort_list = [('Cntrl 2','2L','Adult'),('Cntrl 1','2L','Adult'),('Cntrl 6','2L','Adult'),('Cntrl 3','2L','Adult')]
        (scale,avg_var) = coestimate_scaling_factors(cohort_list,control_count_data,N_reps)
        # print("Coestimated scaling: ")


    coest = ML_coestimate_cohort_scale_param(expected_2L_freq_dict,control_count_data)
    print(str(coest))

    return


# Run when this script is called directly
if __name__ == '__main__':
    main()


