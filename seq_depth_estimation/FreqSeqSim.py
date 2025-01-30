
# For simulation of inversion frequency sequencing,
# using counts as a readout of genotype, instead of modeling sequence

# 
import numpy as np

# Number of bootstrap reps per read depth combination
N_samples = 10000
freq_gradation = 0.05
max_coverage = 5
read_gradation = .01

N_emb = 2000
N_pat = 600
N_adu = 600

file_name = "bootstrap_freq_2.csv"








# Simultaneous simulation and writing to file of Freq-seq Drosophila bootstraps
def sim_write_two_allele_bootstrap(N_samples, freq_gradation, max_coverage, read_gradation,
    N_pat, N_emb, N_adu, file_name):
    N_read_trials = int(max_coverage//read_gradation)
    N_freq_trials = int(1.0//freq_gradation) - 1
    # Prepare to write a line per simulation
    with open(file_name,'w') as output:
        output.write("Pat_Exp_Freq, Coverage, Pat_Count, Emb_Count, "\
            +"Adu_Count, Pat_Freq, Emb_Freq, Adu_Freq\n")
        for f in np.arange(N_freq_trials):
            for c in np.arange(N_read_trials):
                for i in np.arange(N_samples):
                    # Calculate the correct total reads for each allele sampling depth
                    temp_cover = (c+1)*read_gradation
                    r_p = temp_cover*N_pat
                    r_e = temp_cover*N_emb
                    r_a = temp_cover*N_adu
                    # Calculate the correct expected frequencies
                    temp_freq = freq_gradation + f*freq_gradation
                    pat_null = [temp_freq,1-temp_freq]
                    off_null = [temp_freq/2,1-temp_freq/2]
                    # Resample paternal
                    pat_reads = np.random.multinomial(r_p,pat_null)
                    # Resample embryo offspring
                    emb_alleles = np.random.multinomial(N_emb,off_null)
                    e_allele_freqs = emb_alleles/sum(emb_alleles)
                    emb_reads = np.random.multinomial(r_e,e_allele_freqs)
                    # Resample adult offspring
                    adu_alleles = np.random.multinomial(N_adu,off_null)
                    a_allele_freqs = adu_alleles/sum(adu_alleles)
                    adu_reads = np.random.multinomial(r_a,a_allele_freqs)
                    # Write the output
                    output.write(str(temp_freq)+", "+str(temp_cover)+", "\
                        +str(pat_reads[0])+", "+str(emb_reads[0])+", "\
                        +str(adu_reads[0])+", "+str((pat_reads/r_p)[0])+", "\
                        +str((emb_reads/r_e)[0])+", "+str((adu_reads/r_a)[0])+"\n")
    return

# Running the basic bootstrap
sim_write_two_allele_bootstrap(N_samples, freq_gradation, max_coverage, 
    read_gradation, N_pat, N_emb, N_adu, file_name)





# Scratch graveyard

# # Takes a genotype set (alleles or reads) and resamples for bootstrapping
# # uses a multinomial based on the sample frequencies
# def resample(counts,N):
#     sampledCounts = [0]*len(counts)
#     pvals = counts/sum(counts)
#     return np.random.multinomial(N,pvals)

# def two_allele_bootstrap(N_samples, freq_gradation, max_coverage, read_gradation,
#     N_pat, N_emb, N_adu):
#     N_read_trials = int(max_coverage//read_gradation)
#     N_freq_trials = int(1.0//freq_gradation) - 1
#     # print(N_read_trials)
#     # print(N_samples)
#     # first three of six elements are the counts, second three are the frequencies
#     samples = np.zeros((N_freq_trials,N_read_trials,N_samples,6,2))
#     # freqs = np.zeros((N_read_trials,N_samples,3,2))
#     freqs = np.zeros((N_freq_trials))
#     coverage = np.zeros((N_read_trials))
#     # print (counts)
#     # print (coverage)
#     for f in np.arange(N_freq_trials):
#         for c in np.arange(N_read_trials):
#             for i in np.arange(N_samples):
#                 # Calculate the correct total reads for each allele sampling depth
#                 temp_cover = (c+1)*read_gradation
#                 r_p = temp_cover*N_pat
#                 r_e = temp_cover*N_emb
#                 r_a = temp_cover*N_adu
#                 # Calculate the correct expected frequencies
#                 temp_freq = freq_gradation + f*freq_gradation
#                 pat_null = [temp_freq,1-temp_freq]
#                 off_null = [temp_freq/2,1-temp_freq/2]
#                 # Resample paternal
#                 pat_reads = np.random.multinomial(r_p,pat_null)
#                 # Resample embryo offspring
#                 emb_alleles = np.random.multinomial(N_emb,off_null)
#                 e_allele_freqs = emb_alleles/sum(emb_alleles)
#                 emb_reads = np.random.multinomial(r_e,e_allele_freqs)
#                 # Resample adult offspring
#                 adu_alleles = np.random.multinomial(N_adu,off_null)
#                 a_allele_freqs = adu_alleles/sum(adu_alleles)
#                 adu_reads = np.random.multinomial(r_a,a_allele_freqs)
#                 # Store the output
#                 samples[f,c,i,] = [pat_reads,emb_reads,adu_reads,pat_reads/r_p,emb_reads/r_e,adu_reads/r_a]
#                 # freqs[f,c,i,] = [pat_reads/r_p,emb_reads/r_e,adu_reads/r_a]
#                 freqs[f] = temp_freq
#                 coverage[c] = temp_cover
#     return (samples,freqs,coverage)



# def write_two_allele_bootstrap(samples, freqs, coverage, file_name):
#     # print(samples)
#     # print(freqs)
#     # print(coverage)
#     with open(file_name,'w') as output:
#         output.write("Pat_Exp_Freq, Coverage, Pat_Count, Emb_Count, "\
#             +"Adu_Count, Pat_Freq, Emb_Freq, Adu_Freq\n")
#         for f in np.arange(len(samples)):
#             for c in np.arange(len(samples[0])):
#                 for n in np.arange(len(samples[0,0])):
#                     output.write(str(freqs[f])+", "+str(coverage[c])+", "\
#                         +str(samples[f,c,n,0,0])+", "+str(samples[f,c,n,1,0])+", "\
#                         +str(samples[f,c,n,2,0])+", "+str(samples[f,c,n,3,0])+", "\
#                         +str(samples[f,c,n,4,0])+", "+str(samples[f,c,n,5,0])+"\n")
#     return



# # Running the basic bootstrap
# (samples,freqs,coverage) = two_allele_bootstrap(N_samples, freq_gradation, 
#     max_coverage, read_gradation, N_pat, N_emb, N_adu)

# write_two_allele_bootstrap(samples,freqs,coverage,file_name)