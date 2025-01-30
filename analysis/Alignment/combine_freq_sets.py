
# For generating a combined read count file from several runs of the same library cohort data

# Import packages
import numpy as np
import argparse
import os
# from inv_amplicon_pipeline import write_inv_freq_file, print_inv_freqs


# Define and collect input arguments
def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('--count_file', type=str,
                        help='The .csv file containing the estimated inverted and standard read counts, etc.')
                        # inv_freqs/freq_data.csv
    # parser.add_argument('--library_data', type=str,
    #                     help='The .csv file with the library metadata')
                        # read_metadata.csv
    parser.add_argument('--od', type=str,
                        help='The output directory into which to write combined count/freq data')
                        # bootstrap/
    parser.add_argument('--out_postfix', type=str, default='',
                        help='The output file postfix to use')

    args = parser.parse_args()
    return args



# Extracts the called read count and frequency data
#   returns a dictionary with library name as key and the metadata,
#   and sets of adjusted proportion and counts as value
def parse_count_data_to_dict_list(count_file_path):
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

            count_data[(line,chrom,cohort)] = count_data.get((line,chrom,cohort),[]) + [[inv,sample_ind_N,
                num_inv,num_std,num_recomb,num_mismap,num_multimap]]
    return count_data

# A helper dictionary to make combined read set names
cohort_dict = {
    'Parental ♂':'par',
    'Embryo':'emb',
    'Early ♀':'eaf',
    'Late ♀':'laf',
    'Early ♂':'eam',
    'Late ♂':'lam'
}

def make_comb_cohort_name(line,chrom,cohort):
    name = line+cohort_dict[cohort]+chrom
    return name

def combine_count_data(count_data):
    comb_count_data = []
    for (line,chrom,cohort) in count_data:
        replicate_data = count_data[(line,chrom,cohort)]
        comb_cohort_name = line+cohort_dict[cohort]+chrom
        inv = replicate_data[0][0]
        ind_N = replicate_data[0][1]
        total_inv_count = replicate_data[0][2]
        total_std_count = replicate_data[0][3]
        total_recomb_count = replicate_data[0][4]
        total_mismap_count = replicate_data[0][5]
        total_multimap_count = replicate_data[0][6]
        if len(replicate_data) > 1:
            for rep in replicate_data[1:]:
                total_inv_count += rep[2]
                total_std_count += rep[3]
                total_recomb_count += rep[4]
                total_mismap_count += rep[5]
                total_multimap_count += rep[6]
        inv_prop = total_inv_count/(total_std_count+total_inv_count)
        std_prop = 1-inv_prop
        cohort_data = [comb_cohort_name,line,cohort,chrom,inv,ind_N,inv_prop,std_prop,
            total_inv_count,total_std_count,total_recomb_count,total_mismap_count,total_multimap_count]
        comb_count_data += [cohort_data]
    return comb_count_data


# For writing a csv file containing the minimal metadata and read counts of
#    a combined set of library sequencings
def rewrite_inv_freq_file(comb_count_data,output_directory,postfix=''):
    file_name = output_directory+"freq_data"+postfix+".csv"
    header = "Library,Line,Cohort,Chrom,Inv,Sample_Ind_N,Inv_F,Std_F,N_Inv,N_Std,N_Recomb,N_Mismap,N_Multmap\n"
    with open(file_name,'w') as out_file:
        out_file.write(header)
        for i in range(len(comb_count_data)):
            line = ''
            for datum in comb_count_data[i]:
                line += str(datum)+','
            line = line[:-1]+'\n'
            out_file.write(line)
    return


# Main function, for running the pipeline when the script is called
def main():

    # Parse arguments
    args = parse_args()

    # Assign input arguments to variables
    called_read_count_file = args.count_file
    output_postfix = args.out_postfix
    output_directory = args.od

    # Ensure the output directory exists and the string is sanitized
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    if output_directory[-1] != '/':
        output_directory += '/'

    # library_metadata = parse_library_metadata_to_dict(library_metadata_file_path)
    count_data = parse_count_data_to_dict_list(called_read_count_file)
    comb_count_data = combine_count_data(count_data)
    rewrite_inv_freq_file(comb_count_data,output_directory,output_postfix)

    return


# Run when this script is called directly
if __name__ == '__main__':
    main()


