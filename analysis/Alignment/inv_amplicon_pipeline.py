
# Purpose: For generating estimates of inversion frequency in pools of amplicon sequences

# Notes: 
## 

# Run with the following command:
## python inv_amplicon_pipeline.py --read_meta <line names file> --aln_script <path to the alignment/qc pipeline file> --d_mel_ref <reference genome fasta with path> --rd <read dir> --id <alignment intermediate file directory> --ad <aligned reads dir> --pd <picards dir> --gd <GATK dir> --amp_ref_d <amplicon reference sequence etc dir> --od <output dir for frequency estimations>
## python inv_amplicon_pipeline.py --read_meta read_metadata.csv --aln_script amplicon_alignment_pipeline.sh --d_mel_ref ../../alignment_software/dmel_ref/DmelRef.fasta --rd reads/ --id intermediates/ --ad aligned_bam/ --pd ../../alignment_software/picard-tools-1.79/picard-tools-1.79/ --gd ../../alignment_software/GenomeAnalysisTK-3.4-46/ --amp_ref_d amp_refs/ --od inv_freqs/
# Import libraries
import argparse
import os
import subprocess

# Define functions

# Define and collect input arguments
def parse_args():
	parser = argparse.ArgumentParser()

	parser.add_argument('--read_meta', required=True, type=str,
						help='Provide the path to a .csv with the names and metadata of sequenced read pools')
	parser.add_argument('--aln_script', type=str,
						help='Provide the path to the alignment pipeline script to be used in analysis')
	parser.add_argument('--d_mel_ref', type=str,
						help='Provide the path to the D. mel genomic reference file to be used in analysis')
	parser.add_argument('--rd', type=str,
						help='Provide the path to the sequencing read file directory')
	parser.add_argument('--id', type=str,
						help='Provide the path to the alignment intermediate file directory')
	parser.add_argument('--ad', type=str,
						help='Provide the path to the alignment output file directory')
	parser.add_argument('--hd', type=str,
						help='Provide the path to the aligned amplicon haplotypes file directory')
	parser.add_argument('--pd', type=str,
						help='Provide the path to the Picard directory')
	parser.add_argument('--gd', type=str,
						help='Provide the path to the GATK directory')
	parser.add_argument('--amp_ref_d', type=str,
						help='Provide the path to the amplicon reference sequence directory')
	parser.add_argument('--od', type=str,
						help='Provide the path to the output directory for inversion frequency data')
	parser.add_argument('--from_aln_bam', dest='pipe_stage', default=0,
						action='store_const', const=1,
						help='Proceed to calculate inversions frequencies from already aligned bam files')
	parser.add_argument('--from_read_haps', dest='pipe_stage', default=0,
						action='store_const', const=2,
						help='Proceed to calculate inversions frequencies from extracted read haplotypes')
	parser.add_argument('--aln_only', dest='pipe_stage', default=0,
						action='store_const', const=3,
						help='Perform alignment and QC on raw paired read fasta files, then stop')
	parser.add_argument('--count_reads_only', dest='pipe_stage', default=0,
						action='store_const', const=4,
						help='Perform read haplotype extraction on aligned bam files, then stop')
	parser.add_argument('--use_info_haps', dest='freq_calc_mode', default=0,
						action='store_const', const=2,
						help='Perform frequency calculation using informative position haplotypes. \
						Default is to use fixed difference SNPs')
	parser.add_argument('--use_haps', dest='freq_calc_mode', default=0,
						action='store_const', const=1,
						help='Perform frequency calculation using complete reference haplotypes. \
						Default is to use fixed difference SNPs')

	args = parser.parse_args()
	return args



# Extracts data about the paired-end read pools to be used in analysis
# Includes file names to analyze
def parse_read_metadata(read_metadata_file_path):
	read_metadata = []
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
			
			# paired_read_file_list = []
			# i = 15
			# while i < len(data):
			# 	paired_read_file_list += [(data[i],data[i+1])]
			# 	i+=2

			inv_name = "In("+chrom+")"+inv
			std_name = "In("+chrom+")Std"

			read_metadata += [[replicate,age_stage,chrom,inv,inv_name,std_name,
				amp_name,library_ID,sequencing_name,UDI_ID,i5_barcode,i7_barcode_rc,
				sample_size,DNA_concentration,paired_read_file_1,paired_read_file_2,amp_start,amp_end]]
	print('\nRead pool data to analyze:')
	for pool in read_metadata:
		print(pool)
	print('\n')
	return read_metadata

# For running the alignment script for a single read pool
def align_read_pool(read_alignment_pipeline_script,pool_name,paired_read_file_1,
		paired_read_file_2,read_directory,d_mel_ref,picard_directory,
		GATK_directory,alignment_intermediates_directory,aligned_output_directory):
	args = ['./'+read_alignment_pipeline_script,pool_name,paired_read_file_1,
		paired_read_file_2,read_directory,d_mel_ref,picard_directory,
		GATK_directory,alignment_intermediates_directory,aligned_output_directory]
	print("Trying alignment pipeline command:\n"+'./'+read_alignment_pipeline_script+' '+pool_name+' '+\
		paired_read_file_1+' '+paired_read_file_2+' '+read_directory+' '+d_mel_ref,picard_directory+' '+\
		GATK_directory+' '+alignment_intermediates_directory+' '+aligned_output_directory+'\n')
	try:
		subprocess.check_call(args)
	except subprocess.CalledProcessError as e:
		print('Alignment pipeline for '+pool_name+' threw an exception:')
		print(e.message+'\n')
	return

# For founding a multiprocessing pool to asynchronously call the alignment script on all read pools,
#   and to return the unique read (haplotype) count data in correct order
def align_read_pools_async_pool(read_metadata,read_alignment_pipeline_script,d_mel_ref,
		read_directory,alignment_intermediates_directory,aligned_output_directory,
		picard_directory,GATK_directory):

	import multiprocessing as mp

	process_returns = []

	# Prepare the multiprocessing pool, leave at least 2?
	pool = mp.Pool(mp.cpu_count()-2)

	# Collect and write the unique read counts with apply_async()
	async_objects = []
	print('Performing read QC and alignment, in asynchronous parallel:')
	for metadata in read_metadata:
		pool_name = metadata[8]
		paired_read_file_1 = metadata[14]
		paired_read_file_2 = metadata[15]
		async_objects += [pool.apply_async(align_read_pool,args=(read_alignment_pipeline_script,pool_name,
										paired_read_file_1,paired_read_file_2,read_directory,d_mel_ref,picard_directory,
										GATK_directory,alignment_intermediates_directory,aligned_output_directory))]

	# Wait for the results with .get() -- could use .wait instead, and a for loop, so as to avoid a None list output
	process_returns = [result.get() for result in async_objects]

	# Close the pool and return the read haplotype sets
	pool.close()
	pool.join()
	print('\n')
	return



# Calls the alignment script on all read pools
def align_read_pools(read_metadata,read_alignment_pipeline_script,d_mel_ref,
		read_directory,alignment_intermediates_directory,aligned_output_directory,
		picard_directory,GATK_directory):
	for metadata in read_metadata:
		pool_name = metadata[8]
		paired_read_file_1 = metadata[14]
		paired_read_file_2 = metadata[15]
		align_read_pool(read_alignment_pipeline_script,pool_name,paired_read_file_1,
			paired_read_file_2,read_directory,d_mel_ref,picard_directory,
			GATK_directory,alignment_intermediates_directory,aligned_output_directory)
	return




# For reading a tab separated sequencing read haplotype count data
#    with a single line header and one read pair and count datum per line
def read_unique_read_dict_file(file_name):
	read_haps = {}
	with open(file_name,'r') as seq_file:
		lines = seq_file.readlines()
		for line in lines[1:]:
			data = line.strip().split('\t')
			r1 = data[0].strip()
			r2 = data[1].strip()
			count = int(data[2].strip())
			read_haps[(r1,r2)] = count
	return read_haps

# For reading in read pool unique read count data
def read_sep_unique_read_files(read_metadata, out_dir):
	read_hap_sets = []
	for metadata in read_metadata:
		pool_name = metadata[8]

		file_name = out_dir+pool_name+'_haps.txt'

		print('Reading read pool '+pool_name+\
			' paired read haplotype (unique sequence) count data:')
		print(file_name)

		read_hap_sets += [read_unique_read_dict_file(file_name)]
	return read_hap_sets


# For reading a tab separated sequence count data file
#    with a single line header and one sequence and count datum per line
def read_tab_dict_file(file_name):
	seq_dict = {}
	with open(file_name,'r') as seq_file:
		lines = seq_file.readlines()
		for line in lines[1:]:
			data = line.strip().split('\t')
			seq = data[0].strip()
			count = int(data[1].strip())
			seq_dict[seq] = count
	return seq_dict

def read_sep_unique_read_amp_files(read_metadata, out_dir):
	read_amp_hap_sets = []
	for metadata in read_metadata:
		pool_name = metadata[8]

		file_name = out_dir+pool_name+'_amp_haps.txt'

		print('Reading read pool '+pool_name+\
			' paired read haplotype (unique sequence) count data:')
		print(file_name)

		read_amp_hap_sets += [read_tab_dict_file(file_name)]
	return read_amp_hap_sets


# For reading a dictionary with a single line header and one string element per file line
def read_hap_file(file_name):
	hap_dict = {}
	with open(file_name,'r') as seq_file:
		lines = seq_file.readlines()
		for line in lines[1:]:
			# print(line)
			# print([l.strip() for l in line.split('\t')])
			(seq,count) = [l.strip() for l in line.split('\t')]
			hap_dict[seq] = int(count)
	return hap_dict


# Compares the read haplotypes to similarly-sezed slices of the amplicon haps 
def calc_inv_freqs_from_full_amp_haps(inv_file,std_file,read_haps,read_length):
	inv_haps = read_hap_file(inv_file)
	std_haps = read_hap_file(std_file)

	# Generate the comparable 'read' slices from the amplicon haps
	inv_paired_haps = []
	std_paired_haps = []
	# print(inv_name+' haps')
	for seq in inv_haps:
		# print(seq[:read_length])
		# print(seq[-read_length:]+'\n')
		inv_paired_haps += [(seq[:read_length],seq[-read_length:])]
	# print(std_name+' haps')
	for seq in std_haps:
		# print(seq[:read_length])
		# print(seq[-read_length:]+'\n')
		std_paired_haps += [(seq[:read_length],seq[-read_length:])]

	# Perform some quality control on the reference read sequences to be compared
	shared_read_1 = []
	shared_read_2 = []
	shared_pairs = []
	for (i_r_1,i_r_2) in inv_paired_haps:
		for (s_r_1,s_r_2) in std_paired_haps:
			if i_r_1 == s_r_1:
				if not i_r_1 in shared_read_1:
					shared_read_1 += [i_r_1]
				# else:
				# 	print('Multiply-shared first read: ')
				# 	print(i_r_1)
			if i_r_2 == s_r_2:
				if not i_r_2 in shared_read_2:
					shared_read_2 += [i_r_2]
				# else:
				# 	print('Multiply-shared second read: ')
				# 	print(i_r_2)
			if i_r_1 == s_r_1 and i_r_2 == s_r_2:
				if not (i_r_1,i_r_2) in shared_pairs:
					shared_pairs += [i_r_2]
					print('WARNING - shared read pair: ')
					print(i_r_1)
					print(i_r_2)
				else:
					print('WARNING - multiply-shared read pair: ')
					print(i_r_1)
					print(i_r_2)

	# Compare the reads to reference slices, to estimate a count of the # of inversions sequenced
	inv_count = 0
	std_count = 0
	shared_count = 0
	unmapped_count = 0
	for (read_1,read_2) in read_haps:
		if (read_1,read_2) in inv_paired_haps and (read_1,read_2) in std_paired_haps:
			shared_count += read_haps[(read_1,read_2)]
		elif (read_1,read_2) in inv_paired_haps:
			inv_count += read_haps[(read_1,read_2)]
		elif (read_1,read_2) in std_paired_haps:
			std_count += read_haps[(read_1,read_2)]
		else:
			unmapped_count += read_haps[(read_1,read_2)]
			# print('Unmapped read pair:')
			# print(read_1)
			# print(read_2+'\n')

	# Estimate pool frequency from the counts
	inv_freq = inv_count/(inv_count+std_count)

	return (inv_freq,inv_count,std_count,shared_count,unmapped_count)

# For reading a file with a first line of informative positions, 
def read_info_hap_file(file_name):
	haps = set()
	with open(file_name,'r') as seq_file:
		lines = seq_file.readlines()
		pos_line = lines[0].strip()[1:-1]
		positions = [int(pos.strip()) for pos in pos_line.split(',')]
		for line in lines[1:]:
			haps.add(line.strip())
	return (positions,haps)

# Combine the reads into a single longer haplotype
def combine_reads(read_haps,amp_length):
	print("Amp length: "+str(amp_length))
	combined_reads = {}
	num_overlap_mismatch = 0
	num_overlap = 0
	for (read_1,read_2) in read_haps:
		l_1 = len(read_1)
		l_2 = len(read_2)
		coverage_len = l_1+l_2
		if coverage_len >= amp_length:
			overlap_offset_1 = amp_length - l_2
			overlap_offset_2 = l_1 - overlap_offset_1

			overlap_1 = read_1[overlap_offset_1:]
			overlap_2 = read_2[:overlap_offset_2]

			# print("Amp_L\tR1_L\tR2_L\tR1_O\tR2_O"+str(amp_length))
			# print(str(amp_length)+"\t"+str(l_1)+"\t"+str(l_2)+"\t"+str(overlap_offset_1)+\
			# 	"\t"+str(overlap_offset_2)+"\n")

			overlap_seg = ''
			num_overlap += 1
			mismatch = False
			for i in range(len(overlap_1)):
				if overlap_1[i] == overlap_2[i]:
					overlap_seg += overlap_1[i]
				else:
					mismatch = True
					overlap_seg += 'N'
			if mismatch:
				num_overlap_mismatch += 1
				# print("Overlap mismatch detected:")
				# print(read_1)
				# print(overlap_1)
				# print(overlap_2)
				# print(read_2)

			combination = read_1[:overlap_offset_1] + overlap_seg + read_2[overlap_offset_2:]
			if len(combination) != amp_length:
				print("Read combination length failure - "+str(len(combination))+" not "+str(amp_length))
			combined_reads[combination] = read_haps[(read_1,read_2)]

		else:
			combination = read_1+'N'*(amp_length-coverage_len)+read_2
			if len(combination) != amp_length:
				print("Read combination length failure - "+str(len(combination))+" not "+str(amp_length))
			combined_reads[combination] = read_haps[(read_1,read_2)]

	# Print a brief summary of mismatched overlapping segments
	if num_overlap:
		print(str(num_overlap_mismatch/num_overlap)+" of overlapping segments in read pairs mismatch: "+str(num_overlap_mismatch)+" in "+str(num_overlap))
	else:
		print("No overlapping segments in read pairs")

	return combined_reads



# Estimate a count of the # of inversions sequenced and their frequency in the pool,
#   by comparing the read pair to reference sequences at informative positions
def calc_inv_freqs_from_info_amp_haps(inv_file,std_file,read_haps,amp_length):
	info_pos,inv_haps = read_info_hap_file(inv_file)
	info_pos,std_haps = read_info_hap_file(std_file)

	# Combine the separate paired reads into a single read sequence
	combined_reads = combine_reads(read_haps,amp_length)


	# Prepare to count the reads based on informative reference positions 
	inv_count = 0
	std_count = 0
	shared_count = 0
	unmapped_count = 0

	# read_print_goal = 40
	# read_checks_printed = 0

	# Check the read span against the informative position data for the reference amplicon sequences
	for read in combined_reads:
		info_read = ''
		removed_pos = []
		# Filter out N's in the reduced read as uninformative
		for i in range(len(info_pos)):
			if read[info_pos[i]] == 'N':
				removed_pos += [True]
			else:
				info_read += read[info_pos[i]]
				removed_pos += [False]
		test_inv_haps = set(inv_haps)
		test_std_haps = set(std_haps)
		# Filter those same positions out of the haplotypes to compare
		if any(removed_pos):
			reduced_inv_haps = set()
			for hap in inv_haps:
				reduced_hap = ''
				for i in range(len(info_pos)):
					if not removed_pos[i]:
						reduced_hap += hap[i]
				reduced_inv_haps.add(reduced_hap)
			test_inv_haps = reduced_inv_haps
			reduced_std_haps = set()
			for hap in std_haps:
				reduced_hap = ''
				for i in range(len(info_pos)):
					if not removed_pos[i]:
						reduced_hap += hap[i]
				reduced_std_haps.add(reduced_hap)
			test_std_haps = reduced_std_haps
		# Estimate the assignment of the read based on these informative positions
		if info_read in test_inv_haps and info_read in test_std_haps:
			shared_count += combined_reads[read]
		elif info_read in test_inv_haps:
			inv_count += combined_reads[read]
		elif info_read in test_std_haps:
			std_count += combined_reads[read]
		else:
			unmapped_count += combined_reads[read]
			# print('Unmapped read pair:')
			# print(info_read)
			# print(next(iter(test_inv_haps)))
			# print(next(iter(test_std_haps)))
			# print(read)
		# if read_checks_printed < read_print_goal:
		# 	read_checks_printed += 1
		# 	print('Read comparison check '+str(read_checks_printed)+':')
		# 	print(info_read)
		# 	print(next(iter(test_inv_haps)))
		# 	print(next(iter(test_std_haps)))
		# 	print(read)

	# Estimate pool frequency from the counts
	inv_freq = inv_count/(inv_count+std_count)

	return (inv_freq,inv_count,std_count,shared_count,unmapped_count)




# For reading a tab-separated file with one fixed difference per line and a header line,
#   of the form:  amplicon_offset \t genomic_position \t inv_ATGCN_count_list \t std_ATGCN_count_list
def read_fixed_diff_file(file_name):
	diffs = []
	with open(file_name,'r') as seq_file:
		lines = seq_file.readlines()
		for line in lines[1:]:
			data = [d.strip() for d in line.strip().split('\t')]
			offset = int(data[0])
			genom_pos = int(data[1])
			inv_base_dist = [int(b.strip()) for b in data[2][1:-1].split(',')]
			std_base_dist = [int(b.strip()) for b in data[3][1:-1].split(',')]
			inv_bases = str(data[4])
			std_bases = str(data[5])
			diffs += [(offset, genom_pos, inv_base_dist, std_base_dist, inv_bases, std_bases)]
	return diffs


# Calculate the inversion frequencies based solely on fixed differences
def calc_inv_freqs_from_fixed_diffs(fixed_diff_file,read_amp_haps,amp_length):

	# Extract and process the fixed differences
	fixed_diffs = read_fixed_diff_file(fixed_diff_file)

	# Combine the separate paired reads into a single read sequence
	# combined_reads = combine_reads(read_haps,amp_length)read_amp_haps
	combined_reads = read_amp_haps

	print(len(next(iter(combined_reads))))

	# Prepare to count the reads based on fixed differences
	inv_count = 0
	std_count = 0
	recombinant_count = 0
	mismapped_count = 0
	multiply_mapped_count = 0

	read_print_goal = 40
	read_checks_printed = 0

	for read in combined_reads:

		# Check the read span against the reference fixed difference data
		full_match_inv = True
		full_match_std = True
		mismatch = False
		for diff in fixed_diffs:
			match_inv = False
			match_std = False
			offset = diff[0]
			inv_bases = diff[4]
			std_bases = diff[5]
			for base in inv_bases:
				if read[offset] == base:
					match_inv = True
			for base in std_bases:
				if read[offset] == base:
					match_std = True
			full_match_inv = full_match_inv and match_inv
			full_match_std = full_match_std and match_std
			if not match_inv and not match_std:
				mismatch = True


		# Assign the read based on matching to fixed differences
		if full_match_inv and full_match_std:
			print("WARNING - reads should NEVER fully match both sets of fixed differences. Coding error likely.")
			# For fixed differences, this should never happen, as the bases should never be shared
			multiply_mapped_count += combined_reads[read]
		elif full_match_inv:
			inv_count += combined_reads[read]
		elif full_match_std:
			std_count += combined_reads[read]
		elif mismatch:
			# Likely something went wrong - the base observed wasn't expected in either karyotype
			mismapped_count += combined_reads[read]
			# print('WARNING - Unmatched bases at fixed diff sites:')
			# read_bases = ''
			# for diff in fixed_diffs:
			# 	read_bases += read[diff[0]]
			# 	diff_line = str(diff[0])
			# 	for d in diff[1:]:
			# 		diff_line+='\t'+str(d)
			# 	print(diff_line)
			# print("Read bases at offsets: "+read_bases)
			# print(read)
		else:
			# If there are no mismatches, but neither fully matched, it is recombinant:
			recombinant_count += combined_reads[read]

		# if read_checks_printed < read_print_goal:
		# 	read_checks_printed += 1
		# 	print('Read comparison check '+str(read_checks_printed)+':')
		# 	read_bases = ''
		# 	for diff in fixed_diffs:
		# 		read_bases += read[diff[0]]
		# 		diff_line = str(diff[0])
		# 		for d in diff[1:]:
		# 			diff_line+='\t'+str(d)
		# 		print(diff_line)
		# 	print("Read bases at offsets: "+read_bases)
		# 	print(read)

	# Estimate pool frequency from the counts
	inv_freq = inv_count/(inv_count+std_count)

	return (inv_freq,inv_count,std_count,recombinant_count,mismapped_count,multiply_mapped_count)


# Takes the aligned .bam files and calculates the 
def calculate_inversion_frequencies(read_metadata,read_amp_hap_sets,amplicon_ref_seq_directory,
		read_length,mode):

	inv_freqs = []
	for i in range(len(read_metadata)):
		pool_name = read_metadata[i][8]
		print('Calculating inversion frequencies for '+pool_name)

		# Generate a panel of expected reference reads from the amplicon haplotype data
		inv_name = read_metadata[i][4]
		std_name = read_metadata[i][5]
		amp_name = read_metadata[i][6]

		if mode == 0:
			amp_fixed_diff_file = amplicon_ref_seq_directory+amp_name+'_fixdif.txt'

			amp_start = read_metadata[i][16]
			amp_end = read_metadata[i][17]
			amp_length = amp_end - amp_start

			inv_freqs += [calc_inv_freqs_from_fixed_diffs(amp_fixed_diff_file,read_amp_hap_sets[i],amp_length)]

		if mode == 1:
			amp_inv_haps_file = amplicon_ref_seq_directory+amp_name+'_'+inv_name+'_haps.txt'
			amp_std_haps_file = amplicon_ref_seq_directory+amp_name+'_'+std_name+'_haps.txt'

			inv_freqs += [calc_inv_freqs_from_full_amp_haps(amp_inv_haps_file,
				amp_std_haps_file,read_amp_hap_sets[i],read_length)]

		if mode == 2:
			amp_inv_info_haps_file = amplicon_ref_seq_directory+amp_name+'_'+inv_name+'_info_haps.txt'
			amp_std_info_haps_file = amplicon_ref_seq_directory+amp_name+'_'+std_name+'_info_haps.txt'

			amp_start = read_metadata[i][16]
			amp_end = read_metadata[i][17]
			amp_length = amp_end - amp_start

			inv_freqs += [calc_inv_freqs_from_info_amp_haps(amp_inv_info_haps_file,
				amp_std_info_haps_file,read_amp_hap_sets[i],amp_length)]

	return inv_freqs

# Simply print the calculated inversion frequencies to stdout
def print_inv_freqs(read_metadata,inv_freqs):
	print('\n\nLibrary\t\tInv_Freq\t\tStd_Freq\t\tN_Inv\tN_Std\tN_Recom\tN_Mismap\tMismap_F\tN_Multmap\n')
	for i in range(len(read_metadata)):
		library_name = read_metadata[i][8]
		inv_freq = inv_freqs[i][0]
		inv_count = inv_freqs[i][1]
		std_count = inv_freqs[i][2]
		recombinant_count = inv_freqs[i][3]
		mismapped_count = inv_freqs[i][4]
		multiply_mapped_count = inv_freqs[i][5]
		prop_fail = mismapped_count/(inv_count+std_count+recombinant_count+mismapped_count+multiply_mapped_count)
		percent_fail = 100*prop_fail
		print(library_name+'\t'+str(inv_freq)+'\t'+str(1-inv_freq)+\
			'\t'+str(inv_count)+'\t'+str(std_count)+'\t'+str(recombinant_count)+\
			'\t'+str(mismapped_count)+'\t'+str(prop_fail)+\
			'\t'+str(multiply_mapped_count)+'\n')
	return


# For writing a csv file containing some read pool metadata
#    as well as the inversion karyotype frequency estimates and called haplotype counts
def write_inv_freq_file(read_metadata,inv_freqs,output_directory,postfix=''):
	file_name = output_directory+"freq_data"+postfix+".csv"
	header = "Library,Line,Cohort,Chrom,Inv,Sample_Ind_N,Inv_F,Std_F,N_Inv,N_Std,N_Recomb,N_Mismap,N_Multmap\n"
	with open(file_name,'w') as out_file:
		out_file.write(header)
		for i in range(len(read_metadata)):
			sequencing_name = read_metadata[i][8]
			line = read_metadata[i][0]
			cohort = read_metadata[i][1]
			chrom = read_metadata[i][2]
			inv = read_metadata[i][3]
			ind_N = read_metadata[i][12]
			inv_freq = inv_freqs[i][0]
			std_freq = 1-inv_freq
			inv_count = inv_freqs[i][1]
			std_count = inv_freqs[i][2]
			recomb_count = inv_freqs[i][3]
			mismap_count = inv_freqs[i][4]
			mult_map_count = inv_freqs[i][5]
			line = sequencing_name+','+line+','+cohort+','+chrom+','+inv+','+str(ind_N)+','+\
				str(inv_freq)+','+str(std_freq)+','+str(inv_count)+','+str(std_count)+','+\
				str(recomb_count)+','+str(mismap_count)+','+str(mult_map_count)+\
				'\n'
			out_file.write(line)
	return


# Main function, for running the pipeline when the script is called
def main():

	# Parse arguments
	args = parse_args()

	# Assign input arguments to variables
	read_metadata_file_path = args.read_meta
	d_mel_ref = args.d_mel_ref
	read_alignment_pipeline_script = args.aln_script
	read_directory = args.rd
	alignment_intermediates_directory = args.id
	aligned_output_directory = args.ad
	hap_output_directory = args.hd
	picard_directory = args.pd
	GATK_directory = args.gd
	amplicon_ref_seq_directory = args.amp_ref_d
	output_directory = args.od

	mode_of_action = args.pipe_stage
	perform_alignment = mode_of_action == 0 or mode_of_action == 3
	perform_read_pair_count = mode_of_action <= 1 or mode_of_action == 4
	perform_freq_calc = mode_of_action <= 2

	mode_of_freq_calc = args.freq_calc_mode

	# Perform input value and directory checks
	assert (os.path.exists(read_metadata_file_path)),'Library read pool metadata file path (--read_meta) {} \
		does not reference a file'.format(read_metadata_file_path)

	if perform_alignment:
		# Check if the required directories and files exist
		assert (os.path.exists(read_alignment_pipeline_script)),'Alignment pipeline script (--aln_script) {} \
			does not reference a file'.format(read_alignment_pipeline_script)
		assert (os.path.exists(d_mel_ref)),'D. mel genomic reference file (--d_mel_ref) {} \
			does not reference a file'.format(d_mel_ref)
		assert (os.path.exists(read_directory)),'Raw paired read sequence directory (--rd) {} \
			does not reference a directory'.format(read_directory)
		assert (os.path.exists(picard_directory)),'Picard tools directory (--pd) {} \
			does not reference a directory'.format(picard_directory)
		assert (os.path.exists(GATK_directory)),'GATK tools directory (--gd) {} \
			does not reference a directory'.format(GATK_directory)

		# Make aligned output and intermediate directories if they don't exist yet
		assert (alignment_intermediates_directory!=''),'Alignment pipeline intermediate working directory \
			(--id) not given'
		if not os.path.exists(alignment_intermediates_directory):
			os.makedirs(alignment_intermediates_directory)
		assert (aligned_output_directory!=''),'Aligned paired read bam directory (--ad) \
			not given'
		if not os.path.exists(aligned_output_directory):
			os.makedirs(aligned_output_directory)

		# Sanitize the directory inputs
		if read_directory[-1] != '/':
			read_directory += '/'
		if alignment_intermediates_directory[-1] != '/':
			alignment_intermediates_directory += '/'
		if aligned_output_directory[-1] != '/':
			aligned_output_directory += '/'
		if picard_directory[-1] != '/':
			picard_directory += '/'
		if GATK_directory[-1] != '/':
			GATK_directory += '/'

	if perform_read_pair_count:
		assert (hap_output_directory!=''),'Aligned amplicon haplotype and haplotype count output directory (--hd) \
			not given'
		if not os.path.exists(hap_output_directory):
			os.makedirs(hap_output_directory)
		if hap_output_directory[-1] != '/':
			hap_output_directory += '/'

		if not perform_alignment:
			assert (os.path.exists(aligned_output_directory)),'Aligned paired read bam directory (--ad) {} \
				does not reference a directory'.format(aligned_output_directory)
			if aligned_output_directory[-1] != '/':
				aligned_output_directory += '/'

	if perform_freq_calc or perform_read_pair_count:
		assert (os.path.exists(amplicon_ref_seq_directory)),'Amplicon reference sequence directory (--amp_ref_d) {} \
			does not reference a directory'.format(amplicon_ref_seq_directory)
		if amplicon_ref_seq_directory[-1] != '/':
			amplicon_ref_seq_directory += '/'

	if perform_freq_calc:
		assert (output_directory!=''),'Estimated inversion frequency and read counts output \
			directory (--od) not given'
		if not os.path.exists(output_directory):
			os.makedirs(output_directory)
		if output_directory[-1] != '/':
			output_directory += '/'

		if not perform_read_pair_count:
			assert (os.path.exists(hap_output_directory)),'Aligned amplicon haplotype and haplotype count \
				output directory (--hd) {} does not reference a directory'.format(output_directory)
			if hap_output_directory[-1] != '/':
				hap_output_directory += '/'


	read_metadata = parse_read_metadata(read_metadata_file_path)

	if perform_alignment:
		align_read_pools_async_pool(read_metadata,read_alignment_pipeline_script,d_mel_ref,
			read_directory,alignment_intermediates_directory,aligned_output_directory,
			picard_directory,GATK_directory)
		# align_read_pools(read_metadata,read_alignment_pipeline_script,d_mel_ref,
		# 	read_directory,alignment_intermediates_directory,aligned_output_directory,
		# 	picard_directory,GATK_directory)

	if perform_read_pair_count:
		# Import the asynchronous extraction method from extrct_read_haps
		from extract_read_haps import extract_and_write_async_pool

		read_amp_hap_sets = extract_and_write_async_pool(read_metadata,
			aligned_output_directory,hap_output_directory,amplicon_ref_seq_directory)
		# read_hap_sets = extract_read_haps(read_metadata,aligned_output_directory)
		# write_sep_read_hap_files(read_hap_sets,read_metadata,output_directory)

	if perform_freq_calc:
		if not perform_read_pair_count:
			read_amp_hap_sets = read_sep_unique_read_amp_files(read_metadata,hap_output_directory)

		read_length = 151 #!!!!!! may need to modify
		inv_freqs = calculate_inversion_frequencies(read_metadata,read_amp_hap_sets,amplicon_ref_seq_directory,
			read_length,mode_of_freq_calc)

		print_inv_freqs(read_metadata,inv_freqs)

		write_inv_freq_file(read_metadata,inv_freqs,output_directory)

	return


# Run when this script is called directly
if __name__ == '__main__':
	main()


