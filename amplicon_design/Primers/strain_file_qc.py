

# A script for general troubleshooting and quality control of input genomic sequences
#   for inv_primer_gen.py

import inv_primer_gen as ipg
import numpy as np
import argparse
import os


# To try and find the cause of universal masking in the count files
def strain_mask_counts(genome_dir,strains,chrom,inv_name,out_dir):

	# Calculate the chromosome length
	first_strain = strains[0]
	# file_name = genome_dir+'Chr'+chrom+'/'+first_strain+'_Chr'+chrom+'.fas1k'
	prefix = genome_dir + 'Chr' + chrom + '/' + first_strain + '_Chr' + chrom
	file_name = ipg.pick_strain_file(prefix)
	first_file = open(file_name,'r')
	lines = first_file.readlines()
	first_file.close()
	line_len = len(lines[0].strip())
	num_lines = len(lines)
	remainder_len = len(lines[num_lines-1].strip())
	chrom_len = line_len*(num_lines-1)+remainder_len
	print(inv_name+' has total length '+str(chrom_len)+' bases')


	mask_counts = []
	mask_percents = []
	print('Counting masking for '+inv_name)
	for strain in strains:
		prefix = genome_dir+'Chr'+chrom+'/'+strain+'_Chr'+chrom
		file_name = ipg.pick_strain_file(prefix)
		strain_file = open(file_name,'r')
		lines = strain_file.read()
		# print(lines)
		strain_file.close()
		mask_count = lines.count('N')
		mask_percent = float(mask_count)/chrom_len
		print('Strain '+strain+':\t'+str(mask_count)+'\t'+str(mask_percent))
		mask_counts += [mask_count]
		mask_percents += [mask_percent]
	return((mask_counts,mask_percents))

# Generates the site-base counts array for each inversion
# Assumes only two haplotypes, the specific inversion and standard
def gen_strain_mask_counts(genome_dir,inv_data,inv_strain_dict,out_dir):
	mask_counts = []
	mask_percents = []
	strain_sets = []
	chroms = []
	for inv_datum in inv_data:
		(chrom,inv,br_1,br_2,inv_name,std_name) = inv_datum
		# Inverted
		strains = inv_strain_dict[inv_name]
		strain_sets += [strains]
		(counts,percents) = strain_mask_counts(genome_dir,
								strains,chrom,inv_name,out_dir)
		mask_counts += [counts]
		mask_percents += [percents]
		chroms += [chrom]
		# Standard
		strains = inv_strain_dict[std_name]
		strain_sets += [strains]
		(counts,percents) = strain_mask_counts(genome_dir,
								strains,chrom,std_name,out_dir)
		mask_counts += [counts]
		mask_percents += [percents]
		chroms += [chrom]
	return((mask_counts,mask_percents,strain_sets,chroms))

# Prints the set of strains with masking over the filter_threshold
def filter_strain_masking(filter_threshold,mask_counts,mask_percents,strain_sets,chroms):
	print('Top masked strains:')
	strain_mask_dict = {}
	strain_chrom_dict = {}
	filtered_strains = []
	for i in np.arange(len(mask_counts)):
		for j in np.arange(len(mask_counts[i])):
			percent_masked = mask_percents[i][j]
			strain = strain_sets[i][j]
			if strain in strain_mask_dict:
				if strain_mask_dict[strain] < percent_masked:
					strain_mask_dict[strain] = percent_masked
					strain_chrom_dict[strain] = chroms[i]
			else:
				strain_mask_dict[strain] = percent_masked
				strain_chrom_dict[strain] = chroms[i]
			if percent_masked > filter_threshold:
				# print(strain+':\t'+str(percent_masked))
				filtered_strains += [strain]
	unique_strains = sorted(list(set(filtered_strains)))
	unique_strains.sort(key = lambda x: strain_mask_dict[x], reverse = True)
	print("Sorted list of over-masked strains")
	print("Strain\tMax Masked\t\tChrom")
	for strain in unique_strains:
		print(strain+'\t'+str(strain_mask_dict[strain])+\
				'\t'+str(strain_chrom_dict[strain]))
	# print(unique_strains)
	return

# For generating a dictionary that returns the largest masking percent over the chromosomes
#   for each strain
def gen_strain_mask_dict(mask_counts,mask_percents,strain_sets):
	strain_mask_dict = {}
	for i in np.arange(len(mask_counts)):
		for j in np.arange(len(mask_counts[i])):
			percent_masked = mask_percents[i][j]
			strain = strain_sets[i][j]
			if strain in strain_mask_dict:
				if strain_mask_dict[strain] < percent_masked:
					strain_mask_dict[strain] = percent_masked
			else:
				strain_mask_dict[strain] = percent_masked
	return(strain_mask_dict)

# Uses the filter_threshold to remove strains from the original isp file in the new output
def write_filtered_isp(in_isp_file,out_isp_file,filter_threshold,strain_mask_dict):
	with open(in_isp_file,'r') as in_isp:
		with open(out_isp_file,'w') as out_isp:
			lines = in_isp.readlines()
			out_isp.write(lines[0])
			print('Filtering out:')
			for line in lines[1:]:
				in_strain = line.split(',')[0]
				if strain_mask_dict[in_strain] < filter_threshold:
					out_isp.write(line)
				else:
					print(in_strain)
	return



# Checks for bases that are present in > threshold haplotypes in one pool,
#   but less than threshold haplotypes in the other, then returns the 
# Helper function for check_potential_misclassification
def check_orphan_SNP(pos_1,pos_2,):
	if (pos_1[4] != 0) or (pos_2[4] != 0):
		return 3
	else:
		shared_bases = np.logical_and(pos_1[0:4],pos_2[0:4])
		num_shared = np.sum(shared_bases)
		base_presence = np.logical_or(pos_1[0:4],pos_2[0:4])
		num_present = np.sum(base_presence)
		if num_present == 1:
			return 0
		elif num_shared == 0:
			return 1
		# elif num_shared > 1 or num_present > num_shared:
		else:
			return 2

# Return an array with base counts by position for one line in the strain file set
# Helper function for check_potential_misclassification
def count_line(strains,):
	counts = 

# Access chromosomes simulataneously to count the instances of 
#   position frequency file, used for inversion-specific strain lists
def check_potential_misclassification(genome_dir,inv_strains,std_strains,
		chrom,inv_name,std_name,out_dir):
	# Check the first file for length data
	# SHOULDN'T SIMULTANEOUSLY OPEN THE WHOLE THING :(
	print('Starting '+inv_name+' parameter check')
	first_strain = strains[0]
	print('Strain used: '+first_strain)
	prefix = genome_dir + 'Chr' + chrom + '/' + first_strain + '_Chr' + chrom
	file_name = pick_strain_file(prefix)
	first_file = open(file_name,'r')
	lines = first_file.readlines()
	first_file.close()
	# print(lines[0].strip())
	line_len = len(lines[0].strip())
	# print(line_len)
	num_lines = len(lines)
	# print(num_lines)
	# print(lines[num_lines-1].strip())
	remainder_len = len(lines[num_lines-1].strip())
	# print(remainder_len)
	chrom_len = line_len*(num_lines-1)+remainder_len
	# print(chrom_len)
	print("Chrom Len: "+str(chrom_len))


	# Open all of the strain files and start the potential mislabel dict
	pot_mislabel = {}
	print('Starting '+inv_name+' strain count checks')
	inv_strain_files = []
	print('Inv:')
	for strain in inv_strains:
		print(strain)
		prefix = genome_dir + 'Chr' + chrom + '/' + strain + '_Chr' + chrom
		file_name = ipg.pick_strain_file(prefix)
		inv_strain_files += [open(file_name,'r')]
		pot_mislabel[strain] = 0
	std_strain_files = []
	print('Std:')
	for strain in std_strains:
		print(strain)
		prefix = genome_dir + 'Chr' + chrom + '/' + strain + '_Chr' + chrom
		file_name = ipg.pick_strain_file(prefix)
		std_strain_files += [open(file_name,'r')]
		pot_mislabel[strain] = 0

	# Begin counting bases for inv/std pools
	i_start = 0
	for l in np.arange(num_lines-1):
		inv_counts_line = np.zeros((line_len,4))
		inv_seqs = {}
		for strain_file in inv_strain_files:
			seq = list(map(ipg.base_dict.get,strain_file.readline().strip()))
			if seq in inv_seqs:
				inv_seqs[seq]
			for i in np.arange(line_len):
				inv_counts_line[i,seq[i]] += 1
		std_counts_line = np.zeros((line_len,4))
		for strain_file in std_strain_files:
			seq = list(map(ipg.base_dict.get,strain_file.readline().strip()))
			for i in np.arange(line_len):
				std_counts_line[i,seq[i]] += 1
		# Check for sites that are off from shared status by 1 or 2 bases
		for i in np.arange(line_len): 
			if check_fixed_poly_base(pos_count_1[i,:],pos_count_2[i,:],m_tol)





	# Finish the final line
	for strain_file in strain_files:
		seq = list(map(base_dict.get,strain_file.readline().strip()))
		for i in np.arange(remainder_len):
			freq_dat[i_start+i,seq[i]] += 1
	freq_dat.flush()
	# Close the strain files
	for strain_file in strain_files:
		strain_file.close()
	# Return the memory-reduced array
	return(freq_dat)

# Generates the site-base counts array for each inversion
# Assumes only two haplotypes, the specific inversion and standard
def gen_inv_site_freqs(genome_dir,inv_data,inv_strain_dict,out_dir):
	inv_site_freqs = []
	for inv_datum in inv_data:
		(chrom,inv,br_1,br_2,inv_name,std_name) = inv_datum
		strains = inv_strain_dict[inv_name]
		inv_site_freqs += [concurrent_pos_freq_gen(genome_dir,
								strains,chrom,inv_name,out_dir)]
		strains = inv_strain_dict[std_name]
		inv_site_freqs += [concurrent_pos_freq_gen(genome_dir,
								strains,chrom,std_name,out_dir)]
	return inv_site_freqs


# For handling direct calls to InvPrimerGen.py,
#   and to allow specifying different input files
# Assumes only two haplotypes, the specific inversion and standard
def main(args):
	strain_file_path = args.s
	inversion_breakpoint_file_path = args.i
	barcode_file_path = args.b
	parameter_file_path = args.p
	primer3_core_directory = args.p3
	genome_dir = args.d
	out_dir = args.od
	out_isp = args.os
	filter_threshold = args.ft

	assert (strain_file_path!=''),'Strain segregation file path (--s) {} \
		not given'.format(strain_file_path)
	assert (inversion_breakpoint_file_path!=''),'Region of interest file path (--i) {} \
		not given'.format(inversion_breakpoint_file_path)
	assert (barcode_file_path!=''),'Barcode sequences file path (--b) {} \
		not given'.format(barcode_file_path)
	assert (parameter_file_path!=''),'Paramter file path (--p) {} \
		not given'.format(parameter_file_path)

	if out_dir[-1] != '/':
		out_dir += '/'
	if primer3_core_directory[-1] != '/':
		primer3_core_directory += '/'
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)
	if not os.path.exists(out_dir+'data/'):
		os.makedirs(out_dir+'data/')
	if not os.path.exists(out_dir+'primer3/'):
		os.makedirs(out_dir+'primer3/')

	inv_strain_dict = ipg.parse_inv_strain_presence(strain_file_path)
	inv_data = ipg.parse_inv_data(inversion_breakpoint_file_path)
	# (illum_A,illum_B,prim_len,min_gap,max_gap,buff) = ipg.parse_parameters(parameter_file_path)

	(mask_counts,mask_percents,strain_sets,chroms) = gen_strain_mask_counts(genome_dir,
													inv_data,inv_strain_dict,out_dir)

	# filter_strain_masking(filter_threshold,mask_counts,mask_percents,strain_sets,chroms)

	strain_mask_dict = gen_strain_mask_dict(mask_counts,mask_percents,strain_sets)

	write_filtered_isp(strain_file_path,out_isp,filter_threshold,strain_mask_dict)

	return







# Runs if InvPrimerGen.py is the main script
if __name__ == "__main__":
	# Collect input arguments
	# CHECK argparse documentation for other input ideas (positional arguments? other?)
	parser = argparse.ArgumentParser(
		description=__doc__,
		formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument('--s',
						help='the input csv list of Strain haploid genomes and inversion presences',
						type=str,
						default='inv_strain_presence.csv')
	parser.add_argument('--i',
						help='the input csv list of Inversions and their breakpoints',
						type=str,
						default='inv_breaks.csv')
	parser.add_argument('--b',
						help='the Barcode list to pull from for FREQ-seq primers',
						type=str,
						default='barcodes.csv')
	parser.add_argument('--p',
						help='the Parameters:\n\
						the illumina adaptor A to append to a forward primer,\n\
						the Illumina adaptor B to append to a reverse primer,\n\
						the primer length,\n\
						the gap length between primers\n\
						the buffer distance from breakpoints to ignore',
						type=str,
						default='parameters.txt')
	parser.add_argument('--d',
						help='the Directory in which the genomes can be found',
						type=str,
						default='')
						# '/Volumes/4_TB_RAID_Set_1/DPGP2plus/wrap1kb/'
	parser.add_argument('--p3',
						help='the Directory in which primer3_core can be found: \
						generally ending with primer3/src/; defaults to primer3/src/',
						type=str,
						default='primer3/src/')
						# '/Volumes/4_TB_RAID_Set_1/DPGP2plus/wrap1kb/'
	parser.add_argument('--od',
						help='candidate region and position frequency files output directory',
						type=str,
						default='primer_output/')
	parser.add_argument('--ft',
						help='the percent masking threshold above which to flag strains',
						type=float,
						default=0.20)
	parser.add_argument('--os',
						help='output file for filtered strain inversion presence file',
						type=str,
						default='filtered_isp.csv')

	args = parser.parse_args()

	main(args)

