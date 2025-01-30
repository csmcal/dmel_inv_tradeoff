

# A script for generating primers targeting SNP-dense regions within inversions,
# as well as oligomers appending Illumina adaptors to the primers,
# for testing polymorphism and inversion frequency in pooled data
# following the form of Illumina metagenomic protocols.
#
# Originally run for Inversions in Zambian D. melanogaster, Pool Lab, 2018
#
# Uses np.memmap to try and reduce RAM costs, but could use pytables or h5py

import numpy as np
import argparse
import os


# A dictionary translating bases to integers in ATGCN
base_dict = {
	'A': 0,
	'T': 1,
	'G': 2,
	'C': 3,
	'N': 4
}

# A dictionary translating integers to bases in ATGCN
int_dict = {
	0: 'A',
	1: 'T',
	2: 'G',
	3: 'C',
	4: 'N'
}

# A dictionary translating poly integer bases to characters
poly_trans = {
	0: 'A',
	1: 'T',
	2: 'G',
	3: 'C',
	4: 'X',
	5: 'Y',
	6: 'N'
}


def parse_inv_strain_presence(strain_file_path):
	inv_strain_dict = {}
	with open(strain_file_path) as strain_file:
		lines = strain_file.readlines()
		inv_names = lines[0].strip().split(',')[1:]
		num_inv = len(inv_names)
		inv_names_std = []
		for inv_name in inv_names:
			(chrom,inv) = inv_name.split(')')
			inv_names_std += [inv_name]
			inv_names_std += [chrom+')Std']
		inv_strain_presences = [[] for i in np.arange(2*num_inv)]
		for line in lines[1:]:
			data = line.strip().split(',')
			strain_name = str(data[0])
			for i in np.arange(num_inv):
				if data[i+1] == 'T':
					inv_strain_presences[2*i] += [strain_name]
				else:
					inv_strain_presences[2*i+1] += [strain_name]
		inv_strain_dict = dict(zip(inv_names_std,inv_strain_presences))
	# print(inv_strain_dict)
	print("Preparing analysis of:")
	for i in np.arange(len(inv_names_std)):
		print(inv_names_std[i]+': '+str(len(inv_strain_presences[i]))+' strains')
	return inv_strain_dict
		

def parse_inv_data(breakpoint_file_path):
	inv_data = []
	with open(breakpoint_file_path) as breakpoint_file:
		lines = breakpoint_file.readlines()
		for line in lines[1:]:
			data = line.strip().split(',')
			chrom = str(data[0])
			inv = str(data[1])
			bp_1 = int(data[2])
			bp_2 = int(data[3])
			inv_name = "In("+chrom+")"+inv
			std_name = "In("+chrom+")Std"
			inv_data += [[chrom,inv,bp_1,bp_2,inv_name,std_name]]
	# print(inv_data)
	return inv_data

def parse_barcodes(barcode_file_path):
	barcodes = []
	with open(barcode_file_path) as barcode_file:
		for line in barcode_file:
			barcode = line.strip()
			barcodes += [barcode]
	# print(barcodes)
	return barcodes
		
# Parses the parameter file ordered as 'parameter name: value \n'
def parse_parameters(parameter_file_path):
	params = []
	with open(parameter_file_path) as parameter_file:
		lines = parameter_file.readlines()
		params = []
		for line in lines:
			params += [line.split(':')[1]]
		# Extract the stubby adaptor sequences
		params[0] = params[0].strip()
		params[1] = params[1].strip()
		# Cast the length parameters to int
		params[2] = int(params[2])
		params[3] = int(params[3])
		params[4] = int(params[4])
		params[5] = int(params[5])
		params[6] = int(params[6])
		params[7] = int(params[7])
		params[8] = int(params[8])
		# Extract the Illumina index adaptor sequences
		params[9] = params[9].strip()
		params[10] = params[10].strip()
		params[11] = params[11].strip()
		params[12] = params[12].strip()
	# print(params)
	return params

# To access the unmasked file, when IBD and admixture are present,
# with .unmask and .NoAdFilt preferred in that order
def pick_strain_file(prefix):
	unmasked_name = prefix + '.unmask'
	non_admix_name = prefix + '.NoAdFilt'
	default_name = prefix + '.fas1k'
	if os.path.isfile(unmasked_name):
		print('using .unmask')
		return(unmasked_name)
	elif os.path.isfile(non_admix_name):
		print('using .NoAdFilt')
		return(non_admix_name)
	else:
		return(default_name)


# Access chromosomes simulataneously to build a chromosome
# 	position frequency file, used for inversion-specific strain lists
def concurrent_pos_count_gen(genome_dir,strains,chrom,inv_name,out_dir):
	# Check the first file for length data
	# SHOULDN'T SIMULTANEOUSLY OPEN THE WHOLE THING :(
	print('Starting '+inv_name+' first strain check')
	first_strain = strains[0]
	print('Strain '+first_strain)
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
	# Create the position-count array
	dat_file = out_dir+'data/'+inv_name+'_counts.dat'
	count_dat = np.memmap(dat_file,mode='w+',shape=(chrom_len,5))
	# Add the first strain sequence to the site-base frequency array
	print('recording first strain counts')
	i = 0
	while i < num_lines-1:
		seq = list(map(base_dict.get,lines[i].strip()))
		# print(i)
		for j in np.arange(line_len):
			# print(i*line_len+j)
			count_dat[i*line_len+j,seq[j]] += 1
		i += 1
	seq = list(map(base_dict.get,lines[i].strip()))
	# print(i)
	for j in np.arange(remainder_len):
		count_dat[i*line_len+j,seq[j]] += 1
	count_dat.flush()
	# Open the rest of the strain files
	print('Starting '+inv_name+' general strains')
	strain_files = []
	for strain in strains[1:]:
		print('Strain '+strain)
		prefix = genome_dir + 'Chr' + chrom + '/' + strain + '_Chr' + chrom
		file_name = pick_strain_file(prefix)
		# file_name = genome_dir + 'Chr' + chrom + '/' + \
		# 				strain + '_Chr' + chrom + '.fas1k'
		strain_files += [open(file_name,'r')]
	i_start = 0
	for l in np.arange(num_lines-1): # faster to generate the ranges first?
		for strain_file in strain_files:
			seq = list(map(base_dict.get,strain_file.readline().strip()))
			for i in np.arange(line_len):
				count_dat[i_start+i,seq[i]] += 1

		count_dat.flush()
		i_start += line_len
	# Finish the final line
	for strain_file in strain_files:
		seq = list(map(base_dict.get,strain_file.readline().strip()))
		for i in np.arange(remainder_len):
			count_dat[i_start+i,seq[i]] += 1
	count_dat.flush()
	# Close the strain files
	for strain_file in strain_files:
		strain_file.close()
	# Return the memory-reduced array
	return(count_dat)

# Generates the site-base counts array for each inversion
# Assumes only two haplotypes, the specific inversion and standard
def gen_inv_site_counts(genome_dir,inv_data,inv_strain_dict,out_dir):
	inv_site_counts = []
	for inv_datum in inv_data:
		(chrom,inv,br_1,br_2,inv_name,std_name) = inv_datum
		strains = inv_strain_dict[inv_name]
		inv_site_counts += [concurrent_pos_count_gen(genome_dir,
								strains,chrom,inv_name,out_dir)]
		strains = inv_strain_dict[std_name]
		inv_site_counts += [concurrent_pos_count_gen(genome_dir,
								strains,chrom,std_name,out_dir)]
	return inv_site_counts

# For writing individual position count files
def write_pos_count_file(file_name,inv_pos_count):
	with open(file_name,'w') as output:
		for pos in inv_pos_count:
			line = str(pos[0])
			for base_count in pos[1:]:
				line += '\t'+str(base_count)
			line += '\n'
			output.write(line)
	return

# For writing the site-base counts in human readable format
def write_pos_counts(inv_pos_counts,inv_data,out_dir):
	for i in np.arange(len(inv_data)):
		(chrom,inv,br_1,br_2,inv_name,std_name) = inv_data[i]
		print('Writing '+inv_name+' human readable position count data')
		file_name = out_dir+inv_name+'_pos_count.txt'
		write_pos_count_file(file_name,inv_pos_counts[2*i])
		print('Writing '+std_name+' human readable position count data')
		file_name = out_dir+std_name+'_pos_count.txt'
		write_pos_count_file(file_name,inv_pos_counts[2*i+1])
	return

def parse_pos_count_file(file_name):
	pos_count = np.genfromtxt(file_name,dtype=np.uint8) #delimiter = '\t' is unnecessary
	return(pos_count)

def read_pos_counts(inv_data,out_dir):
	pos_counts = []
	for (chrom,inv,br_1,br_2,inv_name,std_name) in inv_data:
		print('Reading '+inv_name+' position-base counts')
		pos_counts += [parse_pos_count_file(out_dir+inv_name+'_pos_count.txt')]
		print('Reading '+std_name+' position-base counts')
		pos_counts += [parse_pos_count_file(out_dir+std_name+'_pos_count.txt')]
	return(pos_counts)

def parse_pos_count_mem_file(dat_file):
	pos_count_dat = np.memmap(dat_file,mode='r+')
	pos_count_dat = np.memmap.reshape(pos_count_dat,(-1,5))
	return(pos_count_dat)

def read_pos_count_mems(inv_data,out_dir):
	pos_counts = []
	for (chrom,inv,br_1,br_2,inv_name,std_name) in inv_data:
		print('Reading '+inv_name+' position-base counts dat file')
		pos_counts += [parse_pos_count_mem_file(out_dir+'data/'+inv_name+'_counts.dat')]
		print('Reading '+std_name+' position-base counts dat file')
		pos_counts += [parse_pos_count_mem_file(out_dir+'data/'+std_name+'_counts.dat')]
	return(pos_counts)



# Takes base counts at two positions
#  and gives a boolean if there are no shared bases
#  and not too many strains with masking at the position
def check_fixed(pos_1,pos_2):
	# if num_masked > m_tol:
	# 	return False
	# else:
	shared_bases = np.logical_and(pos_1[0:4],pos_2[0:4])
	num_shared = np.sum(shared_bases)
	base_presence = np.logical_or(pos_1[0:4],pos_2[0:4])
	num_present = np.sum(base_presence)
	if (num_present == 0) or (num_present == 1 and num_shared == 0):
		print('check: no recorded base counts after ignoring masking')
	elif num_shared == 0:
		return True
	else:
		return False

# Inefficient function to get the bases and their counts at fixed difference SNPs
def extract_bases(pos_counts):
	bases = []
	counts = []
	for i in np.arange(4):
		if pos_counts[i] > 0:
			bases += [int_dict[i]]
			counts += [pos_counts[i]]
	return bases, counts

# Returns a list object with every fixed difference SNP location for this inversion
def find_fixed_SNPs(pos_count_1,pos_count_2,mask_tolerance):
	chrom_len = pos_count_1.shape[0]
	total_strains = sum(pos_count_1[0])+sum(pos_count_2[0])
	# m_tol = total_strains * mask_tolerance
	m_tol = mask_tolerance
	SNP_list = []
	i = 0
	while i < chrom_len:
		pos_1 = pos_count_1[i,:]
		pos_2 = pos_count_2[i,:]
		num_masked = pos_1[4] + pos_2[4]
		if num_masked <= m_tol and check_fixed(pos_1,pos_2):
			print("\nFound difference at "+str(i)+":\n"+str(pos_1)+"\n"+str(pos_2)+"\n")
			inv_bases, inv_base_counts = extract_bases(pos_1[0:4])
			std_bases, std_base_counts = extract_bases(pos_2[0:4])
			SNP_list += [[i,inv_bases,std_bases,inv_base_counts,std_base_counts,num_masked]]
		i += 1
	return SNP_list


# Returns a list object with every fixed difference SNP location
#   for each inversion in inv_data, using the position count files
def gen_fixed_SNPs(pos_counts,inv_data,out_dir,mask_tolerance):
	fixed_SNPs = []
	for i in np.arange(len(inv_data)):
		(chrom,inv,br_1,br_2,inv_name,std_name) = inv_data[i]
		print("Generating fixed difference SNP set for "+str(inv_name))
		# dat_file = out_dir+'data/'+inv_name+'_poly.dat'
		fixed_SNPs += [find_fixed_SNPs(pos_counts[2*i],pos_counts[2*i+1],mask_tolerance)]
	return fixed_SNPs


# For writing individual position count files
def write_fixed_SNPs_file(file_name,SNP_list):
	with open(file_name,'w') as output:
		header = 'Position\tInv_Bases\tStd_Bases\tInv_Base_Count\tStd_Base_Count\tMasking_Count\n'
		output.write(header)
		for site in SNP_list:
			inv_bases_str = ''
			inv_base_count_str = ''
			for i in np.arange(len(site[1])):
				inv_bases_str += site[1][i]+','
				inv_base_count_str += str(site[3][i])+','
			inv_bases_str = inv_bases_str[:-1]
			inv_base_count_str = inv_base_count_str[:-1]

			std_bases_str = ''
			std_base_count_str = ''
			for i in np.arange(len(site[2])):
				std_bases_str += site[2][i]+','
				std_base_count_str += str(site[4][i])+','
			std_bases_str = std_bases_str[:-1]
			std_base_count_str = std_base_count_str[:-1]

			line = str(site[0])+'\t'+inv_bases_str+'\t'+std_bases_str+'\t'+\
				inv_base_count_str+'\t'+std_base_count_str+'\t'+str(site[5])+'\n'
			output.write(line)
	return

# For writing the site-base counts in human readable format
def write_fixed_SNPs(fixed_SNPs,inv_data,out_dir):
	for i in np.arange(len(inv_data)):
		(chrom,inv,br_1,br_2,inv_name,std_name) = inv_data[i]
		print('Writing '+inv_name+' fixed SNP position count data')
		file_name = out_dir+inv_name+'_fixed_SNPs.txt'
		write_fixed_SNPs_file(file_name,fixed_SNPs[i])
	return




# Translates base frequency counts at two positions into an integer,
# preserving the identity of fixed bases:
#   0: fixed in both for A
#   1: fixed in both for T
#   2: fixed in both for G
#   3: fixed in both for C
#   4: no shared bases
#   5: some shared bases
#   6: masking (N) present in either
def check_fixed_poly_base(pos_1,pos_2,m_tol):
	if (pos_1[4] + pos_2[4]) > m_tol:
		return 6
	else:
		shared_bases = np.logical_and(pos_1[0:4],pos_2[0:4])
		num_shared = np.sum(shared_bases)
		base_presence = np.logical_or(pos_1[0:4],pos_2[0:4])
		num_present = np.sum(base_presence)
		if (num_present == 0) or (num_present == 1 and num_shared == 0):
			print('check: no recorded base counts after ignoring masking')
		if num_present == 1:
			return np.nonzero(base_presence)[0][0]
		elif num_shared == 0:
			return 4
		else:
			return 5

# Translates base frequency counts at two positions into an integer:
#   0: fixed in both for the same base
#   1: no shared bases
#   2: some shared bases
#   3: masking (N) present in either
def check_fixed_poly(pos_1,pos_2):
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

# Returns a mapping of two polymorphism sets to 
#   shared fixed base at the position (0)
#   no shared bases at the position (1)
#	some shared bases at the position (2)
#   masking in either - an N (3)
def map_fixed_poly(pos_count_1,pos_count_2,dat_file,mask_tolerance):
	chrom_len = pos_count_1.shape[0]
	total_strains = sum(pos_count_1[0])+sum(pos_count_2[0])
	# m_tol = total_strains * mask_tolerance
	m_tol = mask_tolerance
	poly_map = np.memmap(dat_file,mode='w+',shape=(chrom_len)) # default dtype fine
	for i in np.arange(chrom_len):
		# poly_map[i] = check_fixed_poly(pos_count_1[i,:],pos_count_2[i,:])
		poly_map[i] = check_fixed_poly_base(pos_count_1[i,:],pos_count_2[i,:],m_tol)
	poly_map.flush()
	return poly_map


# Returns mappings for each inversion in a set,
#   with alternating inv and std in the pos_count array list
def gen_poly_maps(pos_counts,inv_data,out_dir,mask_tolerance):
	poly_maps = []
	for i in np.arange(len(inv_data)):
		(chrom,inv,br_1,br_2,inv_name,std_name) = inv_data[i]
		print("Generating polymorphism map for "+str(inv_name))
		dat_file = out_dir+'data/'+inv_name+'_poly.dat'
		poly_maps += [map_fixed_poly(pos_counts[2*i],pos_counts[2*i+1],
				dat_file,mask_tolerance)]
	return poly_maps

# For writing an individual poly map to fas1k file
def write_poly_map(file_name,poly_map):
	line_len = 1000
	chrom_len = poly_map.shape[0]
	reg_line_num = chrom_len//line_len
	remainder_len = chrom_len%line_len
	with open(file_name,'w') as output:
		for i in np.arange(reg_line_num):
			line = ''
			for j in np.arange(line_len):
				line += str(poly_map[i*line_len+j])
			line += '\n'
			output.write(line)
		line = ''
		for i in np.arange(remainder_len):
			line += str(poly_map[reg_line_num*line_len+i])
		line += '\n'
		output.write(line)
	return

# For writing polymorphism maps to human-readable .fas1k files
def write_poly_maps(poly_maps,inv_data,out_dir):
	for i in np.arange(len(inv_data)):
		(chrom,inv,br_1,br_2,inv_name,std_name) = inv_data[i]
		file_name = out_dir+inv_name+'_poly_map.fas1k'
		write_poly_map(file_name,poly_maps[i])
	return

def parse_poly_map(poly_file):
	# SHOULDN'T SIMULTANEOUSLY OPEN THE WHOLE THING :(
	map_file = open(poly_file,'r')
	lines = map_file.readlines()
	map_file.close()
	line_len = len(lines[0].strip())
	num_lines = len(lines)
	remainder_len = len(lines[num_lines-1].strip())
	chrom_len = line_len*(num_lines-1)+remainder_len
	# Create the polymorphism array
	poly_map = np.empty(chrom_len,dtype=np.uint8)
	# Add the data by line
	i = 0
	while i < num_lines-1:
		seq = lines[i].strip()
		for j in np.arange(line_len):
			poly_map[i*line_len+j] = int(seq[j])
		i += 1
	seq = lines[i].strip()
	for j in np.arange(remainder_len):
		poly_map[i*line_len+j] = int(seq[j])
	return(poly_map)

def read_poly_maps(inv_data,out_dir):
	poly_maps = []
	for (chrom,inv,br_1,br_2,inv_name,std_name) in inv_data:
		print('Reading '+inv_name+' polymorphism map')
		poly_maps += [parse_poly_map(out_dir+inv_name+'_poly.fas1k')]
	return(poly_maps)

def load_poly_map_mem(poly_dat):
	poly_map = np.memmap(poly_dat,mode='r+')
	return(poly_map)

def load_poly_map_mems(inv_data,out_dir):
	poly_maps = []
	for (chrom,inv,br_1,br_2,inv_name,std_name) in inv_data:
		print('Reading '+inv_name+' polymorphism map dat')
		poly_maps += [load_poly_map_mem(out_dir+'data/'+inv_name+'_poly.dat')]
	return(poly_maps)

# Maps the polymorphism data to regions:
#   0: Conserved regions with fixed shared sequence > min_prim_len
#   1: Regions between, containing no N masking and some fixed differences
#   2: Regions between, containing no N masking, but also no fixed differences
#   3: Regions between, containing at least one N masking position
# Gives output as a list of region identities, 
#   a list of start positions by genomic sequence
#   a list of number of fixed differences
def map_regions(poly_calls, inv_start, inv_end, min_prim_len, buffer_len, prim_buff):
	region_types = []
	region_starts = []
	num_fixed_diffs = []
	i = inv_start + buffer_len
	start = i
	cons_run = 0
	reg_type = 0
	num_diff = 0
	min_cons_reg = min_prim_len + prim_buff
	# Start the region scan
	while i < inv_end-buffer_len:
		poly_state = poly_calls[i]
		# If the base is conserved:
		# keep going unless the conserved run = min len
		if poly_state == 0:
			if reg_type != 0:
				cons_run += 1
				if cons_run == min_cons_reg:
					# print('starting conserved region')
					region_types += [reg_type]
					region_starts += [start]
					num_fixed_diffs += [num_diff]
					start = i-cons_run+1
					reg_type = 0
					num_diff = 0
					cons_run = 0
		# If the base is fixed different:
		elif poly_state == 1:
			if reg_type == 0:
				# print('starting difference region')
				region_types += [reg_type]
				region_starts += [start]
				num_fixed_diffs += [num_diff] # or += [0]
				start = i
				reg_type = 1
				num_diff = 1
			# elif reg_type == 1: # Redundant to final else
			# 	num_diff += 1
			elif reg_type == 2:
				# print('switching to difference region')
				reg_type = 1
				num_diff = 1
			else: # Might want to allow N's in amplicon? Consider
				num_diff += 1
			cons_run = 0
		# If the base is partially shared:
		# ignore it unless breaking a conserved region
		elif poly_state == 2:
			if reg_type == 0:
				# print('starting polymorphic region')
				region_types += [reg_type]
				region_starts += [start]
				num_fixed_diffs += [num_diff] # or += [0]
				start = i
				reg_type = 2
				num_diff = 0 # redundant
			cons_run = 0
		# If the base is masked somewhere:
		# elif poly_state == 3:
		else:
			if reg_type == 0:
				# print('starting masked region')
				region_types += [reg_type]
				region_starts += [start]
				num_fixed_diffs += [num_diff] # or += [0]
				start = i
				reg_type = 3
				num_diff = 0 # redundant
			# All other cases, change region to masked
			else:
				# print('switching to masked region')
				reg_type = 3
			cons_run = 0
		i += 1
	if reg_type == 0:
		# print('ending on conserved region')
		region_types += [reg_type]
		region_starts += [start]
		num_fixed_diffs += [num_diff] # or += [0]
		# Add a masked tail, ending the sequence
		region_types += [3]
		region_starts += [i]
		num_fixed_diffs += [0] # Hmmmm
	else:
		# print('ending on non-conserved region')
		region_types += [3]
		region_starts += [start]
		num_fixed_diffs += [num_diff] # Hmmmm
	return((region_types,region_starts,num_fixed_diffs))

# Maps the polymorphism data to regions:
#   0: Conserved regions with fixed shared sequence > min_prim_len
#   1: Regions between, containing no N masking and some fixed differences
#   2: Regions between, containing no N masking, but also no fixed differences
#   3: Regions between, containing at least one N masking position
# Gives output as a list of region identities, 
#   a list of start positions by genomic sequence
#   a list of number of fixed differences
def map_regions_exp(poly_calls, inv_start, inv_end, min_prim_len, buffer_len, prim_buff):
	region_types = []
	region_starts = []
	num_fixed_diffs = []
	i = inv_start + buffer_len
	start = i
	cons_run = 0
	reg_type = 0
	num_diff = 0
	min_cons_reg = min_prim_len + prim_buff
	# Check for parameter mismatch to chrom/inversion data
	if i >= inv_end-buffer_len:
	    raise Exception('buffer larger than inversion: inv_start+buffer_len >= inv_end-buffer_len')
	# Start the region scan
	while i < inv_end-buffer_len:
		poly_state = poly_calls[i]
		# If the base is conserved:
		# keep going unless the conserved run = min len
		if poly_state < 4:
			if reg_type != 0:
				cons_run += 1
				if cons_run == min_cons_reg:
					# print('starting conserved region')
					region_types += [reg_type]
					region_starts += [start]
					num_fixed_diffs += [num_diff]
					start = i-cons_run+1
					reg_type = 0
					num_diff = 0
					cons_run = 0
		# If the base is fixed different:
		elif poly_state == 4:
			if reg_type == 0 or reg_type == 3: # Must break masking regs if ignoring some
				# print('starting difference region')
				region_types += [reg_type]
				region_starts += [start]
				num_fixed_diffs += [num_diff] # or += [0]
				start = i
				reg_type = 1
				num_diff = 1
			# elif reg_type == 1: # Redundant to final else
			# 	num_diff += 1
			elif reg_type == 2:
				# print('switching to difference region')
				reg_type = 1
				num_diff = 1
			else: # Might want to allow N's in amplicon? Consider
				num_diff += 1
			cons_run = 0
		# If the base is partially shared:
		# ignore it unless breaking a conserved region
		elif poly_state == 5:
			if reg_type == 0:
				# print('starting polymorphic region')
				region_types += [reg_type]
				region_starts += [start]
				num_fixed_diffs += [num_diff] # or += [0]
				start = i
				reg_type = 2
				num_diff = 0 # redundant
			cons_run = 0
		# If the base is masked somewhere:
		else: # if poly_state == 6:
			if reg_type == 0 or reg_type == 1: # Must break masking regs if ignoring some
				# print('starting masked region')
				region_types += [reg_type]
				region_starts += [start]
				num_fixed_diffs += [num_diff] # or += [0]
				start = i
				reg_type = 3
				num_diff = 0 # redundant
			# All other cases, change region to masked
			else:
				# print('switching to masked region')
				reg_type = 3
			cons_run = 0
		i += 1
	if reg_type == 0:
		# print('ending on conserved region')
		region_types += [reg_type]
		region_starts += [start]
		num_fixed_diffs += [num_diff] # or += [0]
		# Add a masked tail, ending the sequence
		region_types += [3]
		region_starts += [i]
		num_fixed_diffs += [0] # Hmmmm
	else:
		# print('ending on non-conserved region')
		region_types += [3]
		region_starts += [start]
		num_fixed_diffs += [num_diff] # Hmmmm
	# print(region_types)
	# print(region_starts)
	# print(num_fixed_diffs)
	return((region_types,region_starts,num_fixed_diffs))

def gen_region_maps(poly_maps,inv_data,prim_len,buff,prim_buff):
	region_maps = []
	for i in np.arange(len(inv_data)):
		(chrom,inv,br_1,br_2,inv_name,std_name) = inv_data[i]
		print("Generating region map for "+str(inv_name))
		# region_maps += [map_regions(poly_maps[i],br_1,br_2,prim_len,buff,prim_buff)]
		region_maps += [map_regions_exp(poly_maps[i],br_1,br_2,prim_len,buff,prim_buff)]
	return region_maps


# Writes a .csv of the region id, start and fixed difference count
def write_region_csv(out_file,region_types,region_starts,num_fixed_diffs):
	with open(out_file,'w') as output:
		for i in np.arange(len(region_types)):
			output.write(str(region_types[i])+', '+str(region_starts[i])\
				+', '+str(num_fixed_diffs[i])+'\n')
	return

def write_region_outputs(region_maps,inv_data,out_dir):
	for i in np.arange(len(inv_data)):
		(chrom,inv,br_1,br_2,inv_name,std_name) = inv_data[i]
		file_name = out_dir+inv_name+'_regions.csv'
		(region_types,region_starts,num_fixed_diffs) = region_maps[i]
		print('Writing polymorphic region .csv for '+inv_name)
		write_region_csv(file_name,region_types,region_starts,num_fixed_diffs)
	return









# Proposes amplicons for PCR from region data
def propose_amplicons_no_mask(region_data,min_gap,max_gap,prim_len,prim_buff):
	(region_types,region_starts,num_fixed_diffs) = region_data
	amplicon_starts = []
	amplicon_ends = []
	targets = []
	target_lengths = []
	spanned_diffs = []
	min_cons_reg = prim_len + prim_buff
	i = 0
	while i < len(region_types):
		if region_types[i] == 0:
			j = i+1
			start = region_starts[i]
			var_start = region_starts[i+1]
			done = False
			spans_diff = False
			fixed_diff = 0
			while not done:
				reg_type = region_types[j]
				reg_start = region_starts[j]
				if (reg_start - var_start > max_gap) or reg_type == 3:
					done = True
				elif reg_type == 1:
					spans_diff = True
					fixed_diff += num_fixed_diffs[j]
				elif reg_type == 0 and (region_starts[j+1] - start - 2*min_cons_reg > min_gap) \
							and spans_diff:
					middle = (var_start+reg_start)//2
					amp_start = max(start,middle-max_gap//2-min_cons_reg)
					amp_end = min(region_starts[j+1],middle+max_gap//2+min_cons_reg)
					amplicon_starts += [amp_start]
					amplicon_ends += [amp_end]
					targets += [var_start]
					target_lengths += [reg_start - var_start]
					spanned_diffs += [fixed_diff]
				j += 1
		i += 1
	return((amplicon_starts,amplicon_ends,targets,target_lengths,spanned_diffs))

# Proposes amplicons for PCR from region data,
#   takes sequence data from the polymorphism map
#   returns 
def propose_amplicons(region_data,min_gap,max_gap,prim_len,prim_buff):
	(region_types,region_starts,num_fixed_diffs) = region_data
	amplicon_starts = []
	amplicon_ends = []
	targets = []
	target_lengths = []
	spanned_diffs = []
	# spanned_masks = []
	min_cons_reg = prim_len + prim_buff
	num_regions = len(region_types)
	i = 0
	while i < len(region_types):
		if region_types[i] == 0:
			j = i+1
			start = region_starts[i]
			var_start = region_starts[i+1]
			done = False
			spans_diff = False
			fixed_diff = 0
			# num_masked = 0
			while not done:
				reg_type = region_types[j]
				reg_start = region_starts[j]
				if (reg_start - var_start > max_gap) or (j == num_regions - 1):
					done = True
				elif reg_type == 1:
					spans_diff = True
					fixed_diff += num_fixed_diffs[j]
				# elif reg_type == 3:
				# 	num_masked += 1
				elif reg_type == 0 and (region_starts[j+1] - start - 2*min_cons_reg > min_gap)\
							and spans_diff:
					middle = (var_start+reg_start)//2
					amp_start = max(start,middle-max_gap//2-min_cons_reg)
					amp_end = min(region_starts[j+1],middle+max_gap//2+min_cons_reg)
					amplicon_starts += [amp_start]
					amplicon_ends += [amp_end]
					targets += [var_start]
					target_lengths += [reg_start - var_start]
					spanned_diffs += [fixed_diff]
					# spanned_masks += [num_masked]
				j += 1
		i += 1
	# print(amplicon_starts)
	# print(amplicon_ends)
	# return((amplicon_starts,amplicon_ends,targets,\
	# 		target_lengths,spanned_diffs,spanned_masks))
	return((amplicon_starts,amplicon_ends,targets,\
			target_lengths,spanned_diffs))

# For retrieving a sequence from lines extracted from .fas1k
def subseq_from_lines(start,end,line_len,lines):
	# Retrieve the sequence
	curr_line = start // line_len
	start_index = start %  line_len
	end_line = end // line_len
	end_index = end % line_len
	seq = ''
	if curr_line == end_line:
		seq += lines[curr_line].strip()[start_index:end_index]
	else:
		seq += lines[curr_line].strip()[start_index:]
		curr_line += 1
		while curr_line < end_line:
			seq += lines[curr_line].strip()
			curr_line += 1
		seq += lines[curr_line].strip()[:end_index]
	return(seq)


# For extracting many strain subsequences simultaneously
def extract_amplicon_seqs(amp_starts,amp_ends,strain_file):
	# Prep relevant values
	strain = open(strain_file,'r')
	lines = strain.readlines()
	strain.close()
	line_len = len(lines[0].strip())
	num_lines = len(lines)
	remainder_len = len(lines[num_lines-1].strip())
	chrom_len = line_len*(num_lines-1)+remainder_len
	# Retrieve the sequences
	seqs = []
	for i in np.arange(len(amp_starts)):
		# Retrieve the sequence

		# start = amp_starts[i]
		# curr_line = start // line_len
		# start_index = start %  line_len
		# end = amp_ends[i]
		# end_line = end // line_len
		# end_index = end % line_len
		# seq = ''
		# if curr_line == end_line:
		# 	seq += lines[curr_line].strip()[start_index:end_index]
		# else:
		# 	seq += lines[curr_line].strip()[start_index:]
		# 	curr_line += 1
		# 	while curr_line < end_line:
		# 		seq += lines[curr_line].strip()
		# 		curr_line += 1
		# 	seq += lines[curr_line].strip()[:end_index]
		seq = subseq_from_lines(start,end,line_len,lines)
		seqs += [seq]
	return(seqs)

def extract_amplicon_seqs_from_poly(amp_starts,amp_ends,poly_calls):
	# Uses the poly translation dictionary defined initially ^
	# could just use poly_trans = ['A','T','G','C','X','Y','N']
	seqs = []
	for i in np.arange(len(amp_starts)):
		# Retrieve the sequence
		start = amp_starts[i]
		end = amp_ends[i]
		seq = ''
		for j in np.arange(start,end):
			seq += poly_trans[poly_calls[j]]
		# print(seq)
		seqs += [seq]
	return(seqs)

# def propose_amplicon_sets(region_maps, inv_data, min_gap, max_gap, prim_len, prim_buff, 
# 		inv_strain_dict, genome_dir, out_dir):
# 	proposed_amplicon_sets = []
# 	for i in np.arange(len(inv_data)):
# 		(chrom,inv,br_1,br_2,inv_name,std_name) = inv_data[i]
# 		print("Generating amplicon proposals for "+str(inv_name))
# 		(amp_starts,amp_ends,targets,tar_lens,diffs,masks) = \
# 			propose_amplicons(region_maps[i],min_gap,max_gap,prim_len,prim_buff)
# 		# Retrieve the sequence of each amplicon
# 		strain = inv_strain_dict[inv_name][0]
# 		prefix = genome_dir + 'Chr' + chrom + '/' + strain + '_Chr' + chrom
# 		strain_file = pick_strain_file(prefix)
# 		# strain_file = genome_dir+'Chr'+chrom+'/'+strain+'_Chr'+chrom+'.fas1k'
# 		file_name = out_dir+inv_name+'_amp_seqs.txt'
# 		seqs = extract_amplicon_seqs(amp_starts,amp_ends,strain_file)
# 		proposed_amplicon_sets += [(seqs,amp_starts,amp_ends,targets,tar_lens,diffs,masks)]
# 	return(proposed_amplicon_sets)

# To count the indels and failed base calls of an arbitrary sequence,
# Assumes runs of 'N' represent indels, and single 'N's represent failed base calls
def count_indels_unassigned_bases(seq):
	i = 0
	indel_count = 0
	failed_base_count = 0
	# print(seq)
	while i < len(seq):
		if seq[i] == 'N':
			run = 0
			while seq[i] == 'N' and i < len(seq):
				# print('run: '+str(i))
				run += 1
				i += 1
			if run > 1:
				indel_count += 1
			else:
				failed_base_count += 1
		else:
			i += 1
	# print('indels: '+str(indel_count)+'unknown: '+str(failed_base_count))
	return(indel_count,failed_base_count)

# To count the indels, failed base calls,
#   and polymorphic bases shared between the two inversion states
def count_indel_fail_poly(seq):
	(indel_count,fail_count) = count_indels_unassigned_bases(seq)
	poly_count = seq.count('Y')
	return((poly_count,fail_count,indel_count))

# Wrapping indel, poly, and failed base counts for multiple sequences
def get_indel_fail_poly_sets(seqs):
	indels = []
	fails = []
	polys = []
	for seq in seqs:
		(poly,fail,indel) = count_indel_fail_poly(seq)
		indels += [indel]
		fails += [fail]
		polys += [poly]
	return((polys,fails,indels))


def propose_amplicon_sets_poly(region_maps, inv_data, min_gap, max_gap, prim_len, prim_buff, 
		poly_maps, genome_dir, out_dir):
	proposed_amplicon_sets = []
	for i in np.arange(len(inv_data)):
		(chrom,inv,br_1,br_2,inv_name,std_name) = inv_data[i]
		print("Generating amplicon proposals for "+str(inv_name))
		(amp_starts,amp_ends,targets,tar_lens,diffs) = \
			propose_amplicons(region_maps[i],min_gap,max_gap,prim_len,prim_buff)
		# Retrieve the sequence of each amplicon
		seqs = extract_amplicon_seqs_from_poly(amp_starts,amp_ends,poly_maps[i])
		(polys,fails,indels) = get_indel_fail_poly_sets(seqs)
		proposed_amplicon_sets += [(seqs,amp_starts,amp_ends,targets,tar_lens,
			diffs,polys,fails,indels)]
	return(proposed_amplicon_sets)


# For retrieving a subsequence from a .fas1k strain file
def extract_strain_subseq(start,end,strain_file):
	# Prep relevant values
	strain = open(strain_file,'r')
	lines = strain.readlines()
	strain.close()
	line_len = len(lines[0].strip())
	seq = subseq_from_lines(start,end,line_len,lines)
	return(seq)

def hap_mask_screen(inv_haps,std_haps):
	i_seqs = []
	i_freqs = []
	for seq, freq in inv_haps:
		i_seqs += [key]
		i_freqs += [freq]
	s_seqs = []
	s_freqs = []
	for seq, freq in std_haps:
		s_seqs += [key]
		s_freqs += [freq]
	seq_len = len(i_seqs[0])
	for i in np.arange(seq_len):
		pass
		



# # Runs the haplotype analysis for the amplicon generation
# # Requires strain file access for haplotype data
# def check_haplotypes_and_misclass(inv_strain_files,std_strain_files,start,end):
# 	# inv_haps = []
# 	# std_haps = []
# 	inv_haps = {}
# 	std_haps = {}
# 	i_hap_strain_dict = {}
# 	s_hap_strain_dict = {}
# 	shared_haps = False
# 	for strain_file in inv_strain_files:
# 		seq = extract_strain_subseq(start,end,strain_file)
# 		# inv_haps.add(seq)
# 		# inv_haps+=[seq]
# 		if not seq in inv_haps:
# 			inv_haps[seq] = 0
# 		else:
# 			inv_haps[seq] = inv_haps[seq]+1
# 		if not seq in i_hap_strain_dict:
# 			i_hap_strain_dict[seq] = [strain_file]
# 		else:
# 			i_hap_strain_dict[seq] = i_hap_strain_dict[seq]+[strain_file]
# 	for strain_file in std_strain_files:
# 		seq = extract_strain_subseq(start,end,strain_file)
# 		if seq in inv_haps:
# 			shared_haps = True
# 		# std_haps.add(seq)
# 		# std_haps+=[seq]
# 		if not seq in std_haps:
# 			std_haps[seq] = 0
# 		else:
# 			std_haps[seq] = std_haps[seq]+1
# 		if not seq in s_hap_strain_dict:
# 			s_hap_strain_dict[seq] = [strain_file]
# 		else:
# 			s_hap_strain_dict[seq] = s_hap_strain_dict[seq]+[strain_file]
# 	if not shared_haps:
# 		print("Unique, identifying haplotypes at "+str(start)+':'+str(end))
# 		print("Inv:")
# 		for seq in set(inv_haps):
# 			print(seq)
# 		print("Std:")
# 		for seq in set(std_haps):
# 			print(seq)
# 	else:
# 		for seq in inv_haps:
# 			print(inv_haps)
# 			print(std_haps)
# 			if inv_haps[seq] < 3 and std_haps[seq] > 2:
# 				print('possible misclassification as an inversion: '+\
# 					str(i_hap_strain_dict[seq]))
# 			elif std_haps[seq] < 3 and inv_haps[seq] > 2:
# 				print('possible misclassification as a standard: '+\
# 					str(s_hap_strain_dict[seq]))
# 	return((shared_haps,inv_haps,std_haps))

# Runs the haplotype analysis for the amplicon generation
# Requires strain file access for haplotype data
def check_haplotypes(inv_strain_files,std_strain_files,start,end):
	inv_haps = set({})
	std_haps = set({})
	shared_haps = False
	for strain_file in inv_strain_files:
		seq = extract_strain_subseq(start,end,strain_file)
		inv_haps.add(seq)
	for strain_file in std_strain_files:
		seq = extract_strain_subseq(start,end,strain_file)
		if seq in inv_haps:
			shared_haps = True
		std_haps.add(seq)
	# if not shared_haps:
	# 	print("Unique, identifying haplotypes at "+str(start)+':'+str(end))
	# 	print("Inv:")
	# 	for seq in set(inv_haps):
	# 		print(seq)
	# 	print("Std:")
	# 	for seq in set(std_haps):
	# 		print(seq)
	return((shared_haps,inv_haps,std_haps))


# Proposes amplicons for PCR based on fixed haplotype differences from region data
def propose_haplotype_amplicons(inv_strains,std_strains,region_data,min_gap,max_gap,
				prim_len,prim_buff):
	(region_types,region_starts,num_fixed_diffs) = region_data
	amplicon_starts = []
	amplicon_ends = []
	targets = []
	target_lengths = []
	spanned_diffs = []
	spanned_masks = []
	num_regions = len(region_types)
	min_cons_reg = prim_len + prim_buff
	i = 0
	while i < len(region_types):
		if region_types[i] == 0:
			j = i+1
			start = region_starts[i]
			var_start = region_starts[i+1]
			done = False
			fixed_diff = 0
			num_masked = 0
			while not done:
				reg_type = region_types[j]
				reg_start = region_starts[j]
				if (reg_start - var_start > max_gap) or (j == num_regions - 1):
					done = True
				elif reg_type == 1:
					fixed_diff += num_fixed_diffs[j]
				elif reg_type == 3:
					num_masked += 1
				elif reg_type == 0 and (region_starts[j+1] - start - 2*min_cons_reg > min_gap)\
						and (num_masked == 0):
					hap_data = check_haplotypes(inv_strains,std_strains,var_start,reg_start)
					if not hap_data[0]:
						# CHANGE this to center the amplicon around the differences?
						middle = (var_start+reg_start)//2
						amp_start = max(start,middle-max_gap//2-min_cons_reg)
						amp_end = min(region_starts[j+1],middle+max_gap//2+min_cons_reg)
						amplicon_starts += [amp_start]
						amplicon_ends += [amp_end]
						targets += [var_start]
						target_lengths += [reg_start - var_start]
						spanned_diffs += [fixed_diff]
						spanned_masks += [num_masked]
				j += 1
		i += 1
	return((amplicon_starts,amplicon_ends,targets,\
			target_lengths,spanned_diffs,spanned_masks))

def propose_hap_amplicon_sets(region_maps, inv_data, min_gap, max_gap, prim_len, prim_buff, 
		poly_maps, inv_strain_dict, genome_dir, out_dir):
	proposed_amplicon_sets = []
	for i in np.arange(len(inv_data)):
		(chrom,inv,br_1,br_2,inv_name,std_name) = inv_data[i]
		inv_strain_files = []
		std_strain_files = []
		for strain in inv_strain_dict[inv_name]:
			prefix = genome_dir + 'Chr' + chrom + '/' + strain + '_Chr' + chrom
			file_name = pick_strain_file(prefix)
			inv_strain_files += [file_name]
		for strain in inv_strain_dict[std_name]:
			prefix = genome_dir + 'Chr' + chrom + '/' + strain + '_Chr' + chrom
			file_name = pick_strain_file(prefix)
			std_strain_files += [file_name]
		print("Generating haplotype-based amplicon proposals for "+str(inv_name))
		(amp_starts,amp_ends,targets,tar_lens,diffs,masks) = \
			propose_haplotype_amplicons(inv_strain_files,std_strain_files,
				region_maps[i],min_gap,max_gap,prim_len,prim_buff)
		# Retrieve the sequence of each amplicon
		seqs = extract_amplicon_seqs_from_poly(amp_starts,amp_ends,poly_maps[i])
		(polys,fails,indels) = get_indel_fail_poly_sets(seqs)
		proposed_amplicon_sets += [(seqs,amp_starts,amp_ends,targets,tar_lens,
			diffs,polys,fails,indels)]
	return(proposed_amplicon_sets)

# Writes the amplicon data in csv: start, end, target, target length, fixed diffs, seqs
def write_amplicons(out_file,amp_data):
	(seqs,amp_starts,amp_ends,targets,tar_lens,diffs,polys,fails,indels) = amp_data
	# Write the sequences
	with open(out_file,'w') as output:
		header = "Amp_ID, Start, End, Target_Start, Target_Length, "+\
			"N_Diff, N_Poly, N_Uncalled_Bases, N_Indels, Sequence\n"
		output.write(header)
		for i in np.arange(len(amp_starts)):
			line = str(i)+', '+str(amp_starts[i])+', '+str(amp_ends[i])+', '+str(targets[i])+\
			', '+str(tar_lens[i])+', '+str(diffs[i])+', '+str(polys[i])+', '+str(fails[i])+\
			', '+str(indels[i])+', '+seqs[i]+'\n'
			output.write(line)
	return

def write_amplicon_sets(amp_sets, inv_data, out_dir):
	for i in np.arange(len(inv_data)):
		(chrom,inv,br_1,br_2,inv_name,std_name) = inv_data[i]
		file_name = out_dir+inv_name+'_amps.txt'
		amp_data = amp_sets[i]
		print("Writing amplicon sequences for "+str(inv_name))
		write_amplicons(file_name,amp_data)
	return

# Tries to find small target regions with highest fixed difference counts,
#   per amplicon set
# There is a LOT more that could be done - other forms of polymorphism,
#   distance to the breakpoints, other genomic features
def filter_amplicons(amp_data,num_returned):
	(seqs,amp_starts,amp_ends,targets,tar_lens,diffs,polys,fails,indels) = amp_data
	num_amps = len(amp_starts)
	# a value of 0 prioritizes differences, 1 prioritizes difference density
	length_factor = 0.35
	if num_amps > num_returned:
		scores = []
		for i in np.arange(num_amps):
			score = diffs[i]/(tar_lens[i]**length_factor)
			scores += [score]
		indices = np.argsort(scores)[-num_returned:]
		filtered_amps = ([seqs[i] for i in indices],
			[amp_starts[i] for i in indices],
			[amp_ends[i] for i in indices],
			[targets[i] for i in indices],
			[tar_lens[i] for i in indices],
			[diffs[i] for i in indices],
			[polys[i] for i in indices],
			[fails[i] for i in indices],
			[indels[i] for i in indices])
		return(filtered_amps)
	else:
		return(amp_data)

# Returns a subset of amplicons based on the criteria used in 'filter_amplicons'
def filter_amplicon_sets(amp_sets, num_returned):
	# num_returned = 15
	filtered_amp_sets = []
	for amp_data in amp_sets:
		filtered_amp_sets += [filter_amplicons(amp_data,num_returned)]
	return(filtered_amp_sets)








# For correctly formatting amplicon suggestions to primer3
def write_primer3_input(out_file,amp_data,inv_name):
	(seqs,amp_starts,amp_ends,targets,tar_lens,diffs,polys,fails,indels) = amp_data
	with open(out_file,'w') as BoulderIO_file:
		for i in np.arange(len(amp_starts)):
			# Remove 'X' and 'Y' characters from the sequence
			seq = seqs[i]
			seq = seq.replace('X','N')
			seq = seq.replace('Y','N')
			# Prepare the amplicon entry
			target = targets[i] - amp_starts[i]
			seq_ID = inv_name+'_amp_'+str(i)
			prim3_input_record = 'SEQUENCE_ID='+seq_ID+\
				'\nSEQUENCE_TEMPLATE='+seq+\
				'\nSEQUENCE_TARGET='+str(target)+','+str(tar_lens[i])+\
				'\nSEQUENCE_INTERNAL_EXCLUDED_REGION='+str(target)+','\
				+str(tar_lens[i])+\
				'\n=\n'
			# Write the record
			BoulderIO_file.write(prim3_input_record)
	return

def write_primer3_inputs(amp_sets, inv_data, out_dir):
	for i in np.arange(len(inv_data)):
		(chrom,inv,br_1,br_2,inv_name,std_name) = inv_data[i]
		file_name = out_dir+'primer3/'+inv_name+'_prim3_in.txt'
		amp_data = amp_sets[i]
		print("Writing primer3 input file for "+str(inv_name)+\
			" proposed amplicons in Boulder-IO format")
		write_primer3_input(file_name,amp_data,inv_name)
	return

# Mostly a placeholder for writing settings from the script input
#   to the primer3 input for batch primer generation
#   Sets min and max primer size to primer length +/- 2
def write_primer3_settings(prim_len, min_gap, max_gap, mask_lib, set_file):
	prim_len_range = 2
	settings = ('Primer3 File - http://primer3.org\n'
	'P3_FILE_TYPE=settings\n\n'
	'P3_FILE_ID=A standard set of settings for'
	'polymorphism primer generation\n'
	'PRIMER_OPT_SIZE='+str(prim_len)+'\n'
	'PRIMER_MIN_SIZE='+str(prim_len-prim_len_range)+'\n'
	'PRIMER_MAX_SIZE='+str(prim_len+prim_len_range)+'\n'
	'PRIMER_PRODUCT_SIZE_RANGE='+str(min_gap)+'-'+str(max_gap)+'\n')
	if (mask_lib != '') and os.path.exists(mask_lib):
		settings += 'PRIMER_MISPRIMING_LIBRARY='+mask_lib+'\n='
	else:
		settings += '='
	with open(set_file,'w') as set_out:
		set_out.write(settings)
	return

# A silly helper function to add qoutes to unusual path names
def str_q(directory_string):
	return('\''+directory_string+'\'')


# A wrapper for a command line call to primer3
#   to process the amplicon files 
def call_primer3(in_name,out_name,err_name,settings_name,prim3_path):
	primer3_call = './'+prim3_path+'primer3_core --p3_settings_file='+settings_name+\
	' --output='+str_q(out_name)+' --error='+str_q(err_name)+' '+str_q(in_name)
	# primer3_call = 'sudo ./'+prim3_path+'primer3_core --output='+str_q(out_name)+\
	# ' --error='+str_q(err_name)+' '+str_q(in_name)
	# primer3_call = './'+prim3_path+'primer3_core --p3_settings_file='+settings_name+\
	# ' --output='+str_q(out_name)+' '+str_q(in_name)
	os.system(primer3_call)
	return

# Attempts to pick primer pairs from the region map
def call_primer3_all(inv_data, out_dir, settings_name, prim3_path):
	for i in np.arange(len(inv_data)):
		(chrom,inv,br_1,br_2,inv_name,std_name) = inv_data[i]
		in_name = out_dir+'primer3/'+inv_name+'_prim3_in.txt'
		out_name = out_dir+'primer3/'+inv_name+'_prim3_out.txt'
		err_name = out_dir+'primer3/'+inv_name+'_prim3_err.txt'
		print("Calling primer3 for "+str(inv_name))
		call_primer3(in_name,out_name,err_name,settings_name,prim3_path)
	return

# Manages all primer3 input, output, and calls at ounce
def gen_primer3_proposals(amp_sets, inv_data, prim_len, 
		min_gap, max_gap, out_dir, prim3_path, mask_lib):
	write_primer3_inputs(amp_sets, inv_data, out_dir)
	set_file = out_dir+'primer3/prim3_gen_settings.txt'
	write_primer3_settings(prim_len, min_gap, max_gap, mask_lib, set_file)
	call_primer3_all(inv_data, out_dir, set_file, prim3_path)
	return

# Parses Boulder-IO primer3 output files
#  MODIFY to collect mispriming, self, pair-dimer, hairpin th data
def parse_primer3_output(p3_file):
	primers = []
	text = ''
	with open(p3_file,'r') as primer_file:
		text = primer_file.read()
	sections = text.split('\n=\n')[:-1]
	num_amps = len(sections)
	for i in np.arange(num_amps):
		lines = sections[i].strip().split('\n')
		num_pairs = int(lines[4][-1])
		primer_set = set({})
		pair_line_size = 0
		if num_pairs > 0:
			pair_line_size = (len(lines)-8)//num_pairs
		# print(pair_line_size)
		for j in np.arange(num_pairs):
			p_start = 8+(j*pair_line_size)
			l_seq = lines[p_start+3].split('=')[1] #or [23:] for .split...
			r_seq = lines[p_start+4].split('=')[1]
			if not l_seq in primer_set:
				primer_set.add(l_seq)
				l_penal = float(lines[p_start+1].split('=')[1])
				l_pos = int(lines[p_start+5].split('=')[1].split(',')[0])
				l_len = int(lines[p_start+5].split('=')[1].split(',')[1])
				l_tm = float(lines[p_start+7].split('=')[1])
				l_gc = float(lines[p_start+9].split('=')[1])
				primers += [(i,l_seq,'f',l_penal,l_pos,l_len,l_tm,l_gc)]
			if not r_seq in primer_set:
				primer_set.add(r_seq)
				r_penal = float(lines[p_start+2].split('=')[1])
				r_pos = int(lines[p_start+6].split('=')[1].split(',')[0])
				r_len = int(lines[p_start+6].split('=')[1].split(',')[1])
				r_tm = float(lines[p_start+8].split('=')[1])
				r_gc = float(lines[p_start+10].split('=')[1])
				primers += [(i,r_seq,'r',r_penal,r_pos,r_len,r_tm,r_gc)]
	return(primers)

# For calling the primer3 output parser across inversions
def read_primer3_proposals(inv_data, out_dir):
	prim3_proposals = []
	for i in np.arange(len(inv_data)):
		(chrom,inv,br_1,br_2,inv_name,std_name) = inv_data[i]
		print('Reading '+inv_name+' primer3 output')
		prim3_proposals += [parse_primer3_output(out_dir+\
			'primer3/'+inv_name+'_prim3_out.txt')]
	return(prim3_proposals)


# Writes a tab delineated file of the primer data
def write_primer_proposal_tdf(primer_data,out_file):
	with open(out_file,'w') as output:
		header = "Seq\tAmp_ID\tFor_Rev\tStart_Offset\tLen\tT_m\tGC_rat\tPenalty\n"
		output.write(header)
		for (i,seq,direction,penal,start,length,tm,gc) in primer_data:
			primer_line = seq+"\t"+str(i)+"\t"+direction+"\t"+str(start)+\
				"\t"+str(length)+"\t"+str(tm)+"\t"+str(gc)+"\t"+str(penal)+"\n"
			output.write(primer_line)
	return

# Writes a tab delineated file of the primer data
def write_primer_proposals(inv_data, primer_sets, out_dir):
	for i in np.arange(len(inv_data)):
		(chrom,inv,br_1,br_2,inv_name,std_name) = inv_data[i]
		print('Writing '+inv_name+' primer data')
		out_file = out_dir+'primer_output/'+inv_name+'_primers.txt'
		write_primer_proposal_tdf(primer_sets[i],out_file)
	return





# For generating primer data lines with amplicon information included
def combine_amp_primer3_data(amp_data, primer_data):
	(seqs,amp_starts,amp_ends,targets,tar_lens,diffs,polys,fails,indels) = amp_data
	primers_detailed = []
	for (i,seq,direction,penal,start,length,tm,gc) in primer_data:
		amp_seq = seqs[i]
		genom_start = amp_starts[i] + start
		genom_tar = targets[i]
		tar_rel_start = genom_tar - amp_starts[i]
		added_length = max(genom_tar - genom_start,genom_start - genom_tar - tar_lens[i])
		tar_diffs = amp_seq.count('X')
		if tar_diffs != diffs[i]:
			print('amplicon sequence diffs: '+str(tar_diffs)+\
				' do not match recorded diffs: '+str(diffs[i]))
		# tar_poly = amp_seq.count('Y')
		# (indel_count,failed_base_count) = count_indels_unassigned_bases(amp_seq)
		primers_detailed += [(seq,i,direction,genom_start,length,tm,gc,penal,
			tar_diffs,polys[i],indels[i],fails[i],added_length,tar_rel_start,
			tar_lens[i],amp_seq)]
	return(primers_detailed)

# For reading the primer3 output and combining it with the
#   amplicon data
def combine_amp_primer3_sets(primer_sets, amp_sets):
	detailed_sets = []
	for i in np.arange(len(amp_sets)):
		print('Combining '+inv_name+' primer and amp output')
		detailed_primers = combine_amp_primer3_data(amp_sets[i],primer_sets[i])
		detailed_sets += [detailed_primers]
	return(detailed_sets)

# For reading the primer3 output and combining it with the
#   amplicon data
def process_primer3_output(inv_data, amp_sets, out_dir):
	primer_sets = []
	for i in np.arange(len(inv_data)):
		(chrom,inv,br_1,br_2,inv_name,std_name) = inv_data[i]
		print('Reading '+inv_name+' primer3 output')
		prim3_proposals = parse_primer3_output(out_dir+\
			'primer3/'+inv_name+'_prim3_out.txt')
		print('Combining '+inv_name+' primer and amp output')
		primers = combine_amp_primer3_data(amp_sets[i],prim3_proposals)
		primer_sets += [primers]
	return(primer_sets)

# For selecting specific primers to focus on,
#   based on the target length, primer data, and target diffs and indels
#   and adding the FREQ-seq oligos for those primers selected
def filter_primers(primers_detailed,gc_bias):
	pass


# Writes a tab delineated file of the primer data, in full detail with amp data
def write_detailed_primer_proposal_tdf(detailed_data,out_file):
	with open(out_file,'w') as output:
		header = "Seq\tAmp_ID\tFor_Rev\tLen\tT_m\tGC_rat\tPenalty\tStart\t"+\
			"N_Diff\tN_Shared_Poly\tN_Indel\tN_Fail_Base\tAdded_Amp_Len\t"+\
			"Target_Len\tTarget_Rel_Pos\tAmplicon_Region\n"
		output.write(header)
		for (seq,i,direction,genom_start,length,tm,gc,penal,
			tar_diffs,tar_poly,indel_count,fail_count,added_length,
			target_pos,tar_len,amp_seq) in detailed_data:
			primer_line = seq+"\t"+str(i)+"\t"+direction+"\t"+str(length)+\
				"\t"+str(tm)+"\t"+str(gc)+"\t"+str(penal)+"\t"+str(genom_start)+\
				"\t"+str(tar_diffs)+"\t"+str(tar_poly)+"\t"+str(indel_count)+\
				"\t"+str(fail_count)+"\t"+str(added_length)+"\t"+str(tar_len)+\
				"\t"+str(target_pos)+"\t"+amp_seq+"\n"
			output.write(primer_line)
	return

# For writing detailed primer data sets
def write_detailed_primer_proposals(inv_data, detailed_primer_sets, out_dir):
	for i in np.arange(len(inv_data)):
		(chrom,inv,br_1,br_2,inv_name,std_name) = inv_data[i]
		print('Writing '+inv_name+' primer data')
		out_file = out_dir+'primer_output/'+inv_name+'_primers_detail.txt'
		write_detailed_primer_proposal_tdf(detailed_primer_sets[i],out_file)
	return







# Takes the left and right primers and the stubby Illumina adaptors,
#   and returns a set of stubby-adaptor oligos for library prep of pooled samples
# Must have primers paired one-to-one left-right
def gen_oligos(left_primers,right_primers,adapt_F,adapt_R):
	num_primers = len(left_primers)
	left_oligos = []
	right_oligos = []
	for i in np.arange(num_primers):
		left_oligo = adapt_F+left_primers[i]
		right_oligo = adapt_R+right_primers[i]
		left_oligos += [left_oligo]
		right_oligos += [right_oligo]
	return((left_oligos,right_oligos))

# Takes the left and right primers, the barcodes and Illumina bridging sequences,
#   and returns a set of dual-barcoded oligos for direct library prep of pooled samples
# Must have more barcodes than primers, and primers paired one-to-one left-right
def gen_indexed_oligos(left_primers,right_primers,barcodes,adapt_F,adapt_R,illum_A,illum_B):
	i5_adaptor_segments = illum_A.split('[i5]')
	i7_adaptor_segments = illum_B.split('[i7]')
	num_primers = len(left_primers)
	left_oligos = []
	right_oligos = []
	for i in np.arange(num_primers):
		left_oligo = i5_adaptor_segments[0]+barcodes[2*i]+\
			i5_adaptor_segments[1]+adapt_F+left_primers[i]
		right_oligo = i7_adaptor_segments[0]+barcodes[2*i]+\
			i7_adaptor_segments[1]+adapt_R+right_primers[i]
		left_oligos += [left_oligo]
		right_oligos += [right_oligo]
	return((left_oligos,right_oligos))

# Writes an oligo set to file
def write_oligos(oligos,file_prefix,out_dir):
	out_file = out_dir+file_prefix+'_oligos.txt'
	with open(out_file,'w') as output:
		header = "Seq\tAmp_ID\tFor_Rev\tLen\tT_m\tGC_rat\tPenalty\tStart\t"+\
			"N_Diff\tN_Shared_Poly\tN_Indel\tN_Fail_Base\tAdded_Amp_Len\t"+\
			"Target_Len\tTarget_Rel_Pos\tAmplicon_Region\n"
		output.write(header)
		for (seq,i,direction,genom_start,length,tm,gc,penal,
			tar_diffs,tar_poly,indel_count,fail_count,added_length,
			target_pos,tar_len,amp_seq) in detailed_data:
			primer_line = seq+"\t"+str(i)+"\t"+direction+"\t"+str(length)+\
				"\t"+str(tm)+"\t"+str(gc)+"\t"+str(penal)+"\t"+str(genom_start)+\
				"\t"+str(tar_diffs)+"\t"+str(tar_poly)+"\t"+str(indel_count)+\
				"\t"+str(fail_count)+"\t"+str(added_length)+"\t"+str(tar_len)+\
				"\t"+str(target_pos)+"\t"+amp_seq+"\n"
			output.write(primer_line)
	return


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
	# output_file_prefix = args.of
	masking_library = args.ml
	calc_fixed_SNPs = args.fd

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
	if not os.path.exists(out_dir+'primer_output/'):
		os.makedirs(out_dir+'primer_output/')

	inv_strain_dict = parse_inv_strain_presence(strain_file_path)
	inv_data = parse_inv_data(inversion_breakpoint_file_path)
	(adapt_F,adapt_R,prim_len,min_gap,max_gap,buff,mask_tolerance,prim_buff,num_amps_returned,
		illum_A,illum_B,illum_full_A,illum_full_B) = parse_parameters(parameter_file_path)

	if not args.from_polys:
		if not args.from_counts:

			inv_pos_counts = gen_inv_site_counts(genome_dir,inv_data,inv_strain_dict,out_dir)
			if args.hr:
				write_pos_counts(inv_pos_counts,inv_data,out_dir)

		else:

			inv_pos_counts = read_pos_count_mems(inv_data,out_dir)

		if calc_fixed_SNPs:

			fixed_SNPs = gen_fixed_SNPs(inv_pos_counts,inv_data,out_dir,mask_tolerance)
			write_fixed_SNPs(fixed_SNPs,inv_data,out_dir)
			return

		poly_maps = gen_poly_maps(inv_pos_counts,inv_data,out_dir,mask_tolerance)
		if args.hr:
				write_poly_maps(poly_maps,inv_data,out_dir)

	else:

		poly_maps = load_poly_map_mems(inv_data,out_dir)

	# proposed_regions = propose_prim_regions(poly_maps,inv_data,prim_len,buff)
	# print(proposed_regions)

	region_maps = gen_region_maps(poly_maps,inv_data,prim_len,buff,prim_buff)
	if args.hr:
		write_region_outputs(region_maps,inv_data,out_dir)

	# proposed_amplicon_sets = propose_amplicon_sets(region_maps, inv_data, min_gap, 
	# 	max_gap, prim_len, inv_strain_dict, genome_dir, out_dir)
	if args.hap:
		proposed_amplicon_sets = propose_hap_amplicon_sets(region_maps, inv_data, min_gap, 
			max_gap, prim_len, prim_buff, poly_maps, inv_strain_dict, genome_dir, out_dir)
	else:
		proposed_amplicon_sets = propose_amplicon_sets_poly(region_maps, inv_data, min_gap, 
			max_gap, prim_len, prim_buff, poly_maps, genome_dir, out_dir)
	if args.fa:
		proposed_amplicon_sets = filter_amplicon_sets(proposed_amplicon_sets,
			num_amps_returned)
	if args.hr:
		write_amplicon_sets(proposed_amplicon_sets, inv_data, out_dir)

	gen_primer3_proposals(proposed_amplicon_sets, inv_data, prim_len, min_gap, 
		max_gap, out_dir, primer3_core_directory, masking_library)

	primer_sets = process_primer3_output(inv_data, proposed_amplicon_sets, out_dir)

	write_detailed_primer_proposals(inv_data, primer_sets, out_dir)

	# primer_sets = filter_primers(primer_sets)

	# if args.ol:
	# 	oligo_sets = gen_oligo_sets(left_primers,right_primers,adapt_F,adapt_R)
	# 	write_oligo_sets()

	# if args.iol:
	# 	barcodes = parse_barcodes(barcode_file_path)
	# 	indexed_oligos = gen_indexed_oligos(left_primers,right_primers,barcodes,adapt_F,adapt_R,illum_A,illum_B)
	# 	write_oligo_sets()



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
						help='the Parameters text file, formatted as::\n\
						the illumina adaptor A to append to a forward primer:\n\
						the Illumina adaptor B to append to a reverse primer:\n\
						the primer length:\n\
						the minimum gap length between primers:\n\
						the maximum gap length between primers:\n\
						the buffer distance around breakpoints to exclude from analysis:\n\
						the number of masked bases tolerated before marking the base as masked in the polymorphism map:\n\
						the number of bases added to the primer length to give the minimum fully conserved region considered for primers:\n\
						the maximum number of amplicons to consider from those that pass quality control (this goes in order of generation):\n\
						the illumina adaptor A to append to a forward primer, truncated to sequence not overlapping the stubby adaptor:\n\
						the illumina adaptor B to append to a reverse primer, truncated to sequence not overlapping the stubby adaptor:\n\
						the full illumina adaptor A for a forward primer:\n\
						the full illumina adaptor B for a reverse primer:\n\
						',
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
	parser.add_argument('--od',
						help='candidate region and position frequency files output directory',
						type=str,
						default='primer_output/')
	# parser.add_argument('--of',
	# 					help='output file prefix',
	# 					type=str,
	# 					default='')
	parser.add_argument('--ml',
						help='masking library used by primer3 to account for some mispriming',
						type=str,
						default='primer3/masking_libraries/drosophila_w_transposons.txt')
	parser.add_argument('--from_counts',
						help='a flag set when the inverted and non-inverted \
						sequence counts are already calculated',
						# dest='from_counts',
						action='store_true')
	parser.add_argument('--from_polys',
						help='a flag set when the inverted and non-inverted \
						sequence polymorphism maps are already calculated, \
						supercedes --from_counts',
						action='store_true')
	parser.add_argument('--hr',
						help='a flag set to output human-readable copies of \
						site-base counts, sequence polymorphism, regions, and amplicons',
						action='store_true')
	parser.add_argument('--fa',
						help='a flag set to output a filtered set of amplicons to primer3, \
						instead of all amplicons: outputs a filtered subset of locations',
						action='store_true')
	parser.add_argument('--ol',
						help='a flag set to output short oligomer sequences for each primer, \
						containing just the short adaptors without barcodes',
						action='store_true')
	parser.add_argument('--iol',
						help='a flag set to output long indexed oligomer sequences for each primer, \
						using the unique barcodes from the barcodes file. \
						This is an inflexible approach, as it prevents changing the barcode \
						between sequencing sessions.',
						action='store_true')
	parser.add_argument('--hap',
						help='a flag set to generate amplicons based on haplotype identity, \
						not on fixed SNP differences',
						action='store_true')
	parser.add_argument('--fd',
						help='a flag set to collect all sites without shared SNPs and output them in a readable format',
						action='store_true')

	args = parser.parse_args()

	main(args)

