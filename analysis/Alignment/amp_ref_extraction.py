
# A script for extracting target amplicon sequences from each sequenced genome,
# to generate panels of reference haplotypes for inverted and non-inverted karyotypes
# for identifying polymorphism and inversion frequency in pooled data
# generated via Illumina metagenomic protocols.
#
# Originally run for Inversions in Zambian D. melanogaster, Pool Lab, 2021
#

import numpy as np
import linecache
import argparse
import os



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
	print("\n\nPreparing analysis of:")
	for i in np.arange(len(inv_names_std)):
		print(inv_names_std[i]+': '+str(len(inv_strain_presences[i]))+' strains')
	return inv_strain_dict



def parse_amplicon_data(amplicon_file_path):
	amp_data = []
	with open(amplicon_file_path) as amp_file:
		lines = amp_file.readlines()
		for line in lines[1:]:
			data = line.strip().split(',')
			amp_name = str(data[0])
			chrom = str(data[1])
			inv = str(data[2])
			bp_1 = int(data[3])
			bp_2 = int(data[4])
			inv_name = "In("+chrom+")"+inv
			std_name = "In("+chrom+")Std"
			primer_f = str(data[5])
			primer_r = str(data[6])
			prim_start = int(data[7])
			amp_start = int(data[8])
			amp_end = int(data[9])
			prim_end = int(data[10])
			primer_product_seq = str(data[11])
			primer_product_len = int(data[12])
			oligo_product_seq = str(data[13])
			oligo_product_len = int(data[14])
			amp_data += [[amp_name,chrom,inv,inv_name,std_name,primer_f,primer_r,\
				prim_start,amp_start,amp_end,prim_end,bp_1,bp_2]]
	print('\nAmplicon data and locations to analyze:')
	for amp in amp_data:
		print(amp)
	print('\n')
	return amp_data


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

# A helper function to get a subsequence from 'wrap1kb' file format genomes
def seq_exctraction(strain_file,start_line,end_line,start_offset,end_offset):
	seq = ''
	# account for linecache indexing
	start_line += 1
	end_line += 1
	if start_line == end_line:
		seq = linecache.getline(strain_file,start_line).strip()[start_offset:end_offset]
	else:
		seq += linecache.getline(strain_file,start_line).strip()[start_offset:]
		i = start_line + 1
		while i < end_line:
			seq += linecache.getline(strain_file,i).strip()
			i += 1
		seq += linecache.getline(strain_file,i).strip()[:end_offset]
	return seq


# Access chromosomes to extract a list of target position sequences,
#   specific to the inversion karyotype specified.
#   Assumes the indexes of the target sequence will be in bounds of the file
def amp_seq_extraction(genome_dir,strains,chrom,inv_name,target_start,target_end):
	amp_seqs = []
	# Check the first file for length data
	# SHOULDN'T SIMULTANEOUSLY OPEN THE WHOLE THING :(
	print('\n\nStarting '+inv_name+' amplicon extraction from '+str(target_start)+' to '+str(target_end))

	# Get the line and character position of the amplicon start and end by reading one of the files
	first_strain = strains[0]
	print('Checking file format (line length and total number) from first strain ('+first_strain+')')
	prefix = genome_dir + 'Chr' + chrom + '/' + first_strain + '_Chr' + chrom
	strain_file_name = pick_strain_file(prefix)

	first_line = linecache.getline(strain_file_name,1).strip()
	line_len = len(first_line)

	start_line = target_start//line_len
	start_offset = target_start%line_len
	
	end_line = target_end//line_len
	end_offset = target_end%line_len

	# Extract the first amplicon sequence
	first_seq = seq_exctraction(strain_file_name,start_line,end_line,start_offset,end_offset)
	amp_seqs += [first_seq]
	print('Strain '+first_strain+' extracted sequence: '+first_seq+'\n')
	print('Sequence recorded for strain '+first_strain)

	# Extract haplotypes from the rest of the strain files
	print('Starting '+inv_name+' general strain sequence extraction')
	for strain in strains[1:]:
		prefix = genome_dir + 'Chr' + chrom + '/' + strain + '_Chr' + chrom
		strain_file_name = pick_strain_file(prefix)
		strain_seq = seq_exctraction(strain_file_name,start_line,end_line,start_offset,end_offset)
		amp_seqs += [strain_seq]
		print('Sequence recorded for strain '+strain)

	# Return the extracted haplotypes
	return(amp_seqs)

# Generates the haplotype data for each inversion and amplicon target
# Assumes only two karyotypes, the specific inversion and standard
def gen_inv_amp_seqs(genome_dir,amp_data,inv_strain_dict):
	amp_seqs = []
	for amp_datum in amp_data:
		(amp_name,chrom,inv,inv_name,std_name,primer_f,primer_r,prim_start,\
			amp_start,amp_end,prim_end,bp_1,bp_2) = amp_datum
		strains = inv_strain_dict[inv_name]
		# inv_haps = amp_seq_extraction(genome_dir,strains,chrom,inv_name,amp_start,amp_end)
		inv_seqs = amp_seq_extraction(genome_dir,strains,chrom,inv_name,prim_start,prim_end)
		strains = inv_strain_dict[std_name]
		# std_haps = amp_seq_extraction(genome_dir,strains,chrom,inv_name,amp_start,amp_end)
		std_seqs = amp_seq_extraction(genome_dir,strains,chrom,inv_name,prim_start,prim_end)
		amp_seqs += [(inv_seqs,std_seqs)]
	return amp_seqs


# For writing a list to a file by line, with a header
def write_seq_file(seqs,header,file_name):
	with open(file_name,'w') as out_file:
		out_file.write(header)
		for seq in seqs:
			out_file.write(seq+'\n')
	return

# For writing haplotype sets with a file per amplicon and karyotype
def write_sep_amp_seq_files(amp_seqs, amp_data, out_dir):
	for i in np.arange(len(amp_data)):
		(amp_name,chrom,inv,inv_name,std_name,primer_f,primer_r,prim_start,\
			amp_start,amp_end,prim_end,bp_1,bp_2) = amp_data[i]
		print('Writing amplicon '+str(i)+': '+amp_name+\
			' karyotype orientation '+inv_name+' extracted sequence data')
		file_name = out_dir+amp_name+'_'+inv_name+'_seqs.txt'
		header = '# '+amp_name+','+inv_name+','+primer_f+','+primer_r+','+\
			str(prim_start)+','+str(prim_end)+'\n'
		write_seq_file(amp_seqs[i][0],header,file_name)
		print('Writing amplicon '+str(i)+': '+amp_name+\
			' karyotype orientation '+std_name+' extracted sequence data')
		file_name = out_dir+amp_name+'_'+std_name+'_seqs.txt'
		header = '# '+amp_name+','+std_name+','+primer_f+','+primer_r+','+\
			str(prim_start)+','+str(prim_end)+'\n'
		write_seq_file(amp_seqs[i][1],header,file_name)
	return


# For reading a list with a single line header and one string element per file line
def read_seq_file(file_name):
	seqs = []
	with open(file_name,'r') as seq_file:
		lines = seq_file.readlines()
		for line in lines[1:]:
			seqs += [line.strip()]
	return seqs

# For writing haplotype sets with a file per amplicon and karyotype
def read_sep_amp_seq_files(amp_data, out_dir):
	amp_seqs = []
	for i in np.arange(len(amp_data)):
		(amp_name,chrom,inv,inv_name,std_name,primer_f,primer_r,prim_start,\
			amp_start,amp_end,prim_end,bp_1,bp_2) = amp_data[i]
		print('Reading amplicon '+str(i)+': '+amp_name+\
			' karyotype orientation '+inv_name+' extracted sequence data')
		file_name = out_dir+amp_name+'_'+inv_name+'_seqs.txt'
		inv_seqs = read_seq_file(file_name)
		print('Reading amplicon '+str(i)+': '+amp_name+\
			' karyotype orientation '+std_name+' extracted sequence data')
		file_name = out_dir+amp_name+'_'+std_name+'_seqs.txt'
		std_seqs = read_seq_file(file_name)
		amp_seqs += [(inv_seqs,std_seqs)]
	return amp_seqs

# Generates the haplotype count data from a set of passed sequences
def gen_hap_counts(seqs):
	haps = {}
	for seq in seqs:
		haps[seq] = haps.get(seq, 0) + 1
	return haps

# Generates the haplotype count data for each inversion and amplicon target,
#   from passed amplicon sequences
# Assumes only two karyotypes, the specific inversion and standard
def gen_inv_amp_hap_counts(amp_seqs,amp_data):
	amp_haps = []
	for i in np.arange(len(amp_data)):
		(amp_name,chrom,inv,inv_name,std_name,primer_f,primer_r,prim_start,\
			amp_start,amp_end,prim_end,bp_1,bp_2) = amp_data[i]
		inv_haps = gen_hap_counts(amp_seqs[i][0])
		std_haps = gen_hap_counts(amp_seqs[i][1])
		amp_haps += [(inv_haps,std_haps)]
	return amp_haps

# For writing a list to a file by line, with a header
def write_hap_file(haps,header,file_name):
	with open(file_name,'w') as out_file:
		out_file.write(header)
		for hap, count in haps.items():
			out_file.write(hap+'\t'+str(count)+'\n')
	return

# For writing haplotype sets with a file per amplicon and karyotype
def write_sep_amp_hap_files(amp_haps, amp_data, out_dir):
	for i in np.arange(len(amp_data)):
		(amp_name,chrom,inv,inv_name,std_name,primer_f,primer_r,prim_start,\
			amp_start,amp_end,prim_end,bp_1,bp_2) = amp_data[i]
		print('Writing amplicon '+str(i)+': '+amp_name+\
			' karyotype orientation '+inv_name+' haplotype sequence data')
		file_name = out_dir+amp_name+'_'+inv_name+'_haps.txt'
		header = '# '+amp_name+','+inv_name+','+primer_f+','+primer_r+','+\
			str(prim_start)+','+str(prim_end)+'\n'
		write_hap_file(amp_haps[i][0],header,file_name)
		print('Writing amplicon '+str(i)+': '+amp_name+\
			' karyotype orientation '+std_name+' haplotype sequence data')
		file_name = out_dir+amp_name+'_'+std_name+'_haps.txt'
		header = '# '+amp_name+','+std_name+','+primer_f+','+primer_r+','+\
			str(prim_start)+','+str(prim_end)+'\n'
		write_hap_file(amp_haps[i][1],header,file_name)
	return



# A dictionary translating bases to integers in ATGCN
#   ACGT is the standard for PFMs,
#   but ATGCN is the format of the Pool lab wrap1kb genomes
base_dict = {
	'A': 0,
	'T': 1,
	'G': 2,
	'C': 3,
	'N': 4
}

# For generating position frequency representations from sequence sets
# Assumes all provided sequences are the same length
# Outputs a PFM as a list, where each element is a sequence position with 5 counts:
#   [A; T; G; C; N] - as specified in base_dict
def gen_pos_freqs(seqs):
	pos_freqs = []
	n = len(seqs)
	l = len(seqs[0])
	for i in np.arange(l):
		base_count = [0,0,0,0,0]
		for j in np.arange(n):
			base_count[base_dict[seqs[j][i]]] += 1
		pos_freqs += [base_count]
	return pos_freqs

# Generates Position Frequency Matrices (PFMs) for all amplicon sequence sets
def convert_seqs_to_PFM(amp_seqs):
	amp_PFMs = []
	for seqs_pair in amp_seqs:
		amp_PFMs += [(gen_pos_freqs(seqs_pair[0]),gen_pos_freqs(seqs_pair[1]))]
	return amp_PFMs


# For writing individual position frequency files
def write_pos_frequency_file(pos_freq,file_name):
	with open(file_name,'w') as out_file:
		for pos in pos_freq:
			line = str(pos[0])
			for base_count in pos[1:]:
				line += '\t'+str(base_count)
			line += '\n'
			out_file.write(line)
	return

# For writing the site-base counts in individual files
def write_PFMs(amp_PFMs,amp_data,out_dir):
	for i in np.arange(len(amp_data)):
		(amp_name,chrom,inv,inv_name,std_name,primer_f,primer_r,prim_start,\
			amp_start,amp_end,prim_end,bp_1,bp_2) = amp_data[i]
		print('Writing amplicon '+str(i)+': '+amp_name+' karyotype orientation '\
			+inv_name+' position frequency matrix data')
		file_name = out_dir+amp_name+'_'+inv_name+'_PFM.txt'
		write_pos_frequency_file(amp_PFMs[i][0],file_name)
		print('Writing amplicon '+str(i)+': '+amp_name+' karyotype orientation '\
			+std_name+' position frequency matrix data')
		file_name = out_dir+amp_name+'_'+std_name+'_PFM.txt'
		write_pos_frequency_file(amp_PFMs[i][1],file_name)
	return


# Find and return the offset, genomic position, and base distribution
#   of fixed differences between the two PFMs
def calc_fixed_diffs_no_N(PFM_1,PFM_2,genomic_index):
	# base_order = 'ATGCN'
	l = len(PFM_1)
	b = len(PFM_1[0])
	fixed_diff_data = []
	rev_base_dict = {i: base for base, i in base_dict.items()}
	for i in np.arange(l):
		# print('Base comparison '+str(i)+':')
		# print(PFM_1[i])
		# print(PFM_2[i])
		shared_bases = np.logical_and(PFM_1[i],PFM_2[i])
		num_shared = np.sum(shared_bases)
		if num_shared == 0:
			inv_bases = ''
			for j in np.arange(b):
				if PFM_1[i][j]:
					inv_bases += rev_base_dict[j]
			std_bases = ''
			for j in np.arange(b):
				if PFM_2[i][j]:
					std_bases += rev_base_dict[j]
			diff_data = (i,i+genomic_index,PFM_1[i],PFM_2[i],inv_bases,std_bases)
			fixed_diff_data += [diff_data]
	return fixed_diff_data


# For generating a list of fixed differences between karyotypes and their genomic positions,
#   for each amplicon PFM pair
def gen_amp_no_N_fixed_diffs(amp_PFMs,amp_data):
	fixed_diffs = []
	for i in np.arange(len(amp_data)):
		(amp_name,chrom,inv,inv_name,std_name,primer_f,primer_r,prim_start,\
			amp_start,amp_end,prim_end,bp_1,bp_2) = amp_data[i]
		print('Calculating amplicon '+str(i)+': '+amp_name+' fixed difference data, \
			excluding fixed differences with masked positions')
		fixed_diffs += [calc_fixed_diffs_no_N(amp_PFMs[i][0],amp_PFMs[i][1],prim_start)]
	return fixed_diffs



# For writing individual fixed difference files
def write_fixed_diff_file(fixed_diffs,file_name):
	with open(file_name,'w') as out_file:
		header = 'Offset\tGenomic Position\tInv Dist\tStd Dist\tInv Base\tStd Base\n'
		out_file.write(header)
		for diff in fixed_diffs:
			line = str(diff[0])
			for datum in diff[1:]:
				line += '\t'+str(datum)
			line += '\n'
			out_file.write(line)
			print(line)
	return

# For writing a list of fixed differences between karyotypes and their genomic positions,
#   for fixed differences calculated between each amplicon karyotype pair
def write_no_N_fixed_diffs(fixed_diffs,amp_data,out_dir):
	for i in np.arange(len(amp_data)):
		(amp_name,chrom,inv,inv_name,std_name,primer_f,primer_r,prim_start,\
			amp_start,amp_end,prim_end,bp_1,bp_2) = amp_data[i]
		print('Writing amplicon '+str(i)+': '+amp_name+' fixed difference data, \
			excluding fixed differences with masked positions')
		file_name = out_dir+amp_name+'_noNfixdif.txt'
		write_fixed_diff_file(fixed_diffs[i],file_name)
	return

# Find and return the offset, genomic position, and base distribution
#   of fixed differences between the two PFMs
def calc_fixed_diffs(PFM_1,PFM_2,genomic_index):
	# base_order = 'ATGCN'
	l = len(PFM_1)
	b = len(PFM_1[0])
	fixed_diff_data = []
	rev_base_dict = {i: base for base, i in base_dict.items()}
	for i in np.arange(l):
		# print('Base comparison '+str(i)+':')
		# print(PFM_1[i][:4])
		# print(PFM_2[i][:4])
		shared_bases = np.logical_and(PFM_1[i][:4],PFM_2[i][:4])
		num_shared = np.sum(shared_bases)
		if num_shared == 0:
			inv_bases = ''
			for j in np.arange(b):
				if PFM_1[i][j]:
					inv_bases += rev_base_dict[j]
			std_bases = ''
			for j in np.arange(b):
				if PFM_2[i][j]:
					std_bases += rev_base_dict[j]
			# print('Difference identified at '+str(i)+':')
			# print(PFM_1[i])
			# print(PFM_2[i])
			diff_data = (i,i+genomic_index,PFM_1[i],PFM_2[i],inv_bases,std_bases)
			fixed_diff_data += [diff_data]
	return fixed_diff_data


# For generating a list of fixed differences between karyotypes and their genomic positions,
#   for each amplicon PFM pair
def gen_amp_fixed_diffs(amp_PFMs,amp_data):
	fixed_diffs = []
	for i in np.arange(len(amp_data)):
		(amp_name,chrom,inv,inv_name,std_name,primer_f,primer_r,prim_start,\
			amp_start,amp_end,prim_end,bp_1,bp_2) = amp_data[i]
		print('Calculating amplicon '+str(i)+': '+amp_name+' fixed difference data')
		fixed_diffs += [calc_fixed_diffs(amp_PFMs[i][0],amp_PFMs[i][1],prim_start)]
	return fixed_diffs


# For writing a list of fixed differences between karyotypes and their genomic positions,
#   for fixed differences calculated between each amplicon karyotype pair
def write_fixed_diffs(fixed_diffs,amp_data,out_dir):
	for i in np.arange(len(amp_data)):
		(amp_name,chrom,inv,inv_name,std_name,primer_f,primer_r,prim_start,\
			amp_start,amp_end,prim_end,bp_1,bp_2) = amp_data[i]
		print('Writing amplicon '+str(i)+': '+amp_name+' fixed difference data')
		file_name = out_dir+amp_name+'_fixdif.txt'
		write_fixed_diff_file(fixed_diffs[i],file_name)
	return


# Returns a boolean check on whether the two base distributions are not fixed for the same base
#   ignores sites with any masked bases - this perhaps could instead be handled in hap calling later?
def is_informative(base_dist_1,base_dist_2):
	if (base_dist_1[4] == 0) and (base_dist_2[4] == 0):
		base_presence = np.logical_or(base_dist_1[0:4],base_dist_2[0:4])
		num_present = np.sum(base_presence)
		if num_present > 1:
			# print("Found informative site:")
			# print(base_dist_1)
			# print(base_dist_2)
			return True
	return False

# Generates a list of positions over which the two PFMs may harbor haplotype differences
def gen_informative_pos(PFM_1,PFM_2):
	# base_order = 'ATGCN'
	l = len(PFM_1)
	b = len(PFM_1[0])
	informative_positions = []
	for i in np.arange(l):
		# print("Checking pos "+str(i))
		base_dist_1 = PFM_1[i]
		base_dist_2 = PFM_2[i]
		if is_informative(base_dist_1,base_dist_2):
			informative_positions += [i]
	return informative_positions

# Generates a set of haplotypes for amplicon sequences from both karyotypes, 
#   but only considering those positions flagged as informative
def gen_informative_haps(positions,inv_seqs,std_seqs):
	inv_info_haps = set()
	for seq in inv_seqs:
		info_seq = ''
		for pos in positions:
			info_seq += seq[pos]
		inv_info_haps.add(info_seq)

	std_info_haps = set()
	for seq in std_seqs:
		info_seq = ''
		for pos in positions:
			info_seq += seq[pos]
		std_info_haps.add(info_seq)

	return (inv_info_haps,std_info_haps)


# Generates a list of positions over which the two PFMs may harbor haplotype differences,
#   and DNA strings composed of the unique sequences (haplotypes) in the pool only considering those positions
def gen_informative_pos_haps(amp_seqs,amp_PFMs):
	info_pos_haps = []
	for i in np.arange(len(amp_PFMs)):
		# print("Generating informative haps for amp "+str(i))
		informative_positions = gen_informative_pos(amp_PFMs[i][0],amp_PFMs[i][1])
		(inv_info_haps,std_info_haps) = gen_informative_haps(informative_positions,
			amp_seqs[i][0],amp_seqs[i][1])
		info_pos_haps += [(informative_positions,inv_info_haps,std_info_haps)]
	return info_pos_haps


# For writing individual informative-position haplotype files
def write_info_pos_hap_file(info_pos,info_haps,file_name):
	with open(file_name,'w') as out_file:
		# header = '# '+amp_name+'\n'
		out_file.write(str(info_pos) + '\n')
		for hap in info_haps:
			out_file.write(hap + '\n')
	return

# For writing informative-position haplotype files for each amplicon and karyotype
def write_informative_pos_haps(info_pos_haps,amp_data,out_dir):
	for i in np.arange(len(amp_data)):
		(amp_name,chrom,inv,inv_name,std_name,primer_f,primer_r,prim_start,\
			amp_start,amp_end,prim_end,bp_1,bp_2) = amp_data[i]
		print('Writing amplicon '+str(i)+': '+amp_name+' karyotype orientation '\
			+inv_name+' informative position haplotype data')
		file_name = out_dir+amp_name+'_'+inv_name+'_info_haps.txt'
		write_info_pos_hap_file(info_pos_haps[i][0],info_pos_haps[i][1],file_name)
		print('Writing amplicon '+str(i)+': '+amp_name+' karyotype orientation '\
			+std_name+' informative position haplotype data')
		file_name = out_dir+amp_name+'_'+std_name+'_info_haps.txt'
		write_info_pos_hap_file(info_pos_haps[i][0],info_pos_haps[i][2],file_name)
	return

# For handling direct calls to amp_ref_extraction.py,
#   and to allow specifying different input files
# Assumes only two inversion karyotypes, the specific inversion and standard
def main(args):
	strain_file_path = args.s
	amplicon_data_file_path = args.a
	genome_dir = args.d
	out_dir = args.od
	# output_file_prefix = args.of
	read_from_seq_file = args.rf

	assert (strain_file_path!=''),'Strain karyotype detail file path (--s) {} \
		not given'.format(strain_file_path)
	assert (amplicon_data_file_path!=''),'Amplicon detail file path (--i) {} \
		not given'.format(amplicon_data_file_path)

	if out_dir[-1] != '/':
		out_dir += '/'
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)


	amp_data = parse_amplicon_data(amplicon_data_file_path)

	amp_seqs = []
	if read_from_seq_file:
		amp_seqs = read_sep_amp_seq_files(amp_data,out_dir)
	else:
		inv_strain_dict = parse_inv_strain_presence(strain_file_path)
		amp_seqs = gen_inv_amp_seqs(genome_dir,amp_data,inv_strain_dict)
		write_sep_amp_seq_files(amp_seqs, amp_data, out_dir)

	amp_haps = gen_inv_amp_hap_counts(amp_seqs,amp_data)
	write_sep_amp_hap_files(amp_haps, amp_data, out_dir)

	amp_PFMs = convert_seqs_to_PFM(amp_seqs)
	write_PFMs(amp_PFMs, amp_data, out_dir)

	amp_fixed_diffs_no_N = gen_amp_no_N_fixed_diffs(amp_PFMs,amp_data)
	write_no_N_fixed_diffs(amp_fixed_diffs_no_N,amp_data,out_dir)
	
	amp_fixed_diffs = gen_amp_fixed_diffs(amp_PFMs,amp_data)
	write_fixed_diffs(amp_fixed_diffs,amp_data,out_dir)

	info_pos_haps = gen_informative_pos_haps(amp_seqs,amp_PFMs)
	write_informative_pos_haps(info_pos_haps,amp_data,out_dir)

	return


# Runs if amp_ref_extraction.py is the main script
if __name__ == "__main__":
	# Collect input arguments
	parser = argparse.ArgumentParser(
		description=__doc__,
		formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument('--s',
						help='the input csv list of Strain haploid genomes and inversion presences',
						type=str,
						default='inv_strain_presence.csv')
	parser.add_argument('--a',
						help='the input csv list of amplicons, their positions, and associated inversions',
						type=str,
						default='amplicon_data.csv')
	parser.add_argument('--d',
						help='the Directory in which the genomes can be found',
						type=str,
						default='')
						# '/Volumes/4_TB_RAID_Set_1/DPGP2plus/wrap1kb/'
	parser.add_argument('--od',
						help='output directory, for haplotype sequences and position frequency files',
						type=str,
						default='amp_refs/')
	parser.add_argument('--read_seq_file', dest='rf', default=False, action='store_true')
	# parser.add_argument('--of',
	# 					help='output file prefix',
	# 					type=str,
	# 					default='')

	args = parser.parse_args()

	main(args)


		