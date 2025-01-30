
# A scratch python script for testing pysam functionality

import pysam

inv_start = 11276564
inv_end = 16163902
amp_start = 11352014
amp_end = 11352305
chrom = '2R'
contig = 7

read_file_name = "test_aligned/192par2R_aln.bam"


complement_dict = {
	'A': 'T',
	'G': 'C',
	'C': 'G',
	'T': 'A',
	'N': 'N'
}

def reverse_complement(seq):
	rev_comp = ''
	for character in seq:
		rev_comp = complement_dict[character] + rev_comp
	return rev_comp

def print_read_info(read):
	print(read)
	print()

	# print("cigar string:")
	# print(read.cigarstring)

	# print("query name:")
	# print(read.query_name)

	print("aligned pairs:")
	print(read.get_aligned_pairs())

	print("gapless blocks:")
	print(read.get_blocks())

	# print("forward qualities:")
	# print(read.get_forward_qualities())

	print("forward sequence:")
	print(read.get_forward_sequence())

	print("reverse sequence:")
	print(reverse_complement(read.get_forward_sequence()))

	# print("overlap:")
	# print(read.get_overlap(amp_start,amp_end))

	print("reference positions:")
	print(read.get_reference_positions())

	print("reference sequence:")
	print(read.get_reference_sequence())

	print("reference ID:")
	print(read.reference_id)

	print("is reversed:")
	print(read.is_reverse)

	print("is secondary:")
	print(read.is_secondary)

	print("is supplementary:")
	print(read.is_supplementary)

	print("is paired:")
	print(read.is_paired)

	# print("tag(?):")
	# print(read.get_tag('MQ'))

	# print("tags(?):")
	# print(read.get_tags())

	# print("infer query length(?):")
	# print(read.infer_query_length())

	# print("infer read length(?):")
	# print(read.infer_read_length())
	return


# read_file = pysam.AlignmentFile(read_file_name,"rb")
# print("Reference name for id 6: "+read_file.get_reference_name(6))

# print('\n')
# read_file = pysam.AlignmentFile(read_file_name,"rb")
# unspec_read_iter = read_file.fetch()
# first_unspec_read = next(unspec_read_iter)
# print("First aligned read in the bam file for maternal Zi192N embryos:")
# print_read_info(first_unspec_read)


# print('\n')
# second_unspec_read = next(unspec_read_iter)
# print("Second aligned read in the bam file for maternal Zi192N embryos:")
# print_read_info(second_unspec_read)


# print('\n')
# read_file = pysam.AlignmentFile(read_file_name,"rb")
# read_iter = read_file.fetch(contig=str(contig),start=amp_start,end=amp_end)
# first_read = next(read_iter)
# print("First read aligned to the In(2R)NS amplicon region (spec by contig, start, end) of maternal Zi192N embryos:")
# print_read_info(first_read)



# print('\n')
# read_file = pysam.AlignmentFile(read_file_name,"rb")
# region_string_total = str(contig)+':'+str(amp_start)+'-'+str(amp_end)
# reg_read_iter = read_file.fetch(region=region_string_total)
# first_reg_read = next(reg_read_iter)
# print("First read aligned to the In(2R)NS amplicon region (spec by reg string) of maternal Zi192N embryos:")
# print_read_info(first_reg_read)



# print('\n')
# read_file = pysam.AlignmentFile(read_file_name,"rb")
# rev_read_iter = read_file.fetch(contig=str(contig),start=(amp_end-10),stop=amp_end)
# a_rev_read = next(rev_read_iter)
# print("First read aligned to 10bp from the end of the In(2R)NS amplicon region (spec by contig, start, end) of maternal Zi192N embryos:")
# print_read_info(a_rev_read)



# print('\n')
# read_file = pysam.AlignmentFile(read_file_name,"rb")
# region_string_end = str(contig)+':'+str(amp_end-10)+'-'+str(amp_end)
# print(region_string_end)
# reg_rev_read_iter = read_file.fetch(region=region_string_end)
# first_reg_rev_read = next(reg_rev_read_iter)
# print("First read aligned to 10bp from the end of the In(2R)NS amplicon region (spec by reg string) of maternal Zi192N embryos:")
# print_read_info(first_reg_rev_read)



# sorted_file_name = read_file_name[:-4]+"_sort.bam"
# pysam.sort("-n","-o",sorted_file_name,read_file_name)
# name_sorted_file = pysam.AlignmentFile(sorted_file_name,"rb")


# print('\n')
# sort_iter = name_sorted_file.fetch()
# first_sort_read = next(sort_iter)
# print("First read by name in the bam file for maternal Zi192N embryos:")
# print_read_info(first_sort_read)

# print('\n')
# second_sort_read = next(sort_iter)
# print("Second read by name in the bam file for maternal Zi192N embryos:")
# print_read_info(second_sort_read)


def test_read_pair_comb(read_1,read_2,amp_start,amp_end):
	from extract_read_haps import generate_amp_seq_from_pair
	print("Running amplicon imputation on the first read aligning to the amplicon region:")
	print("Read 1 reference sequence:")
	print(read_1.get_reference_sequence())
	# print("Read 1 reference positions:")
	# print(read_1.get_reference_positions())
	print("Read 1 aligned pairs:")
	print(read_1.get_aligned_pairs())
	print("Read 2 reference sequence:")
	print(read_2.get_reference_sequence())
	# print("Read 2 reference positions:")
	# print(read_2.get_reference_positions())
	print("Read 2 aligned pairs:")
	print(read_2.get_aligned_pairs())
	(combined_read, n_mis_b) = generate_amp_seq_from_pair(read_1,read_2,amp_start,amp_end)
	print("\nCombined read - length "+str(len(combined_read))+", expected amplicon length "+str(amp_end-amp_start)+":")
	print(combined_read)
	return combined_read


def test_192par_read_comb():
	print('\n')
	test_file_2L = "test_aligned/192par2L_aln.bam"
	start_2L = 2238742
	end_2L = 2238967
	contig_2L = 3
	test_file_2R = "test_aligned/192par2R_aln.bam"
	start_2R = 11352014
	end_2R = 11352305
	contig_2R = 7
	test_file_3L = "test_aligned/192par3L_aln.bam"
	start_3L = 2375229
	end_3L = 2375544
	contig_3L = 5
	test_file_3R = "test_aligned/192par3R_aln.bam"
	start_3R = 7609417
	end_3R = 7609723
	contig_3R = 8

	print("testing the first 5 reads of 192par2L")
	read_file = pysam.AlignmentFile(test_file_2L,"rb")
	read_iter = read_file.fetch(contig=str(contig_2L),start=start_2L,end=end_2L)
	read_1 = next(read_iter)
	read_2 = read_file.mate(read_1)
	test_read_pair_comb(read_1,read_2,start_2L,end_2L)
	read_1 = next(read_iter)
	read_2 = read_file.mate(read_1)
	test_read_pair_comb(read_1,read_2,start_2L,end_2L)
	read_1 = next(read_iter)
	read_2 = read_file.mate(read_1)
	test_read_pair_comb(read_1,read_2,start_2L,end_2L)
	read_1 = next(read_iter)
	read_2 = read_file.mate(read_1)
	test_read_pair_comb(read_1,read_2,start_2L,end_2L)
	read_1 = next(read_iter)
	read_2 = read_file.mate(read_1)
	test_read_pair_comb(read_1,read_2,start_2L,end_2L)

	print("testing the first 5 reads of 192par2R")
	read_file = pysam.AlignmentFile(test_file_2R,"rb")
	read_iter = read_file.fetch(contig=str(contig_2R),start=start_2R,end=end_2R)
	read_1 = next(read_iter)
	read_2 = read_file.mate(read_1)
	test_read_pair_comb(read_1,read_2,start_2R,end_2R)
	read_1 = next(read_iter)
	read_2 = read_file.mate(read_1)
	test_read_pair_comb(read_1,read_2,start_2R,end_2R)
	read_1 = next(read_iter)
	read_2 = read_file.mate(read_1)
	test_read_pair_comb(read_1,read_2,start_2R,end_2R)
	read_1 = next(read_iter)
	read_2 = read_file.mate(read_1)
	test_read_pair_comb(read_1,read_2,start_2R,end_2R)
	read_1 = next(read_iter)
	read_2 = read_file.mate(read_1)
	test_read_pair_comb(read_1,read_2,start_2R,end_2R)

	print("testing the first 5 reads of 192par3L")
	read_file = pysam.AlignmentFile(test_file_3L,"rb")
	read_iter = read_file.fetch(contig=str(contig_3L),start=start_3L,end=end_3L)
	read_1 = next(read_iter)
	read_2 = read_file.mate(read_1)
	test_read_pair_comb(read_1,read_2,start_3L,end_3L)
	read_1 = next(read_iter)
	read_2 = read_file.mate(read_1)
	test_read_pair_comb(read_1,read_2,start_3L,end_3L)
	read_1 = next(read_iter)
	read_2 = read_file.mate(read_1)
	test_read_pair_comb(read_1,read_2,start_3L,end_3L)
	read_1 = next(read_iter)
	read_2 = read_file.mate(read_1)
	test_read_pair_comb(read_1,read_2,start_3L,end_3L)
	read_1 = next(read_iter)
	read_2 = read_file.mate(read_1)
	test_read_pair_comb(read_1,read_2,start_3L,end_3L)

	print("testing the first 5 reads of 192par3R")
	read_file = pysam.AlignmentFile(test_file_3R,"rb")
	read_iter = read_file.fetch(contig=str(contig_3R),start=start_3R,end=end_3R)
	read_1 = next(read_iter)
	read_2 = read_file.mate(read_1)
	test_read_pair_comb(read_1,read_2,start_3R,end_3R)
	read_1 = next(read_iter)
	read_2 = read_file.mate(read_1)
	test_read_pair_comb(read_1,read_2,start_3R,end_3R)
	read_1 = next(read_iter)
	read_2 = read_file.mate(read_1)
	test_read_pair_comb(read_1,read_2,start_3R,end_3R)
	read_1 = next(read_iter)
	read_2 = read_file.mate(read_1)
	test_read_pair_comb(read_1,read_2,start_3R,end_3R)
	read_1 = next(read_iter)
	read_2 = read_file.mate(read_1)
	test_read_pair_comb(read_1,read_2,start_3R,end_3R)

	return


def print_read_comb_data(r1,r2,comb):
	# Extract the data of the first read:
	r1_seq = r1.get_forward_sequence()
	if r1.is_reverse:
		r1_seq = reverse_complement(r1_seq)
	r1_pos_pairs = r1.get_aligned_pairs()
	# Extract the data of the second read:
	r2_seq = r2.get_forward_sequence()
	if r2.is_reverse:
		r2_seq = reverse_complement(r2_seq)
	r2_pos_pairs = r2.get_aligned_pairs()

	# print(r1)
	# print(r2)

	print("Read sequences:")
	print(r1_seq)
	print(r2_seq)

	print("Combined sequence:")
	print(comb)

	print("\nRead 1 reference sequence:")
	print(r1.get_reference_sequence())
	print("Read 2 reference sequence:")
	print(r2.get_reference_sequence())
	print("Read 1 aligned pairs:")
	print(r1_pos_pairs)
	print("Read 2 aligned pairs:")
	print(r2_pos_pairs)

	print()

	return

def get_mult_align(read_set,fixed_diffs,buff_len):
	positions = [diff[0] for diff in fixed_diffs]
	extracted_seqs = []
	for (r1,r2,comb) in read_set:
		extr_seq = ''
		for i in range(len(positions)):
			if i == 0:
				extr_seq += comb[positions[i]-buff_len:positions[i]].lower()
				extr_seq += comb[positions[i]]
			else:
				extr_seq += comb[positions[i-1]+1:positions[i]].lower()
				extr_seq += comb[positions[i]]
		extr_seq += comb[positions[-1]+1:positions[-1]+buff_len].lower()
		extracted_seqs += [extr_seq]
	return extracted_seqs




def collect_report_mismatch(bam_file,fixed_diffs,contig,start,end,pool,num,report_recomb=False):
	print('Fixed differences expected: Offset in amp; Offset in genome; Base distribution in inv ref; Base distribution in std ref; Base(s) in inv ref; Base(s) in std ref')
	for diff in fixed_diffs:
		diff_line = str(diff[0])
		for d in diff[1:]:
			diff_line+='\t'+str(d)
		print(diff_line)

	from extract_read_haps import generate_amp_seq_from_pair
	from extract_read_haps import read_pair_generator_tally

	inv_count = 0
	std_count = 0
	recombinant_count = 0
	mismapped_count = 0
	multiply_mapped_count = 0

	inv_read_set = []
	std_read_set = []
	recomb_read_set = []
	mismap_read_set = []
	multimap_read_set = []

	bam = pysam.AlignmentFile(bam_file,"rb")
	for (r1,r2) in read_pair_generator_tally(bam, contig, start, end, pool):
		(combined_read, n_mis_b) = generate_amp_seq_from_pair(r1,r2,start,end)

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
				if combined_read[offset] == base:
					match_inv = True
			for base in std_bases:
				if combined_read[offset] == base:
					match_std = True
			full_match_inv = full_match_inv and match_inv
			full_match_std = full_match_std and match_std
			if not match_inv and not match_std:
				mismatch = True


		# Assign the read based on matching to fixed differences
		if full_match_inv and full_match_std:
			# print("WARNING - reads should NEVER fully match both sets of fixed differences. Coding error likely.")
			# For fixed differences, this should never happen, as the bases should never be shared
			if multiply_mapped_count < num:
				if (r1,r2,combined_read) not in multimap_read_set:
					multimap_read_set += [(r1,r2,combined_read)]
					multiply_mapped_count += 1
		elif full_match_inv:
			if inv_count < num:
				if (r1,r2,combined_read) not in inv_read_set:
					inv_read_set += [(r1,r2,combined_read)]
					# print("\n\n\nInv read "+str(inv_count)+':')
					# print_read_comb_data(r1,r2,combined_read)
					inv_count += 1
		elif full_match_std:
			if std_count < num:
				if (r1,r2,combined_read) not in std_read_set:
					std_read_set += [(r1,r2,combined_read)]
					# print("\n\n\nStd read "+str(std_count)+':')
					# print_read_comb_data(r1,r2,combined_read)
					std_count += 1
		elif mismatch:
			# Something may be wrong - the base observed wasn't expected in either karyotype
			if mismapped_count < num:
				if (r1,r2,combined_read) not in mismap_read_set:
					mismap_read_set += [(r1,r2,combined_read)]
					# print("\n\n\nStd read "+str(mismapped_count)+':')
					# print_read_comb_data(r1,r2,combined_read)
					mismapped_count += 1
		else:
			if recombinant_count < num:
				if (r1,r2,combined_read) not in recomb_read_set:
					recomb_read_set += [(r1,r2,combined_read)]
					# print("\n\n\nMixed read "+str(recombinant_count)+':')
					# print_read_comb_data(r1,r2,combined_read)
					recombinant_count += 1

		if inv_count >= num and std_count >= num and mismapped_count >= num:
			if report_recomb:
				if recombinant_count >= num:
					break
			else:
				break

	print('\nExample called segments at the site(s) of fixed difference(s):')
	buff_len = 8
	print('Inverted:')
	for seq in get_mult_align(inv_read_set,fixed_diffs,buff_len):
		print(seq)
	print('Standard:')
	for seq in get_mult_align(std_read_set,fixed_diffs,buff_len):
		print(seq)
	print('Recombinant:')
	for seq in get_mult_align(recomb_read_set,fixed_diffs,buff_len):
		print(seq)
	print('SNP Mismatching:')
	for seq in get_mult_align(mismap_read_set,fixed_diffs,buff_len):
		print(seq)

	print('\n\nExample full reads and their alignments:')
	print('Inverted:')
	for (r1,r2,combined_read) in inv_read_set:
		print_read_comb_data(r1,r2,combined_read)
	print('Standard:')
	for (r1,r2,combined_read) in std_read_set:
		print_read_comb_data(r1,r2,combined_read)
	print('Recombinant:')
	for (r1,r2,combined_read) in recomb_read_set:
		print_read_comb_data(r1,r2,combined_read)
	print('SNP Mismatching:')
	for (r1,r2,combined_read) in mismap_read_set:
		print_read_comb_data(r1,r2,combined_read)
	return



def test_2L_alignments(bam_2L_path):
	print("testing a selection of 2L reads")
	start_2L = 2238742
	end_2L = 2238967
	contig_2L = 3

	fixed_diffs = [(169,2238911,[0, 44, 0, 0, 0],[149, 0, 0, 0, 0],'T','A')]

	collect_report_mismatch(bam_2L_path,fixed_diffs,str(contig_2L),start_2L,end_2L,bam_2L_path,5)

	return


def test_2R_alignments(bam_2R_path):
	print("testing a selection of 2R reads")
	start_2R = 11352014
	end_2R = 11352305
	contig_2R = 7

	fixed_diffs = [(208, 11352222, [0, 24, 0, 0, 0], [0, 0, 0, 169, 0], 'T', 'C'),
		(247, 11352261, [24, 0, 0, 0, 0], [0, 0, 169, 0, 0], 'A', 'G')]

	collect_report_mismatch(bam_2R_path,fixed_diffs,str(contig_2R),start_2R,end_2R,bam_2R_path,5,True)

	return


def test_3L_alignments(bam_3L_path):
	print("testing a selection of 3L reads")
	start_3L = 2375229
	end_3L = 2375544
	contig_3L = 5

	fixed_diffs = [(247, 2375476, [0, 35, 0, 0, 0], [0, 0, 0, 158, 0], 'T', 'C')]

	collect_report_mismatch(bam_3L_path,fixed_diffs,str(contig_3L),start_3L,end_3L,bam_3L_path,5)

	return


def test_3R_alignments(bam_3R_path):
	print("testing a selection of 3R reads")
	start_3R = 7609417
	end_3R = 7609723
	contig_3R = 8

	fixed_diffs = [(158, 7609575, [0, 34, 0, 0, 0], [159, 0, 0, 0, 0], 'T', 'A')]

	collect_report_mismatch(bam_3R_path,fixed_diffs,str(contig_3R),start_3R,end_3R,bam_3R_path,5)

	return


def check_2R_hap_classes(bam_2R_path):
	print("Counting the proportion of calls of In(2R)NS as TA (Inv), CG (Std), TG (Recombinant), CA (Recombinant), and Other at fixed difference bases in reads from "+bam_2R_path)
	start_2R = 11352014
	end_2R = 11352305
	contig_2R = 7

	fixed_diffs = [(208, 11352222, [0, 24, 0, 0, 0], [0, 0, 0, 169, 0], 'T', 'C'),
		(247, 11352261, [24, 0, 0, 0, 0], [0, 0, 169, 0, 0], 'A', 'G')]


	print('\n\nFixed differences expected: Offset in amp; Offset in genome; Base distribution in inv ref; Base distribution in std ref; Base(s) in inv ref; Base(s) in std ref')
	for diff in fixed_diffs:
		diff_line = str(diff[0])
		for d in diff[1:]:
			diff_line+='\t'+str(d)
		print(diff_line)

	start = start_2R
	end = end_2R
	contig = contig_2R
	bam_file = bam_2R_path
	pool = bam_2R_path

	from extract_read_haps import generate_amp_seq_from_pair
	from extract_read_haps import read_pair_generator_tally

	TA_count = 0
	CG_count = 0
	TG_count = 0
	CA_count = 0
	other_count = 0

	bam = pysam.AlignmentFile(bam_file,"rb")
	for (r1,r2) in read_pair_generator_tally(bam, contig, start, end, pool):
		(combined_read, n_mis_b) = generate_amp_seq_from_pair(r1,r2,start,end)

		# Check the read span against the reference fixed difference data
		if combined_read[208] == 'T':
			if combined_read[247] == 'A':
				TA_count += 1
			elif combined_read[247] == 'G':
				TG_count += 1
			else:
				other_count += 1
		elif combined_read[208] == 'C':
			if combined_read[247] == 'A':
				CA_count += 1
			elif combined_read[247] == 'G':
				CG_count += 1
			else:
				other_count += 1
		else:
			other_count += 1

	prop_TA = TA_count/(TA_count+CG_count+TG_count+CA_count+other_count)
	prop_CG = CG_count/(TA_count+CG_count+TG_count+CA_count+other_count)
	prop_TG = TG_count/(TA_count+CG_count+TG_count+CA_count+other_count)
	prop_CA = CA_count/(TA_count+CG_count+TG_count+CA_count+other_count)
	prop_other = other_count/(TA_count+CG_count+TG_count+CA_count+other_count)
	print('\n\nBase Identity\t\tProportion\t\tCount')
	print('TA (Inv)\t\t'+str(prop_TA)+'\t'+str(TA_count))
	print('CG (Std)\t\t'+str(prop_CG)+'\t'+str(CG_count))
	print('TG (Recombinant)\t'+str(prop_TG)+'\t'+str(TG_count))
	print('CA (Recombinant)\t'+str(prop_CA)+'\t'+str(CA_count))
	print('Other\t\t\t'+str(prop_other)+'\t'+str(other_count))
	print()

	return


"""
RUNNING SPECIFIC TESTS
"""

# test_192par_read_comb()


# test_2L_alignments("test_aligned/192par2L_aln.bam")
# test_2R_alignments("test_aligned/192par2R_aln.bam")
# test_3L_alignments("test_aligned/192par3L_aln.bam")
# test_3R_alignments("test_aligned/192par3R_aln.bam")


check_2R_hap_classes("test_aligned/192par2R_aln.bam")


