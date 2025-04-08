
# Purpose: To provide functions for the extraction of unique read pairs and their counts from aligned bam files

# Notes: 
#   Separated from the rest of the pipeline script so that pysam and multiprocessing dependency
#    is only checked and called when needed

# Import libraries
import multiprocessing as mp
from collections import defaultdict
import pysam




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


# A dictionary translating chromosome names to reference file sequence headers
#   based on the Drosophila melanogaster reference genome
#   from the Drosophila genome Nexus, reference genome release 5
ref_chrom_dict = {
	'Yhet': 1,
	'mtDNA': 2,
	'2L': 3,
	'X': 4,
	'3L': 5,
	'4': 6,
	'2R': 7,
	'3R': 8,
	'Uextra': 9,
	'2RHet': 10,
	'2LHet': 11,
	'3LHet': 12,
	'3RHet': 13,
	'U': 14,
	'XHet': 15,
	'Wolbachia': 16
}

# For caching reads while waiting for the paired read in paired-end read sets
#   adapted from a code snippet in a biostars post by ggydush
# NEEDS an update to handling multiply-mapped reads, so >2 entries per read pair
#   currently just throws those out
def read_pair_generator(bam, contig, start, end):
	"""
	Generate read pairs in a BAM file or within a region string.
	Reads are added to read_dict until a pair is found.
	"""
	read_dict = defaultdict(lambda: [None, None])
	for read in bam.fetch(contig, start, end):
		if not read.is_proper_pair or read.is_secondary or read.is_supplementary:
			continue
		qname = read.query_name
		if qname not in read_dict:
			if read.is_read1:
				read_dict[qname][0] = read
			else:
				read_dict[qname][1] = read
		else:
			if read.is_read1:
				yield read, read_dict[qname][1]
			else:
				yield read_dict[qname][0], read
			del read_dict[qname]


# Iterator/generator for caching reads from a BAM while waiting for the paired read in paired-end read sets,
#   throws out secondary alignments and reads without pairs
def read_pair_generator_tally(bam, contig, start, end, pool_name, snp_positions=None):
	num_unpaired = 0
	num_secondary = 0
	num_low_qual = 0
	num_low_mapq = 0
	num_low_snp_qual = 0
	read_dict = defaultdict(lambda: [None, None])
	
	for read in bam.fetch():
		# Count and throw out unpaired and secondary alignments
		if not read.is_proper_pair:
			num_unpaired += 1
			continue
		if read.is_secondary or read.is_supplementary:
			num_secondary += 1
			continue
			
		# Filter reads with mapping quality <= 20
		if read.mapping_quality <= 20:
			num_low_mapq += 1
			continue
			
		# Check base quality at SNP positions if provided
		if snp_positions is not None:
			qualities = read.query_qualities
			positions = read.get_aligned_pairs()
			if qualities is not None and positions is not None:
				low_quality_snp = False
				for pos in snp_positions:
					# Find the read position that maps to this reference position
					for read_pos, ref_pos in positions:
						if ref_pos == pos and read_pos is not None:
							if qualities[read_pos] <= 20:
								low_quality_snp = True
								break
					if low_quality_snp:
						break
				if low_quality_snp:
					num_low_snp_qual += 1
					continue
		
		# Trim bases with quality < 20 from both ends
		qualities = read.query_qualities
		if qualities is not None:
			# Find first position with quality >= 20
			start_trim = 0
			while start_trim < len(qualities) and qualities[start_trim] < 20:
				start_trim += 1
				
			# Find last position with quality >= 20
			end_trim = len(qualities) - 1
			while end_trim >= 0 and qualities[end_trim] < 20:
				end_trim -= 1
				
			# If entire read is low quality, skip it
			if start_trim > end_trim:
				num_low_qual += 1
				continue
				
			# Trim the read sequence and qualities
			read.query_sequence = read.query_sequence[start_trim:end_trim+1]
			read.query_qualities = qualities[start_trim:end_trim+1]
		
		# Reads are added to read_dict until a pair is found
		qname = read.query_name
		if qname not in read_dict:
			if read.is_read1:
				read_dict[qname][0] = read
			else:
				read_dict[qname][1] = read
		else:
			if read.is_read1:
				yield read, read_dict[qname][1]
			else:
				yield read_dict[qname][0], read
			del read_dict[qname]
			
	print("removed "+str(num_unpaired)+" unpaired, "+str(num_secondary)+\
		" secondary alignments, and "+str(num_low_mapq+num_low_qual)+\
		" low mapping quality reads, and "+str(num_low_snp_qual)+\
		" reads with low quality at SNP positions in library "+pool_name)


# Opens the aligned bam file and generates a dictionary
#   with key (read1,read2) and value the int count of that read pair sequence
def extract_haps(metadata,aligned_output_directory,snp_positions=None):
	# Generate the aligned bam file name and access it with pysam
	pool_name = metadata[8]
	print('Extracting '+pool_name)
	aligned_reads_file = aligned_output_directory+pool_name+'_aln.bam'
	bamfile = pysam.AlignmentFile(aligned_reads_file, "rb")

	# Generate the SAM region string to pull correctly-mapped reads
	chrom = metadata[2]
	amp_start = metadata[16] - 1
	amp_end = metadata[17] - 1
	contig = str(ref_chrom_dict[chrom])
	region = str(contig)+':'+str(amp_start)+'-'+str(amp_end)

	# Collect counts for each unique read pair sequence
	read_haps = {}
	pair_count = 0
	for read_1, read_2 in read_pair_generator_tally(bamfile,contig,amp_start,amp_end,pool_name,snp_positions):
		r1_seq = read_1.get_forward_sequence()
		r2_seq = reverse_complement(read_2.get_forward_sequence())
		read_haps[(r1_seq,r2_seq)] = read_haps.get((r1_seq,r2_seq),0) + 1
		pair_count += 1
	bamfile.close()
	print('counted '+str(pair_count)+' read pairs in pool '+pool_name)
	return read_haps

# Opens the aligned bam file and generates a dictionary
#   with key (read1,read2) and value the int count of that read pair sequence
def extract_read_hap_sets(read_metadata,aligned_output_directory,snp_positions=None):
	read_hap_sets = []
	print('Extracting paired end read unique sequences and counts:')
	for metadata in read_metadata:
		read_haps = extract_haps(metadata,aligned_output_directory,snp_positions)
		read_hap_sets += [read_haps]
	print('\n')
	return read_hap_sets


# To take a read pair of pysam.AlignedSegment objects and generate the consensus sequence
def generate_amp_seq_from_pair_2(r1,r2,amp_start,amp_end):
	# print("generating amp by simultaneous iteration")
	amp_length = amp_end-amp_start

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
	# print("Assembling reads:")
	# # print(r1)
	# # print(r1.get_aligned_pairs())
	# print(len(r1_pos_pairs))
	# # print(r2)
	# # print(r2.get_aligned_pairs())
	# print(len(r2_pos_pairs))
	# print(r1_seq)
	# print(r2_seq)
	# print("amp start  - "+str(amp_start))
	# print("amp end    - "+str(amp_end))
	# print("amp length - "+str(amp_end-amp_start))
	# print()

	# Build the consensus span:
	combined_read = ''
	r1_i = 0
	r2_i = 0
	r1_pair_i = 0
	r2_pair_i = 0 # assumes that the aligned pairs will always be in alignment order
	r1_pair_l = len(r1_pos_pairs)
	r2_pair_l = len(r2_pos_pairs)
	num_overlap_mismatched_bases = 0
	for i in range(amp_length):
		# print()
		# print("offset "+str(i))
		# print("r1_ind "+str(r1_pair_i))
		# print("r2_ind "+str(r2_pair_i))
		# print(combined_read)
		# if r1_pair_i < r1_pair_l:
		# 	print("r1_pair "+str(r1_pos_pairs[r1_pair_i]))
		# if r2_pair_i < r2_pair_l:
		# 	print("r2_pair "+str(r2_pos_pairs[r2_pair_i]))
		# print()

		# Handle gaps or late starts in alignment
		while r1_pair_i < r1_pair_l and (r1_pos_pairs[r1_pair_i][1] is None or r1_pos_pairs[r1_pair_i][1] < i+amp_start):
				r1_pair_i += 1
				r1_i += 1
		while r2_pair_i < r2_pair_l and (r2_pos_pairs[r2_pair_i][1] is None or r2_pos_pairs[r2_pair_i][1] < i+amp_start):
				r2_pair_i += 1
				r2_i += 1

		# Handle reaching the end of the alignments
		if r1_pair_i == r1_pair_l:
			if r2_pair_i == r2_pair_l:
				# If past the end of both, just add a masked N
				combined_read += 'N'
			elif r2_pos_pairs[r2_pair_i][1] > i+amp_start:
				# If the remaining 2nd read aligns farther down, mask w/ N
				combined_read += 'N'
			else:
				# Assumes that otherwise r2_pos_pairs[r2_pair_i][1] == i+amp_start
				if r2_pos_pairs[r2_pair_i][0] is None:
					combined_read += 'N'
					r2_pair_i += 1
				else:
					# If there is a mapped base, add it
					combined_read += r2_seq[r2_i]
					r2_pair_i += 1
					r2_i += 1
		elif r2_pair_i == r2_pair_l:
			if r1_pos_pairs[r1_pair_i][1] > i+amp_start:
				# If the remaining 1st read aligns farther down, mask w/ N
				combined_read += 'N'
			else:
				# Assumes that otherwise r1_pos_pairs[r1_pair_i][1] == i+amp_start
				if r1_pos_pairs[r1_pair_i][0] is None:
					combined_read += 'N'
					r1_pair_i += 1
				else:
					# If there is a mapped base, add it
					combined_read += r1_seq[r1_i]
					r1_pair_i += 1
					r1_i += 1
		elif r1_pos_pairs[r1_pair_i][1] > i+amp_start:
			if r2_pos_pairs[r2_pair_i][1] > i+amp_start:
				# If both sequence pairs have no base, mask it
				combined_read += 'N'
			# if r2_pos_pairs[r2_pair_i][1] == i+amp_start:
			else:
				# Assumes that otherwise r2_pos_pairs[r2_pair_i][1] == i+amp_start
				if r2_pos_pairs[r2_pair_i][0] is None:
					combined_read += 'N'
					r2_pair_i += 1
				else:
					# If there is a mapped base, add it
					combined_read += r2_seq[r2_i]
					r2_pair_i += 1
					r2_i += 1
		elif r2_pos_pairs[r2_pair_i][1] > i+amp_start:
			# Assumes that otherwise r1_pos_pairs[r1_pair_i][1] == i+amp_start
			if r1_pos_pairs[r1_pair_i][0] is None:
				combined_read += 'N'
				r1_pair_i += 1
			else:
				# If there is a mapped base, add it
				combined_read += r1_seq[r1_i]
				r1_pair_i += 1
				r1_i += 1
		else:
			# Assumes that otherwise both have an alignment here
			if r1_pos_pairs[r1_pair_i][0] is None:
				if r2_pos_pairs[r2_pair_i][0] is not None:
					# The 2nd has something mapped, but the 1st has an indel
					r2_i += 1
					num_overlap_mismatched_bases += 1
				combined_read += 'N'
			elif r2_pos_pairs[r2_pair_i][0] is None:
				# The 1st has something mapped, but the 2nd has an indel
				r1_i += 1
				num_overlap_mismatched_bases += 1
				combined_read += 'N'
			# Check whether both agree
			elif r1_seq[r1_i] == r2_seq[r2_i]:
				combined_read += r1_seq[r1_i]
				r1_i += 1
				r2_i += 1
			else:
				combined_read += 'N'
				num_overlap_mismatched_bases += 1
				r1_i += 1
				r2_i += 1
			r1_pair_i += 1
			r2_pair_i += 1

	# print()
	# print(combined_read)
	# print(num_overlap_mismatched_bases)
	return (combined_read,num_overlap_mismatched_bases)

# To take a read pair of pysam.AlignedSegment objects and generate the consensus sequence
#    Doesn't 
def generate_amp_seq_from_pair(r1,r2,amp_start,amp_end):
	amp_length = amp_end-amp_start
	# print("generating amp by separate iteration")

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
	# print("Assembling reads:")
	# # print(r1)
	# # print(r1.get_aligned_pairs())
	# print(len(r1_pos_pairs))
	# # print(r2)
	# # print(r2.get_aligned_pairs())
	# print(len(r2_pos_pairs))
	# print(r1_seq)
	# print(r2_seq)
	# print("amp start  - "+str(amp_start))
	# print("amp end    - "+str(amp_end))
	# print("amp length - "+str(amp_end-amp_start))
	# print()

	# Build the consensus span:
	combined_read_chars = ['N' for i in range(amp_length)]
	modified = [False for i in range(amp_length)]
	num_overlap_mismatched_bases = 0

	# Modify the combined read positions based on the first pair
	for pair in r1_pos_pairs:
		if pair[1] is None:
			continue
		elif pair[1] >= amp_start and pair[1] < amp_end:
			amp_pos = pair[1]-amp_start
			if pair[0] is not None:
				combined_read_chars[amp_pos] = r1_seq[pair[0]]
			modified[amp_pos] = True

	# Modify from the second pair and check for mismatch
	for pair in r2_pos_pairs:
		if pair[1] is None:
			continue
		elif pair[1] >= amp_start and pair[1] < amp_end:
			amp_pos = pair[1]-amp_start
			if pair[0] is not None:
				if modified[amp_pos] and (combined_read_chars[amp_pos] != r2_seq[pair[0]]):
					combined_read_chars[amp_pos] = 'N'
					num_overlap_mismatched_bases += 1
				else:
					combined_read_chars[amp_pos] = r2_seq[pair[0]]
			else:
				if modified[amp_pos] and (combined_read_chars[amp_pos] != 'N'):
					combined_read_chars[amp_pos] = 'N'
					num_overlap_mismatched_bases += 1
			modified[amp_pos] = True

	combined_read = ''.join(combined_read_chars)

	# print()
	# print(combined_read)
	# print(num_overlap_mismatched_bases)

	return (combined_read,num_overlap_mismatched_bases)


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

# Opens the aligned bam file and generates a dictionary
#   with key (read1,read2) and value the int count of that read pair sequence
def extract_combined_amp_seqs(metadata,aligned_output_directory,amp_ref_dir):
	# Generate the aligned bam file name and access it with pysam
	pool_name = metadata[8]
	print('Extracting '+pool_name)
	aligned_reads_file = aligned_output_directory+pool_name+'_aln.bam'
	bamfile = pysam.AlignmentFile(aligned_reads_file, "rb")

	# Generate the SAM region string to pull correctly-mapped reads
	chrom = metadata[2]
	amp_start = metadata[16]
	amp_end = metadata[17]
	contig = str(ref_chrom_dict[chrom])
	region = str(contig)+':'+str(amp_start)+'-'+str(amp_end)

	# Get the positions of the SNPs from the amplicon reference sequence file
	amp_name = metadata[6]
	amp_fixed_diff_file = amp_ref_dir+amp_name+'_fixdif.txt'
	snp_positions = [diff[1] for diff in read_fixed_diff_file(amp_fixed_diff_file)]

	# Collect counts for each unique combined read pair sequence
	amp_haps = {}
	pair_count = 0
	num_overlap_mismatched_bases = 0
	num_mismatched_overlaps = 0
	for read_1, read_2 in read_pair_generator_tally(bamfile,contig,amp_start,amp_end,pool_name,snp_positions):
		(combined_amp, n_mis_b) = generate_amp_seq_from_pair(read_1,read_2,amp_start,amp_end)
		if n_mis_b:
			num_mismatched_overlaps += 1
			num_overlap_mismatched_bases += n_mis_b
		amp_haps[combined_amp] = amp_haps.get(combined_amp,0) + 1
		pair_count += 1
		# print('recorded '+str(pair_count)+'th read pair')
		# if pair_count > 5:
		# 	raise KeyboardInterrupt()
	bamfile.close()

	# Print a brief summary of extraction, mismatched overlapping segments
	if num_mismatched_overlaps:
		print(pool_name+":  "+str(num_mismatched_overlaps/pair_count)+" of segments in read pairs have mismatch in overlap: "+\
			str(num_overlap_mismatched_bases/num_mismatched_overlaps)+" avg mismatches in "+str(num_mismatched_overlaps)+" mismatched overlaps in "+\
			str(pair_count)+' read pairs in pool')
	else:
		print(pool_name+":  No mismatched overlapping segments in read pairs, counted "+str(pair_count)+' read pairs in pool '+pool_name)
	return amp_haps


# Opens the aligned bam file and generates a dictionary
#   with key (read1,read2) and value the int count of that read pair sequence
def extract_read_amp_hap_sets(read_metadata,aligned_output_directory):
	read_amp_hap_sets = []
	print('Extracting paired end read unique sequences and counts:')
	for metadata in read_metadata:
		amp_haps = extract_haps(metadata,aligned_output_directory)
		read_amp_hap_sets += [amp_haps]
	print('\n')
	return read_amp_hap_sets



# For writing a paired read count list to a file by line, with a header
def write_read_hap_file(haps,header,file_name):
	with open(file_name,'w') as out_file:
		out_file.write(header)
		for (r1,r2), count in haps.items():
			out_file.write(r1+'\t'+r2+'\t'+str(count)+'\n')
	return


# For writing a paired read count list to a file by line,
#   and composing the header
def write_read_hap_file_with_header(haps,metadata,out_dir):
	pool_name = metadata[8]
	amp_name = metadata[6]
	chrom = metadata[2]
	amp_start = metadata[16]
	amp_end = metadata[17]

	header = '# '+pool_name+','+amp_name+','+chrom+','+\
		str(amp_start)+','+str(amp_end)+': Read_1, Read_2, Count\n'
	file_name = out_dir+pool_name+'_read_seqs.txt'

	print('Writing read pool '+pool_name+\
		' paired read haplotype (unique sequence) count data')
	print(file_name)
	write_read_hap_file(haps,header,file_name)
	return


# For writing a paired read count list to a file by line, with a header
def write_read_amp_hap_file(haps,header,file_name):
	with open(file_name,'w') as out_file:
		out_file.write(header)
		for amp_seq, count in haps.items():
			out_file.write(amp_seq+'\t'+str(count)+'\n')
	return

# For writing a combined amplicon sequence count list to a file by line,
#   and composing the header
def write_read_amp_hap_file_with_header(haps,metadata,out_dir):
	pool_name = metadata[8]
	amp_name = metadata[6]
	chrom = metadata[2]
	amp_start = metadata[16]
	amp_end = metadata[17]

	header = '# '+pool_name+','+amp_name+','+chrom+','+\
		str(amp_start)+','+str(amp_end)+': Amplicon Read, Count\n'
	file_name = out_dir+pool_name+'_amp_haps.txt'

	print('Writing read pool '+pool_name+\
		' paired read extrated amplicon count data')
	print(file_name)
	write_read_amp_hap_file(haps,header,file_name)

# For writing sequencing read haplotype count data with a file per read pool
def write_sep_read_hap_files(read_hap_sets, read_metadata, out_dir):
	for i in range(len(read_metadata)):
		write_read_hap_file_with_header(read_hap_sets[i],read_metadata[i],out_dir)
	print('\n')
	return

# Handles extraction and writing of unique read counts,
#  returns an order variable along with the unique read/hap counts
#  to allow order-preserving sorting during asynchronous execution
def extract_and_write_flag_order(metadata,aligned_output_directory,out_dir,i):
	# Assemble the haplotype/unique read counts
	read_haps = extract_haps(metadata,aligned_output_directory)

	# Write the output to file
	write_read_hap_file_with_header(read_haps,metadata,out_dir)

	# Return the read output with the order flag for sorting
	return (i,read_haps)


# For founding a multiprocessing pool to asynchronously extract and record unique read pair counts,
#   and to return the unique read (haplotype) count data in correct order, using callbacks
def extract_and_write_async_pool_callbacks(read_metadata,aligned_output_directory,out_dir):
	read_amp_hap_sets = []

	# Prepare the multiprocessing pool, leave at least 2?
	pool = mp.Pool(mp.cpu_count()-2)

	# Collect and write the unique read counts with apply_async()
	async_objects = []
	print('Extracting paired end read unique sequences and counts, in asynchronous parallel:')
	for metadata in read_metadata:
		async_objects += [pool.apply_async(extract_and_write,args=(metadata,aligned_output_directory,out_dir))]
	for i in range(len(read_metadata)):
		read_amp_haps = extract_and_write_flag_order(read_metadata[i],aligned_output_directory,i)
		read_amp_hap_sets += [read_amp_haps]

	# Sort and return the unique read/hap counts
	asynchronous_outputs.sort(key=lambda x: x[0])
	read_hap_sets = [haps for i, haps in asynchronous_outputs]
	print('\n')
	return read_hap_sets



# Handles extraction and writing of unique read counts
def extract_and_write(metadata,aligned_output_directory,out_dir,amp_ref_dir):
	# Assemble the haplotype/unique read counts
	read_amp_haps = extract_combined_amp_seqs(metadata,aligned_output_directory,amp_ref_dir)

	# Write the output to file
	write_read_amp_hap_file_with_header(read_amp_haps,metadata,out_dir)

	# Return the read output
	return read_amp_haps


# For founding a multiprocessing pool to asynchronously extract and record unique read pair counts,
#   and to return the unique read (haplotype) count data in correct order
def extract_and_write_async_pool(read_metadata,aligned_output_directory,out_dir,amp_ref_dir):
	read_amp_hap_sets = []

	# Prepare the multiprocessing pool, leave at least 2?
	pool = mp.Pool(mp.cpu_count()-2)

	# Collect and write the unique read counts with apply_async()
	async_objects = []
	print('Extracting paired end read unique sequences and counts, in asynchronous parallel:')
	for metadata in read_metadata:
		async_objects += [pool.apply_async(extract_and_write,args=(metadata,aligned_output_directory,out_dir,amp_ref_dir))]

	# Wait for the results with .get()
	read_amp_hap_sets = [result.get() for result in async_objects]

	# Close the pool and return the read haplotype sets
	pool.close()
	pool.join()
	print('\n')
	# print('Async read hap outputs?')
	# print(read_hap_sets)
	# print('\n')
	return read_amp_hap_sets




