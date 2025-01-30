#!/bin/bash

#Purpose: For running read sets from pooled amplicons through alignment and QC steps, before estimating frequency

#Run as: 
# ./amplicon_alignment_pipeline.sh <read set name> <left paired end file> <right paired end file> <read dir> <ref genome> <Picard dir> <GATK dir>

#Example run:
# ./../amplicon_alignment_pipeline.sh 251par2L 251par2L_S148_L001_R1_001 251par2L_S148_L001_R2_001 ../reads/ ../../../alignment_software/dmel_ref/DmelRef.fasta ../../../alignment_software/picard-tools-1.79/picard-tools-1.79/ ../../../alignment_software/GenomeAnalysisTK-3.4-46/ ../intermediates/
# ./amplicon_alignment_pipeline.sh 251par2L 251par2L_S148_L001_R1_001 251par2L_S148_L001_R2_001 reads/ ../../alignment_software/dmel_ref/DmelRef.fasta ../../alignment_software/picard-tools-1.79/picard-tools-1.79/ ../../alignment_software/GenomeAnalysisTK-3.4-46/ intermediates/ aligned_bam/


# Script contents:

# Check the number of passed arguments
if [ $# -le 6 ] || [ $# -ge 10 ]
then
  2>&1 printf "\nRequires 7-9 arguments.\nUsage as:\n./amplicon_alignment_pipeline.sh <read set name> <left paired end file> <right paired end file> <read dir> <ref genome> <Picard dir> <GATK dir> <intermediate dir> <output dir>\n"
  exit 1
fi

# Expand the command line arguments
args=("$@")
Read_set=${args[0]}
Left_reads=${args[1]}
Right_reads=${args[2]}
Read_directory=${args[3]}
D_mel_ref=${args[4]}
Picard_directory=${args[5]}
GATK_directory=${args[6]}
intermediate_directory=${args[7]}
output_directory=${args[8]}

# Argument check
printf '\nArgument echo:\n'
echo $Read_set
echo $Left_reads
echo $Right_reads
echo $Read_directory
echo $D_mel_ref
echo $Picard_directory
echo $GATK_directory
echo $intermediate_directory
echo $output_directory
printf '\n'

# Space saving variables
# L_reads=$Read_directory$Left_reads'.fastq'
# R_reads=$Read_directory$Right_reads'.fastq'
L_reads=$Read_directory$Left_reads
R_reads=$Read_directory$Right_reads
L_sai=$intermediate_directory$Left_reads.sai
R_sai=$intermediate_directory$Right_reads.sai
Interm_prefix=$intermediate_directory$Read_set


# Take a time-stamp for log recording
TIMESTAMP=$(date +"%Y-%m-%d_%H-%M-%S")
# Generate the log file
LOG=$Interm_prefix'_aln_'$TIMESTAMP'.log'
printf '\nGenerating log file:\n'
printf $LOG'\n'


# Check the number of read clusters in the input .fastq
printf '\nLine count in one paired read file\n' > $LOG 2>&1
wc -l $L_reads >> $LOG 2>&1


# # Using bwa aln to generate .sam files from paired end reads
# # Generate .sai files for both read directions
# echo bwa aln $D_mel_ref $L_reads ">" $L_sai
# bwa aln $D_mel_ref $L_reads > $L_sai
# echo 'Read direction 1 processed to .sai by bwa mem' >> $LOG 2>&1
# echo bwa aln $D_mel_ref $R_reads '>' $R_sai
# bwa aln $D_mel_ref $R_reads > $R_sai
# echo 'Read direction 2 processed to .sai by bwa mem' >> $LOG 2>&1
# # Combine the read directions into single reads in a sam file
# echo bwa sampe -P $D_mel_ref $L_sai $R_sai $L_reads $R_reads '>' $Interm_prefix.sam
# bwa sampe -P $D_mel_ref $L_sai $R_sai $L_reads $R_reads > $Interm_prefix.sam
# echo 'Reads combined into .sam by bwa sampe' >> $LOG 2>&1


# Using bwa mem to generate .sam files from paired end reads
# REQUIRES running "bwa index reference.fasta" to generate index files first
printf '\n\nUsing bwa mem to generate .sam file from paired end reads\n' >> $LOG 2>&1
echo bwa mem $D_mel_ref $L_reads $R_reads '>' $Interm_prefix.sam
bwa mem $D_mel_ref $L_reads $R_reads > $Interm_prefix.sam
# echo 'Reads combined into .sam by bwa sampe' >> $LOG 2>&1

# Convert the sam file to a bam file
printf '\n\nConverting the sam file to a bam file using samtools view\n' >> $LOG 2>&1
echo samtools view -bS -q 20 $Interm_prefix.sam '>' $Interm_prefix.bam
samtools view -bS -q 20 $Interm_prefix.sam > $Interm_prefix.bam
# echo 'Converted into .bam by samtools view' >> $LOG 2>&1
# echo -e '\n\n' >> $LOG 2>&1



# Run
printf '\n\n' >> $LOG 2>&1
echo java -Xmx5g -jar $Picard_directory'CleanSam.jar' INPUT=$Interm_prefix.bam OUTPUT=$Interm_prefix'_clean.bam'
echo java -Xmx5g -jar $Picard_directory'CleanSam.jar' INPUT=$Interm_prefix.bam OUTPUT=$Interm_prefix'_clean.bam' >> $LOG 2>&1
java -Xmx5g -jar $Picard_directory'CleanSam.jar' INPUT=$Interm_prefix.bam OUTPUT=$Interm_prefix'_clean.bam' >> $LOG 2>&1
printf '\n\n' >> $LOG 2>&1

#
printf '\n\n' >> $LOG 2>&1
echo java -Xmx5g -jar $Picard_directory'SortSam.jar' SORT_ORDER=coordinate INPUT=$Interm_prefix'_clean.bam' OUTPUT=$Interm_prefix'_sort.bam'
echo java -Xmx5g -jar $Picard_directory'SortSam.jar' SORT_ORDER=coordinate INPUT=$Interm_prefix'_clean.bam' OUTPUT=$Interm_prefix'_sort.bam' >> $LOG 2>&1
java -Xmx5g -jar $Picard_directory'SortSam.jar' SORT_ORDER=coordinate INPUT=$Interm_prefix'_clean.bam' OUTPUT=$Interm_prefix'_sort.bam' >> $LOG 2>&1

#
printf '\n\n' >> $LOG 2>&1
echo java -Xmx5g -jar $Picard_directory'AddOrReplaceReadGroups.jar' RGLB=Lane1 RGPL=Illumina RGPU=TTAGGC RGSM=$Read_set INPUT=$Interm_prefix'_sort.bam' OUTPUT=$Interm_prefix'_header.bam'
echo java -Xmx5g -jar $Picard_directory'AddOrReplaceReadGroups.jar' RGLB=Lane1 RGPL=Illumina RGPU=TTAGGC RGSM=$Read_set INPUT=$Interm_prefix'_sort.bam' OUTPUT=$Interm_prefix'_header.bam' >> $LOG 2>&1
java -Xmx5g -jar $Picard_directory'AddOrReplaceReadGroups.jar' RGLB=Lane1 RGPL=Illumina RGPU=TTAGGC RGSM=$Read_set INPUT=$Interm_prefix'_sort.bam' OUTPUT=$Interm_prefix'_header.bam' >> $LOG 2>&1

#
printf '\n\n' >> $LOG 2>&1
echo java -Xmx5g -jar $Picard_directory'BuildBamIndex.jar' INPUT=$Interm_prefix'_header.bam'
echo java -Xmx5g -jar $Picard_directory'BuildBamIndex.jar' INPUT=$Interm_prefix'_header.bam' >> $LOG 2>&1
java -Xmx5g -jar $Picard_directory'BuildBamIndex.jar' INPUT=$Interm_prefix'_header.bam' >> $LOG 2>&1




# 
printf '\n\n' >> $LOG 2>&1
echo java -Xmx5g -jar $GATK_directory'GenomeAnalysisTK.jar' -T RealignerTargetCreator -R $D_mel_ref -I $Interm_prefix'_header.bam' -o $Interm_prefix.intervals
echo java -Xmx5g -jar $GATK_directory'GenomeAnalysisTK.jar' -T RealignerTargetCreator -R $D_mel_ref -I $Interm_prefix'_header.bam' -o $Interm_prefix.intervals >> $LOG 2>&1
java -Xmx5g -jar $GATK_directory'GenomeAnalysisTK.jar' -T RealignerTargetCreator -R $D_mel_ref -I $Interm_prefix'_header.bam' -o $Interm_prefix.intervals >> $LOG 2>&1

# 
printf '\n\n' >> $LOG 2>&1
echo java -Xmx124g -jar $GATK_directory'GenomeAnalysisTK.jar' -T IndelRealigner -R $D_mel_ref -targetIntervals $Interm_prefix.intervals -I $Interm_prefix'_header.bam' -o $output_directory$Read_set'_aln.bam'
echo java -Xmx124g -jar $GATK_directory'GenomeAnalysisTK.jar' -T IndelRealigner -R $D_mel_ref -targetIntervals $Interm_prefix.intervals -I $Interm_prefix'_header.bam' -o $output_directory$Read_set'_aln.bam' >> $LOG 2>&1
java -Xmx40g -jar $GATK_directory'GenomeAnalysisTK.jar' -T IndelRealigner -R $D_mel_ref -targetIntervals $Interm_prefix.intervals -I $Interm_prefix'_header.bam' -o $output_directory$Read_set'_aln.bam' >> $LOG 2>&1



