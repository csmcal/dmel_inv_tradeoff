#!/bin/bash

#Purpose: For easy running of the alignment and frequency analysis pipeline on the Pool Lab genepool server,
#   with command output recording and nohup built in

#Run as: 
# ./inv_pipe_nohup.sh

# Generate log file:
TIMESTAMP=$(date +"%Y-%m-%d_%H-%M-%S")
LOG='nohup_inv_freq_pipe_'$TIMESTAMP'.log'
printf '\nGenerating log file:'
printf $LOG'\n'

# script_invocation="$(printf %q "$BASH_SOURCE")$((($#)) && printf ' %q' "$@")"
script_name="$(printf %q "$BASH_SOURCE")"
echo $script_name >> $LOG 2>&1
printf 'An instance of inv_freq_pipe_marula_nohup.sh, running the inv alignment/frequency analysis pipeline on Pool lab servers\n' >> $LOG 2>&1
echo $TIMESTAMP >> $LOG 2>&1

# peripheral variables for use in input variable definitions
genepool_aln_soft_dir='/raid10/genepool_4_TB_RAID_Set_1/DPGP2plus/wrap1kb/alignment_software/'
aln_soft=$genepool_aln_soft_dir

# Defined variables to pass:
read_metadata_file='read_metadata.csv'
# read_metadata_file='read_metadata_truncated.csv'
D_mel_ref=$aln_soft'dmel_ref/DmelRef.fasta'
alignment_bash_script='amplicon_alignment_pipeline.sh'
read_directory='/raid10/chris/inv_tradeoff_exp/reads/'
intermediate_directory='aln_intermediates/'
aligned_output_directory='aligned_bam/'
amp_haps_output_directory='aln_hap_counts/'
Picard_directory=$aln_soft'picard-tools-1.79/picard-tools-1.79/'
GATK_directory=$aln_soft'GenomeAnalysisTK-3.4-46/'
amplicon_reference_directory='amp_refs/'
inversion_frequency_output_directory='inv_freqs/'
inv_pipeline_python_script='inv_amplicon_pipeline.py'


# Argument recording
printf '\n\nArgument echo:\n' >> $LOG 2>&1
echo $read_metadata_file >> $LOG 2>&1
echo $D_mel_ref >> $LOG 2>&1
echo $alignment_bash_script >> $LOG 2>&1
echo $read_directory >> $LOG 2>&1
echo $intermediate_directory >> $LOG 2>&1
echo $aligned_output_directory >> $LOG 2>&1
echo $amp_haps_output_directory >> $LOG 2>&1
echo $Picard_directory >> $LOG 2>&1
echo $GATK_directory >> $LOG 2>&1
echo $amplicon_reference_directory >> $LOG 2>&1
echo $inversion_frequency_output_directory >> $LOG 2>&1
echo $inv_pipeline_python_script >> $LOG 2>&1
printf '\n\n' >> $LOG 2>&1

# Echo the nohup call
echo nohup python $inv_pipeline_python_script --read_meta $read_metadata_file --aln_script $alignment_bash_script --d_mel_ref $D_mel_ref --rd $read_directory --id $intermediate_directory --ad $aligned_output_directory --hd $amp_haps_output_directory --pd $Picard_directory --gd $GATK_directory --amp_ref_d $amplicon_reference_directory --od $inversion_frequency_output_directory '>>' $LOG '2>&1 &'

# Record the nohup call to the Log file
echo nohup python $inv_pipeline_python_script --read_meta $read_metadata_file --aln_script $alignment_bash_script --d_mel_ref $D_mel_ref --rd $read_directory --id $intermediate_directory --ad $aligned_output_directory --hd $amp_haps_output_directory --pd $Picard_directory --gd $GATK_directory --amp_ref_d $amplicon_reference_directory --od $inversion_frequency_output_directory '>>' $LOG '2>&1 &' >> $LOG 2>&1

# Call the pipeline and exit
nohup python $inv_pipeline_python_script --read_meta $read_metadata_file --aln_script $alignment_bash_script --d_mel_ref $D_mel_ref --rd $read_directory --id $intermediate_directory --ad $aligned_output_directory --hd $amp_haps_output_directory --pd $Picard_directory --gd $GATK_directory --amp_ref_d $amplicon_reference_directory --od $inversion_frequency_output_directory >> $LOG 2>&1 &



