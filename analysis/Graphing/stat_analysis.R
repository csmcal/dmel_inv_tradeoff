# 
# For  the invworking with the chi-square based statistics used in 

# Set working directory
setwd("/Users/cmcallester/cMcAl/Pool Lab/Inversion ∆ Experiment/Amplicon Analysis/Graphing")
# setwd(getSrcDirectory(function(x) {x})[1])



load_stat_data <- function(stat_data_file){
  stat_data <- read.csv(file=stat_data_file, header=TRUE, sep=",")
  # stat_data$Cohort <- factor(stat_data$Cohort, levels = c('Parental ♂','Embryo','Early ♀','Late ♀','Early ♂','Late ♂'))
  # stat_data$Line <- factor(stat_data$Line)
  return(stat_data)
}

combine_2way_stat_data <- function(cost_ben_data,ben_cost_data){
  pval_name = "Fishers_pval_across_lines"
  short_cost_ben_data <- cost_ben_data[,c(1:8)]
  colnames(short_cost_ben_data) <- c(names(short_cost_ben_data)[1:7],pval_name)
  # short_cost_ben_data$Statistic <- "Chi-square on two comparisons, conditional that the first is closer to a decrease and the second an increase"
  short_cost_ben_data$Tested_direction <- "decrease then increase"
  
  short_ben_cost_data <- ben_cost_data[,c(1:8)]
  colnames(short_ben_cost_data) <- c(names(short_ben_cost_data)[1:7],pval_name)
  # short_cost_ben_data$Statistic <- "Chi-square on two comparisons, conditional that the first is closer to an increase and the second a decrease"
  short_ben_cost_data$Tested_direction <- "increase then decrease"
    
  stat_2way_data <- rbind(short_cost_ben_data,short_ben_cost_data)
  stat_2way_data$Benjamini_Yekutieli <- p.adjust(stat_2way_data$Fishers_pval_across_lines,method="BY")
  stat_2way_data$BY_null_rej_at_05 <- stat_2way_data$Benjamini_Yekutieli < 0.05
  stat_2way_data$Bonferroni <- p.adjust(stat_2way_data$Fishers_pval_across_lines,method="bonferroni")
  stat_2way_data$Bon_null_rej_at_05 <- stat_2way_data$Bonferroni < 0.05
  return(stat_2way_data)
}

# Read two-way  CSV data
stat_data_dir = "../Alignment & QC Scripts/significance/comb_2way_nohup_22-01-10/"
cost_ben_filename = "allele_chi2_combined_two_way_cost_ben_bootstrap_embryo_comb_off_two_way_line-combined_pval.csv"
ben_cost_filename = "allele_chi2_combined_two_way_ben_cost_bootstrap_embryo_comb_off_two_way_line-combined_pval.csv"
cost_ben_data <- load_stat_data(paste0(stat_data_dir,cost_ben_filename))
ben_cost_data <- load_stat_data(paste0(stat_data_dir,ben_cost_filename))
stat_2way_data <- combine_2way_stat_data(cost_ben_data,ben_cost_data)


