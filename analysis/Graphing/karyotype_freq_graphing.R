# 
# For plotting the inversion/karyotype frequency results

# Set working directory
setwd("/Users/cmcallester/cMcAl/UWM/Pool/dmel_inv_tradeoff/Analysis/Graphing")
# setwd(getSrcDirectory(function(x) {x})[1])

library(tidyverse)


# Function for getting the cutoff values 
cutoff_summary <- function(reps, alpha=0.025) {
  cutoffs = data.frame(Pat_Exp_Freq=double(), Coverage=double(), 
                                   P_E_Low_Cutoff=double(), P_E_High_Cutoff=double(), 
                                   E_A_Low_Cutoff=double(), E_A_High_Cutoff=double())
  for (f in unique(reps$Pat_Exp_Freq)) {
    for (c in unique(reps$Coverage)) {
      temp_reps <- reps[reps$Pat_Exp_Freq == f & reps$Coverage == c, ]
      p_e_diffs <- temp_reps$Pat_Freq/2.0 - temp_reps$Emb_Freq
      e_a_diffs <- temp_reps$Emb_Freq - temp_reps$Adu_Freq
      
      # To avoid losing column names when adding the first row
      cutoffs[nrow(cutoffs)+1,] <- c(f,c,quantile(p_e_diffs,alpha)[[1]], 
                                     quantile(p_e_diffs,1.0-alpha)[[1]],
                                     quantile(e_a_diffs,alpha)[[1]], 
                                     quantile(e_a_diffs,1.0-alpha)[[1]])
      # cutoffs = rbind(cutoffs,c(f,c,quantile(p_e_diffs,alpha)[[1]],
      #                           quantile(p_e_diffs,1.0-alpha)[[1]],
      #                           quantile(e_a_diffs,alpha)[[1]],
      #                           quantile(e_a_diffs,1.0-alpha)[[1]]))
    }
  }
  return(cutoffs)
}


# Plot generator for cutoff frequency differences between,
# each line is a different pat. freq, plotting coverage by cutoff change in freq
plot_bootstrap_cutoffs <- function(cutoffs, alpha=0.025) {
  
  # Load ggplot library
  library(ggplot2)
  
  # Plotting the paternal-embryo differences using ggplot2
  p <- ggplot() + 
    geom_line(data=cutoffs,aes(x=Coverage, y=P_E_Low_Cutoff, group=factor(Pat_Exp_Freq), 
                               linetype = 'solid', colour=Pat_Exp_Freq)) +
    geom_line(data=cutoffs,aes(x=Coverage, y=P_E_High_Cutoff, group=factor(Pat_Exp_Freq), 
                               linetype = 'dashed', colour=Pat_Exp_Freq)) +
    xlab("Coverage (per allele)") + 
    ylab("Frequency Difference") + 
    scale_x_continuous(breaks=seq(0, 3, .5), minor_breaks=seq(0, 3, .1),
                       limits=c(-.1,3.1)) +
    scale_y_continuous(breaks=seq(-.2, .2, .05), minor_breaks=seq(-.2, .2, .01),
                       limits=c(-.21,.21)) +
    scale_linetype_manual(values = c(1,2),
                          labels=c("Upper",
                                   "Lower"),
                          name="Bound") +
    scale_colour_continuous(name="Paternal Frequency", 
                            guide = guide_colourbar(title.position = "top",direction = "horizontal")) + 
    ggtitle(paste("Paternal-Embryo Bootstrap Differences at ",alpha,"quantile extremes")) +
    theme_bw() +
    theme(legend.justification=c(1,0),
          legend.position=c(0.95,0.05),  # Position legend in bottom right
          legend.box="horizontal",  # Position multiple legends horizontally
          legend.background = element_rect(fill="white",size=0.5, linetype="solid", colour ="black"),
          legend.title = element_text(size=10),
          legend.text = element_text(size=7),
          plot.title = element_text(hjust=1.0, size=13, face="bold"))
  ggsave(paste("P-E Bootstrap Diffs at",alpha,".png"), 
         plot=p,units="in", width=6, height=4, dpi=600, device = 'png')
  
  # Plotting the embryo-adult differences using ggplot2
  p <- ggplot() + 
    geom_line(data=cutoffs,aes(x=Coverage, y=E_A_Low_Cutoff, group=factor(Pat_Exp_Freq), 
                               linetype = 'solid', colour=Pat_Exp_Freq)) +
    geom_line(data=cutoffs,aes(x=Coverage, y=E_A_High_Cutoff, group=factor(Pat_Exp_Freq), 
                               linetype = 'dashed', colour=Pat_Exp_Freq)) +
    xlab("Coverage (per allele)") + 
    ylab("Frequency Difference") + 
    scale_x_continuous(breaks=seq(0, 3, .5), minor_breaks=seq(0, 3, .1),
                       limits=c(-.1,3.1)) +
    scale_y_continuous(breaks=seq(-.2, .2, .05), minor_breaks=seq(-.2, .2, .01),
                       limits=c(-.21,.21)) +
    scale_linetype_manual(values = c(1,2),
                          labels=c("Upper",
                                   "Lower"),
                          name="Bound") +
    scale_colour_continuous(name="Paternal Frequency", 
                            guide = guide_colourbar(title.position = "top",direction = "horizontal")) + 
    ggtitle(paste("Embryo-Adult Bootstrap Differences at ",alpha,"quantile extremes")) +
    theme_bw() +
    theme(legend.justification=c(1,0),
          legend.position=c(0.95,0.05),  # Position legend in bottom right
          legend.box="horizontal",  # Position multiple legends horizontally
          legend.background = element_rect(fill="white",size=0.5, linetype="solid", colour ="black"),
          legend.title = element_text(size=10),
          legend.text = element_text(size=7),
          plot.title = element_text(hjust=1.0, size=13, face="bold"))
  ggsave(paste("E-A Bootstrap Diffs at",alpha,".png"), 
         plot=p,units="in", width=6, height=4, dpi=600, device = 'png')
}



# Plot generator for cutoff frequency differences 
plot_bootstrap_freq_sensitivity <- function(cutoffs, alpha=0.025, coverage=0.3) {
  
  # Get the cutoff values for just the desired coverage
  cutoffs = cutoffs[cutoffs$Coverage == coverage, ]
  
  # Load ggplot library
  library(ggplot2)
  
  # Plotting the paternal-embryo differences using ggplot2
  p <- ggplot() + 
    geom_line(data=cutoffs,aes(x=Pat_Exp_Freq, y=P_E_Low_Cutoff, #group=factor(Coverage), 
                               linetype='solid', colour='Pat-Emb')) +
    geom_line(data=cutoffs,aes(x=Pat_Exp_Freq, y=P_E_High_Cutoff, #group=factor(Coverage), 
                               linetype='dashed', colour='Pat-Emb')) +
    geom_line(data=cutoffs,aes(x=Pat_Exp_Freq, y=E_A_Low_Cutoff, #group=factor(Coverage), 
                               linetype='solid', colour='Emb-Adu')) +
    geom_line(data=cutoffs,aes(x=Pat_Exp_Freq, y=E_A_High_Cutoff, #group=factor(Coverage), 
                               linetype='dashed', colour='Emb-Adu')) +
    xlab("Paternal Frequency)") + 
    ylab("Frequency Change") + 
    # scale_x_continuous(breaks=seq(0, 3, .5), minor_breaks=seq(0, 3, .1),
    #                    limits=c(-.1,3.1)) +
    # scale_y_continuous(breaks=seq(-.2, .2, .05), minor_breaks=seq(-.2, .2, .01),
    #                    limits=c(-.21,.21)) +
    scale_linetype_manual(values = c(1,2),
                          labels=c("Upper",
                                   "Lower"),
                          name="Bound") +
    # scale_colour_discrete(values = c('brown3','cyan'),
    #                       # labels=c("Upper",
    #                       #          "Lower"),
    #                       name="Comparison") +
    scale_colour_discrete(name="Comparison") +
    ggtitle(paste("Paternal-Embryo Bootstrap ",alpha,"quantile changes,",coverage,"read depth/allele")) +
    theme_bw() +
    theme(legend.justification=c(1,0),
          legend.position=c(0.95,0.35),  # Position legend in bottom right
          legend.box="horizontal",  # Position multiple legends horizontally
          legend.background = element_rect(fill="white",size=0.5, linetype="solid", colour ="black"),
          legend.title = element_text(size=10),
          legend.text = element_text(size=7),
          plot.title = element_text(hjust=1.0, size=12, face="bold"))
  ggsave(paste("Bootstrap Diffs at",alpha,"quant,",coverage,"depth.png"), 
         plot=p,units="in", width=6, height=4, dpi=600, device = 'png')
}

# Plot generator for cutoff frequency differences,
# each line is a different coverage, plotting pat. freq by cutoff change in freq
plot_bootstrap_cutoffs_alt <- function(cutoffs, alpha=0.025) {
  
  # A dirty bit of code to restrict the coverage to values with impact
  cutoffs = cutoffs[cutoffs$Coverage <= 2.5, ]
  
  # Load ggplot library
  library(ggplot2)
  
  # Plotting the paternal-embryo differences using ggplot2
  p <- ggplot() + 
    geom_line(data=cutoffs,aes(x=Pat_Exp_Freq, y=P_E_Low_Cutoff, group=factor(Coverage), 
                               linetype = 'solid', colour=Coverage)) +
    geom_line(data=cutoffs,aes(x=Pat_Exp_Freq, y=P_E_High_Cutoff, group=factor(Coverage), 
                               linetype = 'dashed', colour=Coverage)) +
    xlab("Paternal Frequency") + 
    ylab("Frequency Difference") + 
    # scale_x_continuous(breaks=seq(0, 3, .5), minor_breaks=seq(0, 3, .1),
    #                    limits=c(-.1,3.1)) +
    scale_y_continuous(breaks=seq(-.25, .25, .05), minor_breaks=seq(-.25, .25, .01),
                       limits=c(-.26,.26)) +
    scale_linetype_manual(values = c(1,2),
                          labels=c("Upper",
                                   "Lower"),
                          name="Bound") +
    scale_colour_continuous(name="Coverage (per allele)", 
                            guide = guide_colourbar(title.position = "top",direction = "horizontal")) + 
    ggtitle(paste("Paternal-Embryo Bootstrap Differences at ",alpha,"quantile extremes")) +
    theme_bw() +
    theme(legend.justification=c(1,0),
          legend.position=c(0.95,0.05),  # Position legend in bottom right
          legend.box="horizontal",  # Position multiple legends horizontally
          legend.background = element_rect(fill="white",size=0.5, linetype="solid", colour ="black"),
          legend.title = element_text(size=10),
          legend.text = element_text(size=7),
          plot.title = element_text(hjust=1.0, size=13, face="bold"))
  ggsave(paste("P-E Bootstrap Diffs at",alpha,"alt.png"), 
         plot=p,units="in", width=6, height=4, dpi=600, device = 'png')
  
  # Plotting the embryo-adult differences using ggplot2
  p <- ggplot() + 
    geom_line(data=cutoffs,aes(x=Pat_Exp_Freq, y=E_A_Low_Cutoff, group=factor(Coverage), 
                               linetype = 'solid', colour=Coverage)) +
    geom_line(data=cutoffs,aes(x=Pat_Exp_Freq, y=E_A_High_Cutoff, group=factor(Coverage), 
                               linetype = 'dashed', colour=Coverage)) +
    xlab("Paternal Frequency") + 
    ylab("Frequency Difference") + 
    # scale_x_continuous(breaks=seq(0, 3, .5), minor_breaks=seq(0, 3, .1),
    #                    limits=c(-.1,3.1)) +
    scale_y_continuous(breaks=seq(-.25, .25, .05), minor_breaks=seq(-.25, .25, .01),
                       limits=c(-.26,.26)) +
    scale_linetype_manual(values = c(1,2),
                          labels=c("Upper",
                                   "Lower"),
                          name="Bound") +
    scale_colour_continuous(name="Coverage (per allele)", 
                            guide = guide_colourbar(title.position = "top",direction = "horizontal")) + 
    ggtitle(paste("Embryo-Adult Bootstrap Differences at ",alpha,"quantile extremes")) +
    theme_bw() +
    theme(legend.justification=c(1,0),
          legend.position=c(0.95,0.05),  # Position legend in bottom right
          legend.box="horizontal",  # Position multiple legends horizontally
          legend.background = element_rect(fill="white",size=0.5, linetype="solid", colour ="black"),
          legend.title = element_text(size=10),
          legend.text = element_text(size=7),
          plot.title = element_text(hjust=1.0, size=13, face="bold"))
  ggsave(paste("E-A Bootstrap Diffs at",alpha,"alt.png"), 
         plot=p,units="in", width=6, height=4, dpi=600, device = 'png')
}


# Running as a script:

# Pick an alpha value (the quantile high or low to asses the significance cutoff at)
# (here we use 2.5% to represent 5% from the mean considering both directions)
# alpha = 0.025

# Save the replicate file
# save(reps,file="bootstrap_replicates.Rdata")

# Load the replicate file
# load("bootstrap_replicates.Rdata")

# Generate the cutoff summary
# cutoffs <- cutoff_summary(reps,alpha)

# Save the cutoff summary
# save(cutoffs,file="bootstrap_cutoffs.Rdata")

# Load the cutoff summary
# load("bootstrap_cutoffs.Rdata")




plot_all_cohorts <- function(exper_data) {
  library(ggplot2)
  
  # Plotting the paternal-embryo differences using ggplot2
  p <- ggplot() + 
    geom_line(data=exper_data,aes(x=Pat_Exp_Freq, y=P_E_Low_Cutoff, group=factor(Coverage), 
                                  linetype = 'solid', colour=Coverage)) +
    geom_line(data=cutoffs,aes(x=Pat_Exp_Freq, y=P_E_High_Cutoff, group=factor(Coverage), 
                               linetype = 'dashed', colour=Coverage)) +
    xlab("Paternal Frequency") + 
    ylab("Frequency Difference") + 
    # scale_x_continuous(breaks=seq(0, 3, .5), minor_breaks=seq(0, 3, .1),
    #                    limits=c(-.1,3.1)) +
    scale_y_continuous(breaks=seq(-.25, .25, .05), minor_breaks=seq(-.25, .25, .01),
                       limits=c(-.26,.26)) +
    scale_linetype_manual(values = c(1,2),
                          labels=c("Upper",
                                   "Lower"),
                          name="Bound") +
    scale_colour_continuous(name="Coverage (per allele)", 
                            guide = guide_colourbar(title.position = "top",direction = "horizontal")) + 
    ggtitle(paste("Paternal-Embryo Bootstrap Differences at ",alpha,"quantile extremes")) +
    theme_bw() +
    theme(legend.justification=c(1,0),
          legend.position=c(0.95,0.05),  # Position legend in bottom right
          legend.box="horizontal",  # Position multiple legends horizontally
          legend.background = element_rect(fill="white",size=0.5, linetype="solid", colour ="black"),
          legend.title = element_text(size=10),
          legend.text = element_text(size=7),
          plot.title = element_text(hjust=1.0, size=13, face="bold"))
  ggsave(paste("P-E Bootstrap Diffs at",alpha,"alt.png"), 
         plot=p,units="in", width=6, height=4, dpi=600, device = 'png')
}

plot_freqs_by_cohort <- function(freq_data){
  exper_lines = unique(freq_data$line)
  exper_chroms = unique(freq_data$chrom)
  for (exp_line in exper_lines) {
    exp_freq_data = subset(freq_data,freq_data$line == exp_line)
    for (exp_chrom in exper_chroms) {
      exp_chrom_freq_data = subset(exp_freq_data,exp_freq_data$chrom == exp_chrom)
      plot_all_cohorts(exp_chrom_freq_data)
    }
  }
}



############

plot_comb_off_sep_chrom <- function(combined_cohort_data,inv_name) {
  library(ggplot2)
  
  # Plotting the paternal-embryo differences using ggplot2
  p <- ggplot() + 
    geom_line(data=combined_cohort_data,aes(x=Cohort, y=Inversion_Frequency, group=1)) +
    # geom_line(data=cutoffs,aes(x=Pat_Exp_Freq, y=P_E_High_Cutoff, group=factor(Coverage), 
    #                            linetype = 'dashed', colour=Coverage)) +
    xlab("Age Cohort") + 
    ylab("Mean Inversion Frequency") + 
    # scale_y_continuous(limits=c(0,0.5)) +
    coord_cartesian(ylim=c(0,0.5)) +
    ggtitle(paste(inv_name," Mean Frequency Across a Generation")) +
    # theme(legend.justification=c(1,0),
    #       legend.position=c(0.95,0.05),  # Position legend in bottom right
    #       legend.box="horizontal",  # Position multiple legends horizontally
    #       legend.background = element_rect(fill="white",size=0.5, linetype="solid", colour ="black"),
    #       legend.title = element_text(size=10),
    #       legend.text = element_text(size=7),
    #       plot.title = element_text(hjust=1.0, size=13, face="bold")) +
    theme_bw()
  ggsave(paste("Pooled Offspring Frequency Plot ",inv_name,".png"), 
         plot=p,units="in", width=7, height=4, dpi=600, device = 'png')
}

plot_comb_off <- function(combined_cohort_data) {
  library(ggplot2)
  
  # Plotting the paternal-embryo and embryo-combined offspring frequency differences using ggplot2
  p <- ggplot() + 
    geom_line(data=combined_cohort_data,aes(x=Cohort, y=Inversion_Frequency, 
                                            group=Inversion, colour=Inversion)) +
    xlab("Age Cohort") + 
    ylab("Mean Inversion Frequency") + 
    # scale_y_continuous(limits=c(0,0.5)) +
    coord_cartesian(ylim=c(0,0.5)) +
    ggtitle(paste("Inversion Frequency Across a Generation")) +
    theme_bw() +
    theme(plot.title = element_text(hjust=0.5, size=20, face="bold"),
          axis.text=element_text(size=12),
          axis.title=element_text(size=15,face="bold"))
  ggsave(paste("Pooled Offspring Frequency Plot.png"), 
         plot=p,units="in", width=7, height=4, dpi=600, device = 'png')
}


combine_freq_sets <- function(freq_data,cohort,inv_name) {
  avg_freq = mean(freq_data$Inv_F)
  sum_inv = sum(freq_data$N_Inv)
  sum_std = sum(freq_data$N_Std)
  combined_cohort_data <- data.frame(inv_name,cohort,avg_freq,sum_inv,sum_std)
  names(combined_cohort_data) <- c('Inversion','Cohort','Inversion_Frequency','Inv_Read_Count','Std_Read_Count')
  return(combined_cohort_data)
}

plot_combined_off_by_chrom <- function(freq_data) {
  exper_chroms = unique(freq_data$Chrom)
  all_chrom_combined <- data.frame()
  for (exp_chrom in exper_chroms) {
    exp_freq_data = subset(freq_data,freq_data$Chrom == exp_chrom)
    inv_name = paste0('In(',exp_chrom,')',exp_freq_data$Inv[1])
    # print(exp_freq_data)
    
    pat_data = subset(exp_freq_data,exp_freq_data$Cohort == 'Parental ♂')
    print(pat_data$Inv_F)
    combined_cohort_data <- combine_freq_sets(pat_data,'Parental',inv_name)
    print(combined_cohort_data$Inversion_Frequency)
    combined_cohort_data$Inversion_Frequency <- combined_cohort_data$Inversion_Frequency/2
    # print(combined_cohort_data)
    
    
    emb_data = subset(exp_freq_data,exp_freq_data$Cohort == 'Embryo')
    combined_cohort_data <- rbind(combined_cohort_data,combine_freq_sets(emb_data,'Embryo',inv_name))
    
    eaf_data = subset(exp_freq_data,exp_freq_data$Cohort == 'Early ♀')
    laf_data = subset(exp_freq_data,exp_freq_data$Cohort == 'Late ♀')
    eam_data = subset(exp_freq_data,exp_freq_data$Cohort == 'Early ♂')
    lam_data = subset(exp_freq_data,exp_freq_data$Cohort == 'Late ♂')
    
    off_data = combine_freq_sets(rbind(eaf_data,laf_data,eam_data,lam_data),'Adult Offspring',inv_name)
    combined_cohort_data <- rbind(combined_cohort_data,off_data)
    
    combined_cohort_data$Cohort <- factor(combined_cohort_data$Cohort, levels = c('Parental','Embryo','Adult Offspring'))
    # print(combined_cohort_data)
    # print(str(combined_cohort_data))
    
    plot_comb_off_sep_chrom(combined_cohort_data,inv_name)
    all_chrom_combined <- rbind(all_chrom_combined,combined_cohort_data)
  }
  print(all_chrom_combined)
  plot_comb_off(all_chrom_combined)
}

scatter_plot_cohorts <- function(freq_data){
  chrom <- unique(freq_data$Chrom)
  inv <- unique(freq_data$Inv)
  colors <- c("192" = "red", "251" = "blue", "254" = "darkgreen", "418" = "orange")
  
  library(ggplot2)
  
  # Plotting the paternal-embryo differences using ggplot2
  p <- ggplot() + 
    geom_point(data=freq_data,aes(x=Cohort, y=Inv_F, colour=Line)) +
    xlab("Age Cohort") + 
    ylab("Inversion Frequency") + 
    # scale_y_continuous(limits=c(0,1)) +
    # scale_linetype_manual(values = c(1,2),
    #                       labels=c("Upper",
    #                                "Lower"),
    #                       name="Bound") +
    scale_colour_manual(name="Maternal Line", values = colors,
                        guide = guide_legend(title.position = "top",direction = "vertical")) +
    ggtitle(paste0("In(",chrom,")",inv," Frequency by Cohort")) +
    # theme(legend.justification=c(1,0),
    #       legend.position=c(0.95,0.05),  # Position legend in bottom right
    #       legend.box="horizontal",  # Position multiple legends horizontally
    #       legend.background = element_rect(fill="white",size=0.5, linetype="solid", colour ="black"),
    #       legend.title = element_text(size=10),
    #       legend.text = element_text(size=7),
    #       plot.title = element_text(hjust=1.0, size=13, face="bold")) +
    theme_bw()
  ggsave(paste0("In(",chrom,")",inv," Frequency by Cohort.png"),
         plot=p, units="in", width=6, height=4, dpi=600, device = 'png')
  
}

halve_paternal <- function(freq_data){
  freq_data_pat <- subset(freq_data,freq_data$Cohort == 'Parental ♂')
  freq_data_pat$Inv_F <- freq_data_pat$Inv_F/2
  freq_data_non_pat <- subset(freq_data,freq_data$Cohort != 'Parental ♂')
  modified_freq_data <- rbind(freq_data_pat,freq_data_non_pat)
  return(modified_freq_data)
}

scatter_plot_by_chrom <- function(freq_data){
  exper_chroms = unique(freq_data$Chrom)
  for (exp_chrom in exper_chroms) {
    exp_freq_data = subset(freq_data,freq_data$Chrom == exp_chrom)
    mod_freq_data = halve_paternal(exp_freq_data)
    scatter_plot_cohorts(mod_freq_data)
  }
}


scatter_plot_diffs <- function(freq_data){
  chrom <- unique(freq_data$Chrom)
  inv <- unique(freq_data$Inv)
  colors <- c("192" = "red", "251" = "blue", "254" = "darkgreen", "418" = "orange")
  
  library(ggplot2)
  
  # Plotting the paternal-embryo differences using ggplot2
  p <- ggplot() + 
    geom_point(data=freq_data,aes(x=Difference_Cohort, y=Inv_Freq_Diff, colour=Line)) +
    geom_hline(yintercept=0) +
    xlab("Age Cohort") + 
    ylab("Inversion Frequency Difference") + 
    ggtitle(paste0("In(",chrom,")",inv," Parental-Offspring Freq. Diffs. by Cohort")) +
    theme_bw() +
    # scale_y_continuous(limits=c(0,1)) +
    # scale_linetype_manual(values = c(1,2),
    #                       labels=c("Upper",
    #                                "Lower"),
    #                       name="Bound") +
    scale_colour_manual(name="Maternal Line", values = colors,
                        guide = guide_legend(title.position = "top",direction = "vertical"))
    # theme(legend.justification=c(1,0),
    #       legend.position=c(0.95,0.05),  # Position legend in bottom right
    #       legend.box="horizontal",  # Position multiple legends horizontally
    #       legend.background = element_rect(fill="white",size=0.5, linetype="solid", colour ="black"),
    #       legend.title = element_text(size=10),
    #       legend.text = element_text(size=7),
    #       plot.title = element_text(hjust=1.0, size=13, face="bold")) +
  ggsave(paste0("In(",chrom,")",inv," Parental Diffs Scatter.png"), 
         plot=p,units="in", width=6, height=4, dpi=600, device = 'png')
  
}

calculate_par_diffs <- function(freq_data,line,chrom) {
  # print(freq_data)
  
  par_data = subset(freq_data,freq_data$Cohort == 'Parental ♂')
  inv <- par_data$Inv
  combined_diff_data <- data.frame()
  
  emb_data = subset(freq_data,freq_data$Cohort == 'Embryo')
  if (!(dim(emb_data)[1] == 0)) {
    par_emb_diff = emb_data$Inv_F - par_data$Inv_F
    this_diff_frame <- data.frame(chrom,inv,line,"Par-Emb",par_emb_diff)
    names(this_diff_frame) <- c('Chrom','Inv','Line','Difference_Cohort','Inv_Freq_Diff')
    combined_diff_data <- rbind(combined_diff_data,this_diff_frame)
  }
  
  eaf_data = subset(freq_data,freq_data$Cohort == 'Early ♀')
  if (!(dim(eaf_data)[1] == 0)) {
    par_eaf_diff = eaf_data$Inv_F - par_data$Inv_F
    this_diff_frame <- data.frame(chrom,inv,line,"Par-EA♀",par_eaf_diff)
    names(this_diff_frame) <- c('Chrom','Inv','Line','Difference_Cohort','Inv_Freq_Diff')
    combined_diff_data <- rbind(combined_diff_data,this_diff_frame)
  }
  
  laf_data = subset(freq_data,freq_data$Cohort == 'Late ♀')
  if (!(dim(laf_data)[1] == 0)) {
    par_laf_diff = laf_data$Inv_F - par_data$Inv_F
    this_diff_frame <- data.frame(chrom,inv,line,"Par-LA♀",par_laf_diff)
    names(this_diff_frame) <- c('Chrom','Inv','Line','Difference_Cohort','Inv_Freq_Diff')
    combined_diff_data <- rbind(combined_diff_data,this_diff_frame)
  }
  
  eam_data = subset(freq_data,freq_data$Cohort == 'Early ♂')
  if (!(dim(eam_data)[1] == 0)) {
    # print(eam_data)
    # print(str(eam_data))
    par_eam_diff = eam_data$Inv_F - par_data$Inv_F
    this_diff_frame <- data.frame(chrom,inv,line,"Par-EA♂",par_eam_diff)
    names(this_diff_frame) <- c('Chrom','Inv','Line','Difference_Cohort','Inv_Freq_Diff')
    combined_diff_data <- rbind(combined_diff_data,this_diff_frame)
  }
  
  lam_data = subset(freq_data,freq_data$Cohort == 'Late ♂')
  if (!(dim(lam_data)[1] == 0)) {
    par_lam_diff = lam_data$Inv_F - par_data$Inv_F
    this_diff_frame <- data.frame(chrom,inv,line,"Par-LA♂",par_lam_diff)
    names(this_diff_frame) <- c('Chrom','Inv','Line','Difference_Cohort','Inv_Freq_Diff')
    combined_diff_data <- rbind(combined_diff_data,this_diff_frame)
  }
  combined_diff_data$Difference_Cohort <- factor(combined_diff_data$Difference_Cohort,levels=c("Par-Emb","Par-EA♀","Par-LA♀","Par-EA♂","Par-LA♂"))
  return(combined_diff_data)
}

scatter_plot_par_diff_by_chrom <- function(freq_data){
  exper_chroms = unique(freq_data$Chrom)
  exper_lines = unique(freq_data$Line)
  mod_freq_data = halve_paternal(freq_data)
  for (exp_chrom in exper_chroms) {
    chrom_data = subset(mod_freq_data,mod_freq_data$Chrom == exp_chrom)
    diff_data <- data.frame()
    for (exp_line in exper_lines) {
      exper_data = subset(chrom_data,chrom_data$Line == exp_line)
      exper_diff_data <- calculate_par_diffs(exper_data,exp_line,exp_chrom)
      diff_data <- rbind(diff_data,exper_diff_data)
    }
    # print(diff_data)
    scatter_plot_diffs(diff_data)
  }
}








scatter_plot_emb_diffs <- function(freq_data){
  chrom <- unique(freq_data$Chrom)
  inv <- unique(freq_data$Inv)
  colors <- c("192" = "red", "251" = "blue", "254" = "darkgreen", "418" = "orange")
  labs <- c("192" = "ZI192N", "251" = "ZI251N", "254" = "ZI254N", "418" = "ZI418N")
  
  library(ggplot2)
  
  # Plotting the paternal-embryo differences using ggplot2
  p <- ggplot() + 
    geom_point(data=freq_data,aes(x=Difference_Cohort, y=Inv_Freq_Diff, colour=Line)) +
    geom_hline(yintercept=0) +
    # xlab("Age Cohort") + 
    theme(axis.title.x=element_blank()) +
    ylab("Inversion Frequency Difference") + 
    ggtitle(paste0("In(",chrom,")",inv," Intra-Generational Frequency ∆ by Cohort")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle=20, hjust=1)) +
    scale_colour_manual(name="Maternal Line", values = colors, labels=labs,
                        guide = guide_legend(title.position = "top",direction = "vertical")) +
    theme(plot.title = element_text(size=15, face="bold"),
          # axis.text=element_text(size=12),
          axis.title=element_text(size=12,face="bold"))
    # theme(legend.justification=c(1,0),
    #       legend.position=c(0.95,0.05),  # Position legend in bottom right
    #       legend.box="horizontal",  # Position multiple legends horizontally
    #       legend.background = element_rect(fill="white",size=0.5, linetype="solid", colour ="black"),
    #       legend.title = element_text(size=10),
    #       legend.text = element_text(size=7),
    #       plot.title = element_text(hjust=1.0, size=13, face="bold")) +
  ggsave(paste0("In(",chrom,")",inv," Embryo Diffs Scatter.png"), 
         plot=p,units="in", width=7, height=4, dpi=600, device = 'png') # Exporting as pdf gives a vector graphic
  
}

calculate_embryo_diffs <- function(freq_data,line,chrom) {
  # print(freq_data)
  
  emb_data = subset(freq_data,freq_data$Cohort == 'Embryo')
  inv <- emb_data$Inv
  combined_diff_data <- data.frame()
  
  par_data = subset(freq_data,freq_data$Cohort == 'Parental ♂')
  if (!(dim(emb_data)[1] == 0)) {
    par_emb_diff = emb_data$Inv_F - par_data$Inv_F
    this_diff_frame <- data.frame(chrom,inv,line,"Embryo-Paternal/2",par_emb_diff)
    names(this_diff_frame) <- c('Chrom','Inv','Line','Difference_Cohort','Inv_Freq_Diff')
    combined_diff_data <- rbind(combined_diff_data,this_diff_frame)
  }
  
  eaf_data = subset(freq_data,freq_data$Cohort == 'Early ♀')
  if (!(dim(eaf_data)[1] == 0)) {
    emb_eaf_diff = eaf_data$Inv_F - emb_data$Inv_F
    this_diff_frame <- data.frame(chrom,inv,line,"Early Eclosing Adult♀-Embryo",emb_eaf_diff)
    names(this_diff_frame) <- c('Chrom','Inv','Line','Difference_Cohort','Inv_Freq_Diff')
    combined_diff_data <- rbind(combined_diff_data,this_diff_frame)
  }
  
  laf_data = subset(freq_data,freq_data$Cohort == 'Late ♀')
  if (!(dim(laf_data)[1] == 0)) {
    emb_laf_diff = laf_data$Inv_F - emb_data$Inv_F
    this_diff_frame <- data.frame(chrom,inv,line,"Late Eclosing Adult♀-Embryo",emb_laf_diff)
    names(this_diff_frame) <- c('Chrom','Inv','Line','Difference_Cohort','Inv_Freq_Diff')
    combined_diff_data <- rbind(combined_diff_data,this_diff_frame)
  }
  
  eam_data = subset(freq_data,freq_data$Cohort == 'Early ♂')
  if (!(dim(eam_data)[1] == 0)) {
    # print(eam_data)
    # print(str(eam_data))
    emb_eam_diff = eam_data$Inv_F - emb_data$Inv_F
    this_diff_frame <- data.frame(chrom,inv,line,"Early Eclosing Adult♂-Embryo",emb_eam_diff)
    names(this_diff_frame) <- c('Chrom','Inv','Line','Difference_Cohort','Inv_Freq_Diff')
    combined_diff_data <- rbind(combined_diff_data,this_diff_frame)
  }
  
  lam_data = subset(freq_data,freq_data$Cohort == 'Late ♂')
  if (!(dim(lam_data)[1] == 0)) {
    emb_lam_diff = lam_data$Inv_F - emb_data$Inv_F
    this_diff_frame <- data.frame(chrom,inv,line,"Late Eclosing Adult♂-Embryo",emb_lam_diff)
    names(this_diff_frame) <- c('Chrom','Inv','Line','Difference_Cohort','Inv_Freq_Diff')
    combined_diff_data <- rbind(combined_diff_data,this_diff_frame)
  }
  combined_diff_data$Difference_Cohort <- factor(combined_diff_data$Difference_Cohort,
                                                 levels=c("Embryo-Paternal/2","Early Eclosing Adult♀-Embryo","Late Eclosing Adult♀-Embryo","Early Eclosing Adult♂-Embryo","Late Eclosing Adult♂-Embryo"))
  return(combined_diff_data)
}

scatter_plot_embryo_diff_by_chrom <- function(freq_data){
  exper_chroms = unique(freq_data$Chrom)
  exper_lines = unique(freq_data$Line)
  mod_freq_data = halve_paternal(freq_data)
  for (exp_chrom in exper_chroms) {
    chrom_data = subset(mod_freq_data,mod_freq_data$Chrom == exp_chrom)
    diff_data <- data.frame()
    for (exp_line in exper_lines) {
      exper_data = subset(chrom_data,chrom_data$Line == exp_line)
      exper_diff_data <- calculate_embryo_diffs(exper_data,exp_line,exp_chrom)
      diff_data <- rbind(diff_data,exper_diff_data)
    }
    # print(diff_data)
    scatter_plot_emb_diffs(diff_data)
  }
}




gen_comb_off_allele_freq_data <- function(freq_data) {
  exper_chroms = unique(freq_data$Chrom)
  exper_lines = unique(freq_data$Line)
  comb_off_allele_freq_data <- data.frame()
  comb_data_names = c('Inversion','Cohort','Line','Inversion_Frequency','Inv_Allele_Count','Std_Allele_Count')
  for (exp_chrom in exper_chroms) {
    chrom_data = subset(freq_data,freq_data$Chrom == exp_chrom)
    inv_name = paste0('In(',exp_chrom,')',chrom_data$Inv[1])
    # print(chrom_data)
    
    for (line in exper_lines) {
      line_data = subset(chrom_data,chrom_data$Line == line)
      # print(line_data)
      line_name = paste0('ZI',line,'N')
      # print(line_name)
      
      pat_data = subset(line_data,line_data$Cohort == 'Parental ♂')
      cohort = 'Parental'
      total_allele_Ns = pat_data$Sample_Ind_N*4 #Here considering the total parental allele N as 2*paternal
      inv_allele_counts = pat_data$Inv_F*pat_data$Sample_Ind_N*2
      std_allele_counts = total_allele_Ns-inv_allele_counts
      inv_allele_freqs = inv_allele_counts/total_allele_Ns
      num_entries = length(inv_allele_freqs)
      # print(pat_data)
      allele_cohort_data <- data.frame(rep(inv_name,num_entries),rep(cohort,num_entries),rep(line_name,num_entries),inv_allele_freqs,inv_allele_counts,std_allele_counts)
      names(allele_cohort_data) <- comb_data_names
      comb_off_allele_freq_data <- rbind(comb_off_allele_freq_data,allele_cohort_data)
      # names(comb_off_allele_freq_data) <- comb_data_names
      # print(comb_off_allele_freq_data)
      
      
      emb_data = subset(line_data,line_data$Cohort == 'Embryo')
      cohort = 'Embryo'
      total_allele_Ns = emb_data$Sample_Ind_N*2
      inv_allele_counts = emb_data$Inv_F*emb_data$Sample_Ind_N*2
      std_allele_counts = total_allele_Ns-inv_allele_counts
      inv_allele_freqs = inv_allele_counts/total_allele_Ns
      num_entries = length(inv_allele_freqs)
      allele_cohort_data <- data.frame(rep(inv_name,num_entries),rep(cohort,num_entries),rep(line_name,num_entries),inv_allele_freqs,inv_allele_counts,std_allele_counts)
      names(allele_cohort_data) <- comb_data_names
      # print(names(allele_cohort_data))
      # print(names(comb_off_allele_freq_data))
      comb_off_allele_freq_data <- rbind(comb_off_allele_freq_data,allele_cohort_data)
      # print(comb_off_allele_freq_data)
      
      eaf_data = subset(line_data,line_data$Cohort == 'Early ♀')
      laf_data = subset(line_data,line_data$Cohort == 'Late ♀')
      eam_data = subset(line_data,line_data$Cohort == 'Early ♂')
      lam_data = subset(line_data,line_data$Cohort == 'Late ♂')
      off_data = rbind(eaf_data,laf_data,eam_data,lam_data)
      # print(off_data)
      cohort = 'Adult Offspring'
      total_allele_N = sum(off_data$Sample_Ind_N*2)
      inv_allele_count = sum(off_data$Inv_F*off_data$Sample_Ind_N*2)
      std_allele_count = total_allele_N-inv_allele_count
      inv_allele_freq = inv_allele_count/total_allele_N
      num_entries = length(inv_allele_freq)
      combined_line_data <- data.frame(rep(inv_name,num_entries),rep(cohort,num_entries),rep(line_name,num_entries),inv_allele_freq,inv_allele_count,std_allele_count)
      names(combined_line_data) <- comb_data_names
      # print(combined_line_data)
      # print(names(combined_line_data))
      # print(names(comb_off_allele_freq_data))
      comb_off_allele_freq_data <- rbind(comb_off_allele_freq_data,combined_line_data)
    }
  }
  
  # names(comb_off_allele_freq_data) <- comb_data_names
  comb_off_allele_freq_data$Cohort <- factor(comb_off_allele_freq_data$Cohort, levels = c('Parental','Embryo','Adult Offspring'))
  # print(comb_off_allele_freq_data)
  # print(str(comb_off_allele_freq_data))
  
  return(comb_off_allele_freq_data)
}


plot_comb_off_panel_chrom_grp_line <- function(comb_data) {
  library(ggplot2)
  
  # Plotting the paternal-embryo differences using ggplot2
  # colors <- c("192" = "red", "251" = "blue", "254" = "darkgreen", "418" = "orange")
  p <- ggplot() + 
    geom_line(data=comb_data,aes(x=Cohort, y=Inversion_Frequency, group=Line, color=Line)) +
    geom_point(data=comb_data,aes(x=Cohort, y=Inversion_Frequency, group=Line, color=Line)) +
    xlab("Demographic Cohort") + 
    ylab("Inversion Frequency") + 
    # scale_y_continuous(limits=c(0,0.5)) + 
    coord_cartesian(ylim=c(0,0.4)) + 
    ggtitle(paste("Inversion Frequency Changes Within a Generation")) +
    # scale_colour_manual(name="Maternal Line", values = colors,
    #                     guide = guide_legend(title.position = "top",direction = "vertical")) +
    theme_bw() +
    facet_wrap(~ Inversion, ncol=4) +
    theme(axis.text.x = element_text(angle=20, hjust=1)) +
    theme(strip.text.x=element_text(face = "bold.italic"))
  # ggsave(paste("Pooled Offspring Frequency Plot by Line.png"), 
  #        plot=p,units="in", width=7, height=3.5, dpi=600, device = 'png')
  # ggsave(paste("Pooled Offspring Frequency Plot panelled by Inv separated by Line.png"), 
  #        plot=p,units="in", width=7, height=3.5, dpi=600, device = 'png')
  return(p)
  # ggsave(paste(fname),plot=p,units="in", width=7, height=3.5, dpi=600, device = 'png')
}

# plot_comb_off_by_chrom_sep_line <- function(freq_data) {
#   comb_off_allele_freq_data <- gen_comb_off_allele_freq_data(freq_data)
#   print(comb_off_allele_freq_data)
#   plot_sep_line_comb_off(comb_off_allele_freq_data)
# }


plot_comb_off_trends_panel_chrom_grp_line <- function(comb_data) {
  library(ggplot2)
  
  # Plotting the paternal-embryo differences using ggplot2
  # colors <- c("192" = "red", "251" = "blue", "254" = "darkgreen", "418" = "orange")
  p <- ggplot() + 
    geom_line(data=comb_data,aes(x=Cohort, y=Inversion_Frequency, group=interaction(Line,Trend), color=Line, linetype=Trend)) + 
    geom_point(data=comb_data,aes(x=Cohort, y=Inversion_Frequency, group=interaction(Line,Trend), color=Line)) + 
    xlab("Demographic Cohort") + 
    ylab("Inversion Frequency") + 
    # scale_y_continuous(limits=c(0,0.5)) + 
    coord_cartesian(ylim=c(0,0.4)) + 
    ggtitle(paste("Inversion Frequency Changes Within a Generation")) +
    # scale_colour_manual(name="Maternal Line", values = colors,
    #                     guide = guide_legend(title.position = "top",direction = "vertical")) +
    theme_bw() +
    facet_wrap(~ Inversion, ncol=4) +
    theme(axis.text.x = element_text(angle=20, hjust=1)) +
    theme(strip.text.x = element_text(face = "bold.italic"))
  # ggsave(paste("Pooled Offspring Frequency Plot by Line.png"), 
  #        plot=p,units="in", width=7, height=3.5, dpi=600, device = 'png')
  # ggsave(paste("Pooled Offspring Frequency Plot panelled by Inv separated by Line.png"), 
  #        plot=p,units="in", width=7, height=3.5, dpi=600, device = 'png')
  return(p)
}




plot_comb_off_panel_line_grp_chrom <- function(comb_data,fname="Pooled Offspring Frequency Plot panelled by Line separated by Inv.png") {
  # chrom <- unique(freq_data$Chrom)
  # inv <- unique(freq_data$Inv)
  # colors <- c("192" = "red", "251" = "blue", "254" = "darkgreen", "418" = "orange")
  # line.labs <- c("192" = "ZI192N", "251" = "ZI251N", "254" = "ZI254N", "418" = "ZI418N")
  library(ggplot2)
  # print(unique(comb_data$Line))
  
  # Plotting the paternal-embryo differences using ggplot2
  p <- ggplot() + 
    geom_line(data=comb_data,aes(x=Cohort, y=Inversion_Frequency, group=Inversion, color=Inversion)) +
    geom_point(data=comb_data,aes(x=Cohort, y=Inversion_Frequency, group=Inversion, color=Inversion)) +
    xlab("Demographic Cohort") + 
    ylab("Inversion Frequency") + 
    # scale_y_continuous(limits=c(0,0.5)) + 
    coord_cartesian(ylim=c(0,0.4)) +
    ggtitle(paste("Inversion Frequency Changes Within a Generation")) +
    # scale_colour_manual(name="Maternal Line", values = colors,
    #                     guide = guide_legend(title.position = "top",direction = "vertical")) +
    theme_bw() +
    facet_wrap(~ Line, ncol=4) +
    # facet_wrap(~ Line, ncol=4, labeller = labeller(Line = line.labs)) +
    theme(axis.text.x = element_text(angle=20, hjust=1)) +
    theme(strip.text.x=element_text(face = "bold.italic"))
  # ggsave(paste("Pooled Offspring Frequency Plot panelled by Line separated by Inv.png"), 
  #        plot=p,units="in", width=7, height=3.5, dpi=600, device = 'png')
  ggsave(fname,plot=p,units="in", width=7, height=3.5, dpi=600, device = 'png')
}






gen_comb_sex_allele_freq_data <- function(freq_data) {
  exper_chroms = unique(freq_data$Chrom)
  exper_lines = unique(freq_data$Line)
  comb_sex_allele_freq_data <- data.frame()
  comb_data_names = c('Inversion','Cohort','Line','Inversion_Frequency','Inv_Allele_Count','Std_Allele_Count')
  for (exp_chrom in exper_chroms) {
    chrom_data = subset(freq_data,freq_data$Chrom == exp_chrom)
    inv_name = paste0('In(',exp_chrom,')',chrom_data$Inv[1])
    # print(chrom_data)
    
    for (line in exper_lines) {
      line_data = subset(chrom_data,chrom_data$Line == line)
      line_name = paste0('ZI',line,'N')
      # print(line_name)
      
      # pat_data = subset(line_data,line_data$Cohort == 'Parental ♂')
      # cohort = 'Parental'
      # total_allele_Ns = pat_data$Sample_Ind_N*4 #Here considering the total parental allele N as 2*paternal
      # inv_allele_counts = pat_data$Inv_F*pat_data$Sample_Ind_N*2
      # std_allele_counts = total_allele_Ns-inv_allele_counts
      # inv_allele_freqs = inv_allele_counts/total_allele_Ns
      # num_entries = length(inv_allele_freqs)
      # # print(pat_data)
      # allele_cohort_data <- data.frame(rep(inv_name,num_entries),rep(cohort,num_entries),rep(line_name,num_entries),inv_allele_freqs,inv_allele_counts,std_allele_counts)
      # names(allele_cohort_data) <- comb_data_names
      # comb_sex_allele_freq_data <- rbind(comb_sex_allele_freq_data,allele_cohort_data)
      # # names(comb_sex_allele_freq_data) <- comb_data_names
      # # print(comb_sex_allele_freq_data)
      
      
      emb_data = subset(line_data,line_data$Cohort == 'Embryo')
      cohort = 'Embryo'
      total_allele_Ns = emb_data$Sample_Ind_N*2
      inv_allele_counts = emb_data$Inv_F*emb_data$Sample_Ind_N*2
      std_allele_counts = total_allele_Ns-inv_allele_counts
      inv_allele_freqs = inv_allele_counts/total_allele_Ns
      num_entries = length(inv_allele_freqs)
      allele_cohort_data <- data.frame(rep(inv_name,num_entries),rep(cohort,num_entries),rep(line_name,num_entries),inv_allele_freqs,inv_allele_counts,std_allele_counts)
      names(allele_cohort_data) <- comb_data_names
      # print(names(allele_cohort_data))
      # print(names(comb_sex_allele_freq_data))
      comb_sex_allele_freq_data <- rbind(comb_sex_allele_freq_data,allele_cohort_data)
      # print(comb_sex_allele_freq_data)
      
      eaf_data = subset(line_data,line_data$Cohort == 'Early ♀')
      laf_data = subset(line_data,line_data$Cohort == 'Late ♀')
      female_data = rbind(eaf_data,laf_data)
      # print(female_data)
      cohort = 'Female'
      total_allele_N = sum(female_data$Sample_Ind_N*2)
      inv_allele_count = sum(female_data$Inv_F*female_data$Sample_Ind_N*2)
      std_allele_count = total_allele_N-inv_allele_count
      inv_allele_freq = inv_allele_count/total_allele_N
      num_entries = length(inv_allele_freq)
      combined_line_data <- data.frame(rep(inv_name,num_entries),rep(cohort,num_entries),rep(line_name,num_entries),inv_allele_freq,inv_allele_count,std_allele_count)
      names(combined_line_data) <- comb_data_names
      # print(combined_line_data)
      comb_sex_allele_freq_data <- rbind(comb_sex_allele_freq_data,combined_line_data)
      
      eam_data = subset(line_data,line_data$Cohort == 'Early ♂')
      lam_data = subset(line_data,line_data$Cohort == 'Late ♂')
      male_data = rbind(eam_data,lam_data)
      # print(male_data)
      cohort = 'Male'
      total_allele_N = sum(male_data$Sample_Ind_N*2)
      inv_allele_count = sum(male_data$Inv_F*male_data$Sample_Ind_N*2)
      std_allele_count = total_allele_N-inv_allele_count
      inv_allele_freq = inv_allele_count/total_allele_N
      num_entries = length(inv_allele_freq)
      combined_line_data <- data.frame(rep(inv_name,num_entries),rep(cohort,num_entries),rep(line_name,num_entries),inv_allele_freq,inv_allele_count,std_allele_count)
      names(combined_line_data) <- comb_data_names
      # print(combined_line_data)
      comb_sex_allele_freq_data <- rbind(comb_sex_allele_freq_data,combined_line_data)
    }
  }
  
  # names(comb_sex_allele_freq_data) <- comb_data_names
  # comb_sex_allele_freq_data$Cohort <- factor(comb_sex_allele_freq_data$Cohort, levels = c('Male Offspring','Embryo','Female Offspring'))
  comb_sex_allele_freq_data$Cohort <- factor(comb_sex_allele_freq_data$Cohort, levels = c('Embryo','Female','Male'))
  # print(comb_sex_allele_freq_data)
  # print(str(comb_sex_allele_freq_data))
  
  return(comb_sex_allele_freq_data)
}



gen_comb_ecl_allele_freq_data <- function(freq_data) {
  exper_chroms = unique(freq_data$Chrom)
  exper_lines = unique(freq_data$Line)
  comb_ecl_allele_freq_data <- data.frame()
  comb_data_names = c('Inversion','Cohort','Line','Inversion_Frequency','Inv_Allele_Count','Std_Allele_Count')
  for (exp_chrom in exper_chroms) {
    chrom_data = subset(freq_data,freq_data$Chrom == exp_chrom)
    inv_name = paste0('In(',exp_chrom,')',chrom_data$Inv[1])
    # print(chrom_data)
    
    for (line in exper_lines) {
      line_data = subset(chrom_data,chrom_data$Line == line)
      line_name = paste0('ZI',line,'N')
      # print(line_name)
      
      # pat_data = subset(line_data,line_data$Cohort == 'Parental ♂')
      # cohort = 'Parental'
      # total_allele_Ns = pat_data$Sample_Ind_N*4 #Here considering the total parental allele N as 2*paternal
      # inv_allele_counts = pat_data$Inv_F*pat_data$Sample_Ind_N*2
      # std_allele_counts = total_allele_Ns-inv_allele_counts
      # inv_allele_freqs = inv_allele_counts/total_allele_Ns
      # num_entries = length(inv_allele_freqs)
      # # print(pat_data)
      # allele_cohort_data <- data.frame(rep(inv_name,num_entries),rep(cohort,num_entries),rep(line_name,num_entries),inv_allele_freqs,inv_allele_counts,std_allele_counts)
      # names(allele_cohort_data) <- comb_data_names
      # comb_ecl_allele_freq_data <- rbind(comb_ecl_allele_freq_data,allele_cohort_data)
      # # names(comb_ecl_allele_freq_data) <- comb_data_names
      # # print(comb_ecl_allele_freq_data)
      
      
      emb_data = subset(line_data,line_data$Cohort == 'Embryo')
      cohort = 'Embryo'
      total_allele_Ns = emb_data$Sample_Ind_N*2
      inv_allele_counts = emb_data$Inv_F*emb_data$Sample_Ind_N*2
      std_allele_counts = total_allele_Ns-inv_allele_counts
      inv_allele_freqs = inv_allele_counts/total_allele_Ns
      num_entries = length(inv_allele_freqs)
      allele_cohort_data <- data.frame(rep(inv_name,num_entries),rep(cohort,num_entries),rep(line_name,num_entries),inv_allele_freqs,inv_allele_counts,std_allele_counts)
      names(allele_cohort_data) <- comb_data_names
      # print(names(allele_cohort_data))
      # print(names(comb_ecl_allele_freq_data))
      comb_ecl_allele_freq_data <- rbind(comb_ecl_allele_freq_data,allele_cohort_data)
      # print(comb_ecl_allele_freq_data)
      
      eam_data = subset(line_data,line_data$Cohort == 'Early ♂')
      eaf_data = subset(line_data,line_data$Cohort == 'Early ♀')
      early_data = rbind(eaf_data,eam_data)
      # print(early_data)
      cohort = 'Early'
      # cohort = 'Early Eclosing Offspring'
      total_allele_N = sum(early_data$Sample_Ind_N*2)
      inv_allele_count = sum(early_data$Inv_F*early_data$Sample_Ind_N*2)
      std_allele_count = total_allele_N-inv_allele_count
      inv_allele_freq = inv_allele_count/total_allele_N
      num_entries = length(inv_allele_freq)
      combined_line_data <- data.frame(rep(inv_name,num_entries),rep(cohort,num_entries),rep(line_name,num_entries),inv_allele_freq,inv_allele_count,std_allele_count)
      names(combined_line_data) <- comb_data_names
      # print(combined_line_data)
      comb_ecl_allele_freq_data <- rbind(comb_ecl_allele_freq_data,combined_line_data)
      
      lam_data = subset(line_data,line_data$Cohort == 'Late ♂')
      laf_data = subset(line_data,line_data$Cohort == 'Late ♀')
      late_data = rbind(laf_data,lam_data)
      # print(late_data)
      cohort = 'Late'
      # cohort = 'Late Eclosing Offspring'
      total_allele_N = sum(late_data$Sample_Ind_N*2)
      inv_allele_count = sum(late_data$Inv_F*late_data$Sample_Ind_N*2)
      std_allele_count = total_allele_N-inv_allele_count
      inv_allele_freq = inv_allele_count/total_allele_N
      num_entries = length(inv_allele_freq)
      combined_line_data <- data.frame(rep(inv_name,num_entries),rep(cohort,num_entries),rep(line_name,num_entries),inv_allele_freq,inv_allele_count,std_allele_count)
      names(combined_line_data) <- comb_data_names
      # print(combined_line_data)
      comb_ecl_allele_freq_data <- rbind(comb_ecl_allele_freq_data,combined_line_data)
    }
  }
  
  # names(comb_ecl_allele_freq_data) <- comb_data_names
  comb_ecl_allele_freq_data$Cohort <- factor(comb_ecl_allele_freq_data$Cohort, levels = c('Embryo','Early','Late'))
  # comb_ecl_allele_freq_data$Cohort <- factor(comb_ecl_allele_freq_data$Cohort, levels = c('Early Eclosing Offspring','Embryo','Late Eclosing Offspring'))
  # print(comb_ecl_allele_freq_data)
  # print(str(comb_ecl_allele_freq_data))
  
  return(comb_ecl_allele_freq_data)
}



gen_all_comb_allele_freq_data <- function(freq_data) {
  exper_chroms = unique(freq_data$Chrom)
  exper_lines = unique(freq_data$Line)
  comb_allele_freq_data <- data.frame()
  comb_data_names = c('Inversion','Cohort','Line','Inversion_Frequency','Inv_Allele_Count','Std_Allele_Count')
  for (exp_chrom in exper_chroms) {
    chrom_data = subset(freq_data,freq_data$Chrom == exp_chrom)
    inv_name = paste0('In(',exp_chrom,')',chrom_data$Inv[1])
    for (line in exper_lines) {
      line_data = subset(chrom_data,chrom_data$Line == line)
      line_name = paste0('ZI',line,'N')
      
      pat_data = subset(line_data,line_data$Cohort == 'Parental ♂')
      cohort = 'Parental'
      total_allele_Ns = pat_data$Sample_Ind_N*4 #Here considering the total parental allele N as 2*paternal
      inv_allele_counts = pat_data$Inv_F*pat_data$Sample_Ind_N*2
      std_allele_counts = total_allele_Ns-inv_allele_counts
      inv_allele_freqs = inv_allele_counts/total_allele_Ns
      num_entries = length(inv_allele_freqs)
      allele_cohort_data <- data.frame(rep(inv_name,num_entries),rep(cohort,num_entries),rep(line_name,num_entries),inv_allele_freqs,inv_allele_counts,std_allele_counts)
      names(allele_cohort_data) <- comb_data_names
      comb_allele_freq_data <- rbind(comb_allele_freq_data,allele_cohort_data)
      
      
      emb_data = subset(line_data,line_data$Cohort == 'Embryo')
      cohort = 'Embryo'
      total_allele_Ns = emb_data$Sample_Ind_N*2
      inv_allele_counts = emb_data$Inv_F*emb_data$Sample_Ind_N*2
      std_allele_counts = total_allele_Ns-inv_allele_counts
      inv_allele_freqs = inv_allele_counts/total_allele_Ns
      num_entries = length(inv_allele_freqs)
      allele_cohort_data <- data.frame(rep(inv_name,num_entries),rep(cohort,num_entries),rep(line_name,num_entries),inv_allele_freqs,inv_allele_counts,std_allele_counts)
      names(allele_cohort_data) <- comb_data_names
      comb_allele_freq_data <- rbind(comb_allele_freq_data,allele_cohort_data)
      
      eaf_data = subset(line_data,line_data$Cohort == 'Early ♀')
      laf_data = subset(line_data,line_data$Cohort == 'Late ♀')
      eam_data = subset(line_data,line_data$Cohort == 'Early ♂')
      lam_data = subset(line_data,line_data$Cohort == 'Late ♂')
      off_data = rbind(eaf_data,laf_data,eam_data,lam_data)
      cohort = 'Adult Offspring'
      total_allele_N = sum(off_data$Sample_Ind_N*2)
      inv_allele_count = sum(off_data$Inv_F*off_data$Sample_Ind_N*2)
      std_allele_count = total_allele_N-inv_allele_count
      inv_allele_freq = inv_allele_count/total_allele_N
      num_entries = length(inv_allele_freq)
      combined_line_data <- data.frame(rep(inv_name,num_entries),rep(cohort,num_entries),rep(line_name,num_entries),inv_allele_freq,inv_allele_count,std_allele_count)
      names(combined_line_data) <- comb_data_names
      comb_allele_freq_data <- rbind(comb_allele_freq_data,combined_line_data)
      
      
      early_data = rbind(eaf_data,eam_data)
      cohort = 'Early'
      # cohort = 'Early Eclosing Offspring'
      total_allele_N = sum(early_data$Sample_Ind_N*2)
      inv_allele_count = sum(early_data$Inv_F*early_data$Sample_Ind_N*2)
      std_allele_count = total_allele_N-inv_allele_count
      inv_allele_freq = inv_allele_count/total_allele_N
      num_entries = length(inv_allele_freq)
      combined_line_data <- data.frame(rep(inv_name,num_entries),rep(cohort,num_entries),rep(line_name,num_entries),inv_allele_freq,inv_allele_count,std_allele_count)
      names(combined_line_data) <- comb_data_names
      comb_allele_freq_data <- rbind(comb_allele_freq_data,combined_line_data)
      
      late_data = rbind(laf_data,lam_data)
      cohort = 'Late'
      # cohort = 'Late Eclosing Offspring'
      total_allele_N = sum(late_data$Sample_Ind_N*2)
      inv_allele_count = sum(late_data$Inv_F*late_data$Sample_Ind_N*2)
      std_allele_count = total_allele_N-inv_allele_count
      inv_allele_freq = inv_allele_count/total_allele_N
      num_entries = length(inv_allele_freq)
      combined_line_data <- data.frame(rep(inv_name,num_entries),rep(cohort,num_entries),rep(line_name,num_entries),inv_allele_freq,inv_allele_count,std_allele_count)
      names(combined_line_data) <- comb_data_names
      comb_allele_freq_data <- rbind(comb_allele_freq_data,combined_line_data)
      
      female_data = rbind(eaf_data,laf_data)
      cohort = 'Female'
      total_allele_N = sum(female_data$Sample_Ind_N*2)
      inv_allele_count = sum(female_data$Inv_F*female_data$Sample_Ind_N*2)
      std_allele_count = total_allele_N-inv_allele_count
      inv_allele_freq = inv_allele_count/total_allele_N
      num_entries = length(inv_allele_freq)
      combined_line_data <- data.frame(rep(inv_name,num_entries),rep(cohort,num_entries),rep(line_name,num_entries),inv_allele_freq,inv_allele_count,std_allele_count)
      names(combined_line_data) <- comb_data_names
      comb_allele_freq_data <- rbind(comb_allele_freq_data,combined_line_data)
      
      male_data = rbind(eam_data,lam_data)
      cohort = 'Male'
      total_allele_N = sum(male_data$Sample_Ind_N*2)
      inv_allele_count = sum(male_data$Inv_F*male_data$Sample_Ind_N*2)
      std_allele_count = total_allele_N-inv_allele_count
      inv_allele_freq = inv_allele_count/total_allele_N
      num_entries = length(inv_allele_freq)
      combined_line_data <- data.frame(rep(inv_name,num_entries),rep(cohort,num_entries),rep(line_name,num_entries),inv_allele_freq,inv_allele_count,std_allele_count)
      names(combined_line_data) <- comb_data_names
      comb_allele_freq_data <- rbind(comb_allele_freq_data,combined_line_data)
    }
  }
  comb_allele_freq_data$Cohort <- factor(comb_allele_freq_data$Cohort, levels = c('Parental','Embryo','Adult Offspring','Early','Late','Female','Male'))
  return(comb_allele_freq_data)
}

# For separating out the two halves of a trend
#  in order to change the order of the trend-lines or cohorts more easily
#  and label specific comparisons as individually significant
separate_two_part_trend <- function(comb_freq_data,shared_cohort){
  cohorts <- unique(comb_freq_data$Cohort)
  non_shared_cohorts <- cohorts[!cohorts %in% shared_cohort]
  if (length(non_shared_cohorts)!=2){
    stop("Trend separation in 'separate_two_part_trend' has only been written to deal with three cohorts in two-part trends")
  }
  shared_cohort_data <- rbind(mutate(subset(comb_freq_data,Cohort==shared_cohort),Trend='A',.before=4),
                              mutate(subset(comb_freq_data,Cohort==shared_cohort),Trend='B',.before=4))
  cohort_A_data <- mutate(subset(comb_freq_data,Cohort==non_shared_cohorts[1]),Trend='A',.before=4)
  cohort_B_data <- mutate(subset(comb_freq_data,Cohort==non_shared_cohorts[2]),Trend='B',.before=4)
  separated_trend_data <- rbind(cohort_A_data,shared_cohort_data,cohort_B_data)
  rownames(separated_trend_data)<-NULL
  return(separated_trend_data)
}


# For separating out an arbitrary number of two-cohort comparisons from a longer trend
#  in order to change the order of the trend-lines or cohorts more easily
#  and label specific comparisons as individually significant
separate_trend <- function(comb_freq_data,first_cohorts,second_cohorts){
  if (length(first_cohorts)!=length(second_cohorts)){
    stop("Both cohort lists must be the same length to define the list of trends between")
  }
  separated_trend_data <- tibble()
  for (i in seq(length(first_cohorts))) {
    trend_name <- paste0(first_cohorts[i],'-',second_cohorts[i])
    spec_trend_cohort_1 <- mutate(subset(comb_freq_data,Cohort==first_cohorts[i]),Trend=LETTERS[i],.before=4)
    spec_trend_cohort_2 <- mutate(subset(comb_freq_data,Cohort==second_cohorts[i]),Trend=LETTERS[i],.before=4)
    separated_trend_data <- rbind(separated_trend_data,spec_trend_cohort_1,spec_trend_cohort_2)
  }
  rownames(separated_trend_data)<-NULL
  return(separated_trend_data)
}

# To 
add_significance <- function(sep_trend_data){
  sig_trend_data <- cbind(sep_trend_data,Significance=rep(c('Yes','No'),length.out=nrow(sep_trend_data)))
  return(sig_trend_data)
}

load_freq_data <- function(freq_data_file){
  freq_data <- read.csv(file=freq_data_file, header=TRUE, sep=",")
  freq_data$Cohort <- factor(freq_data$Cohort, levels = c('Parental ♂','Embryo','Early ♀','Late ♀','Early ♂','Late ♂'))
  freq_data$Line <- factor(freq_data$Line)
  return(freq_data)
}

load_comb_pval_data <- function(comb_pval_data_file){
  comb_pval_data <- read.csv(file=comb_pval_data_file, header=TRUE, sep=",")
  # comb_pval_data$Cohort <- factor(comb_pval_data$Cohort, levels = c('Parental ♂','Embryo','Early ♀','Late ♀','Early ♂','Late ♂'))
  return(comb_pval_data)
}

# frame_manip_test(freq_data)


# Read inversion frequency CSV data
freq_data <- load_freq_data("../Alignment/inv_freqs/freq_data_combined.csv")

# Read combined p-value CSV data
spec_comb_file = "sep_pval_21-11-12/allele_chi2_combined_all_individual_combined_pval.csv"
comb_pval_data <- load_comb_pval_data(paste0("../Significance/",spec_comb_file))

# Generate plots with all adult offspring combined
comb_off_allele_freq_data <- gen_comb_off_allele_freq_data(freq_data)
# Generate a quartered plot separated into panels by inversion, with combined adult cohorts
plot_comb_off_panel_chrom_grp_line(comb_off_allele_freq_data)
# Generate a quartered plot separated into panels by line, with combined adult cohorts
plot_comb_off_panel_line_grp_chrom(comb_off_allele_freq_data)

trend_sep_comb_off <- separate_trend(comb_off_allele_freq_data,
                                     c('Parental','Embryo'),
                                     c('Embryo','Adult Offspring'))

# Generate plots with eclosion offspring combined
comb_ecl_allele_freq_data <- gen_comb_ecl_allele_freq_data(freq_data)
# Generate a quartered plot separated into panels by inversion, with combined eclosion cohorts
plot_comb_off_panel_chrom_grp_line(comb_ecl_allele_freq_data,fname="Pooled Eclosion Frequency Plot panelled by Inv separated by Line.png")
# Generate a quartered plot separated into panels by line, with combined eclosion cohorts
plot_comb_off_panel_line_grp_chrom(comb_ecl_allele_freq_data,fname="Pooled Eclosion Frequency Plot panelled by Line separated by Inv.png")

trend_sep_comb_ecl <- separate_trend(comb_ecl_allele_freq_data,
                                     c('Embryo','Embryo'),
                                     c('Early','Late'))
ecl_trend_p <- plot_comb_off_trends_panel_chrom_grp_line(trend_sep_comb_ecl)
ecl_trend_p <- ecl_trend_p +
  theme(legend.position="none") +
  scale_linetype_manual(values=c(1,1))
ggsave("Pooled Offspring Frequency by Eclosion Trend Plot panelled by Inv separated by Line.png",
       plot=ecl_trend_p,units="in", width=6, height=3.2, dpi=600, device = 'png')



# Generate plots with m/f offspring combined
comb_sex_allele_freq_data <- gen_comb_sex_allele_freq_data(freq_data)
# Generate a quartered plot separated into panels by inversion, with combined sex cohorts
plot_comb_off_panel_chrom_grp_line(comb_sex_allele_freq_data,fname="Pooled Sex Frequency Plot panelled by Inv separated by Line.png")
# Generate a quartered plot separated into panels by line, with combined sex cohorts
plot_comb_off_panel_line_grp_chrom(comb_sex_allele_freq_data,fname="Pooled Sex Frequency Plot panelled by Line separated by Inv.png")

trend_sep_comb_sex <- separate_trend(comb_sex_allele_freq_data,
                                     c('Embryo','Embryo'),
                                     c('Female','Male'))
sex_trend_p <- plot_comb_off_trends_panel_chrom_grp_line(trend_sep_comb_sex)
sex_trend_p <- sex_trend_p +
  theme(legend.position="none") +
  scale_linetype_manual(values=c(1,1))
ggsave("Pooled Offspring Frequency by Sex Trend Plot panelled by Inv separated by Line.png",
       plot=sex_trend_p,units="in", width=6, height=3.2, dpi=600, device = 'png')

sing_trend_comb_sex <- separate_trend(comb_sex_allele_freq_data,
                                     c('Female'),
                                     c('Male'))
sex_sing_trend_p <- plot_comb_off_trends_panel_chrom_grp_line(sing_trend_comb_sex)
sex_sing_trend_p <- sex_sing_trend_p +
  theme(legend.position="none") +
  scale_linetype_manual(values=c(1,1))
ggsave("Pooled Sex Comparison Plot panelled by Inv separated by Line.png",
       plot=sex_sing_trend_p,units="in", width=6, height=3.2, dpi=600, device = 'png')




all_comb_allele_freq_data <- gen_all_comb_allele_freq_data(freq_data)
sex_off_freq_data <- filter(all_comb_allele_freq_data,
                            Cohort%in%c('Female','Adult Offspring','Male'))
sex_off_freq_data$Cohort <- factor(sex_off_freq_data$Cohort,levels=c('Female','Adult Offspring','Male'))
p<-plot_comb_off_panel_chrom_grp_line(sex_off_freq_data)
p<-p+scale_x_discrete(labels = c('Female','All Adult Offspring','Male'))+
  aes(shape=Cohort)+
  scale_shape_manual(values=c(19, 1, 19)) +
  theme(legend.position="none")
ggsave("Test Pooled Sex Frequency Plot panelled by Inv separated by Line.png",plot=p,units="in", width=6.3, height=3.5, dpi=600, device = 'png')


subset_sex_off_freq_data <- filter(all_comb_allele_freq_data,
                            Cohort%in%c('Female','Adult Offspring','Male'),
                            Inversion%in%c('In(3R)K'))
subset_sex_off_freq_data$Cohort <- factor(subset_sex_off_freq_data$Cohort,levels=c('Female','Adult Offspring','Male'))
p<-plot_comb_off_panel_chrom_grp_line(subset_sex_off_freq_data)+
  scale_x_discrete(labels = c('F','A','M'),
                      expand=expansion(mult=.1))+
  aes(shape=Cohort)+
  scale_shape_manual(values=c(19, 1, 19)) +
  labs(title=NULL,x=NULL,y=NULL) +
  theme(legend.position="none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
        # plot.margin = unit(c(1,5,1,1),'cm'))
# ggsave("3MT Adult Sex Comp.png",plot=p,units="in", width=2.5, height=2.5, dpi=600, device = 'png')
ggsave("3MT Adult Sex Comp.png",plot=p,units="in", width=3, height=2, dpi=600, device = 'png')

subset_off_freq_data <- filter(all_comb_allele_freq_data,
                               Cohort%in%c('Parental','Embryo','Adult Offspring'),
                               Inversion%in%c('In(3R)K'))
                               # Inversion%in%c('In(2R)NS','In(3R)K'))
p<-plot_comb_off_panel_chrom_grp_line(subset_off_freq_data) +
  scale_x_discrete(labels = c('P','E','A'),
                   expand=expansion(mult=.1))+
  labs(title=NULL,x=NULL,y=NULL) +
  theme(legend.position="none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
# ggsave("3MT Reversal Comp.png",plot=p,units="in", width=2.5, height=2.5, dpi=600, device = 'png')
ggsave("3MT Reversal Comp.png",plot=p,units="in", width=3, height=2, dpi=600, device = 'png')


# Generate a plot for each inversion, with each line plotted
plot_combined_off_by_chrom(freq_data)

# Generate a plot for each replicate line, with each inversion plotted
# plot_by_line(freq_data)

# Generate scatter plots for each cohort of each inversion 
scatter_plot_by_chrom(freq_data)

# Generate scatter plots of parental-offspring differences
scatter_plot_par_diff_by_chrom(freq_data)

# Generate scatter plots of parental-embryo and embryo-adult differences
scatter_plot_embryo_diff_by_chrom(freq_data)




# # Run the power comparison by pat freq at a specific coverage
# plot_bootstrap_freq_sensitivity(cutoffs,alpha)
# 
# # Run the alternate plot script
# plot_bootstrap_cutoffs_alt(cutoffs,alpha)


