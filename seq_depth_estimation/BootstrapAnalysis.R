# 
# For plotting the bootstrap results on Freq-seq simulations

# Set working directory
setwd("/Users/cmcallester/Documents/Pool Lab/Inversion Freq-seq/Bootstrapping")
# setwd(getSrcDirectory(function(x) {x})[1])


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
          plot.title = element_text(hjust=1.0, size=13, face="bold")) +
    ggsave(paste("P-E Bootstrap Diffs at",alpha,".png"), 
           units="in", width=6, height=4, dpi=600, device = 'png')
  
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
          plot.title = element_text(hjust=1.0, size=13, face="bold")) +
    ggsave(paste("E-A Bootstrap Diffs at",alpha,".png"), 
           units="in", width=6, height=4, dpi=600, device = 'png')
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
          plot.title = element_text(hjust=1.0, size=12, face="bold")) +
    ggsave(paste("Bootstrap Diffs at",alpha,"quant,",coverage,"depth.png"), 
           units="in", width=6, height=4, dpi=600, device = 'png')
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
          plot.title = element_text(hjust=1.0, size=13, face="bold")) +
    ggsave(paste("P-E Bootstrap Diffs at",alpha,"alt.png"), 
           units="in", width=6, height=4, dpi=600, device = 'png')
  
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
          plot.title = element_text(hjust=1.0, size=13, face="bold")) +
    ggsave(paste("E-A Bootstrap Diffs at",alpha,"alt.png"), 
           units="in", width=6, height=4, dpi=600, device = 'png')
}

# Running as a script:

# Pick an alpha value (the quantile high or low to asses the significance cutoff at)
# (here we use 2.5% to represent 5% from the mean considering both directions)
alpha = 0.025

# Read Bootstrap CSV data
# reps <- read.csv(file="bootstrap_freq_2.csv", header=TRUE, sep=",")

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

# Run the plot script
plot_bootstrap_cutoffs(cutoffs,alpha)

# Run the power comparison by pat freq at a specific coverage
plot_bootstrap_freq_sensitivity(cutoffs,alpha)

# Run the alternate plot script
plot_bootstrap_cutoffs_alt(cutoffs,alpha)


