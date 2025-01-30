# 
# Quick script for plotting the raw variance control data

# Set working directory
setwd("/Users/cmcallester/cMcAl/UWm/Pool/Inversion ∆ Experiment/Amplicon Analysis/Graphing")


load_freq_data <- function(freq_data_file){
  freq_data <- read.csv(file=freq_data_file, header=TRUE, sep=",")
  freq_data$Cohort <- factor(freq_data$Cohort, levels = c('Parental ♂','Embryo','Early ♀','Late ♀','Early ♂','Late ♂'))
  freq_data$Line <- factor(freq_data$Line)
  return(freq_data)
}


scatter_plot_replicates <- function(freq_data,title){
  cohort <- unique(freq_data$Line)
  inv <- unique(freq_data$Inv)
  colors <- c("Control 2" = "red", 
              "251" = "blue", 
              "254" = "darkgreen", 
              "418" = "orange")
  
  library(ggplot2)
  
  # Plotting the paternal-embryo differences using ggplot2
  p <- ggplot() + 
    geom_point(data=freq_data,aes(x=Difference_Cohort, y=Inv_Freq_Diff, colour=Line)) +
    geom_hline(yintercept=0) +
    xlab("Inversion ") + 
    ylab("Inversion Frequency Difference") + 
    ggtitle(title) +
    theme_bw() +
    # scale_y_continuous(limits=c(0,1)) +
    # scale_linetype_manual(values = c(1,2),
    #                       labels=c("Upper",
    #                                "Lower"),
    #                       name="Bound") +
    scale_colour_manual(name="Replicate Set", values = colors,
                        guide = guide_legend(title.position = "top",direction = "vertical")) +
    # theme(legend.justification=c(1,0),
    #       legend.position=c(0.95,0.05),  # Position legend in bottom right
    #       legend.box="horizontal",  # Position multiple legends horizontally
    #       legend.background = element_rect(fill="white",size=0.5, linetype="solid", colour ="black"),
    #       legend.title = element_text(size=10),
    #       legend.text = element_text(size=7),
    #       plot.title = element_text(hjust=1.0, size=13, face="bold")) +
    ggsave(paste0("In(",chrom,")",inv," Parental Diffs Scatter.png"), 
           units="in", width=6, height=4, dpi=600, device = 'png')
  
}

# Read inversion frequency CSV data for the 
freq_data <- load_freq_data("../Alignment & QC Scripts/variance_controls/inv_freqs/freq_data.csv")
