# install.packages("ggpubr")
library(ggpubr)
# install.packages("ggplot2")
library(ggplot2)
# install.packages("ggprism")
library(ggprism)
library(gridExtra)
# #install.packages("stringr")  # Install & load stringr
# library("stringr")
library("rstatix")
# # Set working directory
# setwd(dir = "/Users/ines/Documents/SV/PDM/Code/R/")
# library(Metrics)
# # install.packages("pracma")
# library("pracma")
# install.packages("patchwork")
# install.packages("patchwork", repos = "https://cloud.r-project.org")
library("patchwork")

path <- '/home/localadmin/Documents/MCDS_code/MCDS-dFMRI/Analysis/data/21_dir_benchmark_'
branching <- "no_branching"

# Read CSV into DataFrame
read_csv = read.csv(paste(path, branching, sep = ""))
colnames(read_csv) = c('X','loc','N','T', 'Sb', 'adc', 'b')
read_csv$b <- as.numeric(read_csv$b)
data <- subset(read_csv, N!="100" & N!="1000" & T!="1000" & b > 0)
data$b <- sprintf(data$b, fmt = '%#.3f')
data$T <- as.factor(data$T)
data$N <- as.factor(data$N)
data$b <- as.factor(data$b)


p <- list()
i <- 1
if(branching == "branching")
{
  ylims <- list(c(0.86, 0.89), c(0.57, 0.6), c(0.41, 0.44), c(0.315, 0.345), c(0.12, 0.15))
} else
{
  ylims <- list(c(0.855, 0.895), c(0.565, 0.615), c(0.4, 0.44), c(0.295, 0.335), c(0.08, 0.12))
}
  
b_names <- c(
  `0.195` = "b = 0.195 [um²/ms]",
  `1.002` = "b = 1.002 [um²/ms]",
  `1.998` = "b = 1.998 [um²/ms]",
  `3.018` = "b = 3.018 [um²/ms]",
  `9.926` = "b = 9.926 [um²/ms]"
)
for (b_ in unique(data$b))
{
  p[[i]] <- ggplot(subset(data, b==b_), aes(x=N, y=Sb, fill=T)) + 
  geom_boxplot(aes(x=N, y=Sb, fill=T)) + 
  xlab("N * 1000") +
  ylab("Sb/S0") +
  facet_grid(~b, scales="free", labeller = as_labeller(b_names)) +
  scale_x_discrete(breaks=c("5000","10000","15000", "50000","55000", "75000", "100000", "125000", "150000"),
                    labels=c("5","10","15", "50","55", "75", "100", "125", "150")) +
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=16),
          strip.text.x = element_text(size = 16)) +
  ylim(ylims[[i]][[1]], ylims[[i]][[2]])
  if(i != 1)
  {
    p[[i]] <- ggplot(subset(data, b==b_), aes(x=N, y=Sb, fill=T)) + 
      geom_boxplot(aes(x=N, y=Sb, fill=T)) + 
      facet_grid(~b, scales="free", labeller = as_labeller(b_names)) + 
      ylab("") +
      xlab("N * 1000") +
      scale_x_discrete(breaks=c("5000","10000","15000", "50000","55000", "75000", "100000", "125000", "150000"),
                       labels=c("5","10","15", "50","55", "75", "100", "125", "150")) +
      theme(axis.text=element_text(size=14),
            axis.title=element_text(size=16),
            strip.text.x = element_text(size = 16))+
      ylim(ylims[[i]][[1]], ylims[[i]][[2]])
  }
  
  i <- i + 1
}
ggarrange(plotlist=p, ncol=5, common.legend = TRUE, legend="bottom")
# landscape : 18 x 5

p <- list()
i <- 1
# ylims <- list(c(0.865, 0.89), c(0.56, 0.615), c(0.4, 0.44), c(0.305, 0.34), c(0.09, 0.142))
for (b_ in unique(data$b))
{
    p[[i]] <- ggplot(subset(data, b==b_), aes(x=N, y=adc, fill=T)) + 
      geom_boxplot(aes(x=N, y=adc, fill=T)) + 
      facet_grid(~b, scales="free")

  
  i <- i + 1
}
ggarrange(plotlist=p, ncol=5, common.legend = TRUE, legend="bottom")
