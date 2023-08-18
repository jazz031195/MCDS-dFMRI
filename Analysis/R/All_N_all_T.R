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
ylims <- list(c(0.865, 0.89), c(0.56, 0.615), c(0.4, 0.44), c(0.305, 0.34), c(0.09, 0.142))
for (b_ in unique(data$b))
{
  p[[i]] <- ggplot(subset(data, b==b_), aes(x=N, y=Sb, fill=T)) + 
  geom_boxplot(aes(x=N, y=Sb, fill=T)) + 
  facet_grid(~b, scales="free") +
  ylim(ylims[[i]][[1]], ylims[[i]][[2]])
  
  i <- i + 1
}
do.call(grid.arrange,p)


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
do.call(grid.arrange,p)
