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
colnames(read_csv) = c('X','loc','N','T', 'Sb', 'b')
N_ <- '10000'
read_csv <- subset(read_csv, N==N_)

read_csv$b <- sprintf(read_csv$b, fmt = '%#.3f')
read_csv$T <- as.factor(read_csv$T)
read_csv$N <- as.factor(read_csv$N)
read_csv$b <- as.factor(read_csv$b)


p <- list()
i <- 1
for (b_ in unique(read_csv$b))
{
  if(b_ != "0.000" )
  {
    print(subset(read_csv, b==b_))
    
    stats = compare_means(Sb ~ T,  data = subset(read_csv, b==b_), 
                          paired=FALSE, method="wilcox.test", p.adjust.method="fdr")
    df_p_val <- rstatix::wilcox_test(subset(read_csv, b==b_ ), Sb ~ T, p.adjust.method = "fdr") %>% 
      rstatix::add_xy_position() %>% filter(p.adj.signif != "ns")
    
    idx = which(stats$p < 0.05)
    my_comparisons <- c()
    if (length(idx) > 0){
      for (k in 1:length(idx)){
        stats[idx[k],]$group1
        stats[idx[k],]$group2
        my_comparisons <- append(my_comparisons, 
                                 combn(c(as.numeric(stats[idx[k],]$group1), 
                                         as.numeric(stats[idx[k],]$group2)),2, 
                                       simplify = FALSE))
      }
    }
    p[[i]] <- ggplot(subset(read_csv, b==b_), aes(x=T, y=Sb)) + 
      geom_boxplot() +
      ggtitle(paste('b = ', b_, ', N  = ', N_, sep = "")) +
      theme(plot.title = element_text(hjust = 0.5))
    
    if(length(df_p_val$p) > 0)
    {
      p[[i]] <- p[[i]] + add_pvalue(df_p_val, label = "p.adj.signif", label.size=8, inherit.aes = FALSE)
    }
    
    i <- i + 1
  }
}
do.call(grid.arrange,p)

