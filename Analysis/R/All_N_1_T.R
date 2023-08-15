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
data <- subset(read_csv, N!="100" & N!="1000")
data$b <- sprintf(data$b, fmt = '%#.3f')
data$T <- as.factor(data$T)
data$N <- as.factor(data$N)
data$b <- as.factor(data$b)

read_csv$b <- sprintf(read_csv$b, fmt = '%#.3f')
read_csv$T <- as.factor(read_csv$T)
read_csv$N <- as.factor(read_csv$N)
read_csv$b <- as.factor(read_csv$b)

p1 <- ggplot(subset(read_csv, T=='1000'), aes(x=b, y=Sb, fill=factor(N))) + 
geom_boxplot() +
ggtitle('T = 1000') +
theme(plot.title = element_text(hjust = 0.5))

p2 <- ggplot(subset(read_csv, T=='5000'), aes(x=b, y=Sb, fill=factor(N))) + 
geom_boxplot() +
ggtitle('T = 5000') +
theme(plot.title = element_text(hjust = 0.5))

p3 <- ggplot(subset(read_csv, T=='10000'), aes(x=b, y=Sb, fill=factor(N))) + 
geom_boxplot() +
ggtitle('T = 10000') +
theme(plot.title = element_text(hjust = 0.5)) 

p4 <- ggplot(subset(read_csv, T=='15000'), aes(x=b, y=Sb, fill=factor(N))) + 
geom_boxplot() +
ggtitle('T = 15000') +
theme(plot.title = element_text(hjust = 0.5)) 

print(p1 + p2 + p3 + p4)

for (T_ in unique(data$T))
{
  p <- list()
  i <- 1
  for (b_ in unique(data$b))
  {
    if(b_ != "0.000" )
    {
      print(subset(data, T==T_ & b==b_))

      stats = compare_means(Sb ~ N,  data = subset(data, T==T_ & b==b_ ), 
                            paired=FALSE, method="wilcox.test", p.adjust.method="bonf")
      df_p_val <- rstatix::wilcox_test(subset(data, T==T_ & b==b_ ), Sb ~ N, p.adjust.method = "bonf") %>% 
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
      p[[i]] <- ggplot(subset(data, T==T_ & b==b_), aes(x=N, y=Sb)) + 
                geom_boxplot() +
                ggtitle(paste('b = ', b_, ', T  = ', T_ , sep = "")) +
                theme(plot.title = element_text(hjust = 0.5))
      
      if(length(df_p_val$p) > 0)
      {
        p[[i]] <- p[[i]] + add_pvalue(df_p_val, label = "p.adj.signif", label.size=8, inherit.aes = FALSE)
      }

      i <- i + 1
    }
  }
  do.call(grid.arrange,p)
} 

