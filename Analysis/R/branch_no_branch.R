path <- '/home/localadmin/Documents/MCDS_code/MCDS-dFMRI/Analysis/data/21_dir_benchmark_'
branching <- "no_branching"

# Read CSV into DataFrame
read_csv_no_branch = read.csv(paste(path, branching, sep = ""))
colnames(read_csv_no_branch) = c('X','loc','N','T', 'Sb', 'b')
read_csv_no_branch$b <- sprintf(read_csv_no_branch$b, fmt = '%#.3f')
read_csv_no_branch$b <- as.numeric(read_csv_no_branch$b)
# read_csv_no_branch$T <- as.factor(read_csv_no_branch$T)
# read_csv_no_branch$N <- as.factor(read_csv_no_branch$N)
# read_csv_no_branch$b <- as.factor(read_csv_no_branch$b)
read_csv_no_branch$branch <- "no branch"

branching <- "branching"
# Read CSV into DataFrame
read_csv_branch = read.csv(paste(path, branching, sep = ""))
colnames(read_csv_branch) = c('X','loc','N','T', 'Sb', 'b')
read_csv_branch$b <- sprintf(read_csv_branch$b, fmt = '%#.3f')
read_csv_branch$b <- as.numeric(read_csv_branch$b)
# read_csv_branch$T <- as.factor(read_csv_branch$T)
# read_csv_branch$N <- as.factor(read_csv_branch$N)
# read_csv_branch$b <- as.factor(read_csv_branch$b)
read_csv_branch$branch <- "branch"

data_all = rbind(read_csv_no_branch, read_csv_branch)

stat.test <- data_all[data_all$N=="10000" & data_all$T=="10000" & data_all$b > 0.000,] %>%
  group_by(b) %>%
  wilcox_test(Sb ~ branch, p.adjust.method = "fdr") 
stat.test 
stat.test <- stat.test %>% add_xy_position(x = "branch", fun = "min")
stat.test$y.position <- c(0.874, 0.588, 0.51, 0.413, 0.15)
ggboxplot(
  data_all[data_all$N=="10000" & data_all$T=="10000" & data_all$b > 0,], x = "branch", y = "Sb",
  add = c("mean_se", "dotplot"), facet.by = c("b"), scales="free_y", ncol=5, 
  ylab="Sb/S0", xlab="B values [um²/ms]", strip.position="bottom", color="branch") +
  stat_pvalue_manual(stat.test, hide.ns = TRUE) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  theme(strip.text.x = element_text(size=12,
                                    face="bold"),
        strip.background = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) 





stat.test <- data_all[data_all$N=="10000" & data_all$T=="10000" & data_all$b > 0,] %>%
group_by(b) %>%
wilcox_test(Sb ~ branch) 
stat.test 
stat.test <- stat.test %>% add_xy_position(x = "branch", fun = "mean_se")
stat.test$y.position <- c(0.874, 0.588, 0.51, 0.413, 0.15)
ggplot(data_all[data_all$N=="10000" & data_all$T=="10000" & data_all$b > 0,], 
       aes_string(x="as.factor(branch)", y = "Sb", fill="branch")) + 
geom_boxplot() + facet_wrap(~b, scales="free", ncol=5) +
stat_pvalue_manual(stat.test, hide.ns = TRUE, inherit.aes=FALSE) + 
theme(strip.text.x = element_text(size=12,
                                  face="bold"),
      strip.background = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank())
