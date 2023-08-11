path <- '/home/localadmin/Documents/MCDS_code/MCDS-dFMRI/Analysis/data/21_dir_benchmark_'
branching <- "no_branching"

# Read CSV into DataFrame
read_csv_no_branch = read.csv(paste(path, branching, sep = ""))
colnames(read_csv_no_branch) = c('X','loc','N','T', 'Sb', 'b')
read_csv_no_branch$b <- sprintf(read_csv_no_branch$b, fmt = '%#.3f')
read_csv_no_branch$T <- as.factor(read_csv_no_branch$T)
read_csv_no_branch$N <- as.factor(read_csv_no_branch$N)
read_csv_no_branch$b <- as.factor(read_csv_no_branch$b)
read_csv_no_branch$branch <- "no branch"

branching <- "branching"
# Read CSV into DataFrame
read_csv_branch = read.csv(paste(path, branching, sep = ""))
colnames(read_csv_branch) = c('X','loc','N','T', 'Sb', 'b')
read_csv_branch$b <- sprintf(read_csv_branch$b, fmt = '%#.3f')
read_csv_branch$T <- as.factor(read_csv_branch$T)
read_csv_branch$N <- as.factor(read_csv_branch$N)
read_csv_branch$b <- as.factor(read_csv_branch$b)
read_csv_branch$branch <- "branch"

data_all = rbind(read_csv_no_branch, read_csv_branch)

stats = compare_means(Sb ~ branch,  data = data_all[data_all$N=="10000" & data_all$T=="10000",], 
                      paired=FALSE, method="wilcox.test", p.adjust.method="fdr")
df_p_val <- rstatix::wilcox_test(data_all[data_all$N=="10000" & data_all$T=="10000",], Sb ~ branch, p.adjust.method = "fdr") %>% 
  rstatix::add_xy_position() %>% filter(p.adj.signif != "ns")

p <- ggplot(data_all[data_all$N=="10000" & data_all$T=="10000",], aes(x=b, y=Sb)) + 
  geom_boxplot() +
  ggtitle(paste('b = ', b_, ', N  = ', N_, sep = "")) +
  theme(plot.title = element_text(hjust = 0.5))

if(length(df_p_val$p) > 0)
{
  p <- p + add_pvalue(df_p_val, label = "p.adj.signif", label.size=8, inherit.aes = FALSE)
}

data_all$b <- as.numeric(data_all$b)
ggplot(data = data_all[data_all$N=="10000" & data_all$T=="10000" & data_all$b > 0,], aes(x = branch, y = Sb, fill = branch)) + 
  geom_boxplot() +
  facet_wrap(~b, strip.position = "bottom", scales = "free_x") +
  theme(panel.spacing = unit(0, "lines"), 
        strip.background = element_blank(),
        strip.placement = "outside")

stat.test <- data_all[data_all$N=="10000" & data_all$T=="10000" & data_all$b > 0,] %>%
  group_by(b) %>%
  wilcox_test(Sb ~ branch, p.adjust.method = "fdr") 
stat.test 
stat.test <- stat.test %>% add_xy_position(x = "branch", fun = "mean_se")
ggboxplot(
  data_all[data_all$N=="10000" & data_all$T=="10000" & data_all$b > 0,], x = "branch", y = "Sb", fill = "#00AFBB",
  add = c("mean_se", "dotplot"), facet = c("b"), scales="free_y"
) +
  stat_pvalue_manual(stat.test, hide.ns = TRUE, tip.length = 0, step.increase = 0) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+  theme(strip.background = element_blank())





stat.test <- data_all[data_all$N=="10000" & data_all$T=="10000" & data_all$b > 0,] %>%
  group_by(b) %>%
  wilcox_test(Sb ~ branch) 
stat.test 
stat.test <- stat.test %>% add_xy_position(x = "branch", fun = "mean_se")
ggplot(data_all[data_all$N=="10000" & data_all$T=="10000" & data_all$b > 0,], add="mean_se") + 
  geom_boxplot(aes(x="branch", y = "Sb", fill='branch')) + facet_wrap(~b) +
  stat_pvalue_manual(stat.test, hide.ns = TRUE, tip.length = 0, step.increase = 0) +
  theme(strip.background = element_blank())
