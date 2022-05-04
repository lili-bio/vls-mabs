#setwd("/Users/as/Documents/ResearchProjects/vlsE_antibody_paper/rate4site/")

library(genoPlotR)
library(ggplot2)

r4s <- read.table("r4s_vls.tsv", header = T)
head(r4s)
p <- ggplot(r4s, aes(x=POS, y=SCORE)) 
p1 <- p + annotate(geom = "rect", xmin=1, xmax=30, ymin=-1, ymax=3, fill="grey", alpha=0.5) +
  annotate(geom = "rect", xmin=46, xmax=58, ymin=-1, ymax=3, fill="grey", alpha=0.5) +
  annotate(geom = "rect", xmin=78, xmax=84, ymin=-1, ymax=3, fill="grey", alpha=0.5) +
  annotate(geom = "rect", xmin=91, xmax=116, ymin=-1, ymax=3, fill="grey", alpha=0.5) +
  annotate(geom = "rect", xmin=134, xmax=141, ymin=-1, ymax=3, fill="grey", alpha=0.5) +
    annotate(geom = "rect", xmin=154, xmax=178, ymin=-1, ymax=3, fill="grey", alpha=0.5)
p2 <- p1 + geom_point() + geom_line() +
  theme_classic() + labs(x="AA site", y="Variability score")
p2 + geom_hline(yintercept = 0, color="black", linetype="dashed")  
