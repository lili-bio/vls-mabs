#setwd("/Users/as/Documents/ResearchProjects/vlsE_antibody_paper/")

library(ggtree)
library(ggplot2)

tree <- read.tree("trees/vls_codon_fast_cutoff.tree")
tips <- tree$tip.label

#set up groups
A <- tips[grep("64b|B31|PAbe|BOL26|ZS7", tips)]
B <- tips[grep("29805|80a|N40|WI91-23", tips)]
C <- tips[grep("118a|72a|94a", tips)]
D <- tips[grep("156a|JD1", tips)]

group_info <- list("A"=A,"B"=B,"C"=C,"D"=D)
tree2 <- groupOTU(tree, group_info)

p <- ggtree(tree2, aes(color=group), layout = "daylight", size=0.8)
p + geom_treescale(x=0.5, y=0.5, width = 0.1)



