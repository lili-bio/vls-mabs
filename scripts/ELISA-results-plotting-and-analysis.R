#setwd("C:/Users/as/Documents/ResearchProjects/vlsE_antibody_paper/ELISA_data/")
setwd("C:/Users/weiga/Dropbox/ms-vlsE-shared/bb-antigens")
#setwd("Dropbox/ms-vlsE-shared/bb-antigens/")
library(tidyverse)
library(drc)
library(broom)

#############################################
### Test vlsE peptides against host sera ####
#############################################
human <- read_tsv("data/ELISA_data/OD_211123_human_sera.txt")
human <- human %>% mutate(host = "human") %>% mutate(lyme = if_else(lyme == "NonLyme", "Control", lyme))
mouse <- read_tsv("data/ELISA_data/OD_211123_mouse_sera.txt")
mouse <- mouse %>% mutate(lyme = NA, host = "mouse")
rabbit <- read_tsv("data/ELISA_data/OD_211216_rabbit_sera.txt")
rabbit <- rabbit %>% mutate(lyme = NA, host = "rabbit")

sera <- bind_rows(human, mouse, rabbit)
sera %>% 
  ggplot(aes(x = antigen, y = OD)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, shape = 1, aes(color = lyme), alpha = 0.5) + 
  theme_bw() +
  labs(x="Antigen", y="OD450") +
  theme(axis.title = element_text(size=12)) + 
  facet_wrap(~host) + 
  theme(legend.position = "bottom")

human %>% ggplot(aes(x = lyme, y = OD)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(shape = 1, width = 0.2, shape = 1,) + 
  theme_bw()

#####################
# correlation & gression analysis
####################
human.wide <- human %>% 
  filter(antigen %in% c('IR4', 'IR6', 'VlsE')) %>% 
  select(-5) %>% 
  pivot_wider(names_from = "antigen", values_from = "OD")

human.long <- human.wide %>% 
  mutate(ir4 = IR4- BSA, ir6 = IR6 - BSA, vls = VlsE - BSA) %>% 
  select(-c(2:5)) %>% 
  pivot_longer(2:4, names_to = "IR", values_to = "OD")
  
human.long %>% 
  ggplot(aes(x = vls, y = OD)) + 
  geom_point(shape = 1) + 
  theme_bw() +
  labs(x="VlsE", y="IR") +
  geom_smooth(method = "lm") +
  facet_wrap(~IR) +
  theme(axis.title = element_text(size=12))

# regression
human.long %>% 
  group_by(IR) %>% 
  do(tidy(lm(data = ., vls ~ OD))) %>% 
  filter(!str_detect(term, "Intercept"))

# correlation
human.long %>%
  group_by(IR) %>%
  do(tidy(cor.test(data = ., ~vls + OD)))

human.long %>% 
  filter(IR == 'ir4') %>% 
  lm(data = ., vls ~ OD) %>% 
  summary()

human.long %>% 
  filter(IR == 'ir6') %>% 
  lm(data = ., vls ~ OD) %>% 
  summary()
#################
# ANOVA and t-test
###################
lm_human <- lm(data = human, OD~antigen)
summary(lm_human)

t_human <- t.test(data = human %>% filter(lyme != 'Control'), OD~lyme)
t_human

# check correlation between IR4 and vlsE
# vlse <- human %>% filter(antigen=="VlsE")
# ir4 <- human %>% filter(antigen=="IR4")
# od <- data.frame(VlsE=vlse$OD, IR4=ir4$OD)
# ggplot(data = od, aes(x=IR4, y=VlsE)) + geom_point() +
#  geom_smooth(method = "lm") + theme_bw()
# lm <- lm(data=od, VlsE~IR4)
s# ummary(lm)

# check correlation between IR6 and vlsE
# ir6 <- human %>% filter(antigen=="IR6")
# od <- data.frame(VlsE=vlse$OD, IR6=ir6$OD)
# ggplot(data = od, aes(x=IR6, y=VlsE)) + geom_point() +
#  geom_smooth(method = "lm") + theme_bw()
#lm <- lm(data=od, VlsE~IR6)
# summary(lm)

# boxplot for lyme disease stages
# ggplot(human, aes(x=lyme, y=OD)) + geom_boxplot() +
#  geom_jitter(width = 0.2) + theme_bw() + 
#  labs(x="Lyme disease stage", y="OD450") +
#  theme(axis.title = element_text(size=12))

# mouse serum collection
# ANOVA
lm_mouse <- lm(data = mouse , OD~antigen)
summary(lm_mouse)
# rabbit serum collection

# ANOVA
lm_rabbit <- lm(data = rabbit, OD~antigen)
summary(lm_rabbit)

################################################
### Test vlsE peptides against B cell lines ####
################################################
cells <- read_tsv("data/ELISA_data/OD_cell_lines.txt")
head(cells)
names <- as.character(cells$cell_line)[1:20]
cells$cell_line <- factor(cells$cell_line, levels=names)#for sorting purpose

cells %>% group_by(antigen) %>% summarise(mean(OD), sd(OD), mean(OD) + 3*sd(OD))

cells %>% filter(antigen %in% c("IR4", "IR6", "VlsE")) %>% 
  ggplot(aes(x=cell_line, y=OD)) +
  geom_bar(stat="identity", size=1) +
  theme_bw() + facet_wrap(~antigen, nrow=2) +
  labs(x="Cell line", y="OD450") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.8)) + 
  geom_hline(yintercept = 0.0478, linetype = 2)

###############################################
### Titration of 4 epitope-specific rMAbs  ####
###############################################

od <- read.table("data/ELISA_data/OD_220130_mAb_titer.txt", header = T)
head(od)
models<- drm(OD ~ cont, data = od, fct=LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")), curveid = antibody)
plot(models, type="all", broken = TRUE, 
     legend = TRUE, cex.legend = 1,
     legendPos = c(10, 1),
     xlab = "Concentration of antibody (ng/ml)", ylab = "OD450",
#     xt = c(1,4,16,64,256),
#     xtlab = c(1000,250,62.5,15.6,3.9)
)
#abline(a = 0.5, b=0, lty=2)
#text(20, 0.5, "EC50", pos=3)
ed <- ED(models, 50)
new <- data.frame(antibody = str_remove(rownames(ed), "e:") %>% str_remove(":50"), cont = ed[,1])
od.pred <- predict(models, new, se.fit = T)
df.ed <- tibble(antibody = new$antibody, ec50 = new$cont, se = ed[,2], od=od.pred[,1])

points(x = df.ed$ec50, y = df.ed$od , pch = 10, cex = 2)
legend(3.2, 0.65, pch = 10, bty = "n", legend  = "EC50", pt.cex = 2)


