#setwd("C:/Users/as/Documents/ResearchProjects/vlsE_antibody_paper/ELISA_data/")
setwd("C:/Users/weiga/Dropbox/ms-vlsE-shared/vls-mabs/")
#setwd("~/Dropbox/ms-vlsE-shared/vls-mabs/")

library(tidyverse)
library(drc)
library(broom)

rm(list = ls())

#############################################
### Test vlsE peptides against host sera ####
#############################################
human <- read_tsv("data/OD_211123_human_sera.txt")
human <- human %>% mutate(host = "human", immunization = NA) 
#%>% mutate(LymeStage = if_else(lyme == "NonLyme", "Control", lyme))
mouse <- read_tsv("data/OD_211123_mouse_sera.txt")
mouse <- mouse %>% mutate(LymeStage = NA, host = "mouse", immunization = NA)
rabbit <- read_tsv("data/OD_220621_rabbit_sera.txt")
rabbit <- rabbit %>% mutate(LymeStage = NA, host = "rabbit")
rabbit <- rabbit %>% dplyr::select(1:3,5:6,4)

sera <- bind_rows(human, mouse, rabbit)

# Fig 5B, by host
sera %>% 
  ggplot(aes(x = antigen, y = OD)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, shape = 1) + 
  theme_bw() +
  labs(x="Antigen", y="OD450") +
  theme(axis.title = element_text(size=12)) + 
  facet_wrap(~host)  
#  theme(legend.position = "bottom")

# Fig 5A. human, by stage
human %>% 
#  filter(LymeStage != 'Control') %>% 
  ggplot(aes(x = antigen, y = OD)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(shape = 1, width = 0.2) + 
  facet_wrap(~LymeStage) + 
  theme_bw()

#####################
# correlation & regression analysis
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

# Fig 5A, IR4
t.test(data = human %>% filter(LymeStage %in% c('Control', 'Early') & antigen == 'IR4'), OD~LymeStage) # ns

t.test(data = human %>% filter(LymeStage %in% c('Control', 'Late') & antigen == 'IR4'), OD~LymeStage) #ns


# IR6
t.test(data = human %>% filter(LymeStage %in% c('Control', 'Early') & antigen == 'IR6'), OD~LymeStage)

t.test(data = human %>% filter(LymeStage %in% c('Control', 'Late') & antigen == 'IR6'), OD~LymeStage)

# VlsE
t.test(data = human %>% filter(LymeStage %in% c('Control', 'Early') & antigen == 'VlsE'), OD~LymeStage)

t.test(data = human %>% filter(LymeStage %in% c('Control', 'Late') & antigen == 'VlsE'), OD~LymeStage)

# Early vs Control
t.test(data = human %>% filter(LymeStage %in% c('Control', 'Early')), OD~LymeStage)

# Late vs Control
t.test(data = human %>% filter(LymeStage %in% c('Control', 'Late')), OD~LymeStage)

lm.early <- lm(data = human %>% filter(lyme == 'Early'), OD~antigen)

summary(lm.early)

lm.late <- lm(data = human %>% filter(lyme == 'Late'), OD~antigen)

summary(lm.late)

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
cells <- read_tsv("data/OD_cell_lines.txt")
head(cells)
names <- as.character(cells$cell_line)[1:20]
cells$cell_line <- factor(cells$cell_line, levels=names)#for sorting purpose

cells %>% group_by(antigen) %>% summarise(mean(OD), sd(OD), mean(OD) + 3*sd(OD))

cells %>% filter(antigen %in% c("IR4", "IR6", "VlsE")) %>% 
  ggplot(aes(x=cell_line, y=OD, fill = antigen)) +
  geom_bar(stat="identity", size=1, position = "dodge") +
  theme_bw() + 
#  facet_wrap(~antigen, nrow=3) +
  labs(x="Cell lines", y="OD450") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.8), legend.position = "top") + 
  scale_fill_manual(values = c("#cccccc", "black", "#eeeeee")) + 
  # scale_color_manual(values = c("darkgray", "black", "gray")) + 
  geom_hline(yintercept = 0.0478, linetype = 2) 
#  coord_flip()

###############################################
### Titration of 4 epitope-specific rMAbs  ####
###############################################

od <- read.table("data/OD_220130_mAb_titer.txt", header = T)
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

points(x = df.ed$ec50, y = df.ed$od , pch = 16, cex = 2)
legend(3.2, 0.65, pch = 16, bty = "n", legend  = "EC50", pt.cex = 2)


