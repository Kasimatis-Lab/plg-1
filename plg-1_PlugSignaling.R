library(dplyr)
library(ggplot2)

#### C. elegans ####
#Read in data
ce.pad.obs <- read.table("Dropbox/Plg-1/data_files/FileS7_CE_PlugSignaling.txt", header = TRUE, sep = "\t")

ce.pad.obs$crossGeno   <- interaction(ce.pad.obs$m1.geno, ce.pad.obs$m2.geno, sep = ".")
ce.pad.obs$crossStrain <- interaction(ce.pad.obs$Male1, ce.pad.obs$Male2, sep = ".")

#Analyze frequency in JU346 and JT11398
ju.jt <- subset(ce.pad.obs, Male1 == "JU346" | Male1 == "JT11398")
model.jujt <- lm(SpermPresent ~ m1.geno + Replicate, data = ju.jt)
summary(model.jujt)

jj    <- as.data.frame(summarise(group_by(ju.jt, crossGeno, Replicate), length(SpermPresent), sum(SpermPresent)))
names(jj) <- c("cross", "replicate", "scored", "sperm.present")
jj$freq <- signif(jj$sperm.present / jj$scored, digits = 3)

#Analyze frequency in EG4725 and NIC266
eg.nic <- subset(ce.pad.obs, Male1 == "EG4725" | Male1 == "NIC266")
model.egnic <- lm(SpermPresent ~ m1.geno + Replicate, data = eg.nic)
summary(model.egnic)

en <- as.data.frame(summarise(group_by(eg.nic, crossGeno, Replicate), length(SpermPresent), sum(SpermPresent)))
names(en) <- c("cross", "replicate", "scored", "sperm.present")
en$freq <- signif(en$sperm.present / en$scored, digits = 3)

#Analyze frequency in ED3017 and MY920
ed.my <- subset(ce.pad.obs, Male1 == "ED3017" | Male1 == "MY920")
model.edmy <- lm(SpermPresent ~ m1.geno + Replicate, data = ed.my)
summary(model.edmy)

em <- as.data.frame(summarise(group_by(ed.my, crossGeno, Replicate), length(SpermPresent), sum(SpermPresent)))
names(em) <- c("cross", "replicate", "scored", "sperm.present")
em$freq <- signif(em$sperm.present / em$scored, digits = 3)

y <- as.data.frame(summarise(group_by(pad.obs, cross.S, Replicate), length(SpermPresent), sum(SpermPresent), sum(spermathecae), sum(uterus), sum(plug)))
names(y) <- c("cross", "rep", "scored", "sperm.present", "in.spermathecae", "in.uterus", "in.plug")
y$freq <- signif(y$sperm.present / y$scored, digits = 3)

#Analyze all strains pooled: all sperm present, in spermathecae, in uterus, in plug
byGeno <- as.data.frame(summarise(group_by(ce.pad.obs, crossGeno), length(SpermPresent), sum(SpermPresent), sum(spermathecae), sum(uterus), sum(plug)))
names(byGeno) <- c("cross", "scored", "sperm.present", "in.spermathecae", "in.uterus", "in.plug")
byGeno$freq <- signif(byGeno$sperm.present / byGeno$scored, digits = 3)

model1 <- lm(SpermPresent ~ m1.geno + Replicate, data = ce.pad.obs)
summary(model1)

model2 <- lm(spermathecae ~ m1.geno, data = ce.pad.obs)
summary(model2)

model3 <- lm(uterus ~ m1.geno, data = ce.pad.obs)
summary(model3)

model4 <- lm(plug ~ m1.geno, data = ce.pad.obs)
summary(model4)

model5 <- lm(SpermPresent ~ m2.geno + Replicate, data = ce.pad.obs)
summary(model5)

#Is the observed second male mating different from random chance (0.5)?
prop.test(x = sum(ce.pad.obs$SpermPresent), n = length(ce.pad.obs$SpermPresent), p = 0.5, alternative="two.sided", conf.level = 0.95) 

#Plot
byStrain        <- as.data.frame(summarise(group_by(ce.pad.obs, crossStrain, crossGeno), length(SpermPresent), sum(SpermPresent)))
names(byStrain) <- c("strains", "cross", "scored", "sperm.present")
byStrain$freq   <- signif(byStrain$sperm.present / byStrain$scored, digits = 3)
byStrain$pair  <- c(1, 1, 2, 2, 3, 3, 3, 3, 1, 1, 2, 2)

p1 <- ggplot(data = byStrain, (aes(x = cross, y = freq))) + ylim(0.5, 1) +
  theme_classic() +
  theme(panel.grid = element_blank())

all <- ggplot(byStrain, aes(x = cross, y = freq, shape = factor(pair))) + ylim(0.5, 1) +
  geom_jitter(position = position_jitter(0.05), cex = 3) + 
  stat_summary(fun.data = mean_se, geom = "pointrange", pch = 19, col = alpha(c("#74C476", "#9970AB", "#74C476", "#9970AB"), 0.95)) + 
  labs( x = "Cross Genotype", y = "Frequency Second Male Mated") +
  theme_classic()

####


#### C. remanei ####
#Read in data
cre.pad.obs <- read.table("~/Dropbox/Plg-1/data_files/FileS8_Cre_PlugSignaling.txt", header = TRUE, sep = "\t")

#Examine frequency by strain
byStrainCre        <- as.data.frame(summarise(group_by(cre.pad.obs, Male), length(Sperm.Present), sum(Sperm.Present)))
names(byStrainCre) <- c("male", "scored", "sperm.present")
byStrainCre$freq   <- signif(byStrainCre$sperm.present / byStrainCre$scored, digits = 3)
byStrainCre$cross  <- rep("plus.plus", 3)

#Are the strains significantly different?
model6 <- lm(Sperm.Present ~ Male, data = cre.pad.obs)
summary(model6)

#Is the observed second male mating different from random chance (0.5)?
prop.test(x = sum(cre.pad.obs$Sperm.Present), n = length(cre.pad.obs$Sperm.Present), p = 0.5, alternative="two.sided", conf.level = 0.95) 

#Plot
p2 <- ggplot(data = byStrainCre, (aes(x = cross, y = freq))) + ylim(0.5, 1.01) +
  theme_classic() +
  theme(panel.grid = element_blank())

cre <- ggplot(byStrainCre, aes(x = cross, y = freq, shape = factor(male))) + ylim(0.5, 1) +
  geom_jitter(position = position_jitter(0.08), cex = 3) + 
  stat_summary(fun.data = mean_se, geom = "pointrange", pch = 19, col = alpha(c("#9970AB"), 0.95)) + 
  labs( x = "Cross Genotype", y = "Frequency Second Male Mated") +
  theme_classic()

####



