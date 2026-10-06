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

y <- as.data.frame(summarise(group_by(ce.pad.obs, crossStrain), length(SpermPresent), sum(SpermPresent), sum(spermathecae), sum(uterus), sum(plug)))
names(y) <- c("cross", "scored", "sperm.present", "in.spermathecae", "in.uterus", "in.plug")
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

#Is the observed second male mating different from first male mating success?
mated <- read.table("Dropbox/Plg-1/data_files/FileS_CE_PlugSignalingNull.txt", header = TRUE, sep = "\t")
mated$interactionG <- interaction(mated$m1.geno, mated$m2.geno, sep = ".")

D1 <- as.data.frame(summarise(group_by(mated, interactionG), sum(D1eggs), length(FemaleID)))
names(D1) <- c("treatment", "matedF", "total")
D1$freq <- D1$matedF / D1$total

D2 <- as.data.frame(summarise(group_by(ce.pad.obs, crossGeno), sum(SpermPresent), length(PadID)))
names(D2) <- c("treatment", "sperm", "total")
D2$freq <- D2$sperm / D2$total

prop.test(D2$sperm[1], D2$total[1], p = D1$freq[1])
prop.test(D2$sperm[2], D2$total[2], p = D1$freq[2])
prop.test(D2$sperm[3], D2$total[3], p = D1$freq[3])
prop.test(D2$sperm[4], D2$total[4], p = D1$freq[4])

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
  geom_hline(yintercept = sum(D1$matedF)/sum(D1$total)) +
  labs( x = "Cross Genotype", y = "Frequency Second Male Mated") +
  theme_classic()



#Does a mating plug delay a second male from mating?
delayNull <- read.table("Dropbox/Plg-1/data_files/FileS_CE_PlugDelayNull.txt", header = TRUE, sep = "\t")
delayNull$interaction <- interaction(delayNull$m1.geno, delayNull$m2.geno)

delay <- read.table("Dropbox/Plg-1/data_files/FileS_CE_PlugDelay.txt", header = TRUE, sep = "\t")
delay$interaction <- interaction(delay$m1.geno, delay$m2.geno)
delay$cross <- interaction(delay$D1Male, delay$D2Male)

null <- as.data.frame(summarise(group_by(delayNull, interaction), sum(D1eggs, na.rm = TRUE), length(D1eggs)))
names(null) <- c("treatment", "matedF", "total")
null$freq <- null$matedF / null$total

obs <- as.data.frame(summarise(group_by(delay, cross, interaction), sum(SpermPresent, na.rm = TRUE), length(SpermPresent), print(pair[1])))
names(obs) <- c("cross", "treatment", "sperm", "total", "pair")
obs$freq <- obs$sperm / obs$total

obs2 <- as.data.frame(summarise(group_by(delay, interaction), sum(SpermPresent, na.rm = TRUE), length(SpermPresent)))
names(obs2) <- c("treatment", "sperm", "total")
obs2$freq <- obs2$sperm / obs2$total

prop.test(obs2$sperm[1], obs2$total[1], p = null$freq[1])
prop.test(obs2$sperm[2], obs2$total[2], p = null$freq[2])
prop.test(obs2$sperm[3], obs2$total[3], p = null$freq[3])
prop.test(obs2$sperm[4], obs2$total[4], p = null$freq[4])

prop.test(c(obs2$sperm[1], obs2$sperm[3]), c(obs2$total[1], obs2$total[3]), p = c(0.4, 0.2952))

model6 <- lm(SpermPresent ~ m1.geno + m2.geno + Replicate, data = delay)
summary(model6)

model7 <- lm(D1eggs ~ D1Male, data = delayNull)
summary(model7)
summarize(group_by(delayNull, D1Male), sum(D1eggs, na.rm = TRUE), length(D1eggs))

#Plot
p3 <- ggplot(obs, aes(x = treatment, y = freq, shape = factor(pair))) + ylim(-0.01, 0.5) +
  geom_jitter(position = position_jitter(0.05), cex = 2) + 
  stat_summary(fun.data = mean_se, geom = "pointrange", pch = 19, col = alpha(c("#74C476", "#9970AB", "#74C476", "#9970AB"), 0.95)) + 
  geom_point(aes(x = 1, y = null$freq[1]), color = "red") +
  geom_point(aes(x = 2, y = null$freq[2]), color = "red") +
  geom_point(aes(x = 3, y = null$freq[3]), color = "red") +
  geom_point(aes(x = 4, y = null$freq[4]), color = "red") +
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

#Is the observed second male mating different from first male mating success?
#All D1 males mated (N = 36 per strain)
prop.test(x = byStrainCre$sperm.present, n = byStrainCre$scored, p = c(0.999, 0.972, 0.972))

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



