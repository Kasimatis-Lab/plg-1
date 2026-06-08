library(dplyr)
library(reshape)
library(ggplot2)
library(ggbeeswarm)
library(RColorBrewer)
library(scales)

#Read in pseudo-female fecundity data
fem.fecundity <- read.table("~/Dropbox/Plg-1/data_files/FileS2_CaeNDR_Female_Fecundity.txt", header = TRUE, sep = "\t")
fem.fecundity.censored <- subset(fem.fecundity, Censor == 0)

#Examine strain means
summarise(group_by(JKfecundity.censored, Male), mean(Total), length(Total), sd(Total)/length(Total))

#Does male genotype impact total female fecundity?
#First examine each strain pair independently
model.jj <- glm(Total ~ m1.geno + Replicate, family = poisson(link = "log"), data = subset(fem.fecundity.censored, Male == "JT11398" | Male == "JU346"))
summary(model.jj)

model.em <- glm(Total ~ m1.geno + Replicate, family = poisson(link = "log"), data = subset(fem.fecundity.censored, Male == "ED3017" | Male == "MY920"))
summary(model.em)

model.ne <- glm(Total ~ m1.geno + Replicate, family = poisson(link = "log"), data = subset(fem.fecundity.censored, Male == "NIC266" | Male == "EG4725"))
summary(model.ne)
#Male genotype is significant in the same direction for each strain pair. Replicate is significant across JT/JU and ED/MY, but does not qualitatively impact results.
#Pool strains and replicates for more power.

model1 <- glm(Total ~ m1.geno, family = poisson(link = "log"), data = fem.fecundity.censored)
summary(model1)

#Does fecundity differ across days?
x <- melt(fem.fecundity.censored[, 3:10], id.vars = "Male", variable_name = "days")
names(x) <- c("male", "days", "fecundity")
x$time <- as.numeric(substr(x$days, start = 2, stop = 2))
x$fecundity <- as.numeric(x$fecundity)

y <- melt(fem.fecundity.censored[, c(3,13)], id.vars = "Male", variable_name = "geno")
names(y) <- c("male", "m1geno", "geno")

fem.fecundity.days <- cbind(x, y$geno)
names(fem.fecundity.days) <- c("male", "days", "fecundity", "time", "geno")

model2 <- glm(fecundity ~ geno + time + geno*time, family = poisson(link = "log"), data = subset(fem.fecundity.days, male != "JT11398"))
summary(model2)

pchisq(model2$deviance, model2$df.residual)

#Plot total fecundity
p1 <- ggplot(data = fem.fecundity.censored, (aes(x = m1.geno, y = Total))) + ylim(0, 620) +
  theme_classic() +
  theme(panel.grid = element_blank())

total <- p1 + geom_quasirandom(data = fem.fecundity.censored,
                             aes(x = m1.geno, y = Total, shape = factor(Male)),
                             varwidth = TRUE, width = 0.25, cex = 2.5, alpha = 0.65, color = "gray50") +
  scale_shape_manual(values = c(19, 17, 0, 15, 1, 2)) +
  stat_summary(fun.data = mean_se, geom = "pointrange", color = c("#74C476", "#9970AB"), cex = 0.7) +
  xlab("Male Genotype") + ylab("Female Reproductive Success") + 
  theme_classic()


#Plot fecundity per day
p2 <- ggplot(data = fem.fecundity.days, (aes(x = time, y = fecundity))) + ylim(0, 125) +
  scale_x_continuous(breaks = pretty(z$time, n = 7)) +
  scale_y_continuous(breaks = pretty(z$fecundity, n = 8)) +
  theme_classic() +
  theme(panel.grid = element_blank())

perday <- p2 +
  stat_summary(fun = "mean", geom = "pointrange", aes(x = time, y = fecundity, shape = factor(male), color = factor(geno)), alpha = 0.5, position = position_dodge(0.3)) +
  stat_summary(fun.data = mean_se, geom = "pointrange", aes(x = time, y = fecundity, color = factor(geno))) +
  scale_shape_manual(values = c(19, 17, 0, 15, 1, 2)) +
  scale_color_manual(values = c("#74C476", "#9970AB")) +
  xlab("Days of Adulthood") + ylab("Female Reproductive Success") + 
  theme_classic()


