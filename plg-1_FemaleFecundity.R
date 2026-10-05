library(dplyr)
library(MASS)
library(multcomp)
library(reshape)
library(ggplot2)
library(ggbeeswarm)
library(RColorBrewer)
library(scales)
library(fitdistrplus)

#### Analyze Fecundity Data ####
#Read in pseudo-female fecundity data
fem.fecundity <- read.table("~/Dropbox/Plg-1/data_files/FileS6_PseudoFemale_Fecundity.txt", header = TRUE, sep = "\t")
fem.fecundity.censored <- subset(fem.fecundity, Censor == 0)
fem.fecundity.censored$StrainPair <- as.factor(fem.fecundity.censored$StrainPair)

#Examine strain means
summarise(group_by(fem.fecundity.censored, Male), mean(Total), length(Total), sd(Total)/length(Total))

#Examine distribution of total female fecundity
descdist(fem.fecundity.censored$Total, discrete = FALSE, graph = TRUE, boot = 500)

fit_norm <- fitdist(fem.fecundity.censored$Total, "norm")
plot(fit_norm)

fit_nbinom <- fitdist(fem.fecundity.censored$Total, "nbinom")
plot(fit_nbinom)

fit_pois <- fitdist(fem.fecundity.censored$Total, "pois")

fit.metrics <- lapply(ls(pattern = "fit\\_"), function(variable) {
  fit = get(variable, envir = .GlobalEnv)
  with(fit, data.frame(name = variable, aic, loglik))
})
do.call(rbind, fit.metrics)
#The normal distribution is the best fit: smallest AIC and largest log-likelihood.

#Does male genotype impact total female fecundity?
#First examine each strain pair independently
model.jj <- glm(Total ~ m1.geno + Replicate, family = gaussian(link = "identity"), data = subset(fem.fecundity.censored, Male == "JT11398" | Male == "JU346"))
summary(model.jj)

model.em <- glm(Total ~ m1.geno + Replicate, family = gaussian(link = "identity"), data = subset(fem.fecundity.censored, Male == "ED3017" | Male == "MY920"))
summary(model.em)

model.ne <- glm(Total ~ m1.geno + Replicate, family = gaussian(link = "identity"), data = subset(fem.fecundity.censored, Male == "NIC266" | Male == "EG4725"))
summary(model.ne)

pvals <- c(coef(summary(model.jj))[2,4], coef(summary(model.em))[2,4], coef(summary(model.ne))[2,4])
p.adjust(pvals, method = "holm", n = length(pvals))
#Male genotype is significant in the same direction for JT/JU and ED/MY, but not significant for EG/NIC (though trending the same).
#Replicate is not significant.

#Pool strains and replicates for more power.
model1 <- glm(Total ~ m1.geno + StrainPair, family = gaussian(link = "identity"), data = fem.fecundity.censored)
summary(model1)
PWmodel1 <- glht(model1, linfct = mcp(StrainPair = "Tukey"))
summary(PWmodel1)

#Does fecundity differ across days?
x <- melt(fem.fecundity.censored[, 3:10], id.vars = "Male", variable_name = "days")
names(x) <- c("male", "days", "fecundity")
x$time <- as.numeric(substr(x$days, start = 2, stop = 2))
x$fecundity <- as.numeric(x$fecundity)

y <- melt(fem.fecundity.censored[, c(3,13)], id.vars = "Male", variable_name = "geno")
names(y) <- c("male", "m1geno", "geno")

z <- melt(fem.fecundity.censored[, c(3,14)], id.vars = "Male", variable_name = "pair")
names(z) <- c("male", "strainpair", "pair")

fem.fecundity.days <- cbind(x, y$geno, z$pair)
names(fem.fecundity.days) <- c("male", "days", "fecundity", "time", "geno", "pair")

#Examine distribution of female fecundity by day
descdist(fem.fecundity.days$fecundity, discrete = FALSE, graph = TRUE, boot = 500)

fit_norm <- fitdist(fem.fecundity.days$fecundity, "norm")
plot(fit_norm)

fit_nbinom <- fitdist(fem.fecundity.days$fecundity, "nbinom")
plot(fit_nbinom)

fit_gamma <- fitdist(fem.fecundity.days$fecundity, "gamma")
plot(fit_gamma)

fit.metrics <- lapply(ls(pattern = "fit\\_"), function(variable) {
  fit = get(variable, envir = .GlobalEnv)
  with(fit, data.frame(name = variable, aic, loglik))
})
do.call(rbind, fit.metrics)
#The normal distribution and negative binomial distributions are equally good in fit: smallest AICs and largest log-likelihoods.

model2 <- glm(fecundity ~ geno + time + geno*time + pair, family = gaussian(link = "identity"), data = fem.fecundity.days)
summary(model2)
PWmodel2 <- glht(model2, linfct = mcp(pair = "Tukey"))
summary(PWmodel2)

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
####



#### Analyze Sperm Transferred per Ejaculate ####
sperm <- read.table("Dropbox/plg-1/data_files/SpermCounts.txt", header = TRUE, sep = "\t")

summarise(group_by(sperm, Strain), mean(SpermCount), sd(SpermCount)/length(SpermCount))
summarise(group_by(subset(sperm, Strain != "EG4725" & Strain != "NIC266"), Geno), mean(SpermCount), sd(SpermCount)/length(SpermCount))

power.t.test(n = 35, delta = 6.3, sd = sd(sperm$SpermCount), sig.level = 0.05, power = NULL)

t.test(subset(sperm$SpermCount, sperm$Strain=="ED3017"), subset(sperm$SpermCount, sperm$Strain=="MY920"))
t.test(subset(sperm$SpermCount, sperm$Strain=="JU346"), subset(sperm$SpermCount, sperm$Strain=="JT11398"))
t.test(subset(sperm$SpermCount, sperm$Strain=="EG4725"), subset(sperm$SpermCount, sperm$Strain=="NIC266"))

model3 <- lm(SpermCount ~ Geno, data = sperm)
summary(model3)

model4 <- lm(SpermCount ~ Geno, data = subset(sperm, Strain != "EG4725" & Strain != "NIC266"))
summary(model4)

p1 <- ggplot(data = sperm, (aes(x = SpermCount, y = Strain))) +
  scale_x_continuous(breaks = pretty(sperm$SpermCount, n = 10)) +
  theme_classic() +
  theme(panel.grid = element_blank())

counts <- p1 + geom_quasirandom(data = sperm, aes(x = SpermCount, y = Strain, color = factor(Geno)), 
                                varwidth = TRUE, width = 0.15, cex = 2.5, alpha = 0.65) +
  scale_color_manual(values = c("#74C476", "#9970AB")) +
  stat_summary(fun.data = mean_se, geom = "pointrange", color = "black", cex = 0.7) +
  xlab("Male Genotype") + ylab("Sperm Transferred in Single Ejaculate") + 
  theme_classic()

p2 <- ggplot(data = sperm, aes(x = SpermCount)) +
  geom_histogram(color = "gray50", alpha = 0.5, binwidth = 5) +
  geom_vline(xintercept = mean(subset(sperm$SpermCount, sperm$Geno=="plugging")), color = "#9970AB") +
  geom_vline(xintercept = mean(subset(sperm$SpermCount, sperm$Geno=="nonplugging")), color = "#74C476") +
  scale_x_continuous(breaks = pretty(sperm$SpermCount, n = 10)) +
  theme_classic()
####
