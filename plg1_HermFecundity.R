library(dplyr)
library(lme4)
library(ape)
library(nlme)
library(geiger)
library(phytools)
library(ggplot2)
library(ggbeeswarm)
library(RColorBrewer)
library(scales)

#Read in data
fecundity <- read.table("FileS5_Hermaphrodite_Total_Fecundity.txt", header = TRUE, sep = "\t")
fecundity.censored <- subset(fecundity, Censor == 0 & Strain != "JU2526" & Strain != "QG4006")


#Examine strain means
strain.mean <- as.data.frame(summarise(group_by(fecundity.censored, Strain), mean(Total), (sd(Total) / sqrt(length(Plate.ID))), length(Plate.ID),
                                       levels(as.factor(Plugging))))
names(strain.mean) <- c("strain", "mean", "se", "n.herm", "plugging")
strain.mean <- strain.mean[order(strain.mean$plugging, decreasing = TRUE),]


#Linear model with replicate as a fixed effect
#Does replicate qualitatively change results?
model1 <- lm(Total ~ Plugging + Replicate, data = fecundity.censored)
summary(model1)
#Significant effect of rep3 -- remove and examine model
model2 <- lm(Total ~ Plugging + Replicate, data = subset(fecundity.censored, Replicate != "rep3"))
summary(model2)
#Removing rep3 does not qualitatively change results, so merge replicates.


#PGLS model
support_tree <- read.tree("WI.20250625.hard-filter.min4.tree")

obj <- name.check(support_tree, levels(as.factor(fecundity.censored$Strain)), data.names = levels(as.factor(fecundity.censored$Strain)))
tree.cut <- drop.tip(support_tree, obj$tree_not_data, keep.root.edgue = TRUE)
name.check(tree.cut, levels(as.factor(fecundity.censored$Strain)), data.names = levels(as.factor(fecundity.censored$Strain)))

tree.rooted <- root.phylo(tree.cut, outgroup = "CB4856", resolve.root = TRUE)

fecundityPGLS <- fecundity.censored

fecundityPGLS$Plugging[fecundityPGLS$Plugging == "no"] <- 0
fecundityPGLS$Plugging[fecundityPGLS$Plugging == "yes"] <- 1
fecundityPGLS$Plugging <- as.numeric(fecundityPGLS$Plugging)

grouped <- split(fecundityPGLS, fecundityPGLS$Strain)
total <- lapply(grouped, "[[","Total")
plug <- lapply(grouped, "[[","Plugging")

pic.Tot <- pic.ortho(total, tree.rooted, intra = TRUE)
pic.Plg <- pic.ortho(plug, tree.rooted, intra = TRUE)

picModel <- lm(pic.Tot ~ pic.Plg - 1)
summary(picModel)


#Plot phylogeny
plot_ord <- tree.rooted$tip.label
sum_data <- as.data.frame(summarise(group_by(fecundity.censored, Strain), levels(as.factor(Plugging))))
names(sum_data) <- c("strain", "plugging")
sum_data$color <- ifelse(sum_data$plugging=="yes", "#9970AB", "#74C476")
sumdata_ord <- sum_data[match(plot_ord, sum_data$strain),]

plot(tree.rooted, font = 2, align.tip.label = FALSE, tip.color = sumdata_ord$color)


#Examine male frequency
fecundity.censored$males <- substr(fecundity.censored$Notes, start = 4, stop = 4)
m <- subset(fecundity.censored, males != "d" & males != "")
m$males <- as.numeric(m$males)
x <- as.data.frame(summarise(group_by(m, Plugging), sum(males), length(males), sum(Total)))
names(x) <- c("plugging", "males", "length_males", "total")
x$freq <- signif(x$males / x$total, digits = 3)
prop.test(x = c(x[1,2], x[2,2]), n = c(x[1,4], x[2,4]), alternative = "two.sided")


#Plot fecundity
y <- as.geom_label()
y <- as.data.frame(summarise(group_by(m, Strain), sum(males), sum(Total), levels(as.factor(Plugging))))
names(y) <- c("strain", "males", "total", "plugging")
y$freq <- signif(y$males / y$total, digits = 1)

p1 <- ggplot(data = subset(fecundity.censored, Plugging == "yes"), (aes(x = Strain, y = Total))) + ylim(-5, 410) +
  theme_classic() +
  theme(panel.grid = element_blank())

plg <- p1 + geom_quasirandom(data = subset(fecundity.censored, Plugging == "yes"),
                             aes(x = Strain, y = Total),
                             varwidth = TRUE, width = 0.25, cex = 2.5, alpha = 0.5, color = "gray60") +
  stat_summary(fun.data = mean_se, geom = "pointrange", color = "#9970AB", cex = 0.8) +
  geom_text(data = subset(y, plugging == "yes"), aes(x = strain, y = -5, label = freq), cex = 3) +
  xlab("Strain") + ylab("Total Hermaphrodite Reproductive Success") + 
  theme_classic()


p2 <- ggplot(data = subset(fecundity.censored, Plugging == "no"), (aes(x = Strain, y = Total))) + ylim(-5, 410) +
  theme_bw() +
  theme(panel.grid = element_blank())

non <- p2 + geom_quasirandom(data = subset(fecundity.censored, Plugging == "no"),
                             aes(x = Strain, y = Total),
                             varwidth = TRUE, width = 0.25, cex = 2.5, alpha = 0.5, color = "gray60") +
  stat_summary(fun.data = mean_se, geom = "pointrange", color = "#74C476", cex = 0.8) +
  geom_text(data = subset(y, plugging == "no"), aes(x = strain, y = -5, label = freq), cex = 3) +
  xlab("Strain") + ylab("Total Hermaphrodite Reproductive Success") + 
  theme_classic()

p3 <- ggplot(data = fecundity.censored, (aes(x = Plugging, y = Total))) + ylim(-5, 410) +
  theme_bw() +
  theme(panel.grid = element_blank())

all <- p3 + geom_quasirandom(data = fecundity.censored,
                             aes(x = Plugging, y = Total),
                             varwidth = TRUE, width = 0.25, cex = 2.5, alpha = 0.5, color = "gray60") +
  stat_summary(fun.data = mean_se, geom = "pointrange", color = c("#74C476", "#9970AB"), cex = 0.8) +
  xlab("plg-1 genotype") + ylab("Total Hermaphrodite Reproductive Success") + 
  theme_classic()


