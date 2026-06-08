library(dplyr)
library(ggplot2)
library(ggbeeswarm)
library(RColorBrewer)
library(scales)
library(reshape)

#Read in plug size data
size <- read.table("Dropbox/Plg-1/data_files/FileS5_PlugSize.txt", sep = "\t", header = TRUE)
size$Area <- as.numeric(size$Area)

size.censor <- subset(size, Censor == 0 & is.na(Area) == FALSE)

#Turn days post-mating into a continuous variable
size.censor$Day[size.censor$Day == "W"] <- 24
size.censor$Day[size.censor$Day == "T"] <- 48
size.censor$Day[size.censor$Day == "F"] <- 72

size.censor$Day <- as.numeric(size.censor$Day)

#Examine mean plug area by strain/species and day
summarise(group_by(size.censor, Strain, Day), mean(Area), sd(Area)/length(Area), length(Strain))
summarise(group_by(size.censor, Species, Day), mean(Area), sd(Area)/length(Area), length(Species))

#Does C.elegans plug size change over time?
#Includes observer effects
model1 <- lm(Area ~ Day + Strain*Name, data = subset(size.censor, Species == "elegans"))
summary(model1)

#Does C.remanei plug size change over time?
#Includes observer effects
model2 <- lm(Area ~ Day + Strain*Name, data = subset(size.censor, Species == "remanei"))
summary(model2)

#Pull out C. remanei strain KRK1
model3 <- lm(Area ~ Day, data = subset(size.censor, Strain == "KRK1"))
summary(model3)

#Pull out C. remanei strain QX1539
#Drop 72hr outlier
model4 <- lm(Area ~ Day, data = subset(size.censor, Strain == "QX1539" & Photo != 36))
summary(model4)

#Pull out C. remanei strain QX1557
model5 <- lm(Area ~ Day, data = subset(size.censor, Strain == "QX1557"))
summary(model5)

#Is there a difference in plug size between the species over time
model6 <- lm(Area ~ Species*Day, data = size.censor)
summary(model6)

#Plot C. elegans
#Note: Remove outlier
#Strain species Day Mean    SE
#ED3017 elegans  24 14841.6 55.48937
Ce.grouped <- c();
for (s in 1:3) {
  strain <- c("ED3017", "EG4725", "JU346")
  
  z <- subset(size.censor, Strain == strain[s])
  
  for (i in 1:3) {
    day <- c(24, 48, 72)
    x   <- subset(z, Day == day[i])
    
    for (j in 1:length(levels(as.factor(x$Photo)))) {
      y   <- subset(x, Photo == levels(as.factor(x$Photo))[j])
      ave <- mean(y$Area)
      se  <- sd(y$Area) / length(y$Area)
      
      Ce.grouped <- rbind(Ce.grouped, data.frame(Strain  = as.character(y$Strain)[1],
                                                 Species = as.character(y$Species)[1],
                                                 Day     = as.numeric(y$Day)[1],
                                                 Mean    = ave,
                                                 SE      = se))
    }
  }
}

p1 <- ggplot(data = subset(Ce.grouped, Mean < 10000), (aes(x = Day, y = Mean))) +
  xlim(0, 80) +
  scale_y_continuous(breaks = pretty(Ce.grouped$Mean, n = 12)) +
  theme_classic() +
  theme(panel.grid = element_blank())

ce <- p1 +
  geom_point(data = subset(Ce.grouped, Mean < 10000), aes(x = Day, y = Mean, shape = factor(Strain)), color = "gray60", alpha = 0.35, cex = 2.5) +
  stat_summary(fun.data = mean_se, geom = "pointrange", aes(x = Day, y = Mean, shape = factor(Strain)), color = "#9970AB", cex = 0.5, position = position_dodge(2)) +
  scale_shape_manual(values = c(19, 17, 15)) +
  xlab("Hours Post-Mating") + ylab("Plug Area (um)") + 
  theme_classic()

#Plot C. remanei
Cre.grouped <- c();
for (s in 1:3) {
  strain <- c("KRK1", "QX1539", "QX1557")
  
  z <- subset(size.censor, Strain == strain[s])
  
  for (i in 1:3) {
    day <- c(24, 48, 72)
    x   <- subset(z, Day == day[i])
    
    for (j in 1:length(levels(as.factor(x$Photo)))) {
      y   <- subset(x, Photo == levels(as.factor(x$Photo))[j])
      ave <- mean(y$Area)
      se  <- sd(y$Area) / length(y$Area)
      
      Cre.grouped <- rbind(Cre.grouped, data.frame(Strain  = as.character(y$Strain)[1],
                                                   Species = as.character(y$Species)[1],
                                                   Day     = as.numeric(y$Day)[1],
                                                   Mean    = ave,
                                                   SE      = se))
    }
  }
}

p2 <- ggplot(data = Cre.grouped, (aes(x = Day, y = Mean))) +
  xlim(0, 80) +
  scale_y_continuous(breaks = pretty(Cre.grouped$Mean, n = 8)) +
  theme_classic() +
  theme(panel.grid = element_blank())

cre <- p2 +
  geom_point(data = Cre.grouped, aes(x = Day, y = Mean, shape = factor(Strain)), color = "gray60", alpha = 0.35, cex = 2.5) +
  stat_summary(fun.data = mean_se, geom = "pointrange", aes(x = Day, y = Mean, shape = factor(Strain)), color = "#9970AB", cex = 0.5, position = position_dodge(2)) +
  scale_shape_manual(values = c(19, 17, 15)) +
  xlab("Hours Post-Mating") + ylab("Plug Area (um)") + 
  theme_classic()

dev.copy2pdf(file = "~/Desktop/plugarea_cre.pdf", useDingbats = FALSE, out.type = "pdf", family = "sans", width = 4, height = 3)
