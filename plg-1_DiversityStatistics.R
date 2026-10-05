library(dplyr)
library(scales)
library(ggplot2)
library(RColorBrewer)
library(ggbeeswarm)

#### Analyze Pi ####
#Read in genetic diversity data (genome-wide)
all.chrom <- read.table("Dropbox/Plg-1/data_files/FileS3_GenomeWide_Pi.txt", header = TRUE, sep = "\t")
all.chrom.isotype <- read.table("Dropbox/Plg-1/data_files/FileS_GenomeWide_Isotype_Pi.txt", header = TRUE, sep = "\t")
all.chrom.isotype.censorHW <- read.table("Dropbox/Plg-1/data_files/FileS_GenomeWide_IsotypeNotHW_Pi.txt", header = TRUE, sep = "\t")

#Read in genetic diversity data (chromosome III)
chr3 <- read.table("Dropbox/Plg-1/data_files/FileS4_Chromosome3_Pi.txt", header = TRUE, sep = "\t")
chr3.isotype <- read.table("Dropbox/Plg-1/data_files/FileS_Chromosome3_Isotype_Pi.txt", header = TRUE, sep = "\t")
chr3.isotype.censorHW <- read.table("Dropbox/Plg-1/data_files/FileS_Chromosome3_IsotypeNotHW_Pi.txt", header = TRUE, sep = "\t")

#Summarize pi genome-wide and test if the mean and distribution differs by plg-1 genotype
summarize(group_by(all.chrom, status), mean(PI), sd(PI)/length(PI), min(PI), max(PI), mean(N_VARIANTS))
wilcox.test(subset(all.chrom.isotype$PI, all.chrom.isotype$status=="plugging"), subset(all.chrom.isotype$PI, all.chrom.isotype$status=="nonplugging"))
ks.test(subset(all.chrom.isotype$PI, all.chrom.isotype$status=="plugging"), subset(all.chrom.isotype$PI, all.chrom.isotype$status=="nonplugging"))

#Summarize pi on chromosome III and test if the mean and distribution differs by plg-1 genotype
summarize(group_by(chr3, status), mean(PI), sd(PI)/length(PI), min(PI), max(PI), mean(N_VARIANTS))
wilcox.test(subset(chr3.isotype$PI, chr3.isotype$status=="plugging"), subset(chr3.isotype$PI, chr3.isotype$status=="nonplugging"))
ks.test(subset(chr3.isotype$N_VARIANTS, chr3.isotype$status=="plugging"), subset(chr3.isotype$N_VARIANTS, chr3.isotype$status=="nonplugging"))

#Plot pi across chromosome III by genotype
plot(plg.pi.chr3$BIN_START, plg.pi.chr3$PI, col = alpha("#9970AB", 0.5), pch = 19, las = 1, bty = "l",
     xlab = "Chromosome III Position (bp)", ylab = "pi")
points(non.pi.chr3$BIN_START, non.pi.chr3$PI, col = alpha("#74C476", 0.5), pch = 19)
#add line for plg-1 locus
abline(v = (8852261 + 8853006)/2)

#Plot the number of variants per 10kb window on chromosome III by genotype
p1 <- ggplot(data = chr3, (aes(x = status, y = N_VARIANTS))) +
  theme_classic() +
  theme(panel.grid = element_blank())

variants <- p1 + geom_quasirandom(data = x, aes(x = status, y = N_VARIANTS),
                                  varwidth = TRUE, width = 0.4, cex = 2, alpha = 0.3) +
  stat_summary(fun.data = mean_se, geom = "pointrange", color = c("#74C476", "#9970AB"), cex = 0.7) +
  xlab("plg-1") + ylab("Number Variants per 1kb") + 
  theme_classic()
####


#### Analyze Computationally Inferred Genotype Data ####
#Read in genotypes
geno <- read.table("Dropbox/Plg-1/data_files/FileS2_Computational_Genotype.txt", header = TRUE, sep = "\t")

table(geno$plugging_status)

#Test if the mean number of Cer1 elements differs by plg-1 genotype
wilcox.test(subset(geno$cer1_copies, geno$plugging_status=="plugging"), subset(geno$cer1_copies, geno$plugging_status=="nonplugging"))


#Color phylogeny by genotype
full_tree <- read.tree("Downloads/WI.20250625.hard-filter.min4.tree")
##Order to match tree order
plot_ord <- full_tree$tip.label
geno.plotting <- geno[, 1:3]
geno.plotting$color <- as.character(ifelse(geno.plotting$plugging_status=="nonplugging", "#74C476",
                                          ifelse(geno.plotting$plugging_status=="plugging", "#9970AB", "gray50")))
mdat_ord <- geno.plotting[match(plot_ord, geno.plotting$strain), ]
##Plot it
plot(full_tree, no.margin = TRUE, edge.width = 2, cex = 0.5, tip.color = mdat_ord$color)


#Plot Cer1 copies by plg-1 genotype
p2 <- ggplot(data = subset(geno, plugging_status != "ambiguous"), (aes(x = plugging_status, y = cer1_copies))) +
  theme_classic() +
  theme(panel.grid = element_blank())

cer <- p2 + geom_quasirandom(data = subset(geno, plugging_status != "ambiguous"), aes(x = plugging_status, y = cer1_copies),
                                  varwidth = TRUE, width = 0.4, cex = 2, alpha = 0.3) +
  stat_summary(fun.data = mean_se, geom = "pointrange", color = c("#74C476", "#9970AB"), cex = 0.7) +
  xlab("plg-1") + ylab("Cer1 copies") + 
  theme_classic()


dev.copy2pdf(file="~/Desktop/2026.05.25_Chr3pi.pdf", useDingbats = FALSE, family = "sans", width = 6, height = 4.5)
####


#### Isotype-Strain Matching ####
strain.data <- read.table("20250625_c_elegans_strain_data.csv", header = TRUE, sep = ",")
new.geno <- merge(geno, strain.data, by.x = "strain", by.y = "strain")
new.geno <- subset(new.geno, plugging_status != "ambiguous")

iso.matching <- c();
for (i in 1:length(levels(as.factor(new.geno$isotype)))) {
  x <- subset(new.geno, isotype == levels(as.factor(new.geno$isotype))[i])$plugging_status
  y <- length(unique(x)) == 1
  
  out <- data.frame(isotype = levels(as.factor(new.geno$isotype))[i],
                    matching = y,
                    status = x[1],
                    hawaiian = subset(new.geno, isotype == levels(as.factor(new.geno$isotype))[i])$pacific_island[1],
                    latitude = subset(new.geno, isotype == levels(as.factor(new.geno$isotype))[i])$latitude[1],
                    longitude = subset(new.geno, isotype == levels(as.factor(new.geno$isotype))[i])$longitude[1])
  
  iso.matching <- rbind(iso.matching, out)
}

p1 <- ggplot(iso.matching, aes(x = longitude, y = latitude)) +
  xlim(-160, 180) +
  xlab("Longitude") + ylab("Latitude") + 
  theme_classic() +
  theme(panel.grid = element_blank())

mapped <- p1 + geom_point(data = iso.matching, aes(x = longitude, y = latitude, color = factor(status)), alpha = 0.25, cex = 2.5) +
          scale_color_manual(values = c("#74C476", "#9970AB"))

####
