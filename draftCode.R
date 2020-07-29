#florida keys
setwd("/Users/student/Documents/GitHub/floridaKeysMcavSnp")

if (!require("pacman")) install.packages("pacman")
pacman::p_load("adegenet", "dendextend", "flextable", "ggdendro", "hierfstat", "Imap", "kableExtra", "paletteer", "patchwork", "officer", "poppr", "RColorBrewer", "reshape2", "StAMPP", "tidyverse", "vcfR", "vegan", "WGCNA")

#-------------WITH BAD BC'S, RE'S, TECH REPS, AND CLONES
# clustering / PCoA based on identity by state (IBS) based on single read resampling
# (for low and/or uneven coverage)
bams=read.table("bams")[,1] # list of bam files
goods=c(1:length(bams))

ma = as.matrix(read.table("noSamplesConcat.ibsMat"))
imnames(ma)=list(bams[goods],bams[goods])
hc=hclust(as.dist(ma),"ave")
plot(hc,cex=0.5)  # this shows how similar clones are

ma=ma[goods,goods]
dimnames(ma)=list(bams[goods],bams[goods])
hc=hclust(as.dist(ma),"ave")
plot(hc,cex=0.75) # without clones

#-------------WITH BAD BC'S AND RE'S CONCAT BUT STILL WITH CLONES AND TECH REPS
# clustering / PCoA based on identity by state (IBS) based on single read resampling
# (for low and/or uneven coverage)
bams=read.table("bams")[,1] # list of bam files
goods=c(1:length(bams))

ma = as.matrix(read.table("fkMcavClones.ibsMat"))
dimnames(ma)=list(bams[goods],bams[goods])
hc=hclust(as.dist(ma),"ave")
plot(hc,cex=0.5)  # this shows how similar clones are

#-------------ROUGH DENDROGRAM NO CLONES
bamsNoClones=read.table("bamsNoClones")[,1] # list of bam files
goodsNoClones=c(1:length(bamsNoClones))

ma = as.matrix(read.table("fkMcavNoClones.ibsMat"))
dimnames(ma)=list(bamsNoClones[goodsNoClones],bamsNoClones[goodsNoClones])
hc=hclust(as.dist(ma),"ave")
plot(hc,cex=0.5)  # this shows how similar clones are

#------------------Dendrogram with Clones----------
cloneBams = read.table("sampleList")[,1] # list of bam files
cloneMa = as.matrix(read.table("fkMcavClones.ibsMat")) # reads in IBS matrix produced by ANGSD
dimnames(cloneMa) = list(cloneBams,cloneBams)
clonesHc = hclust(as.dist(cloneMa),"ave")

cloneMeta = read.csv("inds2pops.csv") # list of bams files and their populations
clonePops = cloneMeta$pop
cloneDepth = cloneMeta$depth

cloneDend = cloneMa %>% as.dist() %>% hclust(.,"ave") %>% as.dendrogram()
cloneDData = cloneDend %>% dendro_data()

# Making the branches hang shorter so we can easily see clonal groups
cloneDData$segments$yend2 = cloneDData$segments$yend
for(i in 1:nrow(cloneDData$segments)) {
  if (cloneDData$segments$yend2[i] == 0) {
    cloneDData$segments$yend2[i] = (cloneDData$segments$y[i] - 0.03)}}

cloneDendPoints = cloneDData$labels
cloneDendPoints$pop = clonePops[order.dendrogram(cloneDend)]
cloneDendPoints$depth=cloneDepth[order.dendrogram(cloneDend)]
rownames(cloneDendPoints) = cloneDendPoints$label

# Making points at the leaves to place symbols for populations
point = as.vector(NA)
for(i in 1:nrow(cloneDData$segments)) {
  if (cloneDData$segments$yend[i] == 0) {
    point[i] = cloneDData$segments$y[i] - 0.03
  } else {
    point[i] = NA}}

cloneDendPoints$y = point[!is.na(point)]

techReps = c("009-1", "009-2", "009-3", "159-1", "159-2", "159-3", "171-1", "171-2", "171-3")

cloneDendA = ggplot() +
  geom_segment(data = segment(cloneDData), aes(x = x, y = y, xend = xend, yend = yend2), size = 0.5) +
  geom_point(data = cloneDendPoints, aes(x = x, y = y, fill = pop, shape = depth), size = 4, stroke = 0.25) +
  scale_fill_brewer(palette = "Dark2", name = "Population") +
  scale_shape_manual(values=c(25, 24))+
  geom_hline(yintercept = 0.1, color = "red", lty = 5, size = 0.75) + # creating a dashed line to indicate a clonal distance threshold
  geom_text(data = subset(cloneDendPoints, subset = label %in% techReps), aes(x = x, y = (y - .015), label = label), angle = 90) + # spacing technical replicates further from leaf
  geom_text(data = subset(cloneDendPoints, subset = !label %in% techReps), aes(x = x, y = (y - .010), label = label), angle = 90) +
  geom_text(data = subset(cloneDendPoints, subset = label %in% c("218", "223", "221", "226")), aes(x = (x + .5), y = (y - .018), label = "*"), angle = 90, size = 6) + # labeling natural clones with asterisks
  #coord_cartesian(xlim = c(3, 84)) +
  labs(y = "Genetic distance (1 - IBS)") +
  guides(fill = guide_legend(override.aes = list(shape = 21)))+
  theme_classic()

cloneDend = cloneDendA + theme(
  axis.title.x = element_blank(),
  axis.text.x = element_blank(),
  axis.line.x = element_blank(),
  axis.ticks.x = element_blank(),
  axis.title.y = element_text(size = 12, color = "black", angle = 90),
  axis.text.y = element_text(size = 10, color = "black"),
  axis.line.y = element_line(),
  axis.ticks.y = element_line(),
  panel.grid = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  plot.background = element_blank(),
  legend.key = element_blank(),
  legend.title = element_text(size = 12),
  legend.text = element_text(size = 10),
  legend.position = "left")

cloneDend

ggsave("cloneDend.tiff", plot = cloneDend, height = 4.75, width = 30, units = "in", dpi = 300)

#----------------------Dendogram No Clones
#plotPops = c("TER-North", "TER-South", "Lower Keys", "Upper Keys")
bamsNoClones = read.table("bamsNoClones")[,1] # list of bam file
snpMa = as.matrix(read.table("fkMcavNoClones.ibsMat"))
snpI2P = read.csv("inds2popsNoClones.csv") # 2-column tab-delimited table of individual assignments to populations; must be in the same order as samples in the bam list or vcf file.
row.names(snpI2P) = snpI2P[,1]

snpDend = snpMa %>% scale %>% dist %>%
  hclust %>% as.dendrogram

snpDData = dendro_data(snpDend)
snpDendPoints = snpDData$labels
snpDendPoints$site = snpI2P[,2][order.dendrogram(snpDend)]
snpDendPoints$depth = snpI2P[,3][order.dendrogram(snpDend)]

snpDendA = ggplot() + 
  geom_segment(data = segment(snpDData), aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_point(data = snpDendPoints, aes(x = x, y = y, fill = site, shape = depth), size = 5) +
  scale_shape_manual(values=c(25, 24))+
  scale_fill_brewer(palette = "Dark2", name = "Population") +
  ggtitle("   SNP") +
  guides(fill = guide_legend(override.aes = list(shape = 21)))+
  theme_dendro()

snpDend = snpDendA + theme(
  legend.key = element_blank(),
  legend.title = element_text(size = 14),
  legend.text = element_text(size = 12))

snpDend

ggsave("snpDend.tiff", plot = snpDend, height = 4.75, width = 35, units = "in", dpi = 300)

#####-----------Read in vcf file
fkVcf = read.vcfR("fkMcavNoClonesRenamed.vcf.gz")
fkGenlight = vcfR2genlight(fkVcf, n.cores = 2) # Converts the vcf file into a file format that poppr uses the "genlight" format
locNames(fkGenlight) = paste(fkVcf@fix[,1],fkVcf@fix[,2],sep="_") 
popData = read.csv("inds2popsNoClones.csv") # Reads in population data for each sample
strata(fkGenlight) = data.frame(popData)
setPop(fkGenlight) = ~pop

#####-----------PCA

snpPca1 = glPca(fkGenlight, nf = 214, n.cores = 1)

snpPcaValues = as.data.frame(snpPca1$scores)
snpPca = data.frame(row.names = rownames(snpPcaValues), sample = rownames(snpPcaValues),
                    site = popData$pop, depth = popData$depth, siteDepth = popData$popsite, PC1 = snpPcaValues[,1],
                    PC2 = snpPcaValues[,2])
head(snpPca)

# percentage of explained variance by axes
snpPcVar = round(snpPca1$eig/sum(snpPca1$eig)*100, 1)

# merge centroid locations into ggplot dataframe
snpPca = merge(snpPca, aggregate(cbind(mean.x=PC1,mean.y=PC2)~siteDepth,snpPca, mean), by="siteDepth")

# SNP PCA biplot
plotPops = c("TER-North", "TER-South", "Lower Keys", "Upper Keys")

fkSnpPcaPlotA = ggplot(snpPca, aes(x = PC1, y = PC2, color = site, fill = site, shape = depth)) +
  geom_hline(yintercept = 0, color = "gray90", size = 0.5) +
  geom_vline(xintercept = 0, color = "gray90", size = 0.5) +
  stat_ellipse(data = subset(snpPca, type = "t", geom = "polygon", alpha = 0.1)) + #ellipse
  geom_point(aes(x = PC1, y = PC2, shape = depth), size = 3) + #individual's points indicated by circles
  scale_shape_manual(values = c(25,24), name = "Depth") +
  geom_point(aes(x = mean.x, y = mean.y, shape = depth), size = 5, color = "black") + #population centroids indicated by triangles
  #scale_shape_manual(values = c(25,24), name = "Depth") +
  scale_fill_manual(values = brewer.pal(name = "Dark2", n = 8), name = " Population", labels = plotPops) +
  scale_color_manual(values = brewer.pal(name = "Dark2", n = 8), name = " Population", labels = plotPops) +
  xlab(paste ("PC 1 (", snpPcVar[1],"%)", sep = "")) + #Prints percent variation explained by first axis
  ylab(paste ("PC 2 (", snpPcVar[2],"%)", sep = "")) + #Prints percent variation explained by second axis
  guides(fill = guide_legend(override.aes = list(shape = 22, size = 4, color = NA), order = 1))+
  ggtitle("SNP") +
  theme_bw()

fkSnpPcaPlot = fkSnpPcaPlotA +
  theme(axis.title.x = element_text(color = "black", size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.title.y = element_text(color = "black", size = 10),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        legend.position = "left",
        panel.border = element_rect(color = "black", size = 1.2),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())



fkSnpPcaPlot

#####-----------SNP PCoA BASED ON NEI'S

fkVcf = read.vcfR("fkMcavNoClonesRenamed.vcf.gz")
fkGenlightPopDepth = vcfR2genlight(fkVcf, n.cores = 2) # Converts the vcf file into a file format that poppr uses the "genlight" format
locNames(fkGenlightPopDepth) = paste(fkVcf@fix[,1],fkVcf@fix[,2],sep="_") 
popData = read.csv("inds2popsNoClonesPopDepth.csv") # Reads in population data for each sample
strata(fkGenlightPopDepth) = data.frame(popData)
setPop(fkGenlightPopDepth) = ~popdepth

#First calculate Nei's genetic distance using the program STAMPP
fkSnpDpop = stamppNeisD(fkGenlightPopDepth, pop = TRUE) # Nei's 1972 distance between pops 
#This mimics the by population PCoA of GenAlEx 

# This is the actual PCoA step; you can use any distance matrix here (euclidean = PCA), 
# We are using Nei's genetic distance
fkMds = cmdscale(fkSnpDpop, eig = TRUE, x.ret = TRUE)

# Determine percent variation captured on each axis 
# Calculate the eigenvalues so later we can figure out % variation shown on each Principal Coordinate
fkSnpPcoaVar = round(fkMds$eig/sum(fkMds$eig)*100, 1)
fkSnpPcoaVar

# Format data to plot
fkSnpPcoaValues = fkMds$points

fkSnpPcoaValues =as.data.frame(fkSnpPcoaValues, sample = rownames(fkSnpPcoaValues))
fkSnpPcoa = data.frame(site = factor(c("TER-North", "TER-North", "TER-South", "TER-South", "Lower Keys", "Lower Keys", "Upper Keys", "Upper Keys")),
                         depth = factor(c("Mesophotic","Shallow", "Mesophotic","Shallow", "Mesophotic","Shallow", "Mesophotic","Shallow"), 
                                        levels=c("Shallow","Mesophotic")), 
                         PCo1 = fkSnpPcoaValues[,1], PCo2 = fkSnpPcoaValues[,2])

# SNP PCoA biplot
fkSnpPcoaPlotA = ggplot(fkSnpPcoa, aes(x = PCo1, y = PCo2, fill = site, shape = depth)) +
  geom_hline(yintercept = 0, color = "gray90", size = 0.5) +
  geom_vline(xintercept = 0, color = "gray90", size = 0.5) + 
  geom_point(size = 4) +
  scale_y_reverse() +
  scale_shape_manual(values = c(24,25), name = "Depth") +
  scale_fill_manual(values = brewer.pal(name = "Dark2", n = 8), name = "Population", labels = plotPops) +
  xlab(paste ("PCo 1 (", fkSnpPcoaVar[1],"%)", sep = "")) +
  ylab(paste ("PCo 2 (", fkSnpPcoaVar[2],"%)", sep = "")) +
  guides(fill = guide_legend(override.aes = list(shape = 22, size = 4, color = NA), order = 1))+
  ggtitle("SNP") +
  theme_bw()

fkSnpPcoaPlot = fkSnpPcoaPlotA + 
  theme(axis.title.x = element_text(color = "black", size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.title.y = element_text(color = "black", size = 10),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        legend.position = "left",
        legend.title = element_text(color = "black", size = 10),
        legend.text = element_text(color = "black", size = 10),
        legend.key = element_blank(),
        panel.border = element_rect(color = "black", size = 1.2),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

fkSnpPcoaPlot

######-------------MISHA'S WAY

########-----PCoA with IBS
# This is the actual PCoA step; you can use any distance matrix here (euclidean = PCA), 
# We are using Nei's genetic distance
#plotPops = c("TER-North", "TER-South", "Lower Keys", "Upper Keys")
snpMa = as.matrix(read.table("fkMcavNoClones.ibsMat"))
fkMds = cmdscale(snpMa, eig = TRUE, x.ret = TRUE)

# Determine percent variation captured on each axis 
# Calculate the eigenvalues so later we can figure out % variation shown on each Principal Coordinate
fkSnpPcoaVar = round(fkMds$eig/sum(fkMds$eig)*100, 1)
fkSnpPcoaVar

# Format data to plot
fkSnpPcoaValues = fkMds$points
fkSnpPcoaValues

snpI2P = read.csv("inds2popsNoClones.csv") # 2-column tab-delimited table of individual assignments to populations; must be in the same order as samples in the bam list or vcf file.
row.names(snpI2P) = snpI2P[,1]
fkSnpPcoaValues=cbind(snpI2P, fkSnpPcoaValues)
fkSnpPcoaValues =as.data.frame(fkSnpPcoaValues, sample = rownames(fkSnpPcoaValues))
colnames(fkSnpPcoaValues)[5] <- "PCo1"
colnames(fkSnpPcoaValues)[6] <- "PCo2"
fkSnpPcoaValues

snpPCoA = merge(fkSnpPcoaValues, aggregate(cbind(mean.x=PCo1,mean.y=PCo2)~popsite, fkSnpPcoaValues, mean), by="popsite")

# SNP PCoA biplot
#plotPops = c("TER-North", "TER-South", "Lower Keys", "Upper Keys")

fkSnpPcoaPlotA = ggplot(snpPCoA, aes(x = PCo1, y = PCo2, color = pop, fill = pop, shape = depth, linetype = depth)) +
  geom_hline(yintercept = 0, color = "gray90", size = 0.5) +
  geom_vline(xintercept = 0, color = "gray90", size = 0.5) +
  stat_ellipse(data = subset(snpPCoA, type = "t", geom = "polygon", alpha = 0.1)) + #ellipse
  scale_linetype_manual(values=c(1,2))+
  geom_point(aes(x = PCo1, y = PCo2, shape = depth), size = 3, alpha = 0.3) + #individual's points indicated by circles
  scale_shape_manual(values = c(25,24), name = "Depth") +
  geom_point(aes(x = mean.x, y = mean.y, shape = depth), size = 5, color = "black") + #population centroids indicated by triangles
  #scale_shape_manual(values = c(25,24), name = "Depth") +
  scale_fill_manual(values = brewer.pal(name = "Dark2", n = 8), name = " Population") +
  scale_color_manual(values = brewer.pal(name = "Dark2", n = 8), name = " Population") +
  xlab(paste ("PCo 1 (", fkSnpPcoaVar[1],"%)", sep = "")) + #Prints percent variation explained by first axis
  ylab(paste ("PCo 2 (", fkSnpPcoaVar[2],"%)", sep = "")) + #Prints percent variation explained by second axis
  guides(fill = guide_legend(override.aes = list(shape = 22, size = 4, color = NA), order = 1))+
  theme_bw()

fkSnpPcoaPlot = fkSnpPcoaPlotA +
  theme(axis.title.x = element_text(color = "black", size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.title.y = element_text(color = "black", size = 10),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        legend.position = "left",
        panel.border = element_rect(color = "black", size = 1.2),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())



fkSnpPcoaPlot

ggsave("fkSnpPcoaPlot.tiff", plot = fkSnpPcoaPlot, height = 5, width = 7, units = "in", dpi = 300)

#########----------ADMIXTURE K4

df <- read.csv("fkMcavNoClones_k4.csv")
df$sample <- factor(df$sample, levels= df$sample[order(df$cluster4)])
#df$sample <- factor(df$sample, levels= df$sample[order(df$cluster3)])
#df$sample <- factor(df$sample, levels= df$sample[order(df$cluster2)])

mdat = melt(df, id.vars=c("sample", "pop"), variable.name="Ancestry", value.name="Fraction")

p = ggplot(mdat, aes(x=sample, y=Fraction, fill=Ancestry)) +
  geom_bar(stat="identity", position="stack") +
  facet_grid(. ~ pop, drop=TRUE, space="free", scales="free")

#col2 = c("turquoise", "blue", "green", "purple")

#names(col2) = levels(mdat$Ancestry)

p2 = ggplot(mdat, aes(x=sample, y=Fraction, fill=Ancestry, order=sample)) +
  geom_bar(stat="identity", position="fill", width=1, colour="grey25") +
  facet_grid(~fct_inorder(pop), scales = "free", switch = "x", space = "free") +
  labs(x = "Population", y = "Ancestry") +
  ggtitle("K4 NGSAdmixture Plot") +
  theme(plot.title = element_text(hjust = 0.5),
  panel.grid=element_blank(),
  panel.background=element_rect(fill=NA, colour="grey25"),
  panel.spacing.x=grid:::unit(0, "lines"),
  panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"),
  #axis.text.x=element_text(size=12, angle=90)
  axis.text.x = element_blank(),
  axis.ticks.x=element_blank(),
  strip.background=element_blank(),
  strip.text=element_text(size=12, angle=90),
  legend.key=element_blank(),
  legend.position = "none",
  legend.title = element_blank()) +
  scale_x_discrete(expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0)) +
  scale_fill_manual(values = brewer.pal(name = "YlGnBu", n = 4), name = "Cluster") +
  guides(fill=guide_legend(override.aes=list(colour=NULL)))

p2

ggsave("admixturePlot.tiff", plot = p2, width = 30, height = 15, units = "cm", dpi = 300)
#########----------ADMIXTURE K3

df <- read.csv("fkMcavNoClones_k3.csv")
df$sample <- factor(df$sample, levels= df$sample[order(df$cluster3)])
df$sample <- factor(df$sample, levels= df$sample[order(df$cluster2)])

mdat = melt(df, id.vars=c("sample", "pop"), variable.name="Ancestry", value.name="Fraction")

p = ggplot(mdat, aes(x=sample, y=Fraction, fill=Ancestry)) +
  geom_bar(stat="identity", position="stack") +
  facet_grid(. ~ pop, drop=TRUE, space="free", scales="free")

#col2 = c("turquoise", "blue", "green", "purple")

#names(col2) = levels(mdat$Ancestry)

p2 = ggplot(mdat, aes(x=sample, y=Fraction, fill=Ancestry, order=sample)) +
  geom_bar(stat="identity", position="fill", width=1, colour="grey25") +
  facet_grid(~fct_inorder(pop), scales = "free", switch = "x", space = "free") +
  labs(x = "Population", y = "Ancestry") +
  ggtitle("K4 NGSAdmixture Plot") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.grid=element_blank()) +
  theme(panel.background=element_rect(fill=NA, colour="grey25")) +
  theme(panel.spacing.x=grid:::unit(0, "lines")) +
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid")) +
  theme(axis.text.x=element_text(size=12, angle=90)) +
  theme(axis.ticks.x=element_blank()) +
  theme(strip.background=element_blank()) +
  theme(strip.text=element_text(size=12, angle=90)) +
  theme(legend.key=element_blank()) +
  theme(legend.title = element_blank()) +
  scale_x_discrete(expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0)) +
  scale_fill_manual(values = brewer.pal(name = "PRGn", n = 3), name = "Cluster") +
  guides(fill=guide_legend(override.aes=list(colour=NULL)))

p2

##########--------------------AMOVA
fkVcf = read.vcfR("fkMcavNoClonesRenamed.vcf.gz")
fkGenlightPopDepth = vcfR2genlight(fkVcf, n.cores = 2) # Converts the vcf file into a file format that poppr uses the "genlight" format
locNames(fkGenlightPopDepth) = paste(fkVcf@fix[,1],fkVcf@fix[,2],sep="_") 
popData = read.csv("inds2popsNoClonesPopDepth.csv") # Reads in population data for each sample
strata(fkGenlightPopDepth) = data.frame(popData)
setPop(fkGenlightPopDepth) = ~popdepth
amova <- poppr.amova(fkGenlightPopDepth, ~popdepth) #Runs AMOVA
amova
set.seed(1999)
amovasignif   <- randtest(amova, nrepet = 99) #Calculates significance levels of the AMOVA with 999 permutations
amovasignif

##########--------------------Pairwise Fst
fkVcf = read.vcfR("fkMcavNoClonesRenamed.vcf.gz")
fkGenlightPopDepth = vcfR2genlight(fkVcf, n.cores = 2) # Converts the vcf file into a file format that poppr uses the "genlight" format
locNames(fkGenlightPopDepth) = paste(fkVcf@fix[,1],fkVcf@fix[,2],sep="_") 
popData = read.csv("inds2popsNoClonesPopDepth.csv") # Reads in population data for each sample
strata(fkGenlightPopDepth) = data.frame(popData)
setPop(fkGenlightPopDepth) = ~popdepth

set.seed(694)
fk.fst<-stamppFst(fkGenlightPopDepth, nboots = 99, percent = 95, nclusters = 4) #99 permutations
fk.fst$Fsts
fk.fst$Pvalues

pops = c("TER-North-Mesophotic", "TER-North-Shallow", "TER-South-Mesophotic", "TER-South-Shallow", "Lower Keys-Shallow", "Lower Keys-Mesophotic", "Upper Keys-Mesophotic", "Upper Keys-Shallow")

# reads in fst matrix
snpFstMa = as.data.frame(fk.fst$Fsts)
names(snpFstMa) = pops
row.names(snpFstMa) = pops

snpFstMa$Pop = factor(row.names(snpFstMa), levels = unique(pops))
levels(snpFstMa$Pop)

snpQMa = data.frame(fk.fst$Pvalues)
names(snpQMa)= pops
row.names(snpQMa) = pops
snpQMa$Pop = factor(row.names(snpQMa), levels = unique(pops))

snpFst = melt(snpFstMa, id.vars = "Pop", value.name = "Fst", variable.name = "Pop2", na.rm = TRUE)
snpFst = snpFst[snpFst$Pop != snpFst$Pop2,]
snpFst$Fst = round(snpFst$Fst, 3)
snpFst = snpFst %>% mutate(Fst = replace(Fst, Fst < 0, 0))
head(snpFst)

snpQ = melt(snpQMa, id.vars = "Pop", value.name = "Pval", variable.name = "Pop2", na.rm = TRUE)
snpQ = snpQ[snpQ$Pop != snpQ$Pop2,]
snpQ$Qval = p.adjust(snpQ$Pval, method = "BH")
head(snpQ)

snpHeatmapA = ggplot(data = snpFst, aes(Pop, Pop2, fill = Fst))+ 
  geom_tile(color = "white")+ 
  scale_fill_gradient2(low = "white", high = "red", midpoint = 0, limit = c(0, 0.1), 
                       space = "Lab", name = expression(paste(italic("F")[ST])))+ 
  geom_text(data = snpFst, aes(Pop, Pop2, label = Fst), color = ifelse(snpQ$Qval <= 0.05,"black", "darkgrey"), size = ifelse(snpQ$Qval < 0.05, 6, 5)) +
  guides(fill=guide_colorbar(barwidth = 1, barheight = 12, title.position = "top", title.hjust = 0.5))+                     
  scale_y_discrete(position = "right")+
  scale_x_discrete(labels = str_wrap(snpFst$Pop, width = 6)) +
  #ggtitle("   SNP") +
  theme_minimal()

snpHeatmap = snpHeatmapA + theme(
  axis.text.x = element_text(vjust = 1, size = 16, hjust = 0.5, color = "black"), 
  axis.text.y = element_text(size = 16, color = "black"),
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.grid.major = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.ticks = element_blank(),
  legend.justification = c(1, 0),
  legend.position = "right",
  legend.direction = "vertical",
  legend.title = element_text(size = 16),
  legend.text = element_text(size = 14),
  plot.title = element_text(size = 16)
)

snpHeatmap

ggsave("snpHeatMap.tiff", plot = snpHeatmap, width = 34, height = 15, units = "cm", dpi = 300)
