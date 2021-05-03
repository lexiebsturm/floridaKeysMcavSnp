library(plyr)

dat = read.table("fkMcav.baye_fst.txt",header=T)
head(dat)
table(dat[,"qval"]<0.1)
dat$locus = c(1:nrow(dat))
outs=which(dat[,"qval"]<0.1)

mcavGo=read.table("Mcavernosa_gene2go.tab")
names(mcavGo)=c("gene", "GO")

genes = read.table("mcav_gene_regions.tab")
names(genes) = c("chromo","start","end","gene")

# expand gene regions ± 2000 bp
genes$start = genes$start -2000
genes$end = genes$end +2000

genes = full_join(genes, mcavGo, by = "gene")

snpLoci = read.table("fkMcavNoClones.mafs.gz", header = TRUE)
snpLoci$locus = c(1:nrow(snpLoci))

snpsByGene = snpLoci %>% dplyr::select(locus, chromo, position) %>%
  full_join(., genes, by = "chromo") %>%
  filter(position>start, position<end, !is.na(GO))

snpsByGene$outlier = 0
for(i in 1:nrow(snpsByGene)){
  if(snpsByGene$locus[i] %in% outs){
    snpsByGene$outlier[i] = 1
  }else{
    snpsByGene$outlier[i] = 0
  }
}

snpsByGene$outlier

input=snpsByGene %>% select(gene, outlier)

write.csv(input, file="input.csv", row.names = FALSE)

# two columns of comma-separated values: gene id, To perform standard GO enrichment analysis based on Fisher's exact test, use binary measure (0 or 1, i.e., either sgnificant or not).

# First, press command-D on mac or ctrl-shift-H in Rstudio and navigate to the directory containing scripts and input files. Then edit, mark and execute the following bits of code, one after another.


# Edit these to match your data file names: 
input="input.csv"
#input="heats.csv" # two columns of comma-separated values: gene id, continuous measure of significance. To perform standard GO enrichment analysis based on Fisher's exact test, use binary measure (0 or 1, i.e., either sgnificant or not).
goAnnotations="Mcavernosa_gene2go.tab"
#goAnnotations="amil_defog_iso2go.tab" # two-column, tab-delimited, one line per gene, multiple GO terms separated by semicolon. If you have multiple lines per gene, use nrify_GOtable.pl prior to running this script.
goDatabase="go.obo" # download from http://www.geneontology.org/GO.downloads.ontology.shtml
goDivision="MF" # either MF, or BP, or CC
source("gomwu.functions.R")


# ------------- Calculating stats
# It might take a few minutes for MF and BP. Do not rerun it if you just want to replot the data with different cutoffs, go straight to gomwuPlot. If you change any of the numeric values below, delete the files that were generated in previos runs first.

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="/usr/bin/perl", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,   # a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           #	Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead. 
           #Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
           Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
)
# do not continue if the printout shows that no GO terms pass 10% FDR.