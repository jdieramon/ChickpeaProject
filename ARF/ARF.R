
# ------------------------------------------------------------------------------
## Genome-wide identification of the ARF Gene Family in Cicer arietinum
# ------------------------------------------------------------------------------

## -- Global Settings --
rm(list=ls())

## Load functions
source("ARFfunctions.R")

## Load Robjects for this script
load("CaARF.RData")

## Dependencies
library(dplyr)
library(ggplot2)
library(rentrez)
library(Biostrings)
library(GenomicRanges)


## Load Blast hits data 
xp = read.csv("ARFblastHit.csv") # puedo meterlo en CaARF.RData object!

# ------------------------------------------------------------------------------
### Build Table 1
# ------------------------------------------------------------------------------

## Initialize the vectors containing the features that we need for each sequence
LOC <-  vector("character") 
Chr <-  vector("character")
chr_s <-  vector("integer")
chr_e <-  vector("integer")
exon <-  vector("integer") 
AA <-  vector("integer")
mol_wt <-  vector("integer")

## Extract info feature from NCBI
for(xi in as.character(xp$XP)) {
  xpinfo <- entrez_summary(db = "protein", id = xi)
  AA = c(AA, xpinfo$slen)
  
  protein_gb <- entrez_fetch(db = "protein", id = xi, rettype = "gp")
  mol_wt = c(mol_wt, extract_mol.wt_from_xp(protein_gb)/1000) # kDa
  
  xplink <- entrez_link(dbfrom = "protein", id = xi, db = "gene")
  genesummary = entrez_summary(db = "gene", id = xplink$links[1])
  LOC = c(LOC, genesummary$name)
  Chr = c(Chr, genesummary$chromosome)
  
  # defensive programming
  # some XP may not have chrstart, chrstop or exon count
  if(length(genesummary$genomicinfo$chrstart) != 0) {
    chr_s = c(chr_s, genesummary$genomicinfo$chrstart)
    chr_e = c(chr_e, genesummary$genomicinfo$chrstop)
    exon = c(exon, genesummary$genomicinfo$exoncount)  
  }
  
  else{
    chr_s = c(chr_s, 0)
    chr_e = c(chr_e, 0)
    exon = c(exon, 0)  
    
    }
}


## If we check out the output, we will detect 2 errors. 
## Those errors have to be manually corrected.

# mol_wt: 1 missing
# table(mol_wt == 0) 
# as.character(xp$XP)[which(mol_wt==0)] 
# Manually compute the mol_wt using the Compute Mw tool (Expasy)
# mol_wt[which(mol_wt==0)] = 76888 / 1000

# exon : 0
# table(exon == 0)
# as.character(xp$XP)[which(exon == 0)]
# Manually find the exon number by NCBI search
# https://www.ncbi.nlm.nih.gov/gene/?term=101514889%5Buid%5D
# exon[which(exon==0)] = 4


## Build the data set 
t1 <- data.frame(XP = xp$XP, LOC, Chr, chr_s, chr_e, AA, mol_wt, exon)

## Add strand
t1 <- t1 %>% mutate(Strand = ifelse(chr_e > chr_s, "+", "-"))
    
## Sort data set by Chr + start coordinate
tsorted = t1 %>% arrange(Chr, chr_s)

## Subset the data set based on 1 gene model / locus
t2 = df_wo_duplicates(tsorted, 2, 1)

## Add Gene ID
t2 <- t2 %>% mutate(GeneID = paste("CaARF", seq(nrow(t2)), sep=""))

## Get info on isoforms 
isof = names(which(table(tsorted$LOC)>1))
length(isof) # 12 loci with isoforms (>1 protein)

## Add number of isoforms 
# Initialize n.isof vector
n.isof <- vector("integer", nrow(t2))  

# Fill in the body of the for loop
for (i in seq(nrow(t2))) {
    n.isof[i] = as.numeric(table(tsorted$LOC)[which(names(table(tsorted$LOC)) == as.character(t2$LOC[i]))])
}

t1 <- t2 %>% mutate(N.Isof = n.isof)

## Show the head of Table 1
t1 <- t1 %>% select(GeneID, XP, LOC, Chr, chr_s, chr_e, Strand, AA, mol_wt, exon, N.Isof)

head(t1)

### ----------------------------------------------------------------------
# write.csv(t1, file = "t1.csv", row.names=FALSE)
# t1 = read.csv("t1.csv")
### ----------------------------------------------------------------------

## Plot some general data
boxplot(t1$exon, main = "ARF Family", ylab = "Exon number")
stripchart(t1$exon, vertical = TRUE, method = "jitter", 
           add = TRUE, col = "grey70", pch = 19, cex = 0.9)

boxplot(t1$AA, main = "ARF Family", ylab = "Protein (aa)")
stripchart(t1$AA, vertical = TRUE, method = "jitter", 
           add = TRUE, col = "grey70", pch = 19, cex = 0.9)

boxplot(t1$mol_wt, main = "ARF Family", ylab = "Molecular Weight (kDa)")
stripchart(t1$mol_wt, vertical = TRUE, method = "jitter", 
           add = TRUE, col = "grey70", pch = 19, cex = 0.9)



# CaARF.RData object contains the number of HMM domains for each protein
# HMM vector

# HMM profiles vs protein length(aa)
t1hmm = t1 %>% bind_cols(hmm = HMM)
p = ggplot(t1hmm, aes(x = hmm, y = AA))

# -------------------------------------------------------------------------
## Figure Manuscript S2B
# -------------------------------------------------------------------------
set.seed(12349)
p + geom_jitter(alpha = .7, col="darkblue", width = 0.1) + 
    scale_x_continuous(breaks = c(2,3),labels = c(2,3)) + 
    labs(title = "ARF HMM profiles", 
         subtitle = "Auxin Response Factor gene family in chickpea", 
         x = "", y = "Length aa")


# -------------------------------------------------------------------------
## Figure Manuscript S2A
# -------------------------------------------------------------------------
q = ggplot(t1p, aes(x = Chr, y = AA))

q + geom_point(alpha = .8, col= "darkblue") + 
    labs(title = "Distribution of ARF proteins per Chromosome", 
         subtitle = "Auxin Response Factor gene family in chickpea", 
         x = "", y = "Length aa")


# CaARF.RData object contains the %AA Identity between CaARF proteins
# and the model species Medicago and Arabidopsis (data frame = Imodels)

# Base plotting system
with(Imodels, boxplot(Identity~Species))

# ggplot
# -------------------------------------------------------------------------
## Figure Manuscript S3
# -------------------------------------------------------------------------
p <- ggplot(Imodels, aes(Species, Identity))
p + geom_boxplot() + ylab("Amino acid Identity (%)") + 
    ylim(30, 100) + xlab("") + 
    ggtitle("Identity between CaARFs and ARFs from model species")

p + geom_boxplot() + geom_jitter(width = 0.3)


## Phylogenetic clades vs n. of exons
# CaARF.RData object contains for each protein, the clade classification 
# according to the phylogenetic analysys
# clades vector

# HMM profiles vs protein length(aa)
table(clades)

t1clades = t1 %>% bind_cols(Clade = clades)

# Base plotting system
boxplot(t1clades$exon ~ t1clades$Clade, 
        ylab = "N. exons", xlab = "ARF proteins", 
        main = "Phylogenetic clades")



## t.test
# Group 1: clade I
g1 <- t1clades %>% filter(Clade %in% c("Ia","Ib","Ic")) %>% select(exon)

# Group 2: clade II
g2 <- t1clades %>% filter(Clade %in% c("IIa","IIb")) %>% select(exon)

summary(g1)
summary(g2)

t.test(g1, g2)
str(t.test(g1, g2))
t.test(g1, g2)$p.value


# ggplot
clades = c(rep("Clade I", length(g1$exon)),
           rep("Clade II", length(g2$exon)))
n.exon = c(g1$exon, g2$exon)

tmp = data.frame(clades, n.exon)


q <- ggplot(tmp, aes(as.factor(clades), n.exon))
q + geom_boxplot(outlier.colour = "NA") +  geom_jitter(width = 0.3) +
    labs(title = "Phylogenetic clades", 
         subtitle = "ARF proteins", 
         x = "", y = "N. exons")

q + geom_boxplot(fill = c("lightgray","steelblue"), outlier.colour = "NA")+  geom_jitter(width = 0.3)
q + geom_boxplot(outlier.colour = "NA") +  
    geom_jitter(aes(colour = clades), width = 0.3)


# ------------------------------------------------------------------------------
##          MR region: Amino Acid Composition
# ------------------------------------------------------------------------------

## MUSCLE alignment
## After alignment, get the coordinates for MR domain (Table S1)
## Those coordinates have been included into the mr Robject

## load the MR coordinates (mr vector) from the CaARF.RData object
# mr

# Amino accid of interest based on literature
letters=c("Q", "S", "P", "G", "L")

# Compute the frequency of each letter for the protein set
mr_domain <- letterFrequencyProteins(as.character(t1$XP), mr, letters)

# Compute the frequency for different enrichment types
mr_domain <-  mr_domain %>% rowwise() %>% mutate(QSL = sum(Q,S,L), 
                                   SPL = sum(S,P,L), 
                                   SGL = sum(S,G,L), 
                                   SPGL = sum(S,P,G,L), 
                                   Max = max(QSL, SPL, SGL, SPGL))
# Identify the max. enrichment
#loop over rows
for(i in 1:nrow(mr_domain)) {
    #loop over columns QSL : SPGL
    #colnames(mr_domain)
    for(j in 6:9) {
        if(mr_domain[i,j] == mr_domain[i,10]) {
            #add column name (enrichment)
            mr_domain$enrichm[i] = colnames(mr_domain)[j]
        }
    }
}
mr_domain <- mr_domain %>% select(c(1:5,11))


## -----------------------------------------------------------------------------
##                    GRanges object
## -----------------------------------------------------------------------------

# Negative widths are not allowed by IRanges, so let's define the true coordinates
# and build a new data frame, which is the input for the function 'makeGRangesFromDataFrame'

tgr <- t1 %>%  mutate(Cstart = ifelse(t1$chr_s > t1$chr_e, t1$chr_e, t1$chr_s), 
                      Cend = ifelse(t1$chr_s < t1$chr_e, t1$chr_e, t1$chr_s)) %>% 
                      select(GeneID, Chr, Strand, Cstart, Cend)

tgr <- makeGRangesFromDataFrame(tgr, start.field = "Cstart", end.field = "Cend", 
                                strand.field = "Strand", ignore.strand = F, 
                                seqnames.field = "Chr", 
                                keep.extra.columns = TRUE)

# Add genome
genome(tgr) = "C. arietinum (ASM33114v1)"
seqinfo(tgr)

# Change strand of the unmapped ARF
strand(tgr[24]) = "*"


## ------------------------------------------------------------------------
##                    TANDEM DUPLICATIONS
## ------------------------------------------------------------------------

# Only Chrs. with more than 1 ARF could exhibit tandem duplication. 
names(which(table(t1$Chr) > 1))

# Subset a gr object by Chr: Check duplications in Chr1
tgr[seqnames(tgr) == "Ca1"]
chr.gr = tgr[seqnames(tgr) == "Ca1"]

getTandems(chr.gr, th = 250000)

# Check duplications in Chr2
chr.gr = tgr[seqnames(tgr) == "Ca2"]
getTandems(chr.gr, th = 250000)

# Check duplications in Chr4
chr.gr = tgr[seqnames(tgr) == "Ca4"]
getTandems(chr.gr, th = 250000)
# Check duplications in Chr6
chr.gr = tgr[seqnames(tgr) == "Ca6"]
getTandems(chr.gr, th = 250000)
# Check duplications in Chr7
chr.gr = tgr[seqnames(tgr) == "Ca7"]
getTandems(chr.gr, th = 250000)

## -----------------------------------------------------------------------------
##                   ALTERNATIVE TRANSCRIPTS
## -----------------------------------------------------------------------------

# Isoform Distribution per Locus
table(t1 %>% select(LOC, N.Isof))

t.isof = t1 %>% select(LOC, N.Isof)
table(t.isof$N.Isof == 1)
table(t.isof$N.Isof == 2) #LOCUS with 2 proteins
table(t.isof$N.Isof == 3) #LOCUS with 3 proteins
table(t.isof$N.Isof == 4) #LOCUS with 4 proteins
table(t.isof$N.Isof == 5) #LOCUS with 5 proteins
table(t.isof$N.Isof == 7) #LOCUS with 7 proteins

# Number of proteins per Chromosome
sort(table(t1$Chr), decreasing = TRUE)

# ----- Locus distribution vs Chromosome length ---------------------------------

# We may want to check if this distribution is real (there must be a biological
# reason) or just random. The appropriate test is a single sample chi-squared test
# (or goodness of fit test); this test will take the distribution of proteins by 
# chromosome and see if the distribution we have is consistent with a given 
# hypothesis, to within random sampling error. For example, could we have gotten 
# this distribution of proteins just by chance, or do some chrs really have more proteins?
# This kind of test can be done with any kind of expected frequencies, but first we 
# will compare against the hypothesis of equal expected frequencies, i.e. the 
# hypothesis that Chrs contains  ARF proteins at the same rate and we got 
# this distribution of proteins just by chance and random sampling.
# 
# # remove the locus that we could not map to any Chr
# table(t1$Chr) 
# chr_dist <- table(t1$Chr)[-8] 
# 
# chisq.test(chr_dist)
# 
# The chi-squared test indicates that the distribution of ARF proteins is likely 
# to be uniform; we cannot reject the null hypothesis (the hypothesis that proteins 
# proteins distribute all over the chromosomes at the same rate) with a high degree
# of confidence.
# 
# Can the locus distribution in Chrs be explained only as a difference of choromosome length?
# # We can get the FDC Frontier, chromosomes size from the NCBI:
# https://www.ncbi.nlm.nih.gov/genome/?term=chickpea 
# size = c(48.36, 36.63, 39.99, 49.19, 48.17, 59.46, 48.96, 16.48)
# # remove the Chr8 as we could not map any ARF into that Chr
# size <- size[-8]
# size / sum(size)
# 
# chisq.test(chr_dist, p = size/sum(size))
# 
# The p-value here is higher, so we cannot reject this simple hypothesis either. 
# There is a strong association between chromosome lengths and locus distribution patterns.



## -----------------------------------------------------------------------------
##         Reference Genes: stability values for qPCR normalization
## -----------------------------------------------------------------------------

## Load the stability values for reference genes from the CaARF.RData object: `qbaseCVdat`
qbaseCVdat

norm1 <- qbaseCVdat %>% filter(bio_rep == 1)
norm2 <- qbaseCVdat %>% filter(bio_rep == 2)
color = c(rep("gray",2),"gray20")

# par(mfrow=c(1,2))

## -----------------------------------------------------------------------------
## Figure Manuscript S1A
## -----------------------------------------------------------------------------
# setEPS()
# plotting code
norm1 %>% ggplot(aes(x=genes, y=M, fill=genes)) +
    geom_bar(stat="identity")+scale_fill_manual(values=color) + 
    theme(legend.position = "none", axis.title.x = element_blank()) + ylim(0, 2) + 
    labs(y= expression(paste("Expression Stability, ", italic("M")))) +
    geom_text(aes(label=round(CV,3)), position=position_dodge(width=1),
              vjust=-0.8) + ggtitle("Biological rep. 1")

# postscript("qBase1.eps")
# dev.off()

## -----------------------------------------------------------------------------
## Figure MS S1B
## -----------------------------------------------------------------------------
# setEPS()
# plotting code
norm2 %>% ggplot(aes(x=genes, y=M, fill=genes)) +
    geom_bar(stat="identity")+scale_fill_manual(values=color) + 
    theme(legend.position = "none", axis.title.x = element_blank()) + ylim(0, 2) + 
    labs(y= expression(paste("Expression Stability, ", italic("M")))) +
    geom_text(aes(label=round(CV,3)), position=position_dodge(width=1),
              vjust=-0.8) + ggtitle("Biological rep. 2")
# postscript("qBase2.eps")
# dev.off()


# ======================================================



#############################
# save(mr, HMM, clades, Imodels, qbaseCVdat, file = "CaARF.RData")
# load("CaARF.RData")
