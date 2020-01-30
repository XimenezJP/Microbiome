library(DESeq2)
library("knitr")
library("BiocStyle")
.cran_packages <- c("ggplot2", "gridExtra")
.bioc_packages <- c("dada2", "phyloseq", "DECIPHER", "phangorn")
.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {
  install.packages(.cran_packages[!.inst])
}
.inst <- .bioc_packages %in% installed.packages()
if(any(!.inst)) {
  source("http://bioconductor.org/biocLite.R")
  biocLite(.bioc_packages[!.inst], ask = F)
}
# Load packages into session, and print package version
sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)

set.seed(100)

miseq_path <- "/data3/projects/cisbi-0136.2/analise_2/genexpet/16s"

list.files(miseq_path)

fnFs <- sort(list.files(miseq_path, pattern="_R1_001.fastq"))
fnRs <- sort(list.files(miseq_path, pattern="_R2_001.fastq"))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sampleNames <- sapply(strsplit(fnFs, "_"), `[`, 1)
# Specify the full path to the fnFs and fnRs
fnFs <- file.path(miseq_path, fnFs)
fnRs <- file.path(miseq_path, fnRs)
fnFs[1:3]           
fnRs[1:3]



plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])



filt_path <- file.path(miseq_path, "filtered") # Place filtered files in filtered/ subdirectory
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sampleNames, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sampleNames, "_R_filt.fastq.gz"))



out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(0),
                     maxN=0, maxEE=c(5), truncQ=4, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE, minLen = 50) # On Windows set multithread=FALSE
View(out)

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sampleNames
names(derepRs) <- sampleNames

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)



plotErrors(errF)

plotErrors(errR)


dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs)
seqtabAll <- makeSequenceTable(mergers)
table(nchar(getSequences(seqtabAll)))
seqtabNoC <- removeBimeraDenovo(seqtabAll)


fastaRef <- "/data3/projects/cisbi-0136.2/analise_2/microbioma/silva_nr_v132_train_set.fa.gz"
fastaRef2 <- "/data3/projects/cisbi-0136.2/analise_2/microbioma/RefSeq-RDP16S_v2_May2018.fa.gz"
fastaRef4 <- "/data3/projects/cisbi-0136.2/analise_2/microbioma/GTDB_bac-arc_ssu_r86.fa.gz"
taxTab <- assignTaxonomy(seqtabNoC, refFasta = fastaRef, multithread=TRUE)
taxTab2 <- assignTaxonomy(seqtabNoC, refFasta = fastaRef2, multithread=TRUE)
taxTab4 <- assignTaxonomy(seqtabNoC, refFasta = fastaRef4, multithread=TRUE)

taxa.print <- taxTab # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)


ps <- phyloseq(otu_table(seqtabNoC, taxa_are_rows=F), 
               tax_table(taxTab))
sequences <- Biostrings::DNAStringSet(taxa_names(ps))
names(sequences) <- taxa_names(ps)
ps <- merge_phyloseq(ps, sequences)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

rank_names(ps)
table(tax_table(ps)[, "Genus"], exclude = NULL)
ps <- subset_taxa(ps, !is.na(Genus) & !Genus %in% c("", "uncharacterized"))

plot_richness(ps, measures=c("Shannon", "Simpson"))


plot_bar(ps, fill='Genus') + 
  geom_bar(aes(color = Genus, fill=Genus), stat="identity", position="stack") +
  labs(x = "", y = "Absolute Abundance\n") +
  theme(panel.background = element_blank())
 

##top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
##ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
##ps.top20 <- prune_taxa(top20, ps.top20)
##plot_bar(ps.top20,fill="Phylum")

##executar esta parte antes de fazer o top 20 porem parece ser apenas mudar o [1:20] e plotar direto
ps_relabund <- transform_sample_counts(ps, function(x) x / sum(x))


mdf = psmelt(ps_relabund)
mdf_sub <- mdf[, c(3,9)]
mdf_sub$Genus <- as.character(mdf_sub$Genus)
teste <- aggregate(mdf_sub[,1], by = list(mdf_sub$Genus), FUN = sum)

library(dplyr)
label <- paste0(round(teste$x * 100, 2), "%")
count.data <- teste %>%
  arrange(desc(Group.1)) %>%
  mutate(lab.ypos = cumsum(x) - 0.5*x)

png(filename = "plot.png", width=15, height=9, res = 500, units = "in")
plot_bar(ps_relabund, fill="Genus") + 
  geom_bar(aes(color = Genus, fill = Genus), stat="identity", position="stack") +
  theme_bw() +
  coord_polar("y", start=0)
dev.off()
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1), legend.title =element_text ()) +
  annotate(geom = "text", y =count.data$lab.ypos, x = 1, label = label, size=3)

mdf = psmelt(ps_relabund)
mdf_sub <- mdf[, c(3,5)]
mdf_sub$Phylum <- as.character(mdf_sub$Phylum)
teste <- aggregate(mdf_sub[,1], by = list(mdf_sub$Phylum), FUN = sum)

library(dplyr)
label <- paste0(round(teste$x * 100, 2), "%")
count.data <- teste %>%
  arrange(desc(Group.1)) %>%
  mutate(lab.ypos = cumsum(x) - 0.5*x)

plot_bar(ps_relabund, fill="Phylum") + 
  geom_bar(aes(color = Phylum, fill = Phylum), stat="identity", position="stack") +
  theme_bw() +
  coord_polar("y", start=0)+
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1), legend.title =element_text ()) + 
  annotate(geom = "text", y =count.data$lab.ypos, x = 1, label = label, size=3)





plot_bar(ps_relabund, fill="Phylum") + 
  geom_bar(aes(color = Phylum, fill = Phylum), stat="identity", position="stack") +
  theme_bw() +
  coord_polar("y", start=0)+
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1), legend.title=element_text())


plot_bar(ps_relabund, fill="Genus") + 
  geom_bar(aes(color = Genus, fill = Genus), stat="identity", position="stack") +
  theme_bw() +
  coord_polar("y", start=0)+
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1), legend.title =element_text ())



















asv_seqs <- colnames(seqtabNoC)
asv_headers <- vector(dim(seqtabNoC)[2], mode="character")

for (i in 1:dim(seqtabNoC)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ASVs.fa")

# count table:
asv_tab <- t(seqtabNoC)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "ASVs_counts+.tsv", sep="\t", quote=F, col.names=NA)
View(asv_tab)
# tax table:
asv_tax <- taxTab
row.names(asv_tax) <- sub(">", "", asv_headers)
write.table(asv_tax, "ASVs_taxonomy+.tsv", sep="\t", quote=F, col.names=NA)
View(asv_tax)













count_tab <- read.table("ASVs_counts+.tsv", header=T, row.names=1,
                        check.names=F, sep="\t")[ , -c(1:4)]


tax_tab <- as.matrix(read.table("ASVs_taxonomy+.tsv", header=T,
                                row.names=1, check.names=F, sep="\t"))

count_tab_phy <- otu_table(count_tab, taxa_are_rows=T)
tax_tab_phy <- tax_table(tax_tab)

ASV_physeq <- phyloseq(count_tab_phy, tax_tab_phy)
phyla_counts_tab <- otu_table(tax_glom(ASV_physeq, taxrank="Phylum"))

plot_richness(ASV_physeq, measures=c("Chao1", "Shannon")) 
phyla_counts_tab <- otu_table(tax_glom(ASV_physeq, taxrank="Phylum")) 
phyla_tax_vec <- as.vector(tax_table(tax_glom(ASV_physeq, taxrank="Phylum"))[,2]) 
rownames(phyla_counts_tab) <- as.vector(phyla_tax_vec)
unclassified_tax_counts <- colSums(count_tab) - colSums(phyla_counts_tab)
# and we'll add this row to our phylum count table:
phyla_and_unidentified_counts_tab <- rbind(phyla_counts_tab, "Unclassified"=unclassified_tax_counts)

temp_major_taxa_counts_tab <- phyla_and_unidentified_counts_tab[!row.names(phyla_and_unidentified_counts_tab) %in% "Firmicutes", ]

# making count table broken down by class (contains classes beyond the
# Proteobacteria too at this point)
class_counts_tab <- otu_table(tax_glom(ASV_physeq, taxrank="Class")) 

# making a table that holds the phylum and class level info
class_tax_phy_tab <- tax_table(tax_glom(ASV_physeq, taxrank="Class")) 

phy_tmp_vec <- class_tax_phy_tab[,2]
class_tmp_vec <- class_tax_phy_tab[,3]
rows_tmp <- row.names(class_tax_phy_tab)
class_tax_tab <- data.frame("Phylum"=phy_tmp_vec, "Class"=class_tmp_vec, row.names = rows_tmp)

# making a vector of just the Proteobacteria classes
proteo_classes_vec <- as.vector(class_tax_tab[class_tax_tab$Phylum == "Firmicutes", "Class"])

# changing the row names like above so that they correspond to the taxonomy,
# rather than an ASV identifier
rownames(class_counts_tab) <- as.vector(class_tax_tab$Class) 

# making a table of the counts of the Proteobacterial classes
proteo_class_counts_tab <- class_counts_tab[row.names(class_counts_tab) %in% proteo_classes_vec, ] 

# there are also possibly some some sequences that were resolved to the level
# of Proteobacteria, but not any further, and therefore would be missing from
# our class table
# we can find the sum of them by subtracting the proteo class count table
# from just the Proteobacteria row from the original phylum-level count table
proteo_no_class_annotated_counts <- phyla_and_unidentified_counts_tab[row.names(phyla_and_unidentified_counts_tab) %in% "Firmicutes", ] - colSums(proteo_class_counts_tab)

# now combining the tables:
major_taxa_counts_tab <- rbind(temp_major_taxa_counts_tab, proteo_class_counts_tab, "Unresolved_Firmicutes"=proteo_no_class_annotated_counts)

# and to check we didn't miss any other sequences, we can compare the column
# sums to see if they are the same
# if "TRUE", we know nothing fell through the cracks
identical(colSums(major_taxa_counts_tab), colSums(count_tab)) 

# now we'll generate a proportions table for summarizing:
major_taxa_proportions_tab <- apply(major_taxa_counts_tab, 2, function(x) x/sum(x)*100)

# if we check the dimensions of this table at this point
dim(major_taxa_proportions_tab)
# we see there are currently 44 rows, which might be a little busy for a
# summary figure
# many of these taxa make up a very small percentage, so we're going to
# filter some out
# this is a completely arbitrary decision solely to ease visualization and
# intepretation, entirely up to your data and you
# here, we'll only keep rows (taxa) that make up greater than 5% in any
# individual sample
temp_filt_major_taxa_proportions_tab <- data.frame(major_taxa_proportions_tab[apply(major_taxa_proportions_tab, 1, max) > 5, ])
# checking how many we have that were above this threshold
dim(temp_filt_major_taxa_proportions_tab) 
# now we have 13, much more manageable for an overview figure

# though each of the filtered taxa made up less than 5% alone, together they
# may add up and should still be included in the overall summary
# so we're going to add a row called "Other" that keeps track of how much we
# filtered out (which will also keep our totals at 100%)
filtered_proportions <- colSums(major_taxa_proportions_tab) - colSums(temp_filt_major_taxa_proportions_tab)
filt_major_taxa_proportions_tab <- rbind(temp_filt_major_taxa_proportions_tab, "Other"=filtered_proportions)


# first let's make a copy of our table that's safe for manipulating
filt_major_taxa_proportions_tab_for_plot <- filt_major_taxa_proportions_tab

# and add a column of the taxa names so that it is within the table, rather
# than just as row names (this makes working with ggplot easier)
filt_major_taxa_proportions_tab_for_plot$Major_Taxa <- row.names(filt_major_taxa_proportions_tab_for_plot)

# now we'll transform the table into narrow, or long, format (also makes
# plotting easier)
library("phyloseq")
library("vegan")
library("DESeq2")
library("ggplot2")
library("dendextend")
library("tidyr")
library("viridis")
library("reshape")
filt_major_taxa_proportions_tab_for_plot.g <- gather(filt_major_taxa_proportions_tab_for_plot, Sample, Proportion, -Major_Taxa)

# take a look at the new table and compare it with the old one
head(filt_major_taxa_proportions_tab_for_plot.g)
head(filt_major_taxa_proportions_tab_for_plot)


# and now we're ready to make some summary figures with our wonderfully
# constructed table

## a good color scheme can be hard to find, i included the viridis package
## here because it's color-blind friendly and sometimes it's been really
## helpful for me, though this is not demonstrated in all of the following :/ 

# one common way to look at this is with stacked bar charts for each taxon per sample:
ggplot(filt_major_taxa_proportions_tab_for_plot.g, aes(x=Sample, y=Proportion, fill=Major_Taxa)) +
  geom_bar(width=0.6, stat="identity") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1), legend.title=element_blank()) +
  labs(x="Sample", y="% of 16S rRNA gene copies recovered", title="All samples")


######## ITS analise


library(dada2)
packageVersion("dada2")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")


path.cut <- "/data3/projects/cisbi-0136.2/analise_2/genexpet/ITS"  ## CHANGE ME to the directory containing the fastq files.
list.files(path.cut)

###FWD <- "GCATCGATGAAGAACGCAGC"  ## CHANGE ME to your forward primer sequence
###REV <- "CTTTCCTCCGCTTATTGATATGC"  ## CHANGE ME...

cutFs <- sort(list.files(path.cut, pattern = "_R1_001.fastq", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_R2_001.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

plotQualityProfile(cutFs[1:2])
plotQualityProfile(cutRs[1:2])

filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))

out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 2), 
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)  # on windows, set multithread = FALSE
head(out)


errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)


derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

table(nchar(getSequences(seqtab.nochim)))

unite.ref <- "/data3/projects/cisbi-0136.2/analise_2/genexpet/ITS/sh_general_release_dynamic_02.02.2019.fasta"  # CHANGE ME to location on your machine
taxa <- assignTaxonomy(seqtab.nochim, unite.ref, multithread = TRUE, tryRC = TRUE)

taxa.print <- taxa  # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=F), 
               tax_table(taxa))
sequences <- Biostrings::DNAStringSet(taxa_names(ps))
names(sequences) <- taxa_names(ps)
ps <- merge_phyloseq(ps, sequences)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

rank_names(ps)
table(tax_table(ps)[, "Genus"], exclude = NULL)
ps <- subset_taxa(ps, !is.na(Genus) & !Genus %in% c("", "uncharacterized"))


ps_relabund <- transform_sample_counts(ps, function(x) x / sum(x))

png(filename = "F1.png", width=15, height=9, res = 500, units = "in")
plot_bar(ps_relabund, fill="Genus") + 
  geom_bar(aes(color = Genus, fill = Genus), stat="identity", position="stack") +
  theme_bw() +
  coord_polar("y", start=0)
dev.off()

mdf = psmelt(ps_relabund)
mdf_sub <- mdf[, c(3,9)]
mdf_sub$Genus <- as.character(mdf_sub$Genus)
teste <- aggregate(mdf_sub[,1], by = list(mdf_sub$Genus), FUN = sum)

library(dplyr)
label <- paste0(round(teste$x * 100, 2), "%")
count.data <- teste %>%
  arrange(desc(Group.1)) %>%
  mutate(lab.ypos = cumsum(x) - 0.5*x)

table(tax_table(ps)[, "Species"], exclude = NULL)
ps <- subset_taxa(ps, !is.na(Species) & !Species %in% c("", "uncharacterized"))

png(filename = "F2.png", width=15, height=9, res = 500, units = "in")
plot_bar(ps_relabund, fill="Species") + 
  geom_bar(aes(color = Species, fill = Species), stat="identity", position="stack") +
  theme_bw() +
  coord_polar("y", start=0)
dev.off()  

mdf = psmelt(ps_relabund)
mdf_sub <- mdf[, c(3,10)]
mdf_sub$Species <- as.character(mdf_sub$Species)
teste <- aggregate(mdf_sub[,1], by = list(mdf_sub$Species), FUN = sum)

library(dplyr)
label <- paste0(round(teste$x * 100, 2), "%")
count.data <- teste %>%
  arrange(desc(Group.1)) %>%
  mutate(lab.ypos = cumsum(x) - 0.5*x)
