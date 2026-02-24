####new silva tax annotiation - use rarefied phyloseq element
library(readr)
library(dada2)
library(ggplot2)
# Load your rarefied ASV table
# Assuming sequences are in the first row
asv_sequences <- colnames(pslybra@otu_table)
# Path to the reference database
silva_db <- "~/Desktop/Sequencing/SequencingJan24/silva_nr99_v138.2_toGenus_trainset.fa.gz"

# Assign taxonomy
taxonomy <- assignTaxonomy(asv_sequences, "~/Desktop/Sequencing/SequencingJan24/silva_nr99_v138.2_toGenus_trainset.fa.gz", multithread = TRUE)

tax_table_138.2 <- taxonomy

write.csv(taxonomy, "tax_table_138.2.csv")
seqtab_nochim_ra <- pslybra@otu_table
metadata.ra <- pslybra@sam_data

pslybra_132_2 <- phyloseq(otu_table(seqtab_nochim_ra, taxa_are_rows=FALSE),
               sample_data(metadata.ra),
               tax_table(tax_table_138.2),
               sample_names)


#Barplots Taxa
#but what about the less abundant taxa? ->https://github.com/joey711/phyloseq/issues/494
# get abundance in %
phyr_n <- transform_sample_counts(pslybra_132_2, function(x) x/sum(x))
# agglomerate taxa
glomr_n <- tax_glom(phyr_n, taxrank = 'Class')
# create dataframe from phyloseq object
datr_n <- psmelt(glomr_n)
# convert Genus to a character vector from a factor because R
datr_n$Class <- as.character(datr_n$Class)
# group dataframe by Phylum, calculate median rel. abundance
datr_n$Class[datr_n$Abundance < 0.01] <-"ZClasses < 1% abund."
#set color palette to accommodate the number of genera
colourCountr = length(unique(datr_n$Class))
library(RColorBrewer)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

p <- ggplot(data=datr_n, aes(x=id, y=Abundance, fill=Class))
p + geom_bar(aes(), stat="identity", position="stack") + scale_fill_manual(values=getPalette(colourCountr)) + theme(legend.position="bottom", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + guides(fill=guide_legend(nrow=5))

##### redo the enrichments
#make a table that is suitable for maaslin
library(reshape2)
T.seqtabnochim = t(pslybra_132_2@otu_table)
c =merge(taxonomy,T.seqtabnochim,by=0)

cmelt = melt(c)

cmelt$longClass = paste(cmelt$Kingdom,cmelt$Phylum,cmelt$Class,sep=';')
cmelt$longOrder = paste(cmelt$Kingdom,cmelt$Phylum,cmelt$Class,cmelt$Order, sep=';')
#make the table by class and by order
longClassTable = dcast(cmelt,variable~longClass,value.var='value',fun=sum)
longOrderTable = dcast(cmelt,variable~longOrder,value.var='value',fun=sum)
write_delim(longClassTable, "classmaas_new.tsv", delim = ' ')
write_delim(metadata.ra, "metadata_new.tsv", delim = ' ')
#maaslin
# if(!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("Maaslin2")

library(Maaslin2)
library(Matrix)
library(readr)
input_data <- read.csv("classmaas_new.tsv",row.names=1,sep=" ")
input_data <- as.data.frame(input_data)
#add metadata for GHG
input_metadata <- read.csv("metadataghg.csv", row.names = 1)
#retain only samples with respective ghg measurements
input_data = input_data[which(rownames(input_data)%in%rownames(input_metadata)),]

library(lme4)
input_metadata[1:5, ]
input_data[1:5, ]

fit_data <- Maaslin2(
  input_data, input_metadata, 'outputranewgram',
  normalization = "CLR",
  transform = "NONE",
  #analysis_method = "LM",
  fixed_effects = c('graminoid'),
  #fixed_effects = c('graminoid'),
  #random_effects = ('site'),
  standardize = FALSE)

fit_data <- Maaslin2(
  input_data, input_metadata, 'outputranewtoc',
  normalization = "CLR",
  transform = "NONE",
  #analysis_method = "LM",
  fixed_effects = c('TOC'),
  #fixed_effects = c('graminoid'),
  #random_effects = ('site'),
  standardize = FALSE)

fit_data <- Maaslin2(
  input_data, input_metadata, 'outputranewasco',
  normalization = "CLR",
  transform = "NONE",
  #analysis_method = "LM",
  fixed_effects = c('Asco'),
  #fixed_effects = c('graminoid'),
  #random_effects = ('site'),
  standardize = FALSE)

fit_data <- Maaslin2(
  input_data, input_metadata, 'outputranewbasidio',
  normalization = "CLR",
  transform = "NONE",
  #analysis_method = "LM",
  fixed_effects = c('Basidio'),
  #fixed_effects = c('graminoid'),
  #random_effects = ('site'),
  standardize = FALSE)


#barplots of significant enrichments of maaslin
library(DECIPHER)
library(mgcv)
library(ggrepel)
library(vegan)
library(ggplot2)
library(dplyr)
library(gridExtra)
theme_set(theme_bw())


df7 <- read.csv2("tocmaas1.csv", sep = ",")
df7$coef <- as.numeric(df7$coef)
normalize <- function(x) {
  2 * ((x - min(x)) / (max(x) - min(x))) - 1
}
df7$norm <- normalize(df7$coef)

df7$color <- ifelse(df7$norm > 0, "slategray1", "mistyrose1")
p7 <- ggplot(data = df7, aes(x = reorder(feature, norm), y = norm, fill = color)) +
  geom_bar(stat = "identity") +
  #labs(title = "Horizontal Bar Plot", x = "Names", y = "Gr") +
  scale_fill_identity(guide = "legend", labels = c("Negative", "Positive"), name = "Enrichment")+
  coord_flip()
p7

df8 <- read.csv2("ascomaas1.csv", sep = ",")
df8$coef <- as.numeric(df8$coef)
normalize <- function(x) {
  2 * ((x - min(x)) / (max(x) - min(x))) - 1
}
df8$norm <- normalize(df8$coef)

df8$color <- ifelse(df8$norm > 0, "slategray1", "mistyrose1")
p8 <- ggplot(data = df8, aes(x = reorder(feature, norm), y = norm, fill = color)) +
  geom_bar(stat = "identity") +
  #labs(title = "Horizontal Bar Plot", x = "Names", y = "Gr") +
  scale_fill_identity(guide = "legend", labels = c("Negative", "Positive"), name = "Enrichment")+
  coord_flip()
p8

df9 <- read.csv2("basidiomaas1.csv", sep = ",")
df9$coef <- as.numeric(df9$coef)
normalize <- function(x) {
  2 * ((x - min(x)) / (max(x) - min(x))) - 1
}
df9$norm <- normalize(df9$coef)

df9$color <- ifelse(df9$norm > 0, "slategray1", "mistyrose1")
p9 <- ggplot(data = df9, aes(x = reorder(feature, norm), y = norm, fill = color)) +
  geom_bar(stat = "identity") +
  #labs(title = "Horizontal Bar Plot", x = "Names", y = "Gr") +
  scale_fill_identity(guide = "legend", labels = c("Negative", "Positive"), name = "Enrichment")+
  coord_flip()
p9
df2 <- read.csv2("grammaas1.csv", sep = ",")
df2$coef <- as.numeric(df2$coef)
normalize <- function(x) {
  2 * ((x - min(x)) / (max(x) - min(x))) - 1
}
df2$norm <- normalize(df2$coef)

df2$color <- ifelse(df2$norm > 0, "slategray1", "mistyrose1")
p2 <- ggplot(data = df2, aes(x = reorder(feature, norm), y = norm, fill = color)) +
  geom_bar(stat = "identity") +
  #labs(title = "Horizontal Bar Plot", x = "Names", y = "Gr") +
  scale_fill_identity(guide = "legend", labels = c("Negative", "Positive"), name = "Enrichment")+
  coord_flip()
p2

grid.arrange(p2, p7, p8, p9, ncol = 2)
