library(dada2); packageVersion("dada2")
path <- "/Volumes/G-tech_RAID/Mi_seq/SOCRATES/DADA2_180912/Seqs/" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)

#シーケンスの拡張子手前の名前に合わせる
fnFs <- sort(list.files(path, pattern="_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

plotQualityProfile(fnFs[1:10])

plotQualityProfile(fnRs[1:10])


# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,210),
                     maxN=0, maxEE=c(1,3), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE

##おかしい時はここでエラー

head(out)





errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
dadaFs[[1]]
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

#Sequence read length check
dim(seqtab.nochim)
table(nchar(getSequences(seqtab.nochim)))

#オプション  ↑の結果を見て、シーケンスのターゲットレングスを絞る場合
#NIPRのシーケンスの場合だと402、426の二つにダブルのピークが出るが、中身はそろぞれ異なるのでここを下限と上限としてフィルターするのが良いのではないか？？
#テストに使った12QAの結果：/Users/uetakejun/Dropbox/Script/sequence_length_filtering_test18080

seqtab.nochim_trimed <- seqtab.nochim[,nchar(colnames(seqtab.nochim)) %in% seq(366,393)]
dim(seqtab.nochim_trimed)
table(nchar(getSequences(seqtab.nochim_trimed)))

#Negative remove: see more https://github.com/benjjneb/dada2/issues/114
rownames(seqtab.nochim_trimed)
nega.samples <- c("TO-nega", "SI-nega") # CHANGE to names of your negative controls
found.nega.samples <- colSums(seqtab.nochim_trimed[nega.samples,])>0
seqtab.nochim_trimed <- seqtab.nochim_trimed[,!found.nega.samples]
rownames(seqtab.nochim_trimed)

#Check nega ZERO and remove, naga label
rowSums(seqtab.nochim_trimed) # nega is zero
nega <- grep("TO-nega|SI-nega", rownames(seqtab.nochim_trimed), invert = TRUE)
seqtab.nochim_trimed <- seqtab.nochim_trimed[nega,]
rownames(seqtab.nochim_trimed)


#Create FASTA
asv_seqs <- colnames(seqtab.nochim_trimed)
asv_headers <- vector(dim(seqtab.nochim_trimed)[2], mode="character")
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "/Volumes/G-tech_RAID/Mi_seq/SOCRATES/DADA2_180912/rep_set.fasta") #need to add header

##Taxonomy##
taxa <- assignTaxonomy(seqtab.nochim_trimed, "/Volumes/G-tech_RAID/Silva_database/silva_nr_v132_train_set.fa_dada2.gz", multithread=TRUE)
write.table(taxa, "/Volumes/G-tech_RAID/Mi_seq/SOCRATES/DADA2_180912/dada2_taxa.txt", sep='\t', row.names=FALSE, quote=FALSE)

##Seq Table
seqtab.nochim_trimed <- t(seqtab.nochim_trimed)#transpose the table
seqtab.nochim_trimed <- cbind('#OTUID' = rownames(seqtab.nochim_trimed), seqtab.nochim_trimed)#Add '#OTUID' to the header (required by biom)
saveRDS(seqtab.nochim_trimed, "/Volumes/G-tech_RAID/Mi_seq/SOCRATES/DADA2_180912/SOCRATES180912.rds") # CHANGE ME to where you want sequence table saved
write.table(seqtab.nochim_trimed, "/Volumes/G-tech_RAID/Mi_seq/SOCRATES/DADA2_180912/dada2_seq_table.txt", sep='\t', row.names=FALSE, quote=FALSE)





##Phyloseq##


##Import meta data##
Meta <- read.delim("/Volumes/G-tech_RAID/Mi_seq/SOCRATES/DADA2_180912/R/Metadata.txt", row.names = 1)

library(phyloseq)
library(ggplot2)

##MAke phyloseq data table
SOCRATES <- phyloseq(otu_table(seqtab.nochim_trimed, taxa_are_rows=FALSE), 
               sample_data(Meta), 
               tax_table(taxa))

SOCRATES

##info
sample_names(SOCRATES)
otu_table(SOCRATES)[1:5, 1:5]
rank_names(SOCRATES)
tax_table(SOCRATES)
sample_variables(SOCRATES)

#Remove samples by category
SOCRATES_ONLY = subset_samples(SOCRATES, Ref=="SOCRATES")

##Alpha diversity 
#phyloseq
plot_richness(SOCRATES_ONLY)
plot_richness(SOCRATES_ONLY, measures=c("Chao1", "Shannon"))
Richness =plot_richness(GP, x="Type", measures=c("Chao1", "Shannon", "InvSimpson"))
Richness + geom_boxplot(color="black", alpha=0.8) + 
  geom_jitter(position=position_jitter(0.1)) +
  theme_linedraw() + 
  scale_fill_brewer(palette = "Paired") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.55, size=10), legend.position = 'none', 
                           strip.text.y = element_text(size=10, face="bold"))

install.packages("DESeq2")
library(microbiomeSeq)
print(SOCRATES_ONLY)




##NORMARIZATION!!
SOCRATES_ONLY_N <- normalise_data(SOCRATES_ONLY, norm.method = "relative", norm.meta = T)

##Alpha but not so good
micp <- plot_anova_diversity(SOCRATES_ONLY_N, method = c("richness", "simpson", "shannon"), 
    grouping_column = "Type", pValueCutoff = 0.05)
print(micp)


ord.res <- ordination(SOCRATES_ONLY_N, distance = "bray", method = "NMDS", grouping_column = "Depth", 
    pvalue.cutoff = 0.05)
p <- plot_ordination(ord.res, method = "NMDS", pvalue.cutoff = 0.05, show.pvalues = T, 
    num.signi.groups = NULL)
print(p)





##Seq Table
seqtab.nochim_trimed <- t(seqtab.nochim_trimed)#transpose the table
seqtab.nochim_trimed <- cbind('#OTUID' = rownames(seqtab.nochim_trimed), seqtab.nochim_trimed)#Add '#OTUID' to the header (required by biom)
saveRDS(seqtab.nochim_trimed, "/Volumes/G-tech_RAID/Mi_seq/SOCRATES/DADA2_180912/R/SOCRATES180912.rds") # CHANGE ME to where you want sequence table saved
write.table(seqtab.nochim_trimed, "/Volumes/G-tech_RAID/Mi_seq/SOCRATES/DADA2_180912/R/dada2_seq_table.txt", sep='\t', row.names=FALSE, quote=FALSE)


ps <- prune_samples(sample_names(ps) != "Mock", ps) # Remove mock sample
ps

