## This is an attempt, after many rounds of analysis, to come up with a better screen for differences
## introduce a control for starting conditions by dividing by the day 1 abundance
## currently, it excludes day 1 for the comparisons as an initial screen
## the pipeline then should move to the stacked.bar.2 script which is then directed at taxa of interest by site

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
library("DESeq2")

## Deseq does not use direct comparisons, it calculates curves and deviations
## changing to relative abundance not necessary

#agglomerate to genus level
physeq.genus = tax_glom(physeq.f, taxrank = "genus")
physeq.genus = subset_samples(physeq.genus, Study_day != 1)

titleD = "POST days 2-6: Infected / Control"
target.set = subset_samples(physeq.genus, Site == "POST")
diagdds = phyloseq_to_deseq2(target.set, ~ C_vs_I)

# method to deal with too many zeros, skip if not needed
geoMeans = apply(counts(diagdds, normalized = FALSE), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)

# always use this
diagdds = DESeq(diagdds, fitType="local")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(target.set)[rownames(sigtab), ], "matrix"))
head(sigtab)

dim(sigtab)

library("ggplot2")
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$phylum = factor(as.character(sigtab$phylum), levels=names(x))

# family order
x = tapply(sigtab$log2FoldChange, sigtab$genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$genus = factor(as.character(sigtab$genus), levels=names(x))
ggplot(sigtab, aes(x=genus, y=log2FoldChange, color=family)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))  + ggtitle(titleD)

write.csv(sigtab, "TS.days26.diffabund.csv")


####### this is aimed at normalizing to starting abundance, pretty difficult to deal with and NOT COMPLETE

#convert to relative abundance
physeq.g.rel = transform_sample_counts(physeq.genus, function(OTU) OTU/sum(OTU) )

#now strip out day 1 samples and merge samples
physeq.d26 = subset_samples(physeq.g.rel, Study_day != 1)
physeq.d1 = subset_samples(physeq.g.rel, Study_day == 1)

physeq.d1 = merge_samples(physeq.d1, "NewPastedVar")
