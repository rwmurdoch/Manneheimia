library("DESeq2")

# calculate geometric means prior to estimate size factors
# because this data set has many zero-values, a "zero-tolerant" gemetric mean function is added
# from : https://github.com/joey711/phyloseq/issues/387
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
testD = physeq.f
titleD = "All samples"

diagdds = phyloseq_to_deseq2(testD, ~ C_vs_I)
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(testD)[rownames(sigtab), ], "matrix"))
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
# genus order
x = tapply(sigtab$log2FoldChange, sigtab$genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$genus = factor(as.character(sigtab$genus), levels=names(x))
ggplot(sigtab, aes(x=genus, y=log2FoldChange, color=phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + ggtitle(titleD)
