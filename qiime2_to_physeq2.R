#make sure this is set to your desired working directory
setwd("/Users/robertmurdoch/Projects/Sam/phyloseq/")

library("ggplot2")
library("phyloseq")
library("ape")

otu_table = read.csv(file = "otu_table.csv", sep=",", row.names=1)
otu_table = as.matrix(otu_table)

taxonomy = read.csv("taxonomy_new.csv", sep=",", row.names=1)
taxonomy = as.matrix(taxonomy)

metadata = read.table("Sample_metadata2.txt", row.names=1, header = TRUE)

phy_tree = read_tree("phyloseq/tree.nwk")

OTU = otu_table(otu_table, taxa_are_rows = TRUE)
TAX = tax_table(taxonomy)
META = sample_data(metadata)
#	(tree	was	already	imported	as	a	phyloseq	object)

#	check	that	your	OTU	names	are	consistent	across	objects
taxa_names(TAX)
taxa_names(OTU)
taxa_names(phy_tree)

#	make	sure	files	have	the	same	sample	names	
sample_names(OTU)
sample_names(META)

#	merge	into	one	phyloseq	object
physeq	=	phyloseq(OTU,	TAX,	META,	phy_tree)
physeq
