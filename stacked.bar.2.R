### new stacked bar plot script based on https://joey711.github.io/phyloseq/plot_bar-examples.html
# 1800830 script was updated to include all Sites (include POST in the analysis)


library("magrittr")
library(ggplot2)
library(vegan)
library(dplyr)
library(scales)
library(grid)
library(reshape2)
library(phyloseq)

###### Simple plots for exploration

physeq.rel = transform_sample_counts(physeq.f, function(OTU) OTU/sum(OTU) )
subset = subset_taxa(physeq.rel, family == "Pasteurellaceae")

plot_bar(subset)

plot_bar(subset, fill="genus")

plot_bar(subset, x="C_vs_I",fill="genus")

plot_bar(subset, facet_grid = ~NewPastedVar,fill="genus")


plot_bar(subset, x = "Study_day", facet_grid = ~NewPastedVar,fill="genus")

# this produces a nice chart of Haemophilus by site and infection
plot_bar(subset, x = "Study_day", facet_grid = C_vs_I~Site,fill="genus")

### this ends up with a nice plot, but still is summing the relative abundances
### http://joey711.github.io/phyloseq-demo/Restroom-Biogeography this describes, about
### halfway down, how to average the relative abundances
### it revolves around merge_data function
### this has to be merged by infection, site, and day
### the only way to accomplish this it to make a new merged metadata variable
### described here https://github.com/joey711/phyloseq/issues/293

###### main plots, plot all taxa agglomerated and by high level organization

### Data preparation:


physeq.rel = transform_sample_counts(physeq.f, function(OTU) OTU/sum(OTU) )
sample_data(physeq.rel)$Site <- gsub("POST","PO",sample_data(physeq.rel)$Site) # POST has to be shortened for this merge process
#physeq.rel = subset_samples(physeq.rel, Site != "POST") #POST SHOULD BE INCLUDED IN ANALYSES
variable1 = as.character(get_variable(physeq.rel, "Site"))
variable2 = as.character(get_variable(physeq.rel, "C_vs_I"))
variable3 = as.character(get_variable(physeq.rel, "Study_day"))
sample_data(physeq.rel)$merge <- mapply(paste0, variable1, variable2, variable3, 
                                            collapse = "_")

physeq.merge = physeq.rel
merge_samples(physeq.merge, "merge")

#metadata is lost here at the merge
physeq.merge = merge_samples(physeq.rel, "merge")

#recovering the lost metadata
variable1 = row.names(sample_data(physeq.merge))
variable2 = substring(variable1, 1,2)
variable3 = substring(variable1, 3,3)

sample_data(physeq.merge)$Site <- variable2
sample_data(physeq.merge)$C_vs_I <- variable3

                                        

# then I have to get relative abundances again

physeq.merge = transform_sample_counts(physeq.merge, function(OTU) OTU/sum(OTU) )
sample_data(physeq.merge)$Site <- gsub('PO',"POST",sample_data(physeq.merge)$Site) #replace PO with POST

## this produces top phyla/class/order etc plots, but must be customized!
physeq_phylum <- tax_glom(physeq.merge, taxrank = "family")
TopNOTUs = names(sort(taxa_sums(physeq_phylum), TRUE)[1:12]) #find the top taxa
physeq.top.phyla = prune_taxa(TopNOTUs, physeq_phylum) #filter the data set to include only the top taxa
plot_bar(physeq.top.phyla, x = "Study_day", facet_grid = C_vs_I~Site,fill="family") +
  ggtitle("Relative Abundances of Top 12 Most Abundant Families") +scale_fill_manual(values = phylum_colors)


## CUSTOM TARGETTED PLOTS BEGIN HERE!!!! 

## this first formula is for focused analysis
## you have to set which taxon you want to explore
physeq.m.subset = subset_taxa(physeq.merge, genus == "Peptoniphilus")

#consider if you want to look only at one site
physeq.m.subset = subset_samples(physeq.m.subset, Site == "POST")

plot_bar(physeq.m.subset, x = "Study_day", facet_grid = C_vs_I~Site,fill="genus") +
  ggtitle("Relative Abundances of Mannheimia in POST") +scale_fill_manual(values = phylum_colors)


## this section is for a chart with too many members, focus on more common
physeq.m.subset = tax_glom(physeq.merge, taxrank = "genus")
topOTUs = names(sort(taxa_sums(physeq.m.subset), TRUE)[1:16])
physeq.top = prune_taxa(topOTUs, order = "Pseudomonadales")

phylum_colors <- c(
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861")
  
plot_bar(physeq.top, x = "Study_day", facet_grid = C_vs_I~Site,fill="genus") +
    ggtitle("Relative Abundances of top 16 most abundant Genera") +scale_fill_manual(values = phylum_colors) +
  ylab("Relative Abundance")
