### install microbiomeSeq and dependencies

#install.packages("git2r")
#install.packages("devtools")
#devtools::install_github("hadley/devtools")
source("https://bioconductor.org/biocLite.R")
biocLite()
library(devtools)  # Load the devtools package
source("https://bioconductor.org/biocLite.R")
biocLite("impute")
biocLite("preprocessCore")
biocLite("GO.db")
#install.packages("visNetwork")


#install_github("umerijaz/microbiomeSeq")  # Install the package
library(microbiomeSeq)  #load the package
library(igraph)
library(plotly)
library(visNetwork)
library(phyloseq)

## First I will try to do a network

######################################################################################

## this seems to only work when transposed, 
testset = subset_samples(physeq.genus, Site == "NL")

#careful, dont use all of these
testset = normalise_data(subset_samples(testset, Study_day!=1), norm.method = "scale")
#testset = normalise_data(testset, norm.method = "scale")
testset = taxa_level(t(testset), which_level = "genus")
co_occr <- co_occurence_network(testset, grouping_column = "C_vs_I", rhos = 0.35,
                                method = "cor", qval_threshold = 0.05, scale.vertex.size = 4, 
                                scale.edge.width = 15, plotNetwork = T, plotBetweennessEeigenvalue = F)

#rhos = 0.35
#select.condition = "I", 

require(visNetwork)

g.NL <- co_occr$net$graph
data <- toVisNetworkData(g.NL)
visNetwork(nodes = data$nodes, edges = data$edges, width = 900, height = 900) %>% visOptions(highlightNearest = TRUE, 
                                                                               nodesIdSelection = TRUE)

######################################################################################
### repeating above::

testset = subset_samples(physeq.genus, Site == "NS")

#careful, dont use all of these
testset = normalise_data(subset_samples(testset, Study_day!=1), norm.method = "scale")
#testset = normalise_data(testset, norm.method = "scale")
testset = taxa_level(t(testset), which_level = "genus")
co_occr <- co_occurence_network(testset, grouping_column = "C_vs_I", rhos = 0.35,
                                method = "cor", qval_threshold = 0.05, scale.vertex.size = 4, 
                                scale.edge.width = 15, plotNetwork = T, plotBetweennessEeigenvalue = F)


require(visNetwork)

g.NS <- co_occr$net$graph
data <- toVisNetworkData(g.NS)
visNetwork(nodes = data$nodes, edges = data$edges, width = 900, height = 900) %>% visOptions(highlightNearest = TRUE, 
                                                                                             nodesIdSelection = TRUE)

# pulling out vectors of the names of members of communities

taxa.roles <- module.roles(g.NS)
write.table(taxa.roles, "taxa.roles.NS")
taxa.roles = as.data.frame(taxa.roles)

NS.comm.1 = taxa.roles[taxa.roles$module == 6,]
NS.comm.1 = row.names(NS.comm.1)
#this removes anything that includes the string "uncultured"
NS.comm.1 = NS.comm.1[lapply(NS.comm.1,function(x) length(grep("uncultured",x,value=FALSE))) == 0]

# getting the larger data set ready

physeq.rel = transform_sample_counts(physeq.f, function(OTU) OTU/sum(OTU) )
physeq.rel = subset_samples(physeq.rel, Site != "POST")
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

## CUSTOMIZATION BEGINS HERE!!!! 
physeq.merge = subset_samples(physeq.merge, Site == "NS")
physeq.merge = tax_glom(physeq.merge, taxrank = "genus")

phylum_colors <- c(
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "black")

## start new communities here:

physeq.m.comm = subset_taxa(physeq.merge, genus %in% NS.comm.1)

plot_bar(physeq.m.comm, x = "Study_day", facet_grid=~C_vs_I, fill="family") +
  ggtitle("Relative Abundance of NS Community 6") + scale_fill_manual(values = phylum_colors)

######################################################################################
### repeating above::

testset = subset_samples(physeq.genus, Site == "PH")

#careful, dont use all of these
testset = normalise_data(subset_samples(testset, Study_day!=1), norm.method = "scale")
#testset = normalise_data(testset, norm.method = "scale")
testset = taxa_level(t(testset), which_level = "genus")
co_occr <- co_occurence_network(testset, grouping_column = "C_vs_I", rhos = 0.35,
                                method = "cor", qval_threshold = 0.05, scale.vertex.size = 4, 
                                scale.edge.width = 15, plotNetwork = T, plotBetweennessEeigenvalue = F)


require(visNetwork)

g.PH <- co_occr$net$graph
data <- toVisNetworkData(g.PH)
visNetwork(nodes = data$nodes, edges = data$edges, width = 900, height = 900) %>% visOptions(highlightNearest = TRUE, 
                                                                                             nodesIdSelection = TRUE)

######################################################################################
### repeating above::

testset = subset_samples(physeq.genus, Site == "TS")

#careful, dont use all of these
testset = normalise_data(subset_samples(testset, Study_day!=1), norm.method = "scale")
#testset = normalise_data(testset, norm.method = "scale")
testset = taxa_level(t(testset), which_level = "genus")
co_occr <- co_occurence_network(testset, grouping_column = "C_vs_I", rhos = 0.35,
                                method = "cor", qval_threshold = 0.05, scale.vertex.size = 4, 
                                scale.edge.width = 15, plotNetwork = T, plotBetweennessEeigenvalue = F)


require(visNetwork)

g.TS <- co_occr$net$graph
data <- toVisNetworkData(g.TS)
visNetwork(nodes = data$nodes, edges = data$edges, width = 900, height = 900) %>% visOptions(highlightNearest = TRUE, 
                                                                                             nodesIdSelection = TRUE)
