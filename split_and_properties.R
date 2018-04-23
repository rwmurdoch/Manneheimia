library("data.table")

#filter the low read samples
physeq.f = prune_samples(sample_sums(physeq) > 1200, physeq)

#generate read depth plots
sdt = data.table(as(sample_data(physeq.f), "data.frame"),
                 TotalReads = sample_sums(physeq.f), keep.rownames = TRUE)
pSeqDepth = ggplot(sdt, aes(TotalReads)) + geom_histogram() + ggtitle("Sequencing Depth")

overall.depth = pSeqDepth

by.sample.depth = pSeqDepth + facet_wrap(~Site) +ggtitle("Sequencing Depth by Site")

p <- plot_richness(physeq.NL, x= "C_vs_I", "Study_day")
(div.NL <- p + geom_boxplot(data = p$data, aes(y = value, color = NULL), 
                            alpha = 0.1))

#make a normalized dataset from the filtered set
## DONT DO THIS, HIHGLY DEBATABLE ##
#physeq.norm = transform_sample_counts(physeq.f, function(x) x / sum(x))

#subset to only samples from each site normalized
#physeq.n.NL = subset_samples(physeq.norm, Site=="NL")
#physeq.n.NS = subset_samples(physeq.norm, Site=="NS")
#physeq.n.PH = subset_samples(physeq.norm, Site=="PH")
#physeq.n.TS = subset_samples(physeq.norm, Site=="TS")
#physeq.n.POST = subset_samples(physeq.norm, Site=="POST")

#subset to only samples from each site integers
physeq.NL = subset_samples(physeq.f, Site=="NL")
physeq.NS = subset_samples(physeq.f, Site=="NS")
physeq.PH = subset_samples(physeq.f, Site=="PH")
physeq.TS = subset_samples(physeq.f, Site=="TS")
physeq.POST = subset_samples(physeq.f, Site=="POST")

#plot overall diversity as a function of site
p <- plot_richness(physeq.f, x= "Site", "C_vs_I")
(div.all.site <- p + geom_boxplot(data = p$data, aes(y = value, color = NULL), 
                            alpha = 0.1)) + ggtitle("Overall Alpha-diversity")

#plot overall diversity as a function of infection status
p <- plot_richness(physeq.f, x= "C_vs_I", "Site")
(div.all.inf <- p + geom_boxplot(data = p$data, aes(y = value, color = NULL), 
                                  alpha = 0.1)) + ggtitle("Overall Alpha-diversity")



#plot diversity for each site depending on infection status
p <- plot_richness(physeq.NL, x= "C_vs_I", "Study_day")
(div.NL <- p + geom_boxplot(data = p$data, aes(y = value, color = NULL), 
                       alpha = 0.1)) + ggtitle("Alpha-diversity: NL")

p <- plot_richness(physeq.NS, x= "C_vs_I", "Study_day")
(div.NS <- p + geom_boxplot(data = p$data, aes(y = value, color = NULL), 
                       alpha = 0.1))+ ggtitle("Alpha-diversity: NS")

p <- plot_richness(physeq.PH, x= "C_vs_I", "Study_day")
(div.PH <- p + geom_boxplot(data = p$data, aes(y = value, color = NULL), 
                       alpha = 0.1))+ ggtitle("Alpha-diversity: PH")

p <- plot_richness(physeq.TS, x= "C_vs_I", "Study_day")
(div.TS <- p + geom_boxplot(data = p$data, aes(y = value, color = NULL), 
                       alpha = 0.1))+ ggtitle("Alpha-diversity: TS")

p <- plot_richness(physeq.POST, x= "C_vs_I", "Study_day")
(div.POST <- p + geom_boxplot(data = p$data, aes(y = value, color = NULL), 
                          alpha = 0.1))+ ggtitle("Alpha-diversity: POST")

# permanova, test if sites have different centroids
dist.all.wuni <- phyloseq::distance(physeq.f, method = "wunifrac")

## make a data frame from the sample_data
sampledf <- data.frame(sample_data(physeq.f))

## Adonis test
#library("vegan")
adonis(dist.all.wuni ~ Site, data = sampledf)

# permanova, test C vs I at different sites
dist.all.wuni.NS <- phyloseq::distance(physeq.NS, method = "wunifrac")
dist.all.wuni.NL <- phyloseq::distance(physeq.NL, method = "wunifrac")
dist.all.wuni.PH <- phyloseq::distance(physeq.PH, method = "wunifrac")
dist.all.wuni.TS <- phyloseq::distance(physeq.TS, method = "wunifrac")
dist.all.wuni.POST <- phyloseq::distance(physeq.POST, method = "wunifrac")

## make a data frame from the sample_data
sampledf.NS <- data.frame(sample_data(physeq.NS))
sampledf.NL <- data.frame(sample_data(physeq.NL))
sampledf.PH <- data.frame(sample_data(physeq.PH))
sampledf.TS <- data.frame(sample_data(physeq.TS))
sampledf.POST <- data.frame(sample_data(physeq.POST))

## Adonis test
#library("vegan")
adonis(dist.all.wuni.NS ~ C_vs_I, data = sampledf.NS)
adonis(dist.all.wuni.NL ~ C_vs_I, data = sampledf.NL)
adonis(dist.all.wuni.PH ~ C_vs_I, data = sampledf.PH)
adonis(dist.all.wuni.TS ~ C_vs_I, data = sampledf.TS)
adonis(dist.all.wuni.POST ~ C_vs_I, data = sampledf.POST)

#Ordinate the whole data set color by differet variables CCA
ord.d = ordinate(physeq.f, method = "CCA", distance = "wunifrac")
plot_ordination(physeq, ord.d, color = "Site") + 
  geom_point(mapping = aes()) + 
  ggtitle("All Samples by Site, CCA: WUniFrac") +
  scale_colour_brewer(type="qual", palette="Set1") +
  scale_fill_brewer(type="qual", palette="Set1") +
  geom_point(size=5) 

plot_ordination(physeq.f, ord.d, color = "C_vs_I") + 
  geom_point(mapping = aes()) +
  ggtitle("All Samples by Infection-status, CCA: WUniFrac")+
  scale_colour_brewer(type="qual", palette="Set1") +
  scale_fill_brewer(type="qual", palette="Set1") +
  geom_point(size=5) 
#facets
plot_ordination(physeq.f, ord.d, color = "C_vs_I") + 
  geom_point(mapping = aes()) +
  ggtitle("CCA: WUniFrac") + facet_wrap("Site")+
  scale_colour_brewer(type="qual", palette="Set1") +
  scale_fill_brewer(type="qual", palette="Set1") +
  geom_point(size=3) 

#Ordinate the whole data set color by differet variables DCA
ord.d = ordinate(physeq.f, method = "DCA", distance = "wunifrac")
plot_ordination(physeq, ord.d, color = "Site") + 
  geom_point(mapping = aes()) + 
  ggtitle("All Samples by Site, DCA: WUniFrac") +
  scale_colour_brewer(type="qual", palette="Set1") +
  scale_fill_brewer(type="qual", palette="Set1") +
  geom_point(size=5) 



#Ordinate the whole data set color by differet variables CCA
ord.c = ordinate(physeq.norm, method = "CCA", distance = "wunifrac")
plot_ordination(physeq.norm, ord.c, color = "Site") + 
  geom_point(mapping = aes()) +
  ggtitle("All Samples, CCA: WUniFrac") + scale_colour_brewer(type="qual", palette="Set1") +
  scale_fill_brewer(type="qual", palette="Set1") +
  geom_point(size=5) 
plot_ordination(physeq.norm, ord.c, color = "C_vs_I") + 
  geom_point(mapping = aes()) +
  ggtitle("All Samples, CCA: WUniFrac") +  scale_colour_brewer(type="qual", palette="Set1") +
  scale_fill_brewer(type="qual", palette="Set1") +
  geom_point(size=3) 
#facets
plot_ordination(physeq.norm, ord.c, color = "C_vs_I") + 
  geom_point(mapping = aes()) +
  ggtitle("CCA: WUniFrac") + facet_wrap("Site") +  scale_colour_brewer(type="qual", palette="Set1") +
  scale_fill_brewer(type="qual", palette="Set1") +
  geom_point(size=3) 

#plot different axes
plot_ordination(physeq.norm, ord.c, color = "C_vs_I", axes = 5:6) + 
  geom_point(mapping = aes()) +
  ggtitle("CCA: WUniFrac") + facet_wrap("Site") +  scale_colour_brewer(type="qual", palette="Set1") +
  scale_fill_brewer(type="qual", palette="Set1") +
  geom_point(size=3) 

 #Ordinate the whole data set color by differet variables NMDS
ord.n = ordinate(physeq.norm, method = "NMDS", distance = "wunifrac")
plot_ordination(physeq.norm, ord.n, color = "Site") + 
  geom_point(mapping = aes()) +
  ggtitle("All Samples, NMDS: WUniFrac") + scale_colour_brewer(type="qual", palette="Set1") +
  scale_fill_brewer(type="qual", palette="Set1") +
  geom_point(size=5) 
plot_ordination(physeq.norm, ord.n, color = "C_vs_I") + 
  geom_point(mapping = aes()) +
  ggtitle("All Samples, NMDS: WUniFrac") +  scale_colour_brewer(type="qual", palette="Set1") +
  scale_fill_brewer(type="qual", palette="Set1") +
  geom_point(size=3) 
#facets
plot_ordination(physeq.norm, ord.n, color = "C_vs_I") + 
  geom_point(mapping = aes()) +
  ggtitle("NMDS: WUniFrac") + facet_wrap("Site") +  scale_colour_brewer(type="qual", palette="Set1") +
  scale_fill_brewer(type="qual", palette="Set1") +
  geom_point(size=2) 

#ordinate each site using DCA
ordNS = ordinate(physeq.NS, method = "DCA", distance = "wunifrac")
plot_ordination(physeq.NS, ordNS, color = "C_vs_I") + 
  geom_point(mapping = aes()) +
  ggtitle("NS: DCA, WUniFrac")+geom_point(size=4)

ordNL = ordinate(physeq.NL, method = "DCA", distance = "wunifrac")
plot_ordination(physeq.NL, ordNL, color = "C_vs_I") + 
  geom_point(mapping = aes()) +
  ggtitle("NL: DCA: WUniFrac")+geom_point(size=4)

ordPH = ordinate(physeq.PH, method = "DCA", distance = "wunifrac")
plot_ordination(physeq.PH, ordPH, color = "C_vs_I") + 
  geom_point(mapping = aes()) +
  ggtitle("PH: DCA, WUniFrac")+geom_point(size=4)
ordTS = ordinate(physeq.TS, method = "DCA", distance = "wunifrac")
plot_ordination(physeq.TS, ordTS, color = "C_vs_I") + 
  geom_point(mapping = aes()) +
  ggtitle("TS: DCA: WUniFrac")+geom_point(size=4)

#ordinate each site using NMDS
ordNS.NMDS = ordinate(physeq.NS, method = "NMDS", distance = "wunifrac")
plot_ordination(physeq.NS, ordNS.NMDS, color = "C_vs_I") + 
  geom_point(mapping = aes()) +
  ggtitle("NS: NMDS, WUniFrac")

ordNL.NMDS = ordinate(physeq.NL, method = "NMDS", distance = "wunifrac")
plot_ordination(physeq.NL, ordNL.NMDS, color = "C_vs_I") + 
  geom_point(mapping = aes()) +
  ggtitle("NL: NMDS: WUniFrac")

ordPH.NMDS = ordinate(physeq.PH, method = "NMDS", distance = "wunifrac")
plot_ordination(physeq.PH, ordPH.NMDS, color = "C_vs_I") + 
  geom_point(mapping = aes()) +
  ggtitle("PH: NMDS, WUniFrac")
ordTS.NMDS = ordinate(physeq.TS, method = "NMDS", distance = "wunifrac")
plot_ordination(physeq.TS, ordTS.NMDS, color = "C_vs_I") + 
  geom_point(mapping = aes()) +
  ggtitle("TS: NMDS: WUniFrac")

#ordinate each site using CCA
ordNS.CCA = ordinate(physeq.NS, method = "CCA", distance = "wunifrac")
plot_ordination(physeq.NS, ordNS.CCA, color = "C_vs_I") + 
  geom_point(mapping = aes()) +
  ggtitle("NS: CCA, WUniFrac")+  scale_colour_brewer(type="qual", palette="Set1") +
  scale_fill_brewer(type="qual", palette="Set1") +
  geom_point(size=5) +
  ggtitle("NS: CCA, WUniFrac")

ordNL.CCA = ordinate(physeq.NL, method = "CCA", distance = "wunifrac")
plot_ordination(physeq.NL, ordNL.CCA, color = "C_vs_I", axes=c(1,2)) + 
  geom_point(mapping = aes()) +
  ggtitle("NL: CCA: WUniFrac") +  scale_colour_brewer(type="qual", palette="Set1") +
  scale_fill_brewer(type="qual", palette="Set1") +
  geom_point(size=5) +ggtitle("NL: CCA, WUniFrac")

ordPH.CCA = ordinate(physeq.PH, method = "CCA", distance = "wunifrac")
plot_ordination(physeq.PH, ordPH.CCA, color = "C_vs_I") + 
  scale_colour_brewer(type="qual", palette="Set1")+
  geom_point(mapping = aes()) +geom_point(size=5)+
  ggtitle("PH: CCA, WUniFrac")

ordTS.CCA = ordinate(physeq.TS, method = "CCA", distance = "wunifrac")
plot_ordination(physeq.TS, ordTS.CCA, color = "C_vs_I") + 
  geom_point(mapping = aes()) +scale_colour_brewer(type="qual", palette="Set1")+
  ggtitle("TS: CCA: WUniFrac") +geom_point(size=5)

## Try building a clustering for a site: http://www.bioconductor.org/packages/3.7/bioc/vignettes/phyloseq/inst/doc/phyloseq-analysis.html
# not very useful IMO and a lot of work to present clearly
#NS.wuni = UniFrac(physeq.NS, weighted = TRUE)
#NS.hclust = hclust(NS.wuni, method="average")
#plot(NS.hclust)

taxa.pal = palette(c(rgb(97,105,216,maxColorValue=255),
             rgb(90,188,73,maxColorValue=255),
             rgb(148,71,189,maxColorValue=255),
             rgb(157,187,54,maxColorValue=255),
             rgb(196,117,229,maxColorValue=255),
             rgb(69,134,53,maxColorValue=255),
             rgb(211,71,174,maxColorValue=255),
             rgb(99,192,124,maxColorValue=255),
             rgb(230,71,142,maxColorValue=255),
             rgb(79,192,171,maxColorValue=255),
             rgb(224,57,96,maxColorValue=255),
             rgb(55,131,93,maxColorValue=255),
             rgb(192,49,36,maxColorValue=255),
             rgb(78,175,216,maxColorValue=255),
             rgb(229,100,44,maxColorValue=255),
             rgb(98,144,227,maxColorValue=255),
             rgb(202,176,56,maxColorValue=255),
             rgb(135,79,159,maxColorValue=255),
             rgb(128,140,36,maxColorValue=255),
             rgb(217,127,205,maxColorValue=255),
             rgb(221,146,48,maxColorValue=255),
             rgb(81,102,165,maxColorValue=255),
             rgb(182,97,37,maxColorValue=255),
             rgb(176,149,217,maxColorValue=255),
             rgb(99,111,45,maxColorValue=255),
             rgb(176,58,128,maxColorValue=255),
             rgb(160,176,102,maxColorValue=255),
             rgb(183,59,99,maxColorValue=255),
             rgb(201,164,88,maxColorValue=255),
             rgb(155,86,134,maxColorValue=255),
             rgb(138,107,43,maxColorValue=255),
             rgb(226,133,173,maxColorValue=255),
             rgb(148,73,35,maxColorValue=255),
             rgb(225,130,133,maxColorValue=255),
             rgb(213,154,108,maxColorValue=255),
             rgb(154,74,92,maxColorValue=255),
             rgb(234,136,100,maxColorValue=255),
             rgb(179,62,65,maxColorValue=255),
             rgb(235,94,85,maxColorValue=255),
             rgb(179,94,77,maxColorValue=255)))         
