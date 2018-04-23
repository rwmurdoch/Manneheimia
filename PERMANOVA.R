# permanova, test if sites have different centroids
dist.all.wuni <- phyloseq::distance(physeq.f, method = "wunifrac")

## make a data frame from the sample_data
sampledf <- data.frame(sample_data(physeq.f))

## Adonis test
#library("vegan")
adonis(dist.all.wuni ~ Site, data = sampledf)
adonis(dist.all.wuni ~ C_vs_I, data = sampledf)
adonis(dist.all.wuni ~ Study_day, data = sampledf)
adonis(dist.all.wuni ~ Tag_Number, data = sampledf)

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

## Adonis test for CvsI
#library("vegan")
adonis(dist.all.wuni.NS ~ C_vs_I, data = sampledf.NS) #sig
adonis(dist.all.wuni.NL ~ C_vs_I, data = sampledf.NL) #sig
adonis(dist.all.wuni.PH ~ C_vs_I, data = sampledf.PH) #sig
adonis(dist.all.wuni.TS ~ C_vs_I, data = sampledf.TS)
adonis(dist.all.wuni.POST ~ C_vs_I, data = sampledf.POST)

## Adonis test for Sample_day
#library("vegan")
adonis(dist.all.wuni.NS ~ Study_day, data = sampledf.NS)
adonis(dist.all.wuni.NL ~ Study_day, data = sampledf.NL) #sig
adonis(dist.all.wuni.PH ~ Study_day, data = sampledf.PH) #sig
adonis(dist.all.wuni.TS ~ Study_day, data = sampledf.TS)
adonis(dist.all.wuni.POST ~ Study_day, data = sampledf.POST)

## Adonis test for Tag_Number
#library("vegan")
adonis(dist.all.wuni.NS ~ Tag_Number, data = sampledf.NS) #sig
adonis(dist.all.wuni.NL ~ Tag_Number, data = sampledf.NL) #sig
adonis(dist.all.wuni.PH ~ Tag_Number, data = sampledf.PH) #sig
adonis(dist.all.wuni.TS ~ Tag_Number, data = sampledf.TS) #sig
adonis(dist.all.wuni.POST ~ Tag_Number, data = sampledf.POST)
