## this script allows for plotting taxa over time with display of every sample and..
## display of a moving average with variance shadow


# define your physeq or subset here
set <- physeq.f
# define your target taxon here
TAXON <- "Mannheimia"
# define chart title here
TITLE = "Relative abundance of Mannheimia"

## You will have to alter manually the variables you want to facet by in the ggplot call at the end of the script

physeq.pasturellaceae = subset_taxa(set, genus == TAXON)

### now I need to subset by infection status and turn each into relative abundance as derived from sum of 
### counts in each sample

# first turn the taxon-specific OTU tables into dataframes

OTU.P = as(otu_table(physeq.pasturellaceae), "matrix")
OTU.P = as.data.frame(OTU.P)
P.sums = as.data.frame(colSums(OTU.P))

# then get the entire study as an OTU data frame and pull out sums for each sample
OTU.Site = as(otu_table(set), "matrix")
OTU.Site = as.data.frame(OTU.Site)
Site.sums = as.data.frame(colSums(OTU.Site))

#and divide the taxon counts by total sequence sums by sample

OTU.P.all.rel = P.sums / Site.sums

# now I have full Taxon OTU tables in relative abundance format

### moving towards a by-sample plot of sum Pasteurellaceae by day


# just add the rel Pasteurellaceae to sub-setted metadata tables
metadata.Site = sample_data(set)
metadata.Site$Study_day <- as.numeric(sample_data(set)$Study_day) #make sure day is numeric, helps with the plotting here
OTU.P.all.rel.sum.meta = cbind(metadata.Site, OTU.P.all.rel)
colnames(OTU.P.all.rel.sum.meta) [7] = "Relative.Abundance"

# now I have added rel-abundance of Pasteurellaceae to the metadata file
# I have to figure a nice wqy to chart this out

library(ggplot2)
library(scales)

#scatter.smooth doesn't handle missing values. A simple way to proceed is to 
#remove the missing data:

#         good <- complete.cases(Mass, HSI)
#         scatter.smooth(Mass[good], HSI[good])

df = OTU.P.all.rel.sum.meta
options(scipen=999)
theme_set(theme_bw())
gg <- ggplot(df, aes(x=Study_day, y=Relative.Abundance)) + 
  geom_point(aes(col=Tag_Number)) + 
  geom_smooth(method = "auto", se=TRUE, level = 0.66) + 
  xlim(c(1, 6)) + 
  ylim(c(-.2, 1)) +
  scale_y_log10(limits=c(0.0000001,1)) +
    labs(subtitle="66% confidence interval", 
       y="log10 relative abundance", 
       x="Study Day", 
       title=TITLE) +
  facet_grid(rows=vars(C_vs_I), cols=vars(Site))
plot(gg)


### this is some code for adding taxonomy to the OTU table, not using it currenlty

# pull out taxonomy map
#Taxa.P = tax_table(physeq.pasturellaceae)
#Taxa.P = as.data.frame(Taxa.P)

#adding taxonomy

##library(data.table)
#setDT(OTU.P, keep.rownames = TRUE)[]
#setDT(Taxa.P, keep.rownames = TRUE)[]
#OTU.full.P = dplyr::left_join(OTU.P, Taxa.P, by = "rn")

#setDT(OTU.P, keep.rownames = TRUE)[]
#setDT(Taxa.P, keep.rownames = TRUE)[]
#OTU.full.P = dplyr::left_join(OTU.P, Taxa.P, by = "rn")

#setDT(OTU.P, keep.rownames = TRUE)[]
#setDT(Taxa.P, keep.rownames = TRUE)[]
#OTU.full.P = dplyr::left_join(OTU.P, Taxa.P, by = "rn")

# so we now have a table of all Pasturellaceae
