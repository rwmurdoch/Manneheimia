otu	<- read.table(file = "phyloseq/otu_table.txt",	header = TRUE)
head(otu)

tax	<- read.table(file="phyloseq/taxonomy.tsv", sep	=	'\t',	header = TRUE)
head(tax)

merged_file	<-	merge(otu, tax, by.x = c("OTUID"), by.y=c("OTUID"))
head(merged_file)

write.table(merged_file, file	=	"combined_otu_tax",	sep = '\t', col.names = 
  TRUE,	row.names	=	FALSE)

#	It	seems	tedious	but	you	need	to	open	the	merged	.txt	file	in	excel	and	split	
#into	two	files:	one	for	taxonomy	(containing	only	the	columns	OTUID	and	
#taxonomic	info)	and	the	other	for	the	OTU	matrix	(containing	only	OTUID	and	
#abundances	in	each	sample).	Note:	for	the	taxonomy	file,	you	need	to	use	data	
#—>	text-to-columns	in	Excel	and	separate	on	semicolon	to	get	columns	for	
#kingdom,	phylum,	class,	etc…	once	you	make	these	two	separate	files	in	excel,	
#save	each	as	a	.csv	
