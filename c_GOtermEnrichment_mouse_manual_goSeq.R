#biocLite("org.Mm.eg.db")
#biocLite("goseq")
#biocLite("biomaRt")
library(biomaRt)
library(org.Mm.eg.db)
library(goseq)
library(geneLenDataBase)

# gene list from http://amp.pharm.mssm.edu/Enrichr/#
# saved as genesGOtest.txt
gns <- read.table("genesGOtest.txt", stringsAsFactors = FALSE)
# retrieve the ENSEMBL symbols for the gene names
anno.mm <- select(org.Mm.eg.db,
               keys = gns$V1, 
               keytype="SYMBOL", # our rownames are ORF identifiers
               columns=c("ENSEMBL","SYMBOL","GENENAME")) # what to return

# in addition, retrieve all possible mouse gene names
# because goseq wants to know about the universe of all genes
mart <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
mm_all <- getBM(attributes = c("ensembl_gene_id"),  mart = mart)

# get length data using goseq's data base of gene lengths
# first, load the data from the geneLenDataBase package
data(list = "mm9.ensGene.LENGTH", package = "geneLenDataBase")

# split the length values per ENSEMBL ID
len_data = split(get("mm9.ensGene.LENGTH")$Length, 
            get("mm9.ensGene.LENGTH")$Gene)

# get the gene names from the length data base
gn_names = get("mm9.ensGene.LENGTH")$Gene

# calculate the median (and other) lengths for the 
# different transcripts
len = data.frame(Gene = names(len_data), 
                 Median = sapply(len_data,median), 
                 Min = sapply(len_data, min),
                 Max = sapply(len_data, max),
                 Count = sapply(len_data, length))

# goseq wants a named binary vector where 1 represents DE, 0 not DE and the names are gene IDs.
gogns <- as.integer(len$Gene %in% anno.mm$ENSEMBL)
names(gogns) <- len$Gene

# weighting for each gene, depending on its length, given by the PW
# Probability Weighting Function (PWF)
# proportion of DE genes is plotted as a function of the transcript length
pwf = goseq::nullp(gogns, bias.data = len$Median) # probably looks bad because the genes weren't taken from a DE set, but from ChIP-seq, I think

GO.wall=goseq(pwf, "mm9", "ensGene")

# can be summarized in http://revigo.irb.hr/
write.table(GO.wall[c("category","over_represented_pvalue")], "test.txt", quote =FALSE, row.names = FALSE, col.names = FALSE)