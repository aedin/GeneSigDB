### GeneSigDB checking : Release 1 ######


# Define where the raw Sigs are (gene_signatures)
GeneSigDB_ReleaseData= "/Users/aliahmed/Documents/Stats/gene_signatures/data"

# Don't edit.. GeneSigDB repository
GeneSigDBdata = "/Users/aliahmed/Documents/Stats/GeneSigDB/data"
GeneSigDBsrc = "/Users/aliahmed/Documents/Stats/GeneSigDB/R"
GeneSigRDa ="GeneSigDB.rda"

# Load R/Bioc Libs
library(AnnotationDbi)
library(biomaRt)
#library(annotate)
source(file.path(GeneSigDBsrc,"GeneSigDBFunctions.R"))
require(dplyr)
require(org.Hs.eg.db)
###########

# Read GeneSigDB file

if (!file.exists(file.path(GeneSigDBdata, GeneSigRDa))){
  GeneSigDBFileName= "GeneSigDB_Working.xlsx"
  GSdb<-readGeneSigDBFile()
  GSdb<-readGeneSigDBFile()
}



################################################################

# 1. Read in GeneSigDB File.
load(file.path(GeneSigDBdata, GeneSigRDa))
GeneSigIndex = GSdb$GeneSigIndex
GeneSigIndex = GeneSigIndex[GeneSigIndex$Release=="R1",]
## Ali
GeneSigIndex = GeneSigIndex[1:302,]
GeneSigIndex = GeneSigIndex[!grepl("NA", rownames(GeneSigIndex)), ]
## Makena
#GeneSigIndex = GeneSigIndex[303:nrow(GeneSigIndex),]


# 2. Read a Gene Signature File
# Number of identifiers
# table(unlist(GeneSigIndex[,colInd]))
# Clone ID   136
# EnsEMBL ID    84
# EntrezGene ID   623
# GenBank ID  1299
# Gene Symbol  4745
# Probe ID  1757
# Protein ID  29
# RefSeq ID   738
# UniGene ID   562

#mapIDs=c("Clone ID", "EnsEMBL ID", "EntrezGene ID","GenBank ID", "Gene Symbol", "Probe ID", "Protein ID", "RefSeq ID", "UniGene ID")
#mapIDs=c( "EnsEMBL ID", "EntrezGene ID","GenBank ID", "Gene Symbol")


#sigID=GeneSigIndex$SigID[13]


#mapSig("16964388-SuppTable1",GeneSigIndex,GeneSigDB_ReleaseData, attributes=c("ensembl_gene_id", "hgnc_symbol", "entrezgene"), species="human", verbose=TRUE)


## change accordingly
numbers = 193:200
##
res<-lapply(numbers,function(i){
  print(i)
  sigID= GeneSigIndex$SigID[i]
  mapSig(sigID,GeneSigIndex,GeneSigDB_ReleaseData, attributes=c("ensembl_gene_id", "hgnc_symbol", "entrezgene"), species="human")
})


### Modifying files
require(readr)
file <- "/Users/aliahmed/Documents/Stats/gene_signatures/data/R1/15790403/15790403-Table2.txt"
#View(read.table(file, header = T, sep = "\t"))
# If above works skip below
x <- read_tsv(file)
View(x)
#write_tsv(x, file)
#View(read.table(file, header = T, sep = "\t"))


