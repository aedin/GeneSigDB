##############################
## ParseSignature.R
## Script to convert gene identifiers in GeneSigDB tables to Common Gene ID (EntrezGeneID or GeneSymbols)
## Aedin Culhane & Ali Ahmed
## email:ali.ahmed01867@gmail.com
## July 12th 2016
##############################



# Define Variables
GeneSigDB_ReleaseData= file.path("../gene_signatures/data")
GeneSigDBdata = file.path("data")
GeneSigDBsrc = file.path("R")
GeneSigRDa ="GeneSigDB.rda"

if (!file.exists(file.path(GeneSigDBdata, GeneSigRDa))){
  GeneSigDBFileName= "GeneSigDB.xls"
  GSdb<-readGeneSigDBFile()
}

# Load R/Bioc Libs
library(AnnotationDbi)
library(biomaRt)
#library(annotate)
source(file.path(GeneSigDBsrc,"GeneSigDBFunctions.R"))

################################################################

# 1. Read in GeneSigDB File.
load(file.path(GeneSigDBdata, GeneSigRDa))
GeneSigIndex = GSdb$GeneSigIndex

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
mapIDs=c( "EnsEMBL ID", "EntrezGene ID","GenBank ID", "Gene Symbol")


sigID=GeneSigIndex$SigID[3]


#mappedSig<-mapSig(sigID,GeneSigIndex,GeneSigDB_ReleaseData, attributes=c("ensembl_gene_id", "hgnc_symbol", "entrezgene"), species="human", verbose=TRUE)


res<-lapply(GeneSigIndex$SigID ,function(sigID) mapSig(sigID,GeneSigIndex,GeneSigDB_ReleaseData, attributes=c("ensembl_gene_id", "hgnc_symbol", "entrezgene"), species="human"))

save(res, file="testing.rda")

# # Using Annotation dbi
# library(TxDb.Hsapiens.UCSC.hg19.knownGene)
# require(org.Hs.eg.db)
# columns(org.Hs.eg.db)
# #Annotation DBi, you are selecting for Accession numbers, ensembl id, entrezid, symbol
# #keytype is Accession Number
# AnnotationDbi::select(org.Hs.eg.db, keys = ids, columns= c("ACCNUM","ENSEMBL"  ,"ENTREZID", "SYMBOL" ), keytype = "ACCNUM")
#
