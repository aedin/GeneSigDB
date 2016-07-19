#######################
## ParseSignature.R
## Script to convert gene identifiers in GeneSigDB tables to Common Gene ID (EntrezGeneID or GeneSymbols)
## Ali Ahmed  
## email:ali.ahmed01867@gmail.com
## July 12th 2016
###########################



# Define Variables
GeneSigDBPath= "~aliahmed/Documents/Stats/GeneSigDB 2"
GeneSigDBdata = file.path(GeneSigDBPath, "data")
GeneSigDBsrc = file.path(GeneSigDBPath, "src")
GeneSigDBData= "/Volumes/gene_signatures/trunk/data/"

GeneSigDBFileName= "GeneSigDB-Table 1.csv"

# Load R/Bioc Libs
library(AnnotationDbi)
library(biomaRt)
#library(annotate)
source(file.path(GeneSigDBsrc,"GeneSigDBFunctions.R"))

################################################################

# 1. Read in GeneSigDB File.
GeneSigIndex<-readGeneSigDBFile(GeneSigDBdata,GeneSigDBFileName) 

# 2. Read a Gene Signature File
#enter in Sig Id that is in the "GeneSigDB-Table 1.csv
#Sig is 
sig<-getSig("10582678-Table1",GeneSigIndex)
ids<-sig[,2]

#Reading a Gene Signature File that Select for something specific()

# 3. Map Genes

## Biomart:
#homo sapiens
require(biomaRt)
mart<-useMart(dataset="hsapiens_gene_ensembl", biomart="ensembl")

#Rats 
require(biomaRt)
mart<-useMart(dataset="rnorvegicus_gene_ensembl", biomart="ensembl")
listDatasets()

#mouse
require(biomaRt)
mart<-useMart(dataset="musculus_gene_ensembl", biomart="ensembl")
listDatasets()

# List all columns
attributes(packageDescription("biomaRt"))

#selecting for ("ensembl_gene_id","embl","entrezgene", "hgnc_symbol")
#keytype= "embyl"
biomaRt::select(mart, keys=ids, columns=c("ensembl_gene_id","embl","entrezgene", "hgnc_symbol"), keytype="embl")
getBM(attributes= c("ensembl_gene_id", 	"entrezgene", "hgnc_symbol"), values=ids, filters="embl", mart)

# Using Annotation dbi
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
require(org.Hs.eg.db)
columns(org.Hs.eg.db)
#Annotation DBi, you are selecting for Accession numbers, ensembl id, entrezid, symbol
#keytype is Accession Number
AnnotationDbi::select(org.Hs.eg.db, keys = ids, columns= c("ACCNUM","ENSEMBL"  ,"ENTREZID", "SYMBOL" ), keytype = "ACCNUM")

