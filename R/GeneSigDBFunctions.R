#################################
## Parsing GeneSigDB R Functions
#################################


readGeneSigDBFile<-function(GeneSigDBPath=GeneSigDBdata,GeneSigDBFileName="GeneSigDB-Table 1.csv") {
  # Function to Read Gene SigDB File(shows the entire GeneSigDB)
  # To call this function type: readGeneSigDBFile(GeneSigDBdata,GeneSigDBFileName)

   GSdb= NULL
  SigFile= file.path(GeneSigDBdata, GeneSigDBFileName)
   if (file.exists(SigFile)){
    GSdb<-read.csv(SigFile, as.is=TRUE)
  } else print(paste("Can't find",GeneSigDBFileName, "in", GeneSigDBPath))
  # Check Data
  return(GSdb)
}


getSig<- function(SigID,GeneSigIndex,...)
  {
# Example SigID = "10582678-Table1"

# To call this function (executing it)
# getSig(SigID,GeneSigDBData=GeneSigDBData)
#

#GeneSigIndex is the GeneSigdb.xls file read using readGeneSigDBFile(GeneSigDBPath,GeneSigDBFileName)
# It turns a data.frame of the signature

  rowInd= GeneSigIndex$SigID==SigID
  print(table(rowInd))

  SigFilePath=file.path(GeneSigDBData,GeneSigIndex$Release[rowInd], GeneSigIndex$PMID[rowInd], GeneSigIndex$FileAssociated[rowInd])
  print(SigFilePath)


  if(file.exists(SigFilePath)) {
    Sig<- read.table(SigFilePath, header=TRUE, sep="\t", as.is=TRUE)
    return(Sig)
  } else print(paste("Can't read in", SigID, "file in", SigFilePath))
  #check Data
}

######################################################

parseSigCols<-function(SigID, GeneSigIndex,...) {
## Extract Cols from Sig

# GeneSigIndex is the GeneSigdb.xls file read using readGeneSigDBFile(GeneSigDBPath,GeneSigDBFileName)
# It turns a data.frame of the signature

  rowInd= GeneSigIndex$SigID==SigID
# Get cols
  colInd=grep("^Column\\d+", colnames(GeneSigIndex))
  SigCols= GeneSigIndex[rowInd, colInd]
  SigCols<- SigCols[!c(is.na(SigCols)| SigCols=="")]

  sig<-getSig(SigID,GeneSigIndex)
}

#' Function to convert identifiers
#'
#' @param ids is a character vector of gene or protein identifiers
#' @param Identifer is a character indicating the identifer types, eg "EnsEMBL ID", "Gene Symbol", "EntrezGene ID" etc
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#' ids =c("IL6", "TP53", "SOX8")
#' mart=getMart("human")
#' parseIDs(ids, Identifer="Gene Symbol")
#'
parseIDs<-function(ids=Sig[,2], identifer="Gene Symbol", ...) {
    mapIds<-switch(identifer,
       "Clone ID"=parseCloneID(ids),
       "EnsEMBL ID"= parseEnsEnsEMBL(ids),
       "GenBank ID"= parseGenBankID(ids),
       "Gene Symbol"= parseGeneSymbol(ids),
       "EntrezGene ID"= parseEntrezGeneID(ids),
       "UniGene ID"= parseUniGeneID(ids),
       "miRBase"= parsemiRBase(ids),
       "Protein ID"= parseProteinID(ids),
       "RefSeq ID"= parseRedSeqID(ids),
       "Probe ID"= parseProbeID(ids)
    )

    # .parseNull returns NULL. These ids are not searched
    if (x=="Secondary Probe ID") .parseNULL(ids)
    if (x=="Gene Description") .parseNULL(ids)
    if (x=="Other Gene Description") .parseNULL(ids)
    if (x=="Geneset Specific Factor") .parseNULL(ids)
    if (x=="Geneset Specific Statistics") .parseNULL(ids)
    if (x=="Chromosome Map") .parseNULL(ids)

}


parseCloneID<-function() {
   # 1. Check format looks correct
   # 2. search
  return("CloneID")
}

parseEnsEnsEMBL<-function(ids, biomart=TRUE) {
  # validateIDs

  if(biomart){
    #we are selecting for "ensembl_gene_id", "embl", "entrezgene", "hgnc_symbol")
    mart<-getMart()
    #mapping<-biomaRt::select(mart, keys=ids, columsn=c("ensembl_gene_id","embl","entrezgene", "hgnc_symbol"), keytype="embl")
    mapping<-getBM(attributes= c("ensembl_gene_id", 	"entrezgene", "hgnc_symbol"), values=ids,
                   filters="ensembl_gene_id", mart=mart)
    return(mapping)
  }
  }

parseGenBankID<-function(ids, biomart=TRUE) {
  if(biomart){
    #selecting for ("ensembl_gene_id","embl","entrezgene", "hgnc_symbol")
    mart<-getMart()
    #mapping<-biomaRt::select(mart, keys=ids, columns=c("ensembl_gene_id","embl","entrezgene", "hgnc_symbol"), keytype="embl")
    mapping<-getBM(attributes= c("embl", 	"entrezgene", "hgnc_symbol"), values=ids, filters="embl", mart=mart)
    return(mapping)
  }
  }

parseGeneSymbol<- function(ids, biomart=TRUE) {
  if(biomart){
    #selecting for ("ensembl_gene_id","embl","entrezgene", "hgnc_symbol")
    mart<-getMart()
    #mapping<-biomaRt::select(mart, keys=ids, columns=c("ensembl_gene_id","embl","entrezgene", "hgnc_symbol"), keytype="embl")
    mapping<-getBM(attributes= c("hgnc_symbol", "ensembl_gene_id",	"entrezgene"), values=ids, filters="hgnc_symbol", mart=mart)
    return(mapping)

  }
  }
parseEntrezGeneID<- function(ids, biomart=TRUE) {
  if(biomart){
    #selecting for ("ensembl_gene_id","embl","entrezgene", "hgnc_symbol")
    mart<-getMart()
    #mapping<-biomaRt::select(mart, keys=ids, columns=c("ensembl_gene_id","embl","entrezgene", "hgnc_symbol"), keytype="embl")
    mapping<-getBM(attributes= c("embl", 	"entrezgene", "hgnc_symbol"), values=ids, filters="embl", mart=mart)
    return(mapping)
  }
  }


parseUniGeneID<-function(ids, biomart=TRUE) {
  # validateIDs

  if(biomart){
    #we are selecting for "ensembl_gene_id", "embl", "entrezgene", "hgnc_symbol")
    mart<-getMart()
    #mapping<-biomaRt::select(mart, keys=ids, columsn=c("ensembl_gene_id","embl","entrezgene", "hgnc_symbol"), keytype="embl")
    mapping<-getBM(attributes= c("embl", 	"entrezgene", "hgnc_symbol"), values=ids,
                   filters="embl", mart=mart)
    return(mapping)
  }
  }

parsemiRBase<-function(ids, biomart=TRUE) {
  # validateIDs

  if(biomart){
    #we are selecting for "ensembl_gene_id", "embl", "entrezgene", "hgnc_symbol")
    mart<-getMart()
    #mapping<-biomaRt::select(mart, keys=ids, columsn=c("ensembl_gene_id","embl","entrezgene", "hgnc_symbol"), keytype="embl")
    mapping<-getBM(attributes= c("embl", 	"entrezgene", "hgnc_symbol"), values=ids,
                   filters="embl", mart=mart)
    return(mapping)
  }
  }


parseProteinID<-function(ids, biomart=TRUE) {
  # validateIDs

  if(biomart){
    #we are selecting for "ensembl_gene_id", "embl", "entrezgene", "hgnc_symbol")
    mart<-getMart()
    #mapping<-biomaRt::select(mart, keys=ids, columsn=c("ensembl_gene_id","embl","entrezgene", "hgnc_symbol"), keytype="embl")
    mapping<-getBM(attributes= c("embl", 	"entrezgene", "hgnc_symbol"), values=ids,
                   filters="embl", mart=mart)
    return(mapping)
  }
  }

parseRefSeqID<-function(ids, biomart=TRUE) {
  # validateIDs

  if(biomart){
    #we are selecting for "ensembl_gene_id", "embl", "entrezgene", "hgnc_symbol")
    mart<-getMart()
    #mapping<-biomaRt::select(mart, keys=ids, columsn=c("ensembl_gene_id","embl","entrezgene", "hgnc_symbol"), keytype="embl")
    mapping<-getBM(attributes= c("embl", 	"entrezgene", "hgnc_symbol"), values=ids,
                   filters="embl", mart=mart)
    return(mapping)
  }
  }


parseProbeID<-function(ids, biomart=TRUE) {
  # validateIDs

  if(biomart){
    #we are selecting for "ensembl_gene_id", "embl", "entrezgene", "hgnc_symbol")
    mart<-getMart()
    #mapping<-biomaRt::select(mart, keys=ids, columsn=c("ensembl_gene_id","embl","entrezgene", "hgnc_symbol"), keytype="embl")
    mapping<-getBM(attributes= c("embl", 	"entrezgene", "hgnc_symbol"), values=ids,
                   filters="embl", mart=mart)
    return(mapping)
  }
  }


.parseNULL<-function (ids) {
  return(NULL)
}


#' Function to define biomart Mart dataset
#'
#' @param species is character, human, mouse, rat or chicken. Default is human
#' @param force is logical. Default is TRUE. If set to FALSE getMart will not create a new mart, if an object calls mart exists in the workspace
#'
#' @return  object of class biomarRt:::Mart
#' @export
#'
#' @examples
#' mart<- getMart("human")
#' class(mart)
#' mart<- getMart("mouse", force=FALSE)

getMart<-function(species="human", force=TRUE){

    require(biomaRt)
    # Force will create a new mart object even if one exists in the current workspace.  To avoid creating a new mart, set force=FALSE
    if (exists("mart")) {
      print("an object  'mart' exists in your workspace")
   }
  if (force==TRUE | !exists("mart")){
  mart<- switch(species,
       mouse= useMart("ENSEMBL_MART_ENSEMBL",dataset="mmusculus_gene_ensembl", host="www.ensembl.org"),
       human= useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="www.ensembl.org"),
     chicken= useMart("ENSEMBL_MART_ENSEMBL",dataset="ggallus_gene_ensembl", host="www.ensembl.org"),
    rat= useMart("ENSEMBL_MART_ENSEMBL",dataset="rnorvegicus_gene_ensembl", host="www.ensembl.org")
   )}
    return(mart)
  }



# if (search="AnnotationDBI")  {
#   # ALTERNATIVE
#   # Using Annotation dbi
#   require(org.Hs.eg.db)
#   columns(org.Hs.eg.db)
#   #Annotation DBi, you are selecting for Accession numbers, ensembl id, entrezid, symbol
#   #keytype is Accession Number
#   mapping<- AnnotationDbi::select(org.Hs.eg.db, keys = ids, columns= c("ACCNUM","ENSEMBL"  ,"ENTREZID", "SYMBOL" ), keytype = "ACCNUM")
#   return(mapping)
#   }
