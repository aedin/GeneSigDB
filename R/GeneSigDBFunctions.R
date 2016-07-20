#################################
## Parsing GeneSigDB R Functions
#################################


readGeneSigDBFile <- function(GeneSigDBPath="data",GeneSigDBFileName="GeneSigDB.xls") {
  # Function to Read Gene SigDB File(shows the entire GeneSigDB)
  # To call this function type: readGeneSigDBFile(GeneSigDBdata,GeneSigDBFileName)
  SigFile <- file.path(GeneSigDBdata, GeneSigDBFileName)
  GSdb <- list()
  if (file.exists(SigFile)) {
    GSdb$GeneSigIndex <- data.frame(readxl::read_excel(SigFile, "GeneSigDB"))
    GSdb$Identifiers <- data.frame(readxl::read_excel(SigFile, "Identifiers"))
    GSdb$Platforms <- data.frame(readxl::read_excel(SigFile, "Platforms"))
    saveRDS(GSdb, file=file.path(GeneSigDBPath, "GeneSigDB.rds"))
    save(GSdb, file=file.path(GeneSigDBPath, "GeneSigDB.rda"))
    } else print(paste("Can't find",GeneSigDBFileName, "in", GeneSigDBPath))
  # Check Data
  return(GSdb)
}


getSig<- function(SigID,GeneSigIndex, GeneSigDB_ReleaseData,...)
  {
# Example SigID = "10582678-Table1"

# To call this function (executing it)
# getSig(SigID,GeneSigDB_ReleaseData=GeneSigDB_ReleaseData)
#

#GeneSigIndex is the GeneSigdb.xls file read using readGeneSigDBFile(GeneSigDBPath,GeneSigDBFileName)
# It turns a data.frame of the signature

  rowInd= GeneSigIndex$SigID%in%SigID
  print(table(rowInd))

  SigFilePath=file.path(GeneSigDB_ReleaseData,GeneSigIndex$Release[rowInd], GeneSigIndex$PMID[rowInd], GeneSigIndex$FileAssociated[rowInd])
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

  rowInd= GeneSigIndex$SigID%in%SigID
# Get cols
  colInd=grep("^Column\\d+", colnames(GeneSigIndex))

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

  SigCols= GeneSigIndex[rowInd, colInd]
  SigCols<- SigCols[!c(is.na(SigCols)| SigCols=="")]


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
parseIDs<-function(ids=Sig[,2], identifer=SigCols[2], ...) {
    identifer<-as.character(identifer)
    print(identifer)
    mapIds<-switch(identifer,
       "Clone ID"=parseCloneID(ids),
       "EnsEMBL ID"= parseEnsEnsEMBL(ids),
       "GenBank ID"= parseGenBankID(ids),
       "Gene Symbol"= parseSymbol(ids),
       "EntrezGene ID"= parseEntrezGeneID(ids),
       "UniGene ID"= parseUniGeneID(ids),
       "miRBase"= parsemiRBase(ids),
       "Protein ID"= parseProteinID(ids),
       "RefSeq ID"= parseRedSeqID(ids),
       "Probe ID"= parseProbeID(ids),
        "Secondary Probe ID"= .parseNULL(ids),
        "Gene Description"= .parseNULL(ids),
        "Other Gene Description"= .parseNULL(ids),
        "Geneset Specific Factor"= .parseNULL(ids),
        "Geneset Specific Statistics"= .parseNULL(ids),
        "Chromosome Map"= .parseNULL(ids))
   return(mapIds)
}

getBMall <- function(attributes, filters = '', values = '', mart, curl = NULL, checkFilters = TRUE, verbose = FALSE, uniqueRows = TRUE, bmHeader = FALSE) {
  spotty <- getBM(attributes, filters, values, mart, curl, checkFilters, verbose, uniqueRows, bmHeader)
  print(spotty)
  x <- as.data.frame(values, stringsAsFactors =FALSE)
  colnames(x) <- filters
  print(x)
  structure(dplyr::left_join(x = x, y = spotty, by = filters, copy=TRUE))
}
parseCloneID<-function() {
   # 1. Check format looks correct
   # 2. search
  return("CloneID")
}

parseEnsEnsEMBL<-function(ids, biomart=TRUE, species="human") {
  # validateIDs

  if(biomart){
    #we are selecting for "ensembl_gene_id", "embl", "entrezgene", "hgnc_symbol")
    mart<-getMart(species)
    #mapping<-biomaRt::select(mart, keys=ids, columsn=c("ensembl_gene_id","embl","entrezgene", "hgnc_symbol"), keytype="embl")
    mapping<-getBM(attributes= c("ensembl_gene_id", 	"entrezgene", "hgnc_symbol"), values=ids,
                   filters="ensembl_gene_id", mart=mart)
    return(mapping)
  }
  }

parseGenBankID<-function(ids, biomart=TRUE, species="human") {
  print("parseGenBankID")
  if(biomart){
    #selecting for ("ensembl_gene_id","embl","entrezgene", "hgnc_symbol")
    mart<-getMart(verbose=FALSE)
    #mapping<-biomaRt::select(mart, keys=ids, columns=c("ensembl_gene_id","embl","entrezgene", "hgnc_symbol"), keytype="embl")
    mapping<-getBM(attributes= c("embl", 	"entrezgene", "hgnc_symbol"), values=ids, filters="embl", mart=mart)
    return(mapping)
  }
  }

parseSymbol<- function(ids, biomart=TRUE, species="human") {
  if(biomart){
    #selecting for ("ensembl_gene_id","embl","entrezgene", "hgnc_symbol")
    mart <- getMart(species, verbose=FALSE)
    #mapping<-biomaRt::select(mart, keys=ids, columns=c("ensembl_gene_id","embl","entrezgene", "hgnc_symbol"), keytype="embl")
    mapping <- getBM(attributes= c("hgnc_symbol", "ensembl_gene_id",	"entrezgene"), values=ids, filters="hgnc_symbol", mart=mart)
    return(mapping)
    }
}

parseEntrezGeneID<- function(ids, biomart=TRUE,species="human") {
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
#' @param verbose is logical. Default is TRUE. If TRUE will print warning if mart exists
#'
#' @return  object of class biomarRt:::Mart
#' @export
#'
#' @examples
#' mart<- getMart("human")
#' class(mart)
#' mart<- getMart("duck")
#'
getMart <- function(species="human", verbose=TRUE){
  species = tolower(species)
  dataset = switch(species,
          "mouse" = "mmusculus_gene_ensembl",
          "human" = "hsapiens_gene_ensembl",
          "hsapiens" = "hsapiens_gene_ensembl",
          "chicken" = "ggallus_gene_ensembl",
          "chick" = "ggallus_gene_ensembl",
          "rat" = "rnorvegicus_gene_ensembl",
           NA
  )
    require(biomaRt)
    if (is.na(dataset)) stop("Can't intepret 'species' parameter, should be human, mouse, rat, or chicken")

    if (exists("mart") & dataset == mart@dataset & verbose){
      print(paste("mart of", dataset, "exists in your workspace"))
      } else {
      mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset)
    }
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
