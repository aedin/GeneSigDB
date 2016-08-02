#################################
## Parsing GeneSigDB R Functions
#################################


readGeneSigDBFile <- function(GeneSigDBPath="data",GeneSigDBFileName="GeneSigDB_Working.xlsx") {
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


#' Title
#'
#' @param SigID GeneSigDB id eg "10582678-Table1"
#' @param GeneSigIndex GeneSigDB data.frame obtained from load()
#' @param GeneSigDB_ReleaseData Where the raw data files exists, "../gene_signatures/data/"
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
getSig<- function(SigID,GeneSigIndex, GeneSigDB_ReleaseData="/Users/aliahmed/Documents/Stats/gene_signatures/data",...)
  {


# To call this function (executing it)
# getSig(SigID,GeneSigDB_ReleaseData=GeneSigDB_ReleaseData)
#

#GeneSigIndex is the GeneSigdb.xls file read using readGeneSigDBFile(GeneSigDBPath,GeneSigDBFileName)
# It turns a data.frame of the signature

  rowInd= GeneSigIndex$SigID%in%SigID
  #print(table(rowInd))

  SigFilePath=file.path(GeneSigDB_ReleaseData,GeneSigIndex$Release[rowInd], GeneSigIndex$PMID[rowInd], GeneSigIndex$FileAssociated[rowInd])
  print(SigFilePath)


  if(file.exists(SigFilePath)) {
    #Sig<- read.table(SigFilePath, header=TRUE, sep="\t", as.is=TRUE, comment.char="", check.names = FALSE)
    Sig <- as.data.frame(readr::read_tsv(SigFilePath))
    return(Sig)
  } else print(paste("Can't read in", SigID, "file in", SigFilePath))
  #check Data
}



#' Title Extract List of Cols from Sig
# GeneSigIndex is the GeneSigdb.xls file read using readGeneSigDBFile(GeneSigDBPath,GeneSigDBFileName)
#'
#' @param sigID a GeneSigDB ID, eg "10821843-Table1"
#' @param GeneSigIndex
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#' parseSigCols("10821843-Table1",GeneSigIndex)
parseSigCols<-function(sigID, GeneSigIndex, mappingColsOnly=TRUE, verbose=FALSE,...) {
  #mapIDs=c("Clone ID", "EnsEMBL ID", "EntrezGene ID","GenBank ID", "Gene Symbol", "Probe ID", "Protein ID", "RefSeq ID", "UniGene ID")
  mapIDs=c( "EnsEMBL ID", "EntrezGene ID","GenBank ID", "Gene Symbol")

  rowInd= GeneSigIndex$SigID%in%sigID
  colInd=grep("^Column\\d+", colnames(GeneSigIndex))
  sigCols= GeneSigIndex[rowInd, colInd]
  sigCols<- sigCols[!c(is.na(sigCols)| sigCols=="")]

  if(mappingColsOnly) {sigCols=sigCols[sigCols%in%mapIDs]}

  if (length(sigCols)<1) {
    print("No Mappable Cols Found")
    return(NULL)
  } else {
    sigInd= as.numeric(sub("Column", "", names(sigCols)))
    names(sigInd)= sigCols
    return(sigInd)}
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
parseIDs<-function(ids=Sig[,2], identifer=SigCols[2],attributes=c("ensembl_gene_id", "hgnc_symbol", "entrezgene"), species="human",...) {
    identifer<-as.character(identifer)
    print(identifer)
    mapIds<-switch(identifer,
       "Clone ID"=parseCloneID(ids),
       "EnsEMBL ID"= getBMall(values=ids, filters="ensembl_gene_id", attributes,species=species),
       "GenBank ID"= getBMall(values=ids, filters="embl", attributes, species=species),
       "Gene Symbol"= parseGeneSymbol(values=ids, filters="hgnc_symbol", attributes=attributes, species=species),
       "EntrezGene ID"= getBMall(values=ids, filters="entrezgene", attributes, species=species),
       "UniGene ID"= getBMall(values=ids, filters="unigene", attributes, species=species),
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

#' Title
#'
#' @param values Values of the filter, e.g. vector of affy IDs. If multiple filters are specified then the argument should be a list of vectors of which the position of each vector corresponds to the position of the filters in the filters argument.
#' @param filters Filters (one or more) that should be used in the query. A possible list of filters can be retrieved using the function listFilters.
#' @param attributes Attributes you want to retrieve. A possible list of attributes can be retrieved using the function listAttributes.
#' @param mart object of class Mart, created with the useMart function.
#' @param curl An optional 'CURLHandle' object, that can be used to speed up getBM when used in a loop.
#' @param checkFilters Sometimes attributes where a value needs to be specified, for example upstream\_flank with value 20 for obtaining upstream sequence flank regions of length 20bp, are treated as filters in BioMarts. To enable such a query to work, one must specify the attribute as a filter and set checkFilters = FALSE for the query to work.
#' @param verbose When using biomaRt in webservice mode and setting verbose to TRUE, the XML query to the webservice will be printed.
#' @param uniqueRows if the result of a query contains multiple identical rows, setting this argument to TRUE (default) will result in deleting the duplicated rows in the query result at the server side.
#' @param bmHeader Boolean to indicate if the result retrieved from the BioMart server should include the data headers or not, defaults to FALSE. This should only be switched on if the default behavior results in errors, setting to on might still be able to retrieve your data in that case
#' @param species Use to obtain mart, human, mouse, rat
#'
#' @return
#' @export
#'
#' @examples
#' ids=c("AA486072" ,"J02931" , "AA478436")
#' getBMall(values=ids, filters= "embl", attributes=c("ensembl_gene_id", "embl", "entrezgene"))
#'
getBMall <- function(values = ids,filters = '', attributes=c("ensembl_gene_id", "hgnc_symbol","entrezgene"),  mart, curl = NULL, checkFilters = TRUE, verbose = FALSE, uniqueRows = TRUE, bmHeader = FALSE, species="human") {

  mart<-getMart(species)
  attributes <- unique(c(c("ensembl_gene_id", "hgnc_symbol","entrezgene"), filters[1]))
  idMatch <- getBM(attributes, c(filters, "transcript_appris"), list(values,TRUE), mart, curl, checkFilters, verbose, uniqueRows, bmHeader)
  if(verbose) print(idMatch)
  print(idMatch)
  print(values)
  x <- as.data.frame(values, stringsAsFactors =FALSE)
  colnames(x) <- filters
  if (verbose) print(x)
  structure(dplyr::left_join(x = x, y = idMatch, by = filters[1], copy=TRUE))
}

parseCloneID<-function() {
   # 1. Check format looks correct
   # 2. search
  return("CloneID")
}


parseGenBankID<-function(ids, biomart=TRUE, species="human") {
  print("parseGenBankID")
  if(biomart){
    #selecting for ("ensembl_gene_id","embl","entrezgene", "hgnc_symbol")
    mart<-getMart(verbose=FALSE)
    #mapping<-biomaRt::select(mart, keys=ids, columns=c("ensembl_gene_id","embl","entrezgene", "hgnc_symbol"), keytype="embl")
    mapping<-getBMall(attributes= c("embl", 	"entrezgene", "hgnc_symbol"), values=ids, filters="embl", mart=mart)
    return(mapping)
  }
  }


parseGeneSymbol<-function(values=ids, filters="hgnc_symbol", ...){
   sym = data.frame(cbind(orig_sym=values, hgnc_symbol=AnnotationDbi::mapIds(org.Hs.eg.db, keys = values, column= "SYMBOL", keytype = "ALIAS")))
   mm<- getBMall(values=sym$hgnc_symbol, filters="hgnc_symbol", attributes, species)
   structure(dplyr::left_join(x = sym, y = mm, by = filters, copy=TRUE))
}






parseProteinID<-function(ids, biomart=TRUE) {
  # validateIDs  - uniprot_sptrembl","uniprot_swissprot", "uniprot_genename, "protein_id"

  if(biomart){
    #we are selecting for "ensembl_gene_id", "embl", "entrezgene", "hgnc_symbol")
    mart<-getMart()
    #mapping<-biomaRt::select(mart, keys=ids, columsn=c("ensembl_gene_id","embl","entrezgene", "hgnc_symbol"), keytype="embl")
    mapping<-getBMall(attributes= c("embl", 	"entrezgene", "hgnc_symbol"), values=ids,
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
    mapping<-getBMall(attributes= c("embl", 	"entrezgene", "hgnc_symbol"), values=ids,
                   filters="embl", mart=mart)
    return(mapping)
  }
}

parsemiRBase<-function(ids, biomart=TRUE) {
# microRNA
# "mirbase_accession"    ,  "mirbase_id"   ,"mirbase_transcript_name"
  if(biomart){
    #we are selecting for "ensembl_gene_id", "embl", "entrezgene", "hgnc_symbol")
    mart<-getMart()
    #mapping<-biomaRt::select(mart, keys=ids, columsn=c("ensembl_gene_id","embl","entrezgene", "hgnc_symbol"), keytype="embl")
    mapping<-getBMall(attributes= c("embl", 	"entrezgene", "hgnc_symbol"), values=ids,
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
    mapping<-getBMall(attributes= c("embl", 	"entrezgene", "hgnc_symbol"), values=ids,
                   filters="embl", mart=mart)
    return(mapping)
  }
  }


.parseNULL<-function (ids) {
  return(NULL)
}


mergeID <- function(mylist,n=NULL){
    if(!is.null(n)) {
      x<-sapply(mylist, function(x)x[,n])
    }
 #  print(x)

    res <- apply(x, 1, function(rr){
        rr<-unique(rr)
        rr<-rr[!is.na(rr)]
        if (length(rr)>1) paste(rr, collapse = "|")
        return(rr)})


  return(res)
}


mapSig<-function(sigID,GeneSigIndex,GeneSigDB_ReleaseData, attributes=c("ensembl_gene_id", "hgnc_symbol", "entrezgene"), species="human", verbose=FALSE, ...){

  sig <- getSig(sigID,GeneSigIndex,GeneSigDB_ReleaseData)
  if(verbose) print(head(sig))
  sigCols <- parseSigCols(sigID, GeneSigIndex,verbose=verbose)
  if(verbose) print(sigCols)

  if(is.null(sigCols)) warning(paste("No mappable Columns in ", sigID))

  if(!is.null(sigCols)){
  # Map Genes
  xx<-lapply(seq_along(sigCols) ,function(i) parseIDs(ids = sig[, sigCols[i]], identifer = names(sigCols)[i], attributes=attributes,verbose=TRUE))


  cols=unique(unlist(sapply(xx, colnames)))
  mapped_orig<-do.call(cbind,lapply(xx, function(g) g[,colnames(g)[!colnames(g)%in%attributes], drop=FALSE]))
#  print(mapped_orig)

  colnames(mapped_orig) = paste0(colnames(sig)[sigCols], "_orig")
  mapped_IDs<-sapply(attributes, function(i) mergeID( xx, i))
  #print(mapped_IDs)
  mapping<-cbind( mapped_orig,mapped_IDs)

  return(mapping)
  }

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
                     "Error"
    )
    if (dataset=="Error") {
        print(paste("Can't intepret", species,   "parameter, should be human, mouse, rat, or chicken"))
        if (exists("mart"))  return(mart)
    } 
    
    #print(dataset)
    
    require(biomaRt)

    if (exists("mart")) {
        if (inherits(mart, "Mart")){
            if (dataset == mart@dataset & verbose){
                print(paste("mart of", dataset, "exists in your workspace"))
                return(mart)
            }
        }
    } 
    
    if(!dataset=="Error") {
            print(paste("Obtaining mart", dataset))
            mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset)
            return(mart)
        }
    
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
