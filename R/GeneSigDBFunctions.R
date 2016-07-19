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
  
  lapply(seq_along(SigCols)), function(x){
    ids=sig[,x]
    if (x=="Clone ID")  parseCloneID(ids)
    if (x=="EnsEMBL ID") parseEnsEnsEMBL(ids)
    if (x=="GenBank ID") parseGenBankID(ids) 
    if (x=="Gene Symbol") parseGeneSymbol(ids) 
    if (x=="EntrezGene ID") parseEntrezGeneID(ids)
    if (x=="UniGene ID") parseUniGeneID(ids)
    if (x=="miRBase") parsemiRBase(ids)
    if (x=="Protein ID") parseProteinID(ids)
    if (x=="RefSeq ID") parseRedSeqID(ids)
    if (x=="Probe ID") parseProbeID(ids)
    if (x=="Secondary Probe ID") parseSecondaryProbeID(ids)
    if (x=="Gene Description") parseNULL(ids)
    if (x=="Other Gene Description") parseNULL(ids)
    if (x=="Geneset Specific Factor") parseNULL(ids)
    if (x=="Geneset Specific Statistics") parseNULL(ids)
    if (x=="Chromosome Map") parseNULL(ids)                                
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
    mapping<-getBM(attributes= c("embl", 	"entrezgene", "hgnc_symbol"), values=ids, filters="embl", mart=mart)
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
<<<<<<< HEAD


=======


>>>>>>> aedin/master
parseSecondaryProbeID<-function(ids, biomart=TRUE) {
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

parseGeneDescription<-function(ids, biomart=TRUE) {
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

parseOtherGeneDescription<-function(ids, biomart=TRUE) {
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

parseGenesetSpecificfactor<-function(ids, biomart=TRUE) {
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

parseGenesetSpecificStatistics<-function(ids, biomart=TRUE) {
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

parseChromosomeMap<-function(ids, biomart=TRUE) {
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
##############################
getMart<-function(ds="hsapiens_gene_ensembl"){
  #mart<-useMart(dataset="hsapiens_gene_ensembl", biomart="ensembl")
  #mart<-useMart(dataset="rnorvegicus_gene_ensembl", biomart="ensembl")  # maybe wrong
  #mart<-useMart(dataset="musculus_gene_ensembl", biomart="ensembl") # maybe wrong
  #listAttributes()
  require(biomaRt)
  #listDatasets()
  if (!exists("mart")) {
    mart<-useMart(dataset=ds, biomart="ensembl") 
  }
    return(mart)
  }
###############################
if (search="AnnotationDBI")  {
  # ALTERNATIVE
  # Using Annotation dbi
  require(org.Hs.eg.db)
  columns(org.Hs.eg.db)
  #Annotation DBi, you are selecting for Accession numbers, ensembl id, entrezid, symbol
  #keytype is Accession Number
  mapping<- AnnotationDbi::select(org.Hs.eg.db, keys = ids, columns= c("ACCNUM","ENSEMBL"  ,"ENTREZID", "SYMBOL" ), keytype = "ACCNUM")
  return(mapping)
  }
