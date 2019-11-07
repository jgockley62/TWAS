library(readr)
library(BioMart)

synapseclient <- reticulate::import("synapseclient")
syn_temp <- synapseclient$Synapse()
syn_temp$login()

#countz <- readr::read_tsv(syn_temp$get(foo, version = 1))
COUNT_OBJ <- syn_temp$get("syn21075911", version = 1)
countz <- read.table(COUNT_OBJ$path, header=T, sep='\t', check.names = F)

myGetGeneLengthAndGCContent <- function(id){
  
  org = 'hsa'
  id.type = 'ensembl_gene_id'
  host = 'sep2019.archive.ensembl.org'
  
  message("Connecting to BioMart ...")
  ensembl <- useMart("ENSEMBL_MART_ENSEMBL", host = host)
  ds <- listDatasets(ensembl)[, "dataset"]
  ds <- grep(paste0("^", org), ds, value = TRUE)
  
  if (length(ds) == 0){
    stop(paste("Mart not found for:", org))
  } else if (length(ds) > 1) {
    message("Found several marts")
    sapply(ds, function(d) message(paste(which(ds == d), d, sep = ": ")))
    n <- readline(paste0("Choose mart (1-", length(ds), ") : "))
    ds <- ds[as.integer(n)]
  }
  
  ensembl <- useDataset(ds, mart = ensembl)
  message(paste0("Downloading sequence", ifelse(length(id) > 1, "s", ""), " ..."))
  if (length(id) > 100)
    message("This may take a few minutes ...")
  attrs <- c(id.type, "ensembl_exon_id", "chromosome_name", "exon_chrom_start", "exon_chrom_end")
  coords <- getBM(filters = id.type, attributes = attrs, values = id, mart = ensembl)
  id <- unique(coords[, id.type])
  coords <- sapply(id, function(i) {
    i.coords <- coords[coords[, 1] == i, 3:5]
    g <- GRanges(i.coords[, 1], IRanges(i.coords[, 2], i.coords[, 3]))
    return(g)
  })
  coords <- lapply(coords[id], reduce)
  len <- plyr::ldply(coords, function(x) sum(IRanges::width(x)), .id = 'ensembl_gene_id') %>%
    dplyr::rename(gne.length = V1)
  
  gc.cont <- getBM(filters = id.type, attributes = c(id.type, 'version', 'hgnc_symbol', 'percentage_gene_gc_content'), values = id, mart = ensembl)
  
  res <- dplyr::full_join( gc.cont, len )
  return(res)
}

ENSGs <- do.call(rbind, strsplit(as.character(countz$feature)[5:dim(countz)[1]], "[.]"))

temp <- myGetGeneLengthAndGCContent( ENSGs )
Temp <- as.data.frame( cbind( paste0(temp$ensembl_gene_id, '.', temp$version), temp[,1:5]))
colnames(Temp) <- c("Gene.ID", "ensembl_gene_id", "position", "hgnc_symbol", "percentage_gc_content", "gene.length")

parentId = 'syn21136908';
activityName = 'Gene Lengths and GC content hg38.p12';
activityDescription = 'Gene Lengths and GC Content Based off of the Sept. 2019 ENTREZ/ENSEMBL Update';
thisFileName <- 'GeneLengthQuant.R'
# Github link
thisRepo <- githubr::getRepo(repository = "jgockley62/TWAS", ref="branch", refName='master')
thisFile <- githubr::getPermlink(repository = thisRepo, repositoryPath=paste0('code/',thisFileName))

### Store files in synapse
activityName = 'Gene Lengths and GC content hg38.p12';
activityDescription = 'Gene Lengths and GC Content Based off of the Sept. 2019 ENTREZ/ENSEMBL Update'

CODE <- syn_temp$store(synapseclient$Folder(name = "GeneStats", parentId = parentId))

Syns_Used <- NULL
# Set annotations
all.annotations = list(
  dataType = 'mRNA',
  dataSubType = 'refencedata',
  summaryLevel = 'gene',
  organism = 'HomoSapiens',
  genomeAssemblyID = 'GRCh38.p12'
)
# Store Length and GC results
write.table(Temp, file = 'GeneLengthsGCcontent.tsv', sep = '\t', quote=F)
ENRICH_OBJ <-  syn_temp$store( synapseclient$File( path='GeneLengthsGCcontent.tsv', name = 'Gene Lengths and GC content hg38.p12', parentId=CODE$properties$id ), used = Syns_Used, activityName = activityName, executed = thisFile, activityDescription = activityDescription)
all.annotations$dataSubType = 'ReferenceData'
syn_temp$setAnnotations(ENRICH_OBJ, annotations = all.annotations)




###########
temp <- myGetGeneLengthAndGCContent( ENSGs )
Temp <- as.data.frame( cbind( paste0(temp$ensembl_gene_id, '.', temp$version), temp[,1:5]))
colnames(Temp) <- c("Gene.ID", "ensembl_gene_id", "position", "hgnc_symbol", "percentage_gc_content", "gene.length")

table( biomart$ensembl_gene_id %in% Temp$ensembl_gene_id)
table( countz$feature %in% Temp$Gene.ID)

##################
org = 'hsa'
id.type = 'ensembl_gene_id'
host = 'jul2019.archive.ensembl.org'
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", host = host)
NewVersion <- getBM( attributes = c(id.type, 'ensembl_gene_id', 'ensembl_gene_id_version', 'hgnc_symbol'), values = id, mart = ensembl)

table(as.character(countz$feature)[ 5:dim(countz)[1] ] %in% NewVersion$ensembl_gene_id_version)

Ver_ENSGs <- countz$feature[ as.character(countz$feature) %in% c("N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous" ) == F ]
ENSGs <- do.call( rbind, strsplit( as.character(countz$feature), "[.]"))[,1][ do.call( rbind, strsplit( as.character(countz$feature), "[.]"))[,1] %in% c("N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous" ) == F ]
Pos <- do.call( rbind, strsplit( as.character(countz$feature), "[.]"))[,2][ do.call( rbind, strsplit( as.character(countz$feature), "[.]"))[,2] %in% c("N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous" ) == F ]

temp <- myGetGeneLengthAndGCContent( ENSGs )




Tab <- as.data.frame( matrix( NA, length(ENSGs), 6) )
colnames(Tab) <- c( "Gene.ID", "ensembl_gene_id", "position", "hgnc_symbol", "percentage_gc_content", "gene.length" )
Tab$Gene.ID <- Ver_ENSGs
Tab$ensembl_gene_id <- ENSGs
Tab$position <- Pos

tab <- as.data.frame( matrix( NA, 0, 6) )
colnames(tab) <- c( "Gene.ID", "ensembl_gene_id", "position", "hgnc_symbol", "percentage_gc_content", "gene.length" )


temp <- myGetGeneLengthAndGCContent( Tab$ensembl_gene_id )
row.names(temp) <- temp$ensembl_gene_id

Temp <- temp[ !duplicated(temp$ensembl_gene_id), ]
row.names(Temp) <- Temp$ensembl_gene_id


row.names( biomart ) <- as.character( biomart$ensembl_gene_id )
biomart$hgnc_symbol <- as.character(biomart$hgnc_symbol)
biomart$percentage_gc_content <- as.numeric(as.character(biomart$percentage_gc_content))
biomart$gene.length <- as.numeric(as.character(biomart$gene.length))

BIO <- biomart[ row.names(Temp), ]



