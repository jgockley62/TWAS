# Load Libraries and se working directory 
##  devtools::install_github('blogsdon/spike/spike/')
library(spike)
library(pheatmap)
library(parallel)
library(doParallel)
library(spike)
library(bcv)
library(reshape2)
library(biomaRt)
library(UpSetR)
library(ComplexHeatmap)
#install.packages("enrichR")
library(enrichR)
setwd( config$filedir )


reticulate::use_python("/usr/bin/python", required = TRUE)
synapseclient <- reticulate::import("synapseclient")
syn_temp <- synapseclient$Synapse()
syn_temp$login()

setwd("~/TWAS/code/")
source("~/TWAS/utilityFunctions/knitfile2synapseClient.R")
source("~/TWAS/utilityFunctions/hook_synapseMdSyntax_plot.R")
#createAndKnitToFolderEntityClient(file = "CMC_ACC_Exp_DataProcessing_V2.Rmd",
#                                  parentId ="syn21478397",
#                                  folderName = 'ACC')

Genes <- c( "APOC1", "CD2AP", "EED", "CEACAM19", "MTCH2", "TREM2", "CLPTM1", "KNOP1" )
ENSGS <- c( "ENSG00000130208", "ENSG00000198087", "ENSG00000074266", "ENSG00000186567", "ENSG00000109919", "ENSG00000095970", "ENSG00000104853", "ENSG00000103550" )
names(ENSGS) <- Genes
# Set Github Provenance links
thisRepo <- githubr::getRepo(repository = "jgockley62/TWAS", ref="branch", refName='master')

CODE <- syn_temp$store(synapseclient$Folder(name = 'Partial Correlations', parentId = 'syn18936948'))

#SynIDs of expression data to pull from
ExpressionDS <- c('syn21291908','syn21292041','syn21285564','syn21285564','syn21285564','syn21285564','syn21291908')
names( ExpressionDS ) <- c('CBE', 'DLPFC', 'FP', 'IFG', 'PHG', 'STG', 'TCX')

#Study ID Translator
Study <- c( 'RosMap', 'Mayo', 'Mayo', 'MSBB', 'MSBB', 'MSBB', 'MSBB')
names(Study) <- c('DLPFC', 'TCX', 'CBE', 'FP', 'IFG', 'PHG', 'STG')

Syn <- list('DLPFC', 'TCX', 'CBE', 'FP', 'IFG', 'PHG', 'STG')

GeneData <- function( id, type ){
  #Retreives gene name info
  #'@id the gene ids to source a list from
  #'@type either ENSG or GNAME for type of gen ids
  
  org = 'hsa'
  ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
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
  
  if(type == 'ENSG'){
    id.type = 'ensembl_gene_id'
    ensembl <- useDataset(ds, mart = ensembl)
    Genes <- getBM(filters = id.type, attributes = c(id.type, 'hgnc_symbol', 'chromosome_name', 'start_position', 'end_position'), values = id, mart = ensembl)
  }else{
    if(type == 'GNAME'){
      id.type = 'hgnc_symbol'
      ensembl <- useDataset(ds, mart = ensembl)
      Genes <- getBM(filters = id.type, attributes = c('ensembl_gene_id', id.type, 'chromosome_name', 'start_position', 'end_position'), values = id, mart = ensembl)
      ##Cut Decmals
    }else{
      #Error Gene names not specififed correctly
      stop("ERROR: SOURCE=Config.yaml Issue=genelisttype must be ENSG or GNAME")
    }
  }
  return(Genes)
}


#Tab <- GeneData( Genes, 'GNAME' )
Tab <- GeneData( ENSGS, 'ENSG' )
Tab <- Tab[ Tab$chromosome_name %in% c(c(1:22),'X','Y'), ]
#Trans used as raw for module 7
Trans <- Tab
Tab$Coord_hg19 <- paste0("chr", Tab$chromosome_name, ":", Tab$start_position, "-",Tab$end_position)
Tab$Interval_hg19 <- NA

for( i in 1:dim(Tab)[1]){
  if(  Tab$start_position[i] < Tab$ensembl_gene_id[i] ){
    Tab$Interval_hg19[i] <- paste0("chr", Tab$chromosome_name[i], ":", as.numeric(Tab$start_position[i])-1e6, "-", as.numeric(Tab$end_position[i])+1e6)
  }else{
    Tab$Interval_hg19[i] <- paste0("chr", Tab$chromosome_name[i], ":", as.numeric(Tab$start_position[i])+1e6, "-", as.numeric(Tab$end_position[i])-1e6)
  }
}

#Useable Chrs for genelist ENSG  translations and wiki analysis preamble
row.names(Trans) <- Trans$ensembl_gene_id
Trans <- Trans[ Trans$chromosome_name %in% c(1:22,"X","Y"), ]

# Make objects for list of Genes Not Expressed In Tissue and dataframes to plot
Missing <- list()
Plots <- list()

FinalTis <- data.frame()
cores <- detectCores()-2 
#cl <- makePSOCKcluster(cores)
#registerDoParallel(cl)
#stopCluster(cl)

for( Tissue in names(Study) ){
  #Tissue<-'CBE'
  Syns_Used <- NULL
  #Tissue <- 'STG'
  
  #Load expression for tissue
  exp <- read.table(syn_temp$get(as.character(ExpressionDS[Tissue]))$path, header =T, sep ='\t', row.names=1)
  colnames( exp ) <- gsub( "X", "", colnames( exp ) )
  Syns_Used <- c(Syns_Used, as.character(ExpressionDS[Tissue] ) )
  
  #Seperate exp by tissue
  if( Tissue == 'DLPFC'){
  }else{
    if( Tissue == 'TCX' | Tissue == 'CBE' ){
      if( Tissue == 'CBE' ){
        slec <- 'CER'
      }else{ slec <- 'TCX' }
      colnames( exp ) <- gsub( "TC", "TCX", colnames( exp ) )
      exp <- exp[ , grepl( slec, colnames(exp)) ]
    }else{
      if( Tissue %in% c('FP', 'IFG', 'PHG', 'STG') ){
        Syns_Used <- c( Syns_Used, 'syn21285520' )
        Meta <- read.table( syn_temp$get('syn21285520')$path, header =T, sep ='\t', row.names=1 )
        exp <- exp[ colnames(exp) %in% row.names(Meta[grepl( Tissue, Meta$Tissue.Diagnosis),])]
      }else{
        stop(paste0("ERROR: SOURCE=Config.yaml Issue=Tissue: ", Tissue," is improper must be one of: CBE, DLPFC, FP, IFG, PHG, STG, TCX"))
      }
    }
  }
  
  #Impute svalues for given gene-patient NA values
  exp <- exp[ rowSums(is.na(exp)) < ncol(exp), ]
  
  foo <- bcv::impute.svd( t(exp) )
  Exp <- foo$x
  row.names(Exp) <- row.names(t(exp)) 
  colnames(Exp) <- colnames(t(exp)) 
  
  #Record Genes missing in tissue of interest
  if( as.numeric(table(row.names(Trans) %in% colnames(Exp) )["TRUE"] ) == dim(Trans)[1] ){
    Missing[[ Tissue ]] <- paste0("All Listed Genes Expressed in ", Tissue)
  }else{
    Missing[[ Tissue ]] <- paste0( "Queried Genes missing from ", Tissue, " : ",
                                   paste(c(Trans[ (row.names(Trans) %in% colnames(Exp) ) == F, ]$hgnc_symbol), collapse = ', ')
    )
  }

  #Pairwise spearman correlation of genes
  XCor <- cor(Exp, method = "spearman")
  Final <- data.frame()
  
  source("~/TWAS/code/utilityFunctions/Parallel_vbsrBootstrap.R")
  
  #Run Partial correlations for each gene
  RUNNe <- function( i=i, Exp=Exp ){ 
    source("~/TWAS/code/utilityFunctions/Parallel_vbsrBootstrap.R")
    OBS <- i
    y <- as.matrix(Exp[,OBS]) 
    colnames(y) <- i
    X <- Exp[,(colnames(Exp) %in% OBS) == F ]
    att <- pvbsrBootstrap( y=y, x=X, nsamp=100, cores=14 )
    names( att ) <- gsub( "intercept", OBS, names(att) )
    att <- c( 'SeedGene' = OBS, att[ colnames(Exp) ]) 
    return( att[c( 'SeedGene', colnames(Exp) )] )
  }
  
  LIST <- Trans$ensembl_gene_id
  foo=data.frame(matrix( 0, length(LIST), (dim(Exp)[2]+1) ))
  colnames(foo) <- c( 'SeedGene', colnames(Exp) )
  
  mark <- Sys.time()
  for( i in 1:length(LIST)){
    foo[i,] <- t( RUNNe( i=LIST[i], Exp=Exp) )
  }
  Sys.time()-mark
  
  row.names( foo ) <-  foo[,1]
  foo <-  foo[, colnames(foo)[ (colnames(foo) %in% "SeedGene") == F ]]
  write.table( foo, file=paste0(Tissue,"_PartialCorMatrix.tsv"), row.names=T, col.names = T, quote=F, sep='\t' )
  DT <- melt(t(foo))
  colnames(DT) <- c( "Hit Gene", "Target Gene", "PartialCor") 
  write.table(DT, file=paste0(Tissue,"_PartialCorTable.tsv"), row.names=T, col.names = T, quote=F, sep='\t' )
}


#Trans
TabOne <- data.frame(data.frame(matrix(0, 0, 4)))
colnames(TabOne) <- c( "Cutoff", "Value", "Tissue", "Gene")

for( Tissue in names(Study) ){
  PCs <- as.data.frame( read.table(paste0(Tissue,"_PartialCorTable.tsv"), sep = '\t', header = T) )
  PCs$Tissue <- Tissue
  
  for( name in names(table(PCs$Target.Gene)) ){
  
    Perc <- data.frame(matrix(0, 101, 4))
    Perc[,1] <- seq(0,100,1)
    Perc[,3] <- Tissue
    Perc[,4] <- Trans[ name, ]$hgnc_symbol
    colnames(Perc) <- c( "Cutoff", "Value", "Tissue", "Gene")
    for( i in 0:100 ){
      amt <- dim( PCs[ PCs$Target.Gene == name & PCs$PartialCor  >= i/100 , ] )[1]
      Rat <- amt/dim( PCs[ PCs$Target.Gene == name ,])[1]
      Perc[i+1,]$Value <- Rat 
    }
    TabOne <- as.data.frame( rbind(TabOne, Perc))
  }
}

names(table(TabOne$Gene))
Plots <- list()
for( Gene in names(table(TabOne$Gene)) ){
  
  Plots[[Gene]] <- ggplot(TabOne[ TabOne$Gene == Gene & TabOne$Tissue != 'CBE', ], aes(colour=Tissue, x = Cutoff, y = Value)) + 
                      geom_line( ) + ggtitle( paste0( 'Target Gene = ',Gene, " Cuttoff = ", 2,'%') ) + xlab('Partial Correlation BootStrap Percent') + 
                      ylab('Percent Assiations >= BootStrap Percent') + theme(plot.title = element_text(hjust = 0.5)) + 
                      geom_vline( xintercept = 2 )
}

TabTwo <- data.frame(data.frame(matrix(0, 0, 4)))
colnames(TabTwo) <- c( "Cutoff", "Value", "Tissue", "Gene")

Genes <- list()
for( name in names(ENSGS) ){
  TabThree <- data.frame(data.frame(matrix(0, 0, 4)))
  colnames(TabThree) <- c( "Cutoff", "Value", "Tissue", "Gene")
  for( Tissue in names(Study) ){
    PCs <- as.data.frame( read.table(paste0(Tissue,"_PartialCorTable.tsv"), sep = '\t', header = T) )
    PCs$Tissue <- Tissue

    message( paste0( Tissue, ' ', dim(PCs[ PCs$Target.Gene == ENSGS[name] & PCs$PartialCor >= 0.02, ])[1]) )
    TabTwo <- as.data.frame( rbind(TabTwo,PCs[ PCs$Target.Gene == ENSGS[name] & PCs$PartialCor >= 0.02, ] ) )
    
    TabThree <- rbind(TabThree,PCs[ PCs$Target.Gene == ENSGS[name] & PCs$PartialCor >= 0.02, ])
  }
  TabThree <- TabThree[ (as.character( TabThree$Tissue ) %in% 'CBE') == F & as.character(TabThree$Target.Gene) == ENSGS[name],] 
  
  Genes[[name]] <- names( table(as.character(TabThree$Hit.Gene))[ table(as.character(TabThree$Hit.Gene)) > 1 ] )
}

All <- c(Genes$APOC1, Genes$CD2AP, Genes$EED, Genes$CEACAM19, Genes$MTCH2, Genes$TREM2, Genes$CLPTM1, Genes$KNOP1)

table(as.numeric(table(All)))

ALL <- as.data.frame(All)
ALL$Target<- c(rep('APOC1', length(Genes$APOC1)),
               rep('CD2AP', length(Genes$CD2AP)),
               rep('EED', length(Genes$EED)),
               rep('CEACAM19', length(Genes$CEACAM19)),
               rep('MTCH2', length(Genes$MTCH2)),
               rep('TREM2', length(Genes$TREM2)),
               rep('CLPTM1', length(Genes$CLPTM1)),
               rep('KNOP1', length(Genes$KNOP1))
)
colnames(ALL) <- c("Hit", "Target")

m1 = make_comb_mat(Genes)
UpSet(m1)

dbs <- listEnrichrDbs()
if (is.null(dbs)) websiteLive <- FALSE

dbs <- c("GO_Biological_Process_2018",
         "GO_Molecular_Function_2018",
         "KEGG_2019_Human",
         "WikiPathways_2019_Human",
         "GO_Cellular_Component_2018")

lapply(Genes, length)
# APOC1 - 114
# CD2AP - 79
# EED - 67
# CEACAM19 - 73
# MTCH2 - 102
# TREM2 - 104
# CLPTM1 - 122
# KNOP1 - 57

Gname <- list()
Gname$APOC1 <- GeneData( Genes$APOC1, 'ENSG' )$hgnc_symbol
Gname$CD2AP <- GeneData( Genes$CD2AP, 'ENSG' )$hgnc_symbol
Gname$EED <- GeneData( Genes$EED, 'ENSG' )$hgnc_symbol
Gname$CEACAM19 <- GeneData( Genes$CEACAM19, 'ENSG' )$hgnc_symbol
Gname$MTCH2 <- GeneData( Genes$MTCH2, 'ENSG' )$hgnc_symbol
Gname$TREM2 <- GeneData( Genes$TREM2, 'ENSG' )$hgnc_symbol
Gname$CLPTM1 <- GeneData( Genes$CLPTM1, 'ENSG' )$hgnc_symbol
Gname$KNOP1 <- GeneData( Genes$KNOP1, 'ENSG' )$hgnc_symbol

Gname$APOC1 <- Gname$APOC1[ (Gname$APOC1 %in% "")==F]
Gname$CD2AP <- Gname$CD2AP[ (Gname$CD2AP %in% "")==F]
Gname$EED <- Gname$EED[ (Gname$EED %in% "")==F]
Gname$CEACAM19 <- Gname$CEACAM19[ (Gname$CEACAM19 %in% "")==F]
Gname$MTCH2 <- Gname$MTCH2[ (Gname$MTCH2 %in% "")==F]
Gname$TREM2 <- Gname$TREM2[ (Gname$TREM2 %in% "")==F]
Gname$CLPTM1 <- Gname$CLPTM1[ (Gname$CLPTM1 %in% "")==F]
Gname$KNOP1 <- Gname$KNOP1[ (Gname$KNOP1 %in% "")==F]

lapply(Gname, length)
# APOC1 - 106
# CD2AP - 77
# EED - 58
# CEACAM19 - 65
# MTCH2 - 98
# TREM2 - 101
# CLPTM1 - 118
# KNOP1 - 46

enriched_APOC1 <- enrichr(Gname$APOC1 , dbs)
enriched_Trem <- enrichr(Gname$TREM2 , dbs)

comb <- c( Gname$APOC1,Gname$TREM2 )
comb <- comb[!duplicated(comb)]
enriched_comb <- enrichr(comb , dbs)  

enriched_CD2AP <- enrichr(Gname$CD2AP, dbs)  
enriched_EED <- enrichr(Gname$EED, dbs)  
enriched_CEACAM19 <- enrichr(Gname$CEACAM19, dbs)  
enriched_MTCH2 <- enrichr(Gname$MTCH2, dbs)  
enriched_CLPTM1 <- enrichr(Gname$CLPTM1, dbs)  
enriched_KNOP1 <- enrichr(Gname$KNOP1, dbs)  

Total <- c(Gname$APOC1, Gname$TREM2, Gname$CD2AP, Gname$EED,Gname$CEACAM19,Gname$MTCH2,Gname$CLPTM1,Gname$KNOP1)
Total<-Total[!duplicated(Total)]
total <-enrichr(Total, dbs)

Other <- c( Gname$CD2AP, Gname$EED,Gname$CEACAM19,Gname$MTCH2,Gname$CLPTM1,Gname$KNOP1)
Other<-Other[!duplicated(Other)]
other <-enrichr(Other, dbs)


PrintR <- function( MOD){
  print((as.data.frame( MOD$GO_Biological_Process_2018[ MOD$GO_Biological_Process_2018$Adjusted.P.value < 0.05, c("Term", "Overlap", "Adjusted.P.value", "Odds.Ratio") ])))
  print((as.data.frame(MOD$GO_Molecular_Function_2018[ MOD$GO_Molecular_Function_2018$Adjusted.P.value < 0.05, c("Term", "Overlap", "Adjusted.P.value", "Odds.Ratio") ])))
  print((as.data.frame(MOD$KEGG_2019_Human[ MOD$KEGG_2019_Human$Adjusted.P.value < 0.05, c("Term", "Overlap", "Adjusted.P.value", "Odds.Ratio") ])))
  print((as.data.frame(MOD$WikiPathways_2019_Human[ MOD$WikiPathways_2019_Human$Adjusted.P.value < 0.05, c("Term", "Overlap", "Adjusted.P.value", "Odds.Ratio") ])))
  print((as.data.frame(MOD$GO_Cellular_Component_2018[ MOD$GO_Cellular_Component_2018$Adjusted.P.value < 0.05, c("Term", "Overlap", "Adjusted.P.value", "Odds.Ratio") ])))
}

PrintR( enriched_APOC1 )
PrintR( enriched_Trem )

PrintR( enriched_comb )


PrintR( enriched_CD2AP )
PrintR( enriched_EED )
PrintR( enriched_CEACAM19 )
PrintR( enriched_MTCH2 )
PrintR( enriched_CLPTM1 )
PrintR( enriched_KNOP1 )

PrintR( total )
PrintR( other )

##################################################################
# - Look at Nieghbors
KNOP1 <- c("VPS35L", "CCP110", "IQCK", "GPRC5B", "GPR139", "GDE1", "TMC5", "CLEC19A", "SYT17", "ITPRIPL2")
#VPS35L==C16orf62
KNOP1_exp <- GeneData( KNOP1, 'GNAME' )

EED <- c("SYTL2", "CCDC89", "CREBZF", "TMEM126A", "CCDC83", "PICALM", "HIKESHI", "MIR6755", "CCDC81", "ME3", "PRSS23")
#HIKESHI==C11orf73
EED_exp <- GeneData( EED, 'GNAME' )

CD2AP <- c('LOC101926962', 'ADGRF5', 'ADGRF1', 'TNFRSF21', 'ADGRF2', 'ADGRF4', 'OPN5', 'PTCHD4')
#GPR116 - ADGRF5
#GPR110 - ADGRF1
#GPR115 - ADGRF4
#GPR111 - ADGRF2
#LOC101926962 - ADGRF5-AS1 
CD2AP_exp <- GeneData( CD2AP, 'GNAME' )



MTCH2 <- c('LRP4', 'SNORD67', 'MIR5582', 'CKAP5', 'F2', 'LRP4-AS',  'ZNF408', 'ARHGAP1', 'C11orf49', 'LRP4',
           'ARFGAP2', 'PACSIN3', 'MIR6745', 'DDB2', 'ACP2', 'NR1H3', 'SPI1', 'MYBPC3', 'MADD', 'RAPSN', 'SLC39A13',
           'PSMC3', 'CELF1', 'PTPMT1', 'KBTBD4', 'C1QTNF4', 'NDUFS3', 'FAM180B', 'NUP160', 'FNBP4', 'AGBL2', 'PTPRJ',
           'OR4C45', 'OR4C3', 'OR4S1', 'OR4X1', 'OR4X2', 'OR4B1', 'OR4A47')
MTCH2_exp <- GeneData( MTCH2, 'GNAME' )
MTCH2_exp <- MTCH2_exp[ MTCH2_exp$chromosome_name ==  11, ]

Rows <- sum( dim( MTCH2_exp )[1],
  dim( CD2AP_exp )[1],
  dim( EED_exp )[1],
  dim( KNOP1_exp )[1]
)

InTis <- as.data.frame(matrix('NULL', Rows, 6))
rownames( InTis ) <- c( MTCH2_exp$ensembl_gene_id, CD2AP_exp$ensembl_gene_id, EED_exp$ensembl_gene_id, KNOP1_exp$ensembl_gene_id )
colnames( InTis ) <- names(ExpressionDS)[ (names(ExpressionDS) == 'CBE') == F ]


for( Tissue in names(ExpressionDS)[ (names(ExpressionDS)== 'CBE') == F ] ){
  exp <- read.table(syn_temp$get(as.character(ExpressionDS[Tissue]))$path, header =T, sep ='\t', row.names=1)
  colnames( exp ) <- gsub( "X", "", colnames( exp ) )
  
  eval(parse( text=paste0( 'InTis$',Tissue ,' <- row.names(InTis) %in% row.names(exp)' )))
  
}


InTis$Tissue <- "No"

for( i  in 1:dim(InTis)[1] ){
  if( FALSE %in% InTis[i,] ){
    
  }else{
    InTis[i,]$Tissue <- 'Yes'
  }
}

table( InTis$Tissue )






















Genes <- list()
for( name in names(ENSGS) ){
  TabTwo <- as.data.frame( rbind(TabTwo,PCs[ PCs$Target.Gene == ENSGS[name] & PCs$PartialCor >= 0.02, ] ) ) 
  TabTwo <- TabTwo[ (as.character( TabTwo$Tissue ) %in% 'CBE') == F,] 
  Genes[[name]] <- names( table(as.character(TabTwo$Hit.Gene))[ table(as.character(TabTwo$Hit.Gene)) > 1 ] )
}

