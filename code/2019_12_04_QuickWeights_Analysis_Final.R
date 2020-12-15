require(foreign)
require(ggplot2)
require(MASS)
require(Hmisc)
require(reshape2)
library(aod)
library(foreach)
library(doMC)
library(rms)
library(foreign)
library(synapser)
library(rms)
library(aod)
library(car)
library(broom)
library(ordinal)
require(gridExtra)
require(data.table)
library(ggplot2)
library(gridExtra)
library(ggpubr)
#synapseclient <- reticulate::import("synapseclient")
#syn_temp <- synapseclient$Synapse()
#syn_temp$login()

#Chr22 <- read.table( file = "/Users/jgockley/Desktop/Projects/TWAS/", header = T, sep="\t", row.names = 1 )
Imputs <- read.table(file=synapser::synGet('syn21346614')$path, sep="", header=T)
row.names(Imputs) <- as.character(Imputs[,1])

New_Pute <- Imputs[ ,grepl("ENSG", colnames(Imputs) ) ==T] 
New_Pute <- New_Pute[ ,grepl("\\.", colnames(New_Pute) ) ==F] 

TotalImpute <- t(New_Pute)

TotalPheno <- Imputs[,1:5]
TotalExpres <- read.table( file=synapser::synGet('syn21346605')$path, header = T, sep = "\t" )
row.names(TotalExpres) <- TotalExpres$Gene_Symbol
colnames(TotalExpres) <- gsub( "X", "", colnames(TotalExpres))


MatchedExpress <- TotalExpres[ row.names(TotalImpute), colnames(TotalImpute)[ colnames(TotalImpute) %in% colnames(TotalExpres)] ]
MatchedImpute <- TotalImpute[ row.names(MatchedExpress),colnames(MatchedExpress) ]

table( row.names(MatchedExpress) %in% row.names(MatchedImpute) )
table( colnames(MatchedImpute) %in% colnames(MatchedExpress) )

foo <- cor.test( as.numeric(MatchedExpress[1,]), as.numeric(MatchedImpute[1,]) , method = "kendall")

Correlation <- matrix(0, dim(MatchedImpute)[1], 3)
row.names(Correlation) <- row.names(MatchedImpute)
colnames(Correlation) <- c( "Z", "Correlation", "P-Val")

for( Gene in row.names(MatchedImpute) ){
  foo <- cor.test( as.numeric(MatchedExpress[Gene,]), as.numeric(MatchedImpute[Gene,]) , method = "kendall")
  Correlation[ Gene, ] <- c( as.numeric(foo$statistic), as.numeric(foo$estimate), as.numeric(foo$p.value) )
}

Correlation <- cbind( Correlation, "Log10_Pval"= -log10(Correlation[,3]))
Correlation <- cbind( Correlation, "FDR"= p.adjust( Correlation[,3], method = "fdr", n = 13650) )
Correlation <- cbind( Correlation, "BH"= p.adjust( Correlation[,3], method = "BH", n = 13650) )
Correlation <- cbind( Correlation, "BF"= p.adjust( Correlation[,3], method = "bonferroni", n = 13650) )

plot( Correlation[,4], Correlation[,2], las =1, pch = 16, cex = .4, bty ='n',
      xlab = "-Log10( PVal )", ylab = "Correlation"
)
abline( v = -log10(0.05/13650))

plot( Correlation[,2],  jitter(-log10( Correlation[,5]), 5), las =1, pch = 16, cex = .4, bty ='n',
      xlab = "Correlation", ylab = "-Log10( PVal )"
)
abline( h = p.adjust( -log10(0.05/13650), method = "fdr", n = 13650) )
#abline( h = 5.466808 )


Correlation <-as.data.frame(Correlation)
hist(Correlation$Correlation, breaks = 100, las=1, main = "Distribution of Correlation", xlab = "Correlation")

BF_SigGenes <- row.names( Correlation[ Correlation$Log10_Pval > -log10(0.05/13650), ] )

SigCor_Impute <- TotalImpute[ BF_SigGenes, ]
Relative_Pheno <- TotalPheno[ TotalPheno$PHENO != -9, ]
Impute_to_test <- SigCor_Impute[ , row.names(Relative_Pheno) ]

foo <- t.test( Impute_to_test[ 'ENSG00000001461', row.names(Relative_Pheno[ Relative_Pheno$PHENO ==1,]) ], Impute_to_test[ 'ENSG00000001461', row.names(Relative_Pheno[ Relative_Pheno$PHENO ==2,]) ]  )
foo$p.value
foo$statistic
foo$conf.int[1]
foo$conf.int[2]

Temp <- matrix(0, length(BF_SigGenes), 4)
row.names(Temp) <- BF_SigGenes
colnames(Temp) <- c("Statistic", "PVal", "ConfA", "ConfB" )

for( name in BF_SigGenes ){
  foo <- t.test( Impute_to_test[ name, row.names(Relative_Pheno[ Relative_Pheno$PHENO ==1,]) ], Impute_to_test[ name, row.names(Relative_Pheno[ Relative_Pheno$PHENO ==2,]) ]  )
  
  Temp[ name,] <- c(
    foo$statistic,
    foo$p.value,
    foo$conf.int[1],
    foo$conf.int[2]
  )
  
}

Temp <- as.data.frame(Temp)

hist(Impute_to_test[ name, row.names(Relative_Pheno[ Relative_Pheno$PHENO ==2,]) ], breaks = 60, col = "red")
hist(Impute_to_test[ name, row.names(Relative_Pheno[ Relative_Pheno$PHENO ==1,]) ], breaks = 60, col = "grey", add = T)

Temp <- as.data.frame(Temp)
Temp$FDR <- p.adjust( Temp$PVal, method = "fdr", n = dim(Temp)[1] ) 
Temp$BH <- p.adjust( Temp$PVal, method = "BH", n = dim(Temp)[1] ) 
Temp$BF <- p.adjust( Temp$PVal, method = "bonferroni", n = dim(Temp)[1] ) 

foo = as.data.frame(t(rbind(Actual = MatchedExpress[5,], Imputed=MatchedImpute[5,]) ))
ggplot( data =foo, aes(x=Actual, y=Imputed) ) + geom_point()
cor(foo$Actual, foo$Imputed)

#TEST Code for ploting imputed vs actual
SinkR <- as.data.frame(matrix( 0, dim(MatchedExpress)[2], 4 ))
colnames( SinkR ) <- c( "Tissue", "Cohort", "Actual", "Imputed" )
row.names( SinkR ) <- colnames(MatchedExpress)

#ROSMAP ID translation: syn3382527
RsTrans <- synGet("syn3382527")
RsTrans <- read.csv(RsTrans$path) 

RosCov <- fread(synGet('syn20800653')$path, header = T)
RosCov$Tissue <- "DLPFC"
RosCov$MatchID <- paste0( RsTrans[ RosCov$SampleID,]$gwas_id, "_", RsTrans[ RosCov$SampleID,]$gwas_id )

MayoCov <- fread(synGet('syn20801992')$path, header = T)
MayoCov$MatchID <- paste0( "0_", MayoCov$Donor_ID )
MayoCov$Tissue <- do.call( rbind, strsplit(MayoCov$SampleID, "_") )[,2]
table(MayoCov$MatchID %in% colnames(MatchedImpute))

MayoCov <- MayoCov[ (MayoCov$MatchID %in% row.names( SinkR )) == T, ]
row.names(MayoCov) <- MayoCov$MatchID

MSSBCov <- fread(synGet('syn20801799')$path, header = T)
MSSBCov$MatchID <- paste0( paste0( do.call( rbind, strsplit(MSSBCov$SampleID, '_'))[,1], do.call( rbind, strsplit(MSSBCov$SampleID, '_'))[,2]), '_', do.call( rbind, strsplit(MSSBCov$SampleID, '_'))[,3])

SinkR[ (grepl("ROS", row.names(SinkR)) == T | grepl("MAP", row.names(SinkR)) == T), ]$Tissue <- "DLFPC"
SinkR[ (grepl("ROS", row.names(SinkR)) == T | grepl("MAP", row.names(SinkR)) == T), ]$Cohort <- "ROSMAP"
SinkR[ grepl("BM", row.names(SinkR)) == T , ]$Cohort <- "MSBB"
SinkR[ grepl("BM10", row.names(SinkR)) == T , ]$Tissue <- "FP"
SinkR[ grepl("BM22", row.names(SinkR)) == T , ]$Tissue <- "STG"
SinkR[ grepl("BM36", row.names(SinkR)) == T , ]$Tissue <- "PHG"
SinkR[ grepl("BM44", row.names(SinkR)) == T , ]$Tissue <- "IFG"

SinkR[ row.names(SinkR) %in% MayoCov$MatchID, ]$Tissue 
SinkR[ SinkR$Cohort == 0, ]$Cohort <- "Mayo"
SinkR[ SinkR$Cohort == "Mayo", ]$Tissue <- "TCX"


PlotR <- function(ENSG, Target, Act, Imp){
  #'@ENSG the gene of interest eg. "ENSG00000131591"
  #'@Target The sink Matrix eg. SinkR
  #'@Act Matrix of actual expression z-values eg. MatchedExpress
  #'@Imp Matrix of imputed expression z-values eg. MatchedImpute
  #EXAMPLE: PlotR( "ENSG00000131591", SinkR, MatchedExpress, MatchedImpute)
  
  #ENSG <- "ENSG00000131591"
  #Target <- SinkR
  #Act <- MatchedExpress
  #Imp <- MatchedImpute
  
  Target$Actual <- as.numeric(Act[ ENSG, ])
  Target$Imputed <- Imp[ ENSG, ]
  
  lm_fit <- lm(Target$Imputed ~ Target$Actual, data=Target)
  summary(lm_fit)
  predicted_df <- data.frame(Fit = predict(lm_fit, Target))
  predicted_df$Act <- Target$Actual
  
  TEST <- cor.test(Target$Actual, Target$Imputed, method = "kendall")
  
  P1 <- ggplot( data =Target, aes(x=Actual, y=Imputed, col=Cohort) ) + 
    geom_point() + 
    theme(axis.title.x=element_blank() ) +
    theme(axis.title.y=element_blank() ) +
    #ggtitle( paste0( ENSG, ": Tau = ", signif(TEST$estimate, 3), " PVal = ", signif(TEST$p.value, 3) ) ) + 
    #theme(plot.title = element_text(hjust = 0.5)) +
    geom_line(color='black',data = predicted_df, aes(x=Act, y=Fit)) 
  
  P2 <- ggplot( data =Target, aes(x=Actual, y=Imputed, col=Tissue) ) + 
    geom_point() + 
    theme(axis.title.x=element_blank() ) +
    theme(axis.title.y=element_blank() ) +
    #ggtitle( paste0( ENSG, ": Tau = ", signif(TEST$estimate, 3), " PVal = ", signif(TEST$p.value, 3) ) ) + 
    #theme(plot.title = element_text(hjust = 0.5)) +
    geom_line(color='black',data = predicted_df, aes(x=Act, y=Fit))
  
  return( grid.arrange(P2, P1, ncol=2, top = paste0( ENSG, ": Tau = ", signif(TEST$estimate, 3), " PVal = ", signif(TEST$p.value, 3) ),
                       bottom = "Actual Expression", left = "Imputed Expression")
  )
  
}


Look <- sample(row.names(MatchedExpress), 20, replace = FALSE, prob = NULL)
for( name in Look ){
  PlotR( name, SinkR, MatchedExpress, MatchedImpute)
}

#PlotR( "ENSG00000131591", SinkR, MatchedExpress, MatchedImpute)
ohlala <- c( "ENSG00000112139", "ENSG00000175895", "ENSG00000154099", "ENSG00000104368" )
for( name in ohlala ){
  PlotR( name, SinkR, MatchedExpress, MatchedImpute)
}

#foo = as.data.frame(t(rbind(Actual = MatchedExpress[5,], Imputed=MatchedImpute[5,]) ))
#ggplot( data =foo, aes(x=Actual, y=Imputed) ) + geom_point()
#cor(foo$Actual, foo$Imputed)


###################
##Look at Kunkle Summary Stats in trained versus untrained genes
#_# SumStats_B  <- fread( synGet('syn18998950')$path, header = T)
SumStats <- fread( synGet('syn18999109')$path, header = T)
TWAS_SNPs <- fread( synGet('syn20835007')$path, header = F)

table(TWAS_SNPs$V2 %in% SumStats$MarkerName)
# 1066683 out of 1068623 SNPs in summary stats (99.82%)

SumStatsFilt <- SumStats[ (SumStats$MarkerName %in% TWAS_SNPs$V2) == T,  ]
SumStatsFilt <-  as.data.frame(SumStatsFilt)

#Load Gene Position File
GP_File <- fread(synGet('syn21413364')$path, head=T)
GP_File$Start <- GP_File$Coord - 500000
GP_File[ GP_File$Start < 0, ]$Start <- 0
GP_File$End <- GP_File$Coord + 500000
GP_File <- as.data.frame(GP_File)
row.names(GP_File) <- GP_File$Gene_Symbol

TrainedGenes <- GP_File[ row.names(Correlation), ]
MisGenes <- GP_File[ row.names(GP_File)[(row.names(GP_File) %in% row.names(Correlation)) == F], ]

InRange <- as.data.frame( matrix( 0, 0, dim(SumStatsFilt)[2]) )
colnames(InRange) <- colnames(SumStatsFilt)
OutRange <- as.data.frame( matrix( 0, 0, dim(SumStatsFilt)[2]) )
colnames(OutRange) <- colnames(SumStatsFilt)

#Build Bed Regions
row.names(SumStatsFilt) <- paste0("chr", SumStatsFilt$Chromosome, ':', SumStatsFilt$Position, '-', SumStatsFilt$Position+1)
GP_File$Coord <- paste0("chr", GP_File$Chr, ':', GP_File$Start, '-', GP_File$End)

#Validate the Vectors:
library(bedr)
is.a.valid  <- is.valid.region(row.names(SumStatsFilt));
is.b.valid  <- is.valid.region(GP_File$Coord);
table(is.a.valid)
table(is.b.valid)

#Sort the regions
Sorter <- function( IN ){
  if(is.sorted.region( IN ) == TRUE){
    return(IN)
  }else{
    a.sort <- bedr.sort.region(IN)
    return(a.sort)
  }
}
SNPS <- Sorter(row.names(SumStatsFilt))
Genes <- Sorter(GP_File$Coord)
GenesTrained <- Sorter( GP_File[row.names(MatchedImpute),]$Coord )
is.sorted.region(SNPS)
is.sorted.region(Genes)
is.sorted.region(GenesTrained)

#Make Intersections:
Within <- bedr(
  input = list(a = SNPS, b = Genes), 
  method = "intersect", 
  params = "-u"
)

Outside <- bedr(
  input = list(a = SNPS, b = Genes), 
  method = "intersect", 
  params = "-v"
)

WithinTrained <- bedr(
  input = list(a = SNPS, b = GenesTrained), 
  method = "intersect", 
  params = "-u"
)

OutsideTrained <- bedr(
  input = list(a = SNPS, b = GenesTrained), 
  method = "intersect", 
  params = "-v"
)

#Partition SNPS
SumStatsFilt$InTrainedGene <- NA
SumStatsFilt$InAnyGene <- NA

SumStatsFilt[Within,]$InAnyGene <- "Within"
SumStatsFilt[Outside,]$InAnyGene <- "Outside"

SumStatsFilt[WithinTrained,]$InTrainedGene <- "ModelTrained"
SumStatsFilt[OutsideTrained,]$InTrainedGene <- "NotTrained"
SumStatsFilt[Outside,]$InTrainedGene <- "NotNearGene"

SumStatsFilt$Pvalue <- as.numeric(SumStatsFilt$Pvalue)

Test<-wilcox.test(as.numeric(SumStatsFilt[Within,]$Pvalue), as.numeric(SumStatsFilt[Outside,]$Pvalue))
TTest <- t.test(as.numeric(SumStatsFilt[Within,]$Pvalue), as.numeric(SumStatsFilt[Outside,]$Pvalue))

Test2 <- wilcox.test(as.numeric(SumStatsFilt[SumStatsFilt$InTrainedGene == 'ModelTrained',]$Pvalue), as.numeric(SumStatsFilt[SumStatsFilt$InTrainedGene == 'NotTrained',]$Pvalue))
TTest2 <- t.test(as.numeric(SumStatsFilt[SumStatsFilt$InTrainedGene == 'ModelTrained',]$Pvalue), as.numeric(SumStatsFilt[SumStatsFilt$InTrainedGene == 'NotTrained',]$Pvalue))

Test3 <- wilcox.test(as.numeric(SumStatsFilt[SumStatsFilt$InTrainedGene == 'ModelTrained',]$Pvalue), as.numeric(SumStatsFilt[SumStatsFilt$InTrainedGene == 'NotNearGene',]$Pvalue))
TTest3 <- t.test(as.numeric(SumStatsFilt[SumStatsFilt$InTrainedGene == 'ModelTrained',]$Pvalue), as.numeric(SumStatsFilt[SumStatsFilt$InTrainedGene == 'NotNearGene',]$Pvalue))

#Raw <- ggplot( data = SumStatsFilt, aes( x = InAnyGene, y = Pvalue, col = InAnyGene) ) + geom_violin() + geom_boxplot(width=0.4) + theme(legend.position = "none")
#Density <- ggplot(data = SumStatsFilt, aes( x = Pvalue, col = InAnyGene) ) + geom_density()

#grid.arrange(Raw, Density, nrow = 1, top = paste0("Enrichment of Kunkle et.al SNP-PValues Within All TWAS Genes. PVal = ", signif( Test$p.value, digits = 4)) )

###Trained Genes
my_comparisons <- list( c("ModelTrained", "NotNearGene"), c("ModelTrained", "NotTrained") )
Raw <- ggplot( data = SumStatsFilt, aes( x = InTrainedGene, y = Pvalue, col = InTrainedGene) ) + geom_violin() + 
  geom_boxplot(width=0.4) + theme(legend.position = "none") + 
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + stat_compare_means( label.y = 1.28)
Density <- ggplot(data = SumStatsFilt, aes( x = Pvalue, col = InTrainedGene) ) + geom_density()

grid.arrange(Raw, Density, nrow = 1, top = paste0("Enrichment of Kunkle et.al SNP-PValues Near Model Trained Genes"))

################################################################################################################################################
HipGWAS <- data.table::fread( synapser::synGet('syn21421307')$path, header = T)
HipGWAS <- as.data.frame(HipGWAS)

#Filter for SNPs
HipGWAS <- HipGWAS[ (HipGWAS$A0 %in% c("A","T","C","G")) ==T, ]
HipGWAS <- HipGWAS[ (HipGWAS$A1 %in% c("A","T","C","G")) ==T, ]
row.names(HipGWAS) <- paste0( HipGWAS$Chrom, ':', as.integer(as.numeric(HipGWAS$Pos-1)), '-', as.integer(as.numeric(HipGWAS$Pos+1)) )

#Filter for Non-Exonic SNPS
library(biomaRt)
mart <- useMart("ensembl", host="grch37.ensembl.org")
mart <- useDataset("hsapiens_gene_ensembl", mart)
att <- listAttributes(mart)
f <- listFilters(mart)
attributes <- c("ensembl_exon_id","ensembl_gene_id", "chromosome_name", "exon_chrom_start", "exon_chrom_end", "start_position","end_position")
filters <- "ensembl_gene_id"
output <- getBM(attributes=attributes, filters=filters, values=GP_File$Gene_Symbol, mart=mart)

#Convert to zero based Bed
output$exon_chrom_end <- output$exon_chrom_end+1

hg19ExonicCoords <- paste0( "chr", output$chromosome_name, ':', output$exon_chrom_start, '-', as.integer(output$exon_chrom_end) )

is.a.valid  <- is.valid.region(hg19ExonicCoords);
is.b.valid  <- is.valid.region(row.names(HipGWAS));
table(is.a.valid)
table(is.b.valid)

Exons <- Sorter(hg19ExonicCoords[!duplicated(hg19ExonicCoords)])
BoneSNPs <- Sorter( row.names(HipGWAS) )
is.sorted.region(Exons)
is.sorted.region(BoneSNPs)

#Filter for non-exonic as these are obviously going to affect 
Non_Exonic_BoneSNP <- bedr(
  input = list(a = row.names(HipGWAS), b = Exons), 
  method = "intersect", 
  params = "-v"
)

##Filter for exonic and check a few entries (~2.7% of SNPs were Exonic)
#_#Exonic_BoneSNP <- bedr(
#_#  input = list(a = row.names(HipGWAS), b = Exons), 
#_#  method = "intersect", 
#_#  params = "-u"
#_#)

Non_Exon_GWAS <- HipGWAS[ Non_Exonic_BoneSNP, ]
#Make Intersections:
Within <- bedr(
  input = list(a = row.names(Non_Exon_GWAS), b = Genes), 
  method = "intersect", 
  params = "-u"
)

Outside <- bedr(
  input = list(a = row.names(Non_Exon_GWAS), b = Genes), 
  method = "intersect", 
  params = "-v"
)

WithinTrained <- bedr(
  input = list(a = row.names(Non_Exon_GWAS), b = GenesTrained), 
  method = "intersect", 
  params = "-u"
)

OutsideTrained <- bedr(
  input = list(a = row.names(Non_Exon_GWAS), b = GenesTrained), 
  method = "intersect", 
  params = "-v"
)

#Partition SNPS
Non_Exon_GWAS$InTrainedGene <- NA
Non_Exon_GWAS$InAnyGene <- NA

Non_Exon_GWAS[Within,]$InAnyGene <- "Within"
Non_Exon_GWAS[Outside,]$InAnyGene <- "Outside"

Non_Exon_GWAS[WithinTrained,]$InTrainedGene <- "ModelTrained"
Non_Exon_GWAS[OutsideTrained,]$InTrainedGene <- "NotTrained"
Non_Exon_GWAS[Outside,]$InTrainedGene <- "NotNearGene"

table(Non_Exon_GWAS$InAnyGene)
table(Non_Exon_GWAS$InTrainedGene)

my_comparisons <- list( c("ModelTrained", "NotNearGene"), c("ModelTrained", "NotTrained") )
Raw <- ggplot( data = Non_Exon_GWAS, aes( x = InTrainedGene, y = GCcorrP, col = InTrainedGene) ) + geom_violin() + 
  geom_boxplot(width=0.2) + theme(legend.position = "none") + 
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + stat_compare_means( label.y = 1.28)
Density <- ggplot(data = Non_Exon_GWAS, aes( x = GCcorrP, col = InTrainedGene) ) + geom_density()

grid.arrange(Raw, Density, nrow = 1, top = paste0("Enrichment of Styrkarsdottir, et.al. Bone Density GWAS SNP-PValues Near Model Trained Genes"))
