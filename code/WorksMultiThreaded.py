#!/usr/bin/python
import sys
import os
import re
import subprocess32 as subprocess
import multiprocessing as mp
from multiprocessing.pool import ThreadPool
import datetime
import shlex, subprocess

#Runs the WEight calculation in parallel
##Must set all variables
##Must have Ordered IDs file in run location!!!!


#Set All Configurations 
GCTA='/gcta_1.92.4beta/gcta64'
PLINK='/usr/local/bin/plink2'
PLINKa='/usr/bin/plink1'
PLIN='/plink'
GEMMA="/fusion_twas-master/gemma-0.98.1-linux-static"
LDREF='/fusion_twas-master/LDREF'

#PRE_GEXP='Processed_ROSMAP_DLPFC_eQTLResidualExpression_Chr21.tsv'
#PRE_GENO='Processed_ROSMAP_DLPFC_eQTLResidualExpression.tsv'
   
#INDV IDs to keep for analysis:
INDV_ID='/fusion_twas-master/Total_IDs.txt'
#PATH TO OUTPUT DIRECTORY (population-specific subdirs will be made)
OUT_DIR='All_Data'

#Number of rows in the gene expression 
##Will loop through one gene at a time and run in parallel 
BATCH_START='1'
BATCH_END='13651'
NR='1_13651'
#BATCH_END=$(wc -l $PRE_GEXP | awk '{ print $1 }' -)

CMD = 'mkdir Genos'
subprocess.Popen(CMD, shell=True).wait()

#This is how I pull the input data from S3 into EC2
#CMD = ''.join([ 'aws s3 cp s3://jkg-s3-synapseencryptedexternalbucket-zszdd03ghnb2/ROSMAP/TWAS_CEU/FINAL_LDREF_CEU_MatchedEXP.bim . --recursive' ])
#subprocess.Popen(CMD, shell=True).wait()

#CMD = ''.join([ 'aws s3 cp s3://jkg-s3-synapseencryptedexternalbucket-zszdd03ghnb2/ROSMAP/TWAS_CEU/FINAL_LDREF_CEU_MatchedEXP.bed . --recursive' ])
#subprocess.Popen(CMD, shell=True).wait()

#CMD = ''.join([ 'aws s3 cp s3://jkg-s3-synapseencryptedexternalbucket-zszdd03ghnb2/ROSMAP/TWAS_CEU/FINAL_LDREF_CEU_MatchedEXP.fam . --recursive' ])
#subprocess.Popen(CMD, shell=True).wait()

for C in range(1,23):
	#split by chr
	CMD = ''.join([ '/plink --bfile All_CEU_ToTrainForTWAS --keep ', str(INDV_ID), ' --chr ', str(C), ' --threads 72 --make-bed --out Genos/AMP-AD_ROSMAP_Imputed_chr', str(C) ])
	subprocess.Popen(CMD, shell=True).wait()
	#Pull Genos
	#CMD = ''.join([ 'aws s3 cp s3://jkg-s3-synapseencryptedexternalbucket-zszdd03ghnb2/ROSMAP/Genotype_Imputed1000G/Binary_Cleaned/chr', str(C), '/ . --recursive' ])
	#subprocess.Popen(CMD, shell=True).wait()
	#Filter for INDV and LDREF SNPS ##Curently set to run on 72cores may need to change!
	#CMD = ''.join([ './plink --bfile AMP-AD_ROSMAP_Rush-Broad_AffymetrixGenechip6_Imputed_chr', str(C), '_ReNamed --keep ', str(INDV_ID), ' --indiv-sort f ', str(INDV_ID), ' --extract ', str(LDREF), '/1000G.EUR.', str(C), '.bim --threads 72 --make-bed --out Genos/AMP-AD_ROSMAP_Imputed_chr', str(C) ])
	#subprocess.Popen(CMD, shell=True).wait()
	#CleanUp
	#CMD = ''.join([ 'rm AMP-AD_ROSMAP_Rush-Broad_AffymetrixGenechip6_Imputed_chr', str(C), '_ReNamed.*'])
	#subprocess.Popen(CMD, shell=True).wait()



# Make the neccesary file paths below the run directory
CMD = ''.join([ 'mkdir temp_pheno' ])
subprocess.Popen(CMD, shell=True).wait()

CMD = ''.join([ 'mkdir All_Data/Genotype' ])
subprocess.Popen(CMD, shell=True).wait()

CMD = ''.join([ 'mkdir All_Data/WEIGHTS' ])
subprocess.Popen(CMD, shell=True).wait()
	
CMD = ''.join([ 'mkdir tmp/', str(NR), '/netResidual_Chrom_Expression/' ])
subprocess.Popen(CMD, shell=True).wait()
CMD = 'mkdir All_Data/ALL/netResidual_Chrom_Expression'
subprocess.Popen(CMD, shell=True).wait()

CMD = ''.join([ 'mkdir --parents tmp/', str(NR), '/netResidual_Chrom_Expression/' ])
subprocess.Popen(CMD, shell=True).wait()
CMD = ''.join([ 'mkdir --parents hsq/', str(NR), '/netResidual_Chrom_Expression/' ])
subprocess.Popen(CMD, shell=True).wait()
CMD = ''.join([ 'mkdir --parents out/', str(NR), '/netResidual_Chrom_Expression/' ])
subprocess.Popen(CMD, shell=True).wait()

# THIS IS DIRECTORY WHERE THE OUTPUT WILL GO:
CMD = ''.join([ 'mkdir --parents ', str(OUT_DIR) ])
subprocess.Popen(CMD, shell=True).wait()
	
CMD = ''.join(['mkdir ', str(OUT_DIR), '/ALL/' ])
subprocess.Popen(CMD, shell=True).wait()

	
##Each indv application of this function will run on a single core
def process( foo ):
	if( 'TSS' in str(foo) ):
		pass
	else:
		LINE = foo.rstrip('\r\n')
		lst = LINE.split('\t')

		CHR=lst[2]
		C=lst[2]
		
		#Define the genomic region of SNPs around the gene
		if(int(lst[3]) < 500000):
			P0 = 1
		else:
			P0=int(lst[3])-500000

		P1=int(lst[3])+500000
		GNAME=lst[1]
		
		OUT=''.join([ 'tmp/', str(NR), '/netResidual_Chrom_Expression/CEU_AncestralMapped_chr', str(C), '.tsv', '.', GNAME ])
		
		PRE_GENO = ''.join([ 'Genos/AMP-AD_ROSMAP_Imputed_chr', str(C) ])
		#PRE_GEXP= 'Processed_ROSMAP_DLPFC_eQTLResidualExpression.tsv'
		PRE_GEXP=str('/fusion_twas-master/ALL_CEU_RNA_Seq_SVAadusted_DiagnosisRegressed.tsv')
		

		#time log Start
		dt_started = datetime.datetime.utcnow()

		#Pull Gene Expression
		CMD = ''.join([ 'cp ', INDV_ID,' temp_pheno/', str(GNAME), '_foo.txt' ])
		subprocess.Popen(CMD, shell=True).wait()
		
		#Pull Gene Expression
		CMD = ''.join([ 'grep ', str(GNAME), ' ', str(PRE_GEXP), ' | sed \'s/\\t/\\n/g\' | sed 1,4d | paste temp_pheno/', str(GNAME), '_foo.txt - > ', str(OUT), '.pheno' ])
		#print CMD
		subprocess.Popen(CMD, shell=True).wait()
		
		#Pull The Genotype Data
		CMD = ''.join([ str(PLIN), ' --bfile ', str(PRE_GENO), ' --pheno ', str(OUT), '.pheno --make-bed --keep ', str(OUT), '.pheno --threads 8 --silent --chr ', str(C), ' --from-bp ', str(P0), ' --to-bp ', str(P1), ' --extract ', str(LDREF), '/1000G.EUR.', str(C), '.bim --out ', str(OUT), '.filt' ])
		#print CMD
		subprocess.Popen(CMD, shell=True).wait()

		#Make the Geno-typing data
		FINAL_OUT = ''.join([ str(OUT_DIR), '/ALL/', str(GNAME) ])
		
		#Run the Modified Weights Calculations
		CMD = ''.join([ 'Rscript FUSION.compute_weights_Parallel_V2.R --bfile ', str(OUT), '.filt --scale ', str(1),' --tmp ', str(OUT), ' --out ', str(FINAL_OUT), ' --verbose 0 --save_hsq --noclean --hsq_p 0.01 --PATH_plink ', str(PLIN), ' --PATH_gcta ', str(GCTA), ' --PATH_gemma ', str(GEMMA), ' --models blup,bslmm,lasso,top1,enet' ])
		print CMD
		subprocess.Popen(CMD, shell=True).wait()

		#Store some extra output data
		CMD = ''.join([ 'cat ', str(FINAL_OUT), '.hsq >> hsq/', str(NR), '.hsq' ])
		subprocess.Popen(CMD, shell=True).wait()

		#time log End and report progress
		dt_ended = datetime.datetime.utcnow()
		return( '\t'.join([ str(GNAME), str(CHR), str(P0), str(P1),  str( (dt_ended - dt_started).total_seconds() ) ]) )

#Change to fit the number of cores you have available
cores = 14
pool = mp.Pool(cores)
jobs = []

#Specify the input file
#This file is fomated w/ rows = Gene Expression x Columns = IDNV PLUS:
##The first colums are: GeneName	ENSG	Chromosome	TSS
##Also has a header (If your file does not you will skip the first line unless you change BATCH_START='1' to 0)
#INPT = str('Processed_ROSMAP_DLPFC_eQTLResidualExpression.tsv')
INPT = str('/fusion_twas-master/ALL_CEU_RNA_Seq_SVAadusted_DiagnosisRegressed.tsv')

#create jobs
with open( INPT, 'r' ) as f:
	#Skips the first line
	next(f)
	for line in f:
		jobs.append( pool.apply_async(process, [line]  ))
print INPT

#wait for all jobs to finish
for job in jobs:
	val = job.get()
	print val

#clean up
pool.close()
