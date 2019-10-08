#!/usr/bin/python

import sys
import os
import re
#from pprint import pprint as pp
import subprocess32 as subprocess
import multiprocessing as mp
from multiprocessing.pool import ThreadPool
import datetime

#Runs the WEight calculation in parallel
##Must set all variables
##Must have Ordered IDs file in run location!!!!

#'@USAGE: python Run_Weights.py <EXPRESSION Data SET> <Geno Prefix> <wc -l EXPRESSION Data SET>
INPT_Exp = sys.argv[1]
INPT_Geno = sys.argv[2]
NumberOfGenes = sys.argv[3]

#Set All Configurations 
GCTA='/gcta_1.92.4beta/gcta64'
PLINK='/usr/local/bin/plink2'
PLINKa='/usr/bin/plink1'
PLIN='/plink'

GEMMA="/fusion_twas-master/gemma-0.98.1-linux-static"
#GEMMA='/usr/local/bin/gemma'
LDREF='/fusion_twas-master/LDREF'

PRE_GEXP=INPT_Exp
PRE_GENO=INPT_Geno
      
#PATH TO OUTPUT DIRECTORY (population-specific subdirs will be made)
OUT_DIR='Training_CV_Groups/All_Data'

#ROWS IN THE MATRIX TO ANALYZE (FOR BATCHED RUNS)
BATCH_START='1'
BATCH_END=NumberOfGenes
NR= str(BATCH_START) + '_' + str(NumberOfGenes)

CMD = ''.join([ 'mkdir --parents tmp/', NR ])
subprocess.Popen(CMD, shell=True).wait()
CMD = ''.join([ 'mkdir --parents hsq/', NR ])
subprocess.Popen(CMD, shell=True).wait()
CMD = ''.join([ 'mkdir --parents out/', NR ])
subprocess.Popen(CMD, shell=True).wait()

##Runs on a single core!
def process( foo ):
	if( 'TSS' in foo):
		pass
	else:
		LINE = foo.rstrip('\r\n')
		lst = LINE.split('\t')

		CHR=lst[2]
		C=lst[2]
		P0=int(lst[3])-500000
		P1=int(lst[3])+500000
		GNAME=lst[1]
		OUT=''.join([ 'tmp/', str(NR), '/', str(PRE_GEXP), '.', GNAME ])

		#print '\t'.join([ GNAME, CHR, P0, P1, str(datetime.datetime.now()) ])
		dt_started = datetime.datetime.utcnow()

		#Pull Gene Expression
		CMD = ''.join([ 'cp Total_IDs.txt ', str(GNAME), '_foo.txt' ])
		subprocess.Popen(CMD, shell=True).wait()

		#Pull Gene Expression
		CMD = ''.join([ 'grep ', str(GNAME), ' ', str(PRE_GEXP), ' | sed \'s/\\t/\\n/g\' | sed 1,4d | paste ', str(GNAME), '_foo.txt - > ', str(OUT), '.pheno' ])
		#print CMD
		subprocess.Popen(CMD, shell=True).wait()

		CMD = ''.join([ str(PLIN), ' --bfile ', str(PRE_GENO), ' --pheno ', str(OUT), '.pheno --make-bed --keep ', str(OUT), '.pheno --threads 8 --silent --chr ', str(C), ' --from-bp ', str(P0), ' --to-bp ', str(P1), ' --extract ', str(LDREF), '/1000G.EUR.', str(C), '.bim --out ', str(OUT), '.filt' ])
		subprocess.Popen(CMD, shell=True).wait()

		FINAL_OUT = ''.join([ str(OUT_DIR), '/ALL/', str(GNAME) ])

		CMD = ''.join([ 'Rscript FUSION.compute_weights_Parallel.R --bfile ', str(OUT), '.filt --tmp ', str(OUT), ' --out ', str(FINAL_OUT), ' --verbose 0 --save_hsq --noclean --hsq_p 0.01 --PATH_plink ', str(PLIN), ' --PATH_gcta ', str(GCTA), ' --PATH_gemma ', str(GEMMA), ' --models blup,bslmm,lasso,top1,enet' ])
		subprocess.Popen(CMD, shell=True).wait()

		CMD = ''.join([ 'cat ', str(FINAL_OUT), '.hsq >> hsq/', str(NR), '.hsq' ])
		subprocess.Popen(CMD, shell=True).wait()

		dt_ended = datetime.datetime.utcnow()
		return( '\t'.join([ str(GNAME), str(CHR), str(P0), str(P1),  str( (dt_ended - dt_started).total_seconds() ) ]) )

#init objects
cores = 6
pool = mp.Pool(cores)
jobs = []

#create jobs
INPT = INPT_Exp
with open( INPT ) as f:
	#Skips the first line
	#next(f)
	for line in f:
		jobs.append( pool.apply_async(process, [line]  ))

#wait for all jobs to finish
for job in jobs:
	val = job.get()
	print val

#clean up
pool.close()


