#!/usr/bin/python
#Altered code from: Gusev et al. “Integrative approaches for large-scale transcriptome-wide association studies” 2016 Nature Genetics

import sys
import os
import re
import argparse
import warnings

#Make the .POS file from the weights list and the list of gene starts and stops to FUSION.Assoc.R script
#'@USAGE': python ../Make_PositionFile.py -GA /fusion_twas-master/Training_CV_Groups/GenePos_Files/Chr21.txt -WL K9.lst -P WEIGHTS/chr21_Weights/ALL/ALL. -T .wgt.RDat -C 21 -o chr21_K9.pos
def get_args():
	'''This function will parse and return command-line arguments'''
	
	# Assign description to the help doc
	parser = argparse.ArgumentParser(
		description='This Script Procces MPRA Read data')

	#Specify arguments
	parser.add_argument(
		'-WL', '--WeightsList', type=str, help="List of Weights", required=True)
	parser.add_argument(
		'-GA', '--GeneAnnots', type=str, help="List of Gene Annotations", required=True)
	parser.add_argument(
		'-P', '--Preceed', type=str, help="String to trim preceeding the ENSG # of the weights list", required=True)
	parser.add_argument(
		'-T', '--Trail', type=str, help="String to trim Trailing the ENSG # of the weights list", required=True)
	parser.add_argument(
		'-C', '--Chr', type=str, help="Chromosmoe", required=False)
	#Pipeline output file Necessary in all modes			
	parser.add_argument(
		'-o', '--out', type=str, help='Specify the out file to store Pipeline Log', required=True)

	#Array for all arguments passed to script
	args = parser.parse_args()
	
	#Assign args to variables	
	WL = args.WeightsList
	GA = args.GeneAnnots
	
	P = args.Preceed
	T = args.Trail

	Chr = args.Chr
	O = args.out

	# Return all variable values
	return WL, GA, P, T, Chr, O

#Match return values from get_arguments() and assign to their respective variables
WL, GA, P, T, Chr, O = get_args()

if Chr == "ALL":
	Chr = [ str(range(1,23)) ]
else:
	pass

GRef = {}
F1 = open(GA)
for line in F1:
	LINE = line.rstrip('\r\n')
	lst = LINE.split( '\t' )
	if(Chr == "ALL" ):
		GRef[ lst[0] ] = '\t'.join([ lst[0], lst[1], lst[2], lst[3] ])
	else:
		if( int(lst[1]) == int(Chr) ):
			GRef[ lst[0] ] = '\t'.join([ lst[0], lst[1], lst[2], lst[3] ])
		else:
			pass

F1.close()

OUT = open( O, "w" )

print >>OUT, '\t'.join([ "WGT", "ID", "CHR", "P0", "P1" ])
F2 = open( WL )
for line in F2:
	LINE = line.rstrip('\r\n')
	ENSG = LINE.replace(P,'')
	ENSG = ENSG.replace(T,'')
	if ENSG in GRef.keys():
		print >>OUT, '\t'.join([ LINE, GRef[ENSG] ])
	else:
		pass
		#print ENSG

F2.close()
OUT.close()


