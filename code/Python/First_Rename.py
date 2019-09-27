import sys
import os
import re

InFile = 'TotalLDRef.bim'
InFileb = 'GWASTotal.bim' #LDREF_FiltSNPs.bim

RS_IDs = {}
RS_Pos = {}

F1 = open(InFile)
for line in F1:
	LINE = line.rstrip('\r\n')
	lst = LINE.split('\t')

	RS_IDs[ ''.join([lst[0], '_', lst[3]]) ] = lst[1]
	RS_Pos[ ''.join([lst[0], '_', lst[3]]) ] = lst[2]
F1.close()

#Open Output
OUT = open( ''.join([ InFileb.replace(".bim", ""), 'Renamed.bim']), 'w' )

F1 = open(InFileb)
for line in F1:
	LINE = line.rstrip('\r\n')
	lst = LINE.split('\t')

	Name = RS_IDs[ ''.join([ lst[0], "_", lst[3] ])]
	Post = RS_Pos[ ''.join([ lst[0], "_", lst[3] ])]

	Entry = '\t'.join([ lst[0], Name, Post, lst[3], lst[4], lst[5] ])
	print >>OUT, Entry

F1.close()
OUT.close()
