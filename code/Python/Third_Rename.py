import sys
import os
import re

InFile = 'TotalLDRef.bim'
InFileb = 'ROSMAP_temp.bim'

RS_Chr = {}
RS_Pos = {}
RS_BP ={}

F1 = open(InFile)
for line in F1:
	LINE = line.rstrip('\r\n')
	lst = LINE.split('\t')

	RS_Chr[ ''.join([lst[1], '_', lst[4], '_', lst[5]]) ] = lst[0]
	RS_Pos[ ''.join([lst[1], '_', lst[4], '_', lst[5]]) ] = lst[2]
	RS_BP[ ''.join([lst[1], '_', lst[4], '_', lst[5]]) ] = lst[3]
F1.close()

#Open Output
OUT = open( ''.join([ InFileb.replace(".bim", ""), 'Renamed.bim']), 'w' )
OUT2 = open( 'FailedSNPS.txt', 'w' )

F1 = open(InFileb)
for line in F1:
	LINE = line.rstrip('\r\n')
	lst = LINE.split('\t')
	LOC = ''.join([ lst[1], '_', lst[4], '_', lst[5] ])
	chrm = RS_Chr.get( LOC, None )
	pos = RS_Pos.get( LOC, None )
	bp = RS_BP.get( LOC, None )

	if( chrm == None ):
	  print >>OUT, LINE
	  print >>OUT2, lst[1]
	else:
	  Entry = '\t'.join([ chrm, lst[1], pos, bp, lst[4], lst[5] ])
	  print >>OUT, Entry

F1.close()
OUT.close()
OUT2.close()