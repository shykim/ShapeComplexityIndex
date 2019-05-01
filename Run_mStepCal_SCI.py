#!/usr/bin/python
## made by SunHyung (email: sunhyung.john.kim@gmail.com)
## This py has the purpose to submit the jobs in longleaf server. 
######################################################################################################### 

import os
import sys
from optparse import OptionParser

Current_Dir = os.getcwd()

def main(opts, argv):
	

	BIN_PATH = '/nas/longleaf/home/shykim/bin/'
	TEMPLATE_PATH = '/proj/NIRAL/atlas/Surface/CIVET_160K/FuncParcel/UNC/'
	
	if(opts.SIDE == 'L'): 	
		TEMPLATE_SURF = TEMPLATE_PATH + 'IBIS_nonBias_sym_Left.obj'
	if(opts.SIDE == 'R'): 	
		TEMPLATE_SURF = TEMPLATE_PATH + 'IBIS_nonBias_sym_Right.obj'

	IN_SURF = argv[0] ##e.g. mid_surf.obj
	_nDIVIDE = opts.nDIVIDE
	
	for i in range (1, _nDIVIDE+1):
		OUT_SCI = IN_SURF[:-4] + '_part' + str(i) + '.txt'	
		os.system('sbatch -n 1 -t 07-00:00:00 --wrap="%smCal_Complex -i %s -o %s -s 3 -k 4 -L -t %s -d %s -p %s"' %(BIN_PATH, IN_SURF,OUT_SCI,TEMPLATE_SURF, _nDIVIDE, i ) )
	
	
if (__name__ == "__main__"):
	parser = OptionParser(usage="%prog in_surf.obj [options]")
	parser.add_option("-d", "--nDivide", dest="nDIVIDE", default=0, type="int", help="How many divide n Vertex")
	parser.add_option("-s", "--SIDE", dest="SIDE", default=False, type="string", help="Which hemisphere")
	#parser.add_option("-p", dest="nPART", "--nPart", default=1, type=int, help="Which part run")
	#parser.add_option("-b",action="store", dest="T1T2Mask",type="string", help="use t1w and t2w bet mask, -b 'operator'(e.g. and,or)", default="" )
	#parser.add_option("-c","--MaskCombine",action="store_true", dest="Combine", default=False, help="Combine AutoSeg and T1/T2 mask, 2 of 3 would be a mask")
	#parser.add_option("-o","--ReOrientation",action="store_true", dest="ReOrient", default=False, help="If input data has LPI, Set orientation RAI")
	(opts, argv) = parser.parse_args()	

	#if (len(argv)<1):
 	#	parser.print_help()
	#	sys.exit(0)

	main(opts, argv)


