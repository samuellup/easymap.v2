# This script formats information from VCF files to simpler VA files

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-a', action="store", dest = 'input')
parser.add_argument('-b', action="store", dest = 'output')
parser.add_argument('-step', action="store", dest = 'step', default="1")

args = parser.parse_args()

#Input file 
input = args.input
f1 = open(input, 'r')
lines = f1.readlines()	

#Output
output = args.output
f2 = open(output, 'w')
f2.write('#CHR	POS	REF	ALT	QUAL	REF_DP	ALT_DP\n')

if str(args.step) == "1" : 
	for i, line in enumerate(lines):
		if not line.startswith('#'): 
			try:
				sp = line.split('\t')
				sp2 = sp[9].split(':')
				DP = sp2[2]
				ADF = sp2[3]
				ADR = sp2[4]
				ADF_ref= (ADF.split(','))[0]
				ADF_alt= (ADF.split(','))[1]
				ADR_ref= (ADR.split(','))[0]
				ADR_alt= (ADR.split(','))[1]
				alt = int(ADF_alt) + int(ADR_alt)
				ref = int(ADF_ref) + int(ADR_ref)
				DP = int(alt) + int(ref)
				f2.write(sp[0] + '\t' + str(sp[1]) + '\t' +sp[3] + '\t' +sp[4] + '\t' + str(sp[5]) + '\t' + str(ref) + '\t' + str(alt) + '\n')
			except:
				pass
if str(args.step) == "vcf_raw" :
	for i, line in enumerate(lines):
		if not line.startswith('#'): 
			try:
				sp = line.split('\t')
				'''
				try: 
					alt = sp[7].split('DP=')[1].split(';')[0]
				except: 
					alt = 1000
				ref = 0
				'''

				alt, ref = "-", "-"
				f2.write(sp[0] + '\t' + str(sp[1]) + '\t' +sp[3] + '\t' +sp[4] + '\t' + str(sp[5]) + '\t' + str(ref) + '\t' + str(alt) + '\n')
			except: 
				pass
	f2.close()