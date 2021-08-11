# This script formats information from VCF files to simpler VA files

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-gff', action="store", dest = 'gff')
parser.add_argument('-cr', action="store", dest = 'cand_reg')
parser.add_argument('-out', action="store", dest = 'out')
parser.add_argument('-f_an', action="store", dest = 'an')


args = parser.parse_args()

#Input  
gff = args.gff
cand_regs = args.cand_reg
fun_an = args.an
#Output
output = args.out
f2 = open(output, 'w')
f2.write('#CHR	GENE_START	GENE_END 	GENE 	GFF_ANOT	FUNCTIONAL_ANOT\n')


for line in open(cand_regs, "r"):
	if line.startswith("?"):
		sp = line.split()
		chr, start, end = str(sp[1]).lower(), sp[2], sp[3]
		if gff != "user_data/n/p": 
			with open(gff, "r") as an: 
				for line in an: 
					if not line.startswith("#"): 
						sp2 = line.split('\t')
						try: 
							if str(sp2[0]).lower() == chr and sp2[2] == "gene": 
								gene_st, gene_end = sp2[3], sp2[4]
								if (gene_st > start or gene_end > start) and (gene_st < end or gene_end < end): 
									f_an = "-"
									gene = str(sp2[8]).split(';')[0][3:]
									try:
										with open(fun_an, "r") as f_anotation: 
											for line in f_anotation: 
												fan = line.split()
												if str(fan[0]).lower() == gene.lower(): 
													f_an = str(line).replace("\t", ";")
									except: pass
									f2.write(str(chr) + '\t' + gene_st  + '\t' + gene_end  + '\t' + gene  + '\t' + str(sp2[8])  + '\t' + f_an   + '\n')
						except: pass
		else: 
			f2.write('\n' + "No GFF input provided, genome structural information not available."  + '\n')