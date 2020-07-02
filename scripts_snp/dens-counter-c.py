import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-in', action="store", dest = 'in_file')
parser.add_argument('-width', action="store", dest = 'width')
parser.add_argument('-step', action="store", dest = 'step')
parser.add_argument('-out', action="store", dest = 'out')
parser.add_argument('-f_input', action="store", dest = 'f_input')

args = parser.parse_args()

width = int(args.width)
step = int(args.step)

in_file = args.in_file
f_input = args.f_input
out = open(args.out, "w")

chromosomes = list()
#This function enables to obtain data regarding chromosomes and lenght 
def read_fasta(fp):
	name, seq = None, []
	for line in fp:
		line = line.rstrip()
		if line.startswith('>'):
			if name: yield (name, ''.join(seq))
			name, seq = line, []
		else:
			seq.append(line)
	if name: yield (name, ''.join(seq))
#From the data of f_input
with open(f_input) as fp:
	for name_contig in read_fasta(fp):
		chromosomes.append([name_contig[0][1:], len(name_contig[1])])
max_chr = 0
for chr in chromosomes: 
	if int(chr[1]) > max_chr: max_chr = chr[1]

recount = list()
for ch in chromosomes: 
	for center_pos in range(width/2, max_chr, step):
		snp_count = 0 
		with open(in_file, "r") as var: 
			for line in var: 
				sp = line.split()
				left = center_pos - width/2
				right = center_pos + width/2
				if str(sp[0]) == str(ch[0]) and int(sp[1]) > left and int(sp[1]) < right: 
					snp_count = snp_count + 1
		recount.append([ch[0], center_pos, snp_count]) 


# Weighted averages 
for i, p in enumerate(recount):
	if i == 0: 
		w_dens = 0.6*recount[i][2] + 0.25*recount[i+1][2] + 0.15*recount[i+2][2]
	elif i == 1: 
		w_dens =  0.2*recount[i-1][2] + 0.55*recount[i][2] + 0.15*recount[i+1][2] + 0.1*recount[i+2][2]
	elif i == len(recount)-2:
		w_dens = 0.1*recount[i-2][2] + 0.15*recount[i-1][2] + 0.55*recount[i][2] + 0.2*recount[i+1][2] 
	elif i == len(recount)-1:
		w_dens = 0.15*recount[i-2][2] + 0.25*recount[i-1][2] + 0.6*recount[i][2]
	else:
		w_dens = 0.1*recount[i-2][2] + 0.15*recount[i-1][2] + 0.5*recount[i][2] + 0.15*recount[i+1][2] + 0.1*recount[i+2][2]
	recount[i].append(w_dens)
	
for p in recount: 
	out.write(str(str(p[0])) + "\t" + str(p[1]) + "\t" + str(p[2]) + "\t" + str(p[3]) + "\n") 
