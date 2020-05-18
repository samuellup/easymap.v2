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
		out.write(str(str(ch[0])) + "\t" + str(center_pos) + "\t" + str(snp_count) + "\n") 
