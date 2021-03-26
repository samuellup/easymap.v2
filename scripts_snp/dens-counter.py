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
out_list = list()
for ch in chromosomes:
	ch_list = list()
	for p in recount:
		if p[0] == ch[0]: 
			ch_list.append(p)

	if len(ch_list) > 3: 
		for i, p in enumerate(ch_list):
			if i == 0: 
				w_avg = 0.6*ch_list[i][2] + 0.25*ch_list[i+1][2] + 0.15*ch_list[i+2][2]
			elif i == 1: 
				w_avg =  0.2*ch_list[i-1][2] + 0.55*ch_list[i][2] + 0.15*ch_list[i+1][2] + 0.1*ch_list[i+2][2]
			elif i == len(ch_list)-2:
				w_avg = 0.1*ch_list[i-2][2] + 0.15*ch_list[i-1][2] + 0.55*ch_list[i][2] + 0.2*ch_list[i+1][2] 
			elif i == len(ch_list)-1:
				w_avg = 0.15*ch_list[i-2][2] + 0.25*ch_list[i-1][2] + 0.6*ch_list[i][2]
			else:
				w_avg = 0.1*ch_list[i-2][2] + 0.15*ch_list[i-1][2] + 0.5*ch_list[i][2] + 0.15*ch_list[i+1][2] + 0.1*ch_list[i+2][2]
			ch_list[i].append(w_avg)
	elif len(ch_list) == 3:
		for i, p in enumerate(ch_list): 
			if i == 0: 
				w_avg = 0.6*ch_list[i][2] + 0.25*ch_list[i+1][2] + 0.15*ch_list[i+2][2]
			elif i == 1: 
				w_avg =  0.2*ch_list[i-1][2] + 0.6*ch_list[i][2] + 0.2*ch_list[i+1][2] 
			elif i == 2: 
				w_avg = 0.15*ch_list[i-2][2] + 0.25*ch_list[i-1][2] + 0.6*ch_list[i][2] 
			ch_list[i].append(w_avg)
	elif len(ch_list) == 2:
		for i, p in enumerate(ch_list): 
			if i == 0: 
				w_avg = 0.75*ch_list[i][2] + 0.25*ch_list[i+1][2]
			elif i == 1: 
				w_avg =  0.25*ch_list[i-1][2] + 0.75*ch_list[i][2] 
			ch_list[i].append(w_avg)
	elif len(ch_list) == 1:
		for i, p in enumerate(ch_list): 
			w_avg = ch_list[i][2] 
		ch_list[i].append(w_avg)

	for p in ch_list: 
		out_list.append(p)

for p in out_list: 
	out.write(str(str(p[0])) + "\t" + str(p[1]) + "\t" + str(p[2]) + "\t" + str(p[3]) + "\n") 
