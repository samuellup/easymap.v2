# Testing: 	python2 qtl-mapping.py -in F2_control_comparison.va  -out mapping_info.txt -out2 candidate_region.txt -width 100000 -step 100000 -f_input genome.fa	-cand_interval 2000000

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-in', action="store", dest = 'in_file')
parser.add_argument('-width', action="store", dest = 'width')
parser.add_argument('-step', action="store", dest = 'step')
parser.add_argument('-out', action="store", dest = 'out')
parser.add_argument('-f_input', action="store", dest = 'f_input')

parser.add_argument('-mut_type', action="store", dest = 'mut_type') # implementation of different window sizes depending on SNP density
parser.add_argument('-cand_interval', action="store", dest= "cand_interval")
parser.add_argument('-out2', action="store", dest = 'out2')

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
	for center_pos in range(step/2, (int(ch[1])-step/2), step):   # for center_pos in range(width/2, (int(ch[1])-width/2), step):
		avgAF = 0
		AF_sum = 0.0
		snp_count = 0 
		left = center_pos - width/2
		right = center_pos + width/2
		with open(in_file, "r") as var: 
			for line in var: 
				sp = line.split()
				if str(sp[0]) == str(ch[0]) and int(sp[1]) > left and int(sp[1]) < right:
					ind = float(sp[9]) + float(sp[10])
					if float(ind) < 1.9 and float(sp[9]) > 0.2 and float(sp[10]) > 0.2 and float(sp[4]) > 100.0 :    # Criteria for selection of mapping SNP 
						snp_count = snp_count + 1
						AF_sum = AF_sum + float(sp[11])
		try: 
			if int(snp_count) > 1: 
				avgAF = float(AF_sum/float(snp_count))
				recount.append([ch[0], center_pos, avgAF]) 
		except: pass

# Weighted averages 
out_list = list()
w_avg_list = list()
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
			w_avg_list.append(w_avg)
	elif len(ch_list) == 3:
		for i, p in enumerate(ch_list): 
			if i == 0: 
				w_avg = 0.6*ch_list[i][2] + 0.25*ch_list[i+1][2] + 0.15*ch_list[i+2][2]
			elif i == 1: 
				w_avg =  0.25*ch_list[i-1][2] + 0.5*ch_list[i][2] + 0.25*ch_list[i+1][2] 
			elif i == 2: 
				w_avg = 0.15*ch_list[i-2][2] + 0.25*ch_list[i-1][2] + 0.6*ch_list[i][2] 
			ch_list[i].append(w_avg)
			w_avg_list.append(w_avg)
	elif len(ch_list) == 2:
		for i, p in enumerate(ch_list): 
			if i == 0: 
				w_avg = 0.7*ch_list[i][2] + 0.3*ch_list[i+1][2]
			elif i == 1: 
				w_avg =  0.3*ch_list[i-1][2] + 0.7*ch_list[i][2] 
			ch_list[i].append(w_avg)
			w_avg_list.append(w_avg)
	elif len(ch_list) == 1:
		for i, p in enumerate(ch_list): 
			w_avg = ch_list[i][2] 
		ch_list[i].append(w_avg)
		w_avg_list.append(w_avg)

	for p in ch_list: 
		out_list.append(p)

del recount
del ch_list

for p in out_list: 
	out.write(str(str(p[0])) + "\t" + str(p[1]) + "\t" + str(p[2]) + "\t" + str(p[3]) + "\n") 

#set starting limit dAF and dAF correction factor (average dAF across the genome)
dAF_lim = 0.4
sum_v=0.0
for v in w_avg_list: sum_v = sum_v + float(v)
try: dAF_correction = sum_v/(float(len(w_avg_list)))
except: dAF_correction = 0.0
del w_avg_list

dAF_peaks = list()
while dAF_lim >= 0.1:
	for p in out_list: 
		if (abs(float(p[3]) - dAF_correction )) > dAF_lim: 
			dAF_peaks.append(p)
	if len(dAF_peaks) >= 1: break
	else: 
		if dAF_lim >= 0.3: dAF_lim = dAF_lim - 0.05
		else: dAF_lim = dAF_lim - 0.02
del out_list

# CR selection 
regs = list()
cr_size = int(args.cand_interval)/2

for peak in dAF_peaks: 
	start_CR = int(peak[1]) - int(cr_size)/2
	if start_CR < 0: start_CR = 0
	end_CR =  int(peak[1]) + int(cr_size)/2
	regs.append([peak[0], start_CR, end_CR])

# Fix overlapping CRs 
for r, reg in enumerate(regs): 
	chrm, x, y = reg[0], int(reg[1]), int(reg[2])
	try: 
		if chrm == regs[r+1][0]:
			j, k = regs[r+1][1], regs[r+1][2]          
			if j == x and k ==  y:
				reg[1], reg[2] = "-", "-"
			elif j == x and y != k:
				if y > k :  regs[r+1][2] = y
				else: reg[2] = y
			elif y == k and x != j:
				if x > j : reg[1] = j  
				else: regs[r+1][1] = x
			else:
				if int(y) >= int(j) and int(k) >= int(x): 
					if int(x) < int(j):    
						reg[2] = j
					if int(x) > int(j): 
						reg[1] = k 
	except: 
		pass

regs_s = []
for reg in regs: 
	if reg not in regs_s: 
		regs_s.append(reg)
regs = regs_s

# Write to output
for reg in regs: 
	with open(args.out2, "a+") as out: 
		out.write('?'+"\t" + str(reg[0]) +  "\t" + str(reg[1]) +"\t" + str(reg[2]) + "\n")

#Write mapping data
try: 
	with open(args.out2, "a+") as out: 
		out.write('#'+"\t" + str("dAF_correction: ") + str(dAF_correction) +  "\t"  + str("dAF_lim: ") +  str(dAF_lim) + "\n")
except: pass