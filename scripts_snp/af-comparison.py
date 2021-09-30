import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-f2_mut', action="store", dest = 'input_mut')
parser.add_argument('-f2_wt', action="store", dest = 'input_wt')
parser.add_argument('-out', action="store", dest = 'output')
parser.add_argument('-out_qtl', action="store", dest = 'output_qtl')
parser.add_argument('-out_test', action="store", dest = 'output_test')
parser.add_argument('-f_input', action="store", dest = 'f_input')
parser.add_argument('-step', action="store", dest = 'step')	
parser.add_argument('-mode', action="store", dest = 'mode')	

args = parser.parse_args()


#Input 
input1 = args.input_mut
f1 = open(input1, 'r')
mut_lines = f1.readlines()	

input2 = args.input_wt
f2 = open(input2, 'r')
wt_lines = f2.readlines()	
f_input = args.f_input
try: step = int(args.step)
except: step = 1
mode=args.mode

#Output
output = args.output
f3 = open(output, 'w')

#Get chromosome length and sequence
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

ch = list()

# Mode = noref
if mode == "noref":
	if step == 1 or step == 2: 
		#From the data of read_fasta, I create a dictionary with the name of the contigs and its lenght
		with open(f_input) as fp:
			for name_contig in read_fasta(fp):
				ch.append(name_contig[0][1:])

		for chr in ch:
			dic_mut = {}
			dic_wt = {}
			list_mut = list()
			list_wt = list()
			list_mut_af = list()


			for i, line in enumerate(mut_lines):
				if not line.startswith('#'):
					inlist=list()
					sp = line.split('\t')
					AF = float(float(sp[6])/(float(sp[6]) + float(sp[5])))
					inlist.append(sp[1])
					inlist.append(int(AF*100))

					if chr == sp[0] and AF > 0.85:
						dic_mut[sp[1]] = [sp[2], sp[3], sp[4], sp[5], sp[6].strip('\n')]
						list_mut.append(sp[1])
						list_mut_af.append(inlist)						
			
			for i, line in enumerate(wt_lines):
				if not line.startswith('#'):
					sp = line.split('\t')
					AF_wt = float(float(sp[6])/(float(sp[6]) + float(sp[5])))
					AF_wt_100 = int(AF_wt*100)
					AF_mut_100 = AF_wt_100					
					for it in list_mut_af:
						if int(it[0]) == int(sp[1]):
							AF_mut_100 = int(it[1])

					dAF = AF_mut_100 - AF_wt_100

					if chr == sp[0] and AF_wt < 0.5 and dAF > 60:
						dic_wt[sp[1]] = [sp[2], sp[3], sp[4], sp[5], sp[6].strip('\n')]
						list_wt.append(sp[1])

			set_2 = frozenset(list_mut)

			intersection = [x for x in list_wt if x in set_2]

			if len(intersection) >= 3:
				for i in intersection:
					if step == 1: 
						f3.write( str(chr) + '\t' + str(i) + '\t' + str(dic_mut[i][0]) +'\t' +  str(dic_mut[i][1]) +'\t' +  str(dic_mut[i][2]) + '\t' + str(dic_wt[i][3]) +'\t' +  str(dic_wt[i][4]) + '\t' + str(dic_mut[i][3]) + '\t' + str(dic_mut[i][4]) + '\n')
					if step == 2: 
						f3.write( str(chr) + '\t' + str(i) + '\t' + str(dic_mut[i][0]) +'\t' +  str(dic_mut[i][1]) +'\t' +  str(dic_mut[i][2]) + '\t'  + str(dic_mut[i][3]) + '\t' + str(dic_mut[i][4]) + '\t' + str(dic_wt[i][3]) +'\t' +  str(dic_wt[i][4]) + '\n')


# Mode = ref
if mode == "ref":
	if step == 1 or step == 2: 
		#From the data of read_fasta, I create a dictionary with the name of the contigs and its lenght
		with open(f_input) as fp:
			for name_contig in read_fasta(fp):
				ch.append(name_contig[0][1:])

		for chr in ch:
			dic_mut = {}
			dic_wt = {}
			list_mut = list()
			list_wt = list()

			for i, line in enumerate(mut_lines):
				if not line.startswith('#'):
					sp = line.split('\t')
					AF = float(float(sp[6])/(float(sp[6]) + float(sp[5])))
					if chr == sp[0] and AF > 0.85:
						dic_mut[sp[1]] = [sp[2], sp[3], sp[4], sp[5], sp[6].strip('\n')]
						list_mut.append(sp[1])

			for i, line in enumerate(wt_lines):
				if not line.startswith('#'):
					sp = line.split('\t')
					if chr == sp[0]:
						dic_wt[sp[1]] = [sp[2], sp[3], sp[4], sp[5], sp[6].strip('\n')]
						list_wt.append(sp[1])

			set_2 = frozenset(list_mut)

			intersection = [x for x in list_wt if x in set_2]


			if len(intersection) >= 3:
				for i in intersection:
					if step == 1: 
						f3.write( str(chr) + '\t' + str(i) + '\t' + str(dic_mut[i][0]) +'\t' +  str(dic_mut[i][1]) +'\t' +  str(dic_mut[i][2]) + '\t' + str(dic_wt[i][3]) +'\t' +  str(dic_wt[i][4]) + '\t' + str(dic_mut[i][3]) + '\t' + str(dic_mut[i][4]) + '\n')
					if step == 2: 
						f3.write( str(chr) + '\t' + str(i) + '\t' + str(dic_mut[i][0]) +'\t' +  str(dic_mut[i][1]) +'\t' +  str(dic_mut[i][2]) + '\t'  + str(dic_mut[i][3]) + '\t' + str(dic_mut[i][4]) + '\t' + str(dic_wt[i][3]) +'\t' +  str(dic_wt[i][4]) + '\n')

if mode == "qtl":
	#From the data of read_fasta, I create a dictionary with the name of the contigs and its lenght
	with open(f_input) as fp:
		for name_contig in read_fasta(fp):
			ch.append(name_contig[0][1:])
	f4 = open(args.output_qtl, 'w')
	f5 = open(args.output_test, 'w')
	for chr in ch:
		dic_mut = {}
		dic_wt = {}
		list_mut = list()
		list_wt = list()

		for i, line in enumerate(mut_lines):
			if not line.startswith('#'):
				sp = line.split('\t')
				if chr == sp[0]:
					dic_mut[sp[1]] = [sp[2], sp[3], sp[4], sp[5], sp[6].strip('\n')]
					list_mut.append(sp[1])

		for i, line in enumerate(wt_lines):
			if not line.startswith('#'):
				sp = line.split('\t')
				if chr == sp[0]:
					dic_wt[sp[1]] = [sp[2], sp[3], sp[4], sp[5], sp[6].strip('\n')]
					list_wt.append(sp[1])
		set_2 = frozenset(list_mut)

		# Intersection variants
		intersection = [x for x in list_wt if x in set_2]
		if len(intersection) >= 3:
			for i in intersection:
				AF_mut = str(float(dic_mut[i][4])/(float(dic_mut[i][3]) + float(dic_mut[i][4])))
				AF_wt = str(float(dic_wt[i][4])/(float(dic_wt[i][3]) + float(dic_wt[i][4])))
				dAF = str(float(AF_mut) - float(AF_wt))
				# f3 goes to mutation mapping module
				ind = float(AF_mut) + float(AF_wt)
				if float(ind) < 1.9 : 
					f3.write( str(chr) + '\t' + str(i) + '\t' + str(dic_mut[i][0]) +'\t' +  str(dic_mut[i][1]) +'\t' +  str(dic_mut[i][2]) + '\t' + str(dic_wt[i][3]) +'\t' +  str(dic_wt[i][4]) + '\t' + str(dic_mut[i][3]) + '\t' + str(dic_mut[i][4]) + '\t' + AF_mut + '\t' + AF_wt + '\t' + dAF + '\n')
				# f4 goes to candidates processing
				if abs(float(dAF)) > 0.3: 
					f4.write(str(chr) + '\t' + str(i) + '\t' + str(dic_mut[i][0]) +'\t' +  str(dic_mut[i][1]) +'\t' +  str(dic_mut[i][2]) + '\t' + str(dic_wt[i][3]) +'\t' +  str(dic_wt[i][4]) + '\t' + str(dic_mut[i][3]) + '\t' + str(dic_mut[i][4]) + '\t' + AF_mut + '\t' + AF_wt + '\t' + dAF + '\n')
		
		# F2 variants
		for line in mut_lines:
			if not line.startswith("#"): 
				sp = line.split()
				if chr == sp[0]:
					if sp[1] in intersection:
						i = sp[1] 
						AF_mut = str(float(dic_mut[i][4])/(float(dic_mut[i][3]) + float(dic_mut[i][4])))
						AF_wt = str(float(dic_wt[i][4])/(float(dic_wt[i][3]) + float(dic_wt[i][4])))
						dAF = str(float(AF_mut) - float(AF_wt))
						f5.write( str(chr) + '\t' + str(i) + '\t' + str(dic_mut[i][0]) +'\t' +  str(dic_mut[i][1]) +'\t' +  str(dic_mut[i][2]) + '\t' + str(dic_mut[i][3]) +'\t' +  str(dic_mut[i][4]) + '\t' + str(dic_wt[i][3]) + '\t' + str(dic_wt[i][4]) + '\t' + AF_mut + '\t' + AF_wt + '\t' + dAF + '\n')
					else: 
						i = sp[1] 
						AF_mut = str(float(dic_mut[i][4])/(float(dic_mut[i][3]) + float(dic_mut[i][4])))
						AF_wt = "-"
						dAF = "-"
						f5.write( str(chr) + '\t' + str(i) + '\t' + str(dic_mut[i][0]) +'\t' +  str(dic_mut[i][1]) +'\t' +  str(dic_mut[i][2]) + '\t' + str(dic_mut[i][3]) + '\t' + str(dic_mut[i][4]) + '\t' + "-" +'\t' +  "-"  + '\t' + AF_mut + '\t' + AF_wt + '\t' + dAF + '\n')

if step == 3: 
	#From the data of read_fasta, I create a dictionary with the name of the contigs and its lenght
	with open(f_input) as fp:
		for name_contig in read_fasta(fp):
			ch.append(name_contig[0][1:])

	for chr in ch:
		dic_mut = {}
		dic_wt = {}
		list_mut = list()
		list_wt = list()

		for i, line in enumerate(mut_lines):
			if not line.startswith('#'):
				sp = line.split('\t')
				AF = float(float(sp[6])/(float(sp[6]) + float(sp[5])))
				if chr == sp[0] and AF > 0.2 and AF < 0.85:
					dic_mut[sp[1]] = [sp[2], sp[3], sp[4], sp[5], sp[6].strip('\n')]
					list_mut.append(sp[1])

		for i, line in enumerate(wt_lines):
			if not line.startswith('#'):
				sp = line.split('\t')
				AF = float(float(sp[6])/(float(sp[6]) + float(sp[5])))
				if chr == sp[0] and AF > 0.2 and AF < 0.85:			
					dic_wt[sp[1]] = [sp[2], sp[3], sp[4], sp[5], sp[6].strip('\n')]
					list_wt.append(sp[1])

		set_2 = frozenset(list_mut)
		intersection = [x for x in list_wt if x in set_2]
		for i in intersection:
			f3.write( str(chr) + '\t' + str(i) + '\t' + str(dic_mut[i][0]) +'\t' +  str(dic_mut[i][1]) +'\t' +  str(dic_mut[i][2]) + '\t'  + str(dic_mut[i][3]) + '\t' + str(dic_mut[i][4]) + '\t' + str(dic_wt[i][3]) +'\t' +  str(dic_wt[i][4]) + '\n')

f3.close()
