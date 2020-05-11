import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-in', action="store", dest = 'in_file')
parser.add_argument('-width', action="store", dest = 'width')
parser.add_argument('-step', action="store", dest = 'step')

args = parser.parse_args()

width = int(args.width)
step = int(args.step)

in_file = args.in_file
file_name = in_file.split(".")[0]
out_name = file_name + "_" + str(width) + "_" + str(step) + ".txt"
out = open(out_name, "w")


chromosomes = list()
with open(in_file, "r") as variants: 
	for line in variants: 
		if line.split()[0] not in chromosomes: 
			chromosomes.append(line.split()[0])

#variant 2
for ch in chromosomes: 
	for center_pos in range(width/2, 32000000, step):
		snp_count = 0 
		with open(in_file, "r") as var: 
			for line in var: 
				sp = line.split()
				left = center_pos - width/2
				right = center_pos + width/2
				if str(sp[0]) == str(ch) and int(sp[1]) > left and int(sp[1]) < right: 
					snp_count = snp_count + 1
		out.write(str(ch) + "\t" + str(center_pos) + "\t" + str(snp_count) + "\n") 










'''
python dens-counter.py -in GACDUS-GFGd.txt -width 50000 -step 20000

#variant 1
with open(in_file, "r") as variants: 
	for ch in chromosomes: 
		center_pos = width/2
		for center_pos in range(width/2, 32000000, step):
			snp_count = 0
			for line in variants:
				sp = line.split()
				if str(sp[0]) == str(ch):
					left = center_pos - width/2
					right = center_pos + width/2
					if int(sp[1]) > left and int(sp[1]) < right:
						snp_count = snp_count + 1
'''			