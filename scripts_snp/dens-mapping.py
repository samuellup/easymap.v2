import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-in', action="store", dest = 'in_file')
parser.add_argument('-out', action="store", dest = 'out')
args = parser.parse_args()

in_file = args.in_file
out = open(args.out, "w")

max_dens=[0,0,0]
for line in open(in_file, "a+"): 
    sp = line.split()
    #print sp[2], max_dens[2]
    if int(sp[2]) > int(max_dens[2]): 
        max_dens =[sp[0],sp[1],sp[2]]

out.write(str(max_dens[0]) + "\t" + str(max_dens[1]) + "\t" + str(max_dens[2]) + "\t" + str(in_file.split(".")[0]) + "\n")
