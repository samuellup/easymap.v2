import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-in', action="store", dest = 'in_file')
parser.add_argument('-out', action="store", dest = 'out')
parser.add_argument('-step', action="store", dest = 'step')
parser.add_argument('-mut_type', action="store", dest = 'mut_type')
parser.add_argument('-cand_interval', action="store", dest= "cand_interval")
parser.add_argument('-pname', action="store", dest= "project")

args = parser.parse_args()

# Retrieving information from intermediate files
if int(args.step) == 1:
    in_file = args.in_file
    out = open(args.out, "a+")
    max_dens=[0,0,0,0]
    threshold=0.9

    for line in open(in_file, "r"): 
        sp = line.split()
        if float(sp[3]) > float(max_dens[3]): 
            max_dens =[sp[0],sp[1],sp[2],sp[3]]

    for line in open(in_file, "r"): 
        sp = line.split()
        if float(sp[3]) > threshold*float(max_dens[3]) :
            out.write(str(sp[0]) + "\t" + str(sp[1]) + "\t" + str(sp[2]) + "\t" + str(sp[3]) + "\t" + str((in_file.split(".")[-2]).split("/")[-1]) + "\n")
# Selecting center for CR 
if int(args.step) == 2: 
    peaks = list()
    if args.mut_type == "EMS": priorities = ["F2_filtered_EMS_hz_dens", "F2_filtered_EMS_dens"]
    if args.mut_type == "all": priorities = ["F2_hz_dens", "F2_control_comparison_dens"]

    for p in reversed(priorities):
        for line in open(args.in_file, "r"):
            if str(p) == str(line.split()[-1]): 
                if float(line.split()[3]) != 0.0: 
                    selected_pos =  [line.split()[0], line.split()[1], str(p)]                                                    
                    try: 
                        chrs=list()
                        for peak in peaks: 
                            chrs.append(str(peak[0]))
                        if str(selected_pos[0]) not in chrs:
                            peaks.append(selected_pos)
                        else:
                            for peak in peaks: 
                                if str(peak[0]) == str(selected_pos[0]):
                                    if abs(int(selected_pos[1]) - int(peak[1])) > 2000000: 
                                        peaks.append(selected_pos)
                    except:
                        peaks.append(selected_pos)
    # CR limits 
    regs = list()
    for peak in peaks: 
        try: 
            selected_chr = str(peak[0])
            selected_data = str(peak[2])
            selected_pos = int(peak[1])
            lines_d = {}
            project = args.project
            for i, line in enumerate(open(project+"/1_intermediate_files/"+selected_data+".txt", "r")):                          
                if selected_chr.lower() == str(line.split()[0]).lower(): 
                    lines_d[i] = line
                    if selected_pos == int(line.split()[1]):
                        selected_index = i

            # SETTING CR COORDINATES
            steepness = 0.5                                                                                        
            max_dens_value = float(lines_d[selected_index].split()[3])
            # Frontwards
            f = selected_index
            stop = False
            while stop == False: 
                try: 
                    n_dens = float(lines_d[f].split()[3])
                    n_dens_next = float(lines_d[f+1].split()[3])
                    if n_dens_next/n_dens > steepness and n_dens_next > max_dens_value*0.1:                                                                         
                        f = f+1
                    else: 
                        try: end_cr = str(lines_d[f+2].split()[1])
                        except: end_cr = str(lines_d[f+1].split()[1])
                        stop=True
                except: 
                    end_cr = str(lines_d[f].split()[1])
                    stop=True
            # Backwards
            f = selected_index
            stop = False
            while stop == False: 
                try:
                    n_dens = float(lines_d[f].split()[3])
                    n_dens_next = float(lines_d[f-1].split()[3])
                    if n_dens_next/n_dens > steepness and n_dens_next > max_dens_value*0.1:                                                                        
                        f = f-1
                    else: 
                        try: start_cr = str(lines_d[f-2].split()[1])
                        except: start_cr = str(lines_d[f-1].split()[1])
                        stop=True
                except: 
                    start_cr = str(lines_d[f].split()[1])
                    stop=True
            
            regs.append([selected_chr, start_cr, end_cr])       
            #with open(args.out, "a+") as out: 
            #    out.write('?'+"\t" + selected_chr +  "\t" + start_cr +"\t" + end_cr + "\n")

        except: 
            selected_pos = 'none'

    # Fix overlapping CRs 
    for r, reg in enumerate(regs): 
        chrm, x, y = reg[0], reg[1], reg[2]
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
    with open(args.out, "a+") as out: 
        for reg in regs: 
            if reg[1] != "-" and reg[2] != "-":
                out.write('?'+"\t" + reg[0] +  "\t" + reg[1] +"\t" + reg[2] + "\n")