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

    for line in open(in_file, "r"): 
        sp = line.split()
        if float(sp[3]) > float(max_dens[3]): 
            max_dens =[sp[0],sp[1],sp[2],sp[3]]

    for line in open(in_file, "r"): 
        sp = line.split()
        if float(sp[3]) > 0.95*float(max_dens[3]) :
            out.write(str(sp[0]) + "\t" + str(sp[1]) + "\t" + str(sp[2]) + "\t" + str(sp[3]) + "\t" + str((in_file.split(".")[0]).split("/")[-1]) + "\n")
    

# Selecting center for CR 
if int(args.step) == 2: 
    regs = list()
    if args.mut_type == "EMS": priorities = ["F2_filtered_EMS_hz_dens", "F2_filtered_EMS_dens", "F2_control_comparison_dens"]
    if args.mut_type == "all": priorities = ["F2_hz_dens", "F2_control_comparison_dens"]

    for p in reversed(priorities):
        for line in open(args.in_file, "r"):
            if str(p) == str(line.split()[-1]): 
                if float(line.split()[3]) != 0.0: 
                    selected_pos =  [line.split()[0], line.split()[1], str(p)]                        #str(line.split()[0])+'_-_'+str(line.split()[1])                            
                    try: 
                        chrs=list()
                        for reg in regs: 
                            chrs.append(str(reg[0]))
                        if str(selected_pos[0]) not in chrs:
                            regs.append(selected_pos)
                        else:
                            for reg in regs: 
                                if str(reg[0]) == str(selected_pos[0]):
                                    if abs(int(selected_pos[1]) - int(reg[1])) > 10000000: 
                                        regs.append(selected_pos)
                    except:
                        regs.append(selected_pos)

    # CR limits 
    for reg in regs: 
        try: 
            selected_chr = str(reg[0])
            selected_data = str(reg[2])
            selected_pos = int(reg[1])
            lines_d = {}
            project = args.project
            for i, line in enumerate(open(project+"/1_intermediate_files/"+selected_data+".txt", "r")):                          
                if selected_chr.lower() == str(line.split()[0]).lower(): 
                    lines_d[i] = line
                    if selected_pos == int(line.split()[1]):
                        selected_index = i

            # SETTING CR COORDINATES
            steepness = 0.20                                                                                         # ARBITRARY VALUE (0.2), calibrate accordingly
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
                    
            with open(args.out, "a+") as out: 
                out.write('?'+"\t" + selected_chr +  "\t" + start_cr +"\t" + end_cr + "\n")

        except: 
            selected_pos = 'none'
