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
    max_dens=[0,0,0]
    for line in open(in_file, "r"): 
        sp = line.split()
        if int(sp[2]) > int(max_dens[2]): 
            max_dens =[sp[0],sp[1],sp[2]]

    out.write(str(max_dens[0]) + "\t" + str(max_dens[1]) + "\t" + str(max_dens[2]) + "\t" + str((in_file.split(".")[0]).split("/")[-1]) + "\n")

# Selecting center for CR 
if int(args.step) == 2: 
    if args.mut_type == "EMS":
        priorities = ["F2_filtered_EMS_hz_dens", "F2_filtered_EMS_dens", "F2_control_comparison_dens"]
        for p in reversed(priorities):
            for line in open(args.in_file, "r"):
                if str(p) == str(line.split()[-1]): 
                    if int(line.split()[2]) != 0: 
                        selected_pos =  [line.split()[0], line.split()[1], str(p)]                        #str(line.split()[0])+'_-_'+str(line.split()[1])

    if args.mut_type == "all": 
        priorities = ["F2_hz_dens", "F2_control_comparison_dens"]
        for p in reversed(priorities):
            for line in open(args.in_file, "r"):
                if str(p) == str(line.split()[-1]): 
                    if int(line.split()[2]) != 0: 
                        selected_pos = selected_pos =  [line.split()[0], line.split()[1], str(p)]

    try: 
        selected_chr = str(selected_pos[0])
        selected_data = str(selected_pos[2])
        selected_pos = int(selected_pos[1])
        lines_d = {}
        project = args.project
        for i, line in enumerate(open(project+"/1_intermediate_files/"+selected_data+".va", "r")):                                             
            if selected_chr.lower() == str(line.split()[0]).lower(): 
                lines_d[i] = line
                if selected_pos == int(line.split()[1]):
                    selected_index = i

        # SETTING CR COORDINATES
        steepness = 0.20                                                                                         # ARBITRARY VALUE (0.2), calibrate accordingly
        max_dens_value = float(lines_d[selected_index].split()[2])
        # Frontwards
        f = selected_index
        stop = False
        while stop == False: 
            try: 
                n_dens = float(lines_d[f].split()[2])
                n_dens_next = float(lines_d[f+1].split()[2])
                if n_dens_next/n_dens > steepness and n_dens_next > max_dens_value*0.1:                                                                         
                    f = f+1
                else: 
                    end_cr = str(lines_d[f+1].split()[1])
                    stop=True
            except: 
                end_cr = str(lines_d[f].split()[1])
                stop=True

        # Backwards
        f = selected_index
        stop = False
        while stop == False: 
            try:
                n_dens = float(lines_d[f].split()[2])
                n_dens_next = float(lines_d[f-1].split()[2])
                if n_dens_next/n_dens > steepness and n_dens_next > max_dens_value*0.1:                                                                        
                    f = f-1
                else: 
                    start_cr = str(lines_d[f-1].split()[1])
                    stop=True
            except: 
                end_cr = str(lines_d[f].split()[1])
                stop=True
                
        with open(args.out, "w") as out: 
            out.write('?'+"\t" + selected_chr +  "\t" + start_cr +"\t" + end_cr)

    except: 
        selected_pos = 'none'