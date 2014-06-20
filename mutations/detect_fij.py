""" Extract mutational data. Only return info for useable mutations """
mutation_data = np.loadtxt("calculated_ddG.dat",dtype=str)
useable_and_core = []
for i in range(mutation_data.shape[0]):
    if (mutation_data[i,0] == "core") and (mutation_data[i,8] == "True"):
        useable_and_core.append(True)
    else:
        useable_and_core.append(False)

useable_and_core = np.array(useable_and_core)
#useable_and_core = np.array([ all([(mutation_data[i,0] == "core"), bool(mutation_data[i,8])]) for i in range(mutation_data.shape[0]) ])

mut_indx = mutation_data[(useable_and_core == True),1]
wt_res = mutation_data[(useable_and_core == True),2]
mut_res = mutation_data[(useable_and_core == True),3]
 

