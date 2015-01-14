import numpy as np

global SKIP_INTERACTIONS
SKIP_INTERACTIONS = [1,8,9]

def get_contacts_epsilons(pairwise_params_file,epsilons_file):
    ''' Grab pairwise_params from file. '''
    model_param_values = np.loadtxt(epsilons_file)
    epsilons = []

    p_lines = [ x.rstrip("\n") for x in open(pairwise_params_file,"r").readlines() ]

    contacts = []

    for i in range(len(p_lines[1:])):
        data = p_lines[1+i].split()
        if not (int(data[3]) in SKIP_INTERACTIONS):
            contacts.append([int(data[0]),int(data[1])])
            if int(data[3]) in [3,5]:
                epsilons.append(-model_param_values[i])
            else:
                epsilons.append(model_param_values[i])

    contacts = np.array(contacts)
    epsilons = np.array(epsilons)

    return contacts,epsilons
