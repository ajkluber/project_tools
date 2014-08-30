""" Compute Jacobian for contact probability function


"""


import numpy as np


from project_tools.parameter_fitting import get_state_bounds,get_states_Vij


def get_target_feature(model):
    """ Get target features """
    name = model.subdir
    iteration = model.Mut_iteration
    
    ## To Do:
    ## - decide on target feature data format

def calculate_average_Jacobian(model):
    """ Calculate the average feature vector (ddG's) and Jacobian """
    
    name = model.subdir
    iteration = model.Mut_iteration

    cwd = os.getcwd()
    sub = "%s/%s/Mut_%d" % (cwd,name,iteration)
    os.chdir(sub)

    temperatures = [ x.split('_')[0] for x in open("T_array_last.txt","r").readlines() ] 
    directories = [ x.rstrip("\n") for x in open("T_array_last.txt","r").readlines() ] 

    epsilons = model.contact_epsilons
    deltas = model.contact_deltas
    sigmas = model.contact_sigmas

    bounds, state_labels = get_state_bounds()
    bounds = [0] + bounds + [model.n_contacts]

    ## Loop over temperatures in Mut subdir. Calculate ddG vector and 
    ## Jacobian for each directory indpendently then save. Save the average
    ## feature vector and Jacobian in the Mut/newton directory.
    sim_feature_all = []
    Jacobian_all = []
    for n in range(len(directories)):
        T = temperatures[n]
        dir = directories[n]
        beta = 1./(GAS_CONSTANT_KJ_MOL*float(T))
        print "  Calculating Jacobian for Mut_%d/%s" % (model.Mut_iteration,dir)
        os.chdir(dir)
        sim_feature, Jacobian = compute_Jacobian_for_directory(model,traj,dir,beta,bounds,state_labels,epsilons,delta,sigmas)
        sim_feature_all.append(sim_feature)
        Jacobian_all.append(Jacobian)
        os.chdir("..")

    sim_feaure_all = np.array(sim_feaure_all)
    Jacobian_all = np.array(Jacobian_all)

    ## Take avg. and use standard deviation as error bars.
    sim_feature_avg = sum(sim_feature_all)/float(len(directories))
    sim_feature_err = np.std(sim_feature_all,axis=0)
    Jacobian_avg = sum(Jacobian_all)/float(len(directories))
    Jacobian_err = np.std(Jacobian_all,axis=0)

    os.chdir(cwd)

    return sim_feature_avg, sim_feature_err, Jacobian_avg, Jacobian_err

def compute_Jacobian_for_directory(model,traj,dir,beta,bounds,state_labels,epsilons,delta,sigmas):
    """ Calculates the feature vector (ddG's) and Jacobian for one directory """

    ## Get trajectory, state indicators, contact energy
    traj,U,TS,N,Uframes,TSframes,Nframes,Vij = get_states_Vij(model,bounds,epsilons,delta,sigmas)

    ## Get Qi
    Qi = np.loadtxt("qimap.dat",dtype=float)
    sim_feature = sum(Qi[TS,:])/TSframes

    ## Initialize Jacobian
    Jacobian = np.zeros((model.n_contacts,model.n_contacts),float)

    ## Compute rows of the Jacobian which are correlation functions of 
    ## contact formation with contact energy.
    for i in range(model.n_contacts):
        Jacobian[i,:] = -beta*(sum(Qi[TS,i]*Vij[TS,:])/TSframes - sum(Qi[TS,i])*Vij[TS,:]/TSframes)

    return sim_feature, Jacobian

if __name__ == "__main__":

    
    ## To Do: 
    ##  -Module testing

    import model_builder as mdb

