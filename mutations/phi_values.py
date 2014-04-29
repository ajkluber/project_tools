import numpy as np
import os

import mdtraj as md
import cplex

import model_builder.models as models
import model_builder.systems as systems


'''
Alexander Kluber

Seed script to figure out how to calculate Phi-values.
'''

global GAS_CONSTANT_KJ_MOL
GAS_CONSTANT_KJ_MOL = 0.0083144621

def calculate_thermodynamic_perturbation(Model,System,append_log,coord="Q"):
    ''' First task is to calculate the perturbations for each mutation for
        each frame in the trajectory.   May be generalized in the future or 
        moved inside Model to deal with Models with multiple parameters per
        interaction (e.g. desolvation barrier, etc.)
    '''
    
    #append_log(System.subdir,"Starting: Calculating_MC2004")

    cwd = os.getcwd()
    sub = cwd+"/"+System.subdir+"/"+System.mutation_active_directory
    T = get_Tf_choice(sub)
    savedir = sub+"/"+T+"_agg"
    beta = 1./(GAS_CONSTANT_KJ_MOL*float(T))

    if not os.path.exists(savedir+"/mut"):
        os.mkdir(savedir+"/mut")

    files = ["M.dat","ddG.dat","eps.dat"] 
    flag = np.array([ not os.path.exists(savedir+"/mut/"+file) for file in files ])
    if np.any(flag):
        print "  One of the following does not exist: M.dat, ddG.dat, eps.dat. Calculating."
        os.chdir(System.subdir)
        ddG, eps, M = calculate_matrix_ddG_eps_M(Model,System,savedir,beta,coord)
        os.chdir(cwd)
    else:
        print "  Loading M.dat, ddG.dat, eps.dat"
        ddG = np.loadtxt(savedir+"/mut/ddG.dat")
        eps = np.loadtxt(savedir+"/mut/eps.dat")
        M = np.loadtxt(savedir+"/mut/M.dat")
    apply_linear_constraints(Model,System,savedir,ddG,eps,M)

    #append_log(System.subdir,"Finished: Calculating_MC2004")

def apply_linear_constraints(Model,System,savedir,ddG,eps,M):
    ''' Search the nullspace dimensions for a solution that minimizes 
        the change in stability and keeps the contacts attractive. Uses
        cplex linear programming package.'''

    ## The general solution is a sum of the particular solution and an
    ## arbitrary vector from the nullspace of M.
    Mpinv = np.linalg.pinv(M)
    x_particular = np.dot(Mpinv,ddG)
    np.savetxt(savedir+"/mut/x_p.dat",x_particular)
    #print x_particular     ## DEBUGGING

    ## Singular value decomposition. As a test you can recover M by,
    ## S = np.zeros(M.shape)
    ## S[:M.shape[1],:M.shape[1]] = np.diag(s)
    ## np.allclose(M,np.dot(u,np.dot(S,v))) --> should be True
    u,s,v = np.linalg.svd(M)
    rank = len(s)

    ## Nullspace basis vectors are the last n-r columns of the matrix v.T. As a check
    ## all entries of the matrix np.dot(M,N) should be extremely small ~0. Because 
    ## they've been sent to the nullspace.
    N = v.T[:,M.shape[0]:]

    ## Objective coefficients are sum of nullspace vectors. This comes from
    ## requiring the same average contact strength.
    objective_coeff = list(sum(N))

    ## Further constraints come from requiring that native contacts remain
    ## attractive.
    right_hand_side = list(eps + x_particular)
    rows = [ list(-N[i,:])  for i in range(len(N)) ] 
    senses = "L"*len(right_hand_side)
   
    ## Set upper and lower bounds on the solution. Arbitrary. Hopefullly these 
    ## don't matter.
    upper_bounds = list(10000.*np.ones(N.shape[1]))
    lower_bounds = list(-10000.*np.ones(N.shape[1]))

    ## Populate cplex linear programming problem
    LP_problem = cplex.Cplex()
    LP_problem.objective.set_sense(LP_problem.objective.sense.minimize)
    LP_problem.variables.add(obj=objective_coeff, ub=upper_bounds, lb=lower_bounds)
    LP_problem.linear_constraints.add(lin_expr=rows, senses=senses, rhs=right_hand_side)

    ## Let cplex do the hard work.
    LP_problem.solve()
    status = LP_problem.solution.get_status()
    solution = LP_problem.solution.get_values()

    print "Cplex summary:"
    print "status: ",status
    print "solution:",solutions

    ## Save new parameters in a BeadBead.dat file and indicate to use the path
    ## to the file in the System object.
    
def calculate_matrix_ddG_eps_M(Model,System,savedir,beta,coord):

    print "  Getting state bounds for coordinate:",coord
    bounds, states = get_state_bounds(savedir,coord)
    num_states = len(states)

    print "  Loading mutants"
    mutants = [ x.split()[1]+x.split()[0]+x.split()[2] for x in open("mutants/mutations.dat","r").readlines()[1:] ]

    print "  Loading trajectory, epsij, deltaij, sigij"
    sigij,epsij,deltaij,interaction_nums,keep_interactions,pairs,traj,traj_dist = load_eps_delta_sig_traj(savedir)
    Fij = get_mutant_fij(mutants,keep_interactions)
    qij = get_Qij(Model,traj_dist,sigij,deltaij,interaction_nums)
    print "  Loading dH for mutants"
    dH = get_mutant_dH(savedir,mutants)

    np.savetxt(savedir+"/mut/eps.dat",epsij)
    np.savetxt(savedir+"/mut/delta.dat",deltaij)
    np.savetxt(savedir+"/mut/sigma.dat",sigij)

    ## DEBUGGING 
    #print qij.shape, dH.shape, Fij.shape
    #print qij.shape, dH.shape
    #print "Qij shape",qij[states[1],:].shape
    #print "dH shape",dH[:,states[1]].shape
    #print "dH shape",sum(dH[:,states[1]])

    ## Load ddG from experiment and theory. 
    print "  Loading ddG from simulation"
    ddGsim_TS_D, ddGsim_N_D = get_sim_ddG(savedir,coord,bounds)
    print "  Loading ddG from experiment"
    ddGexp_TS_D, ddGexp_N_D = get_exp_ddG()
    ddG_all = np.concatenate(((ddGexp_TS_D - ddGsim_TS_D),(ddGexp_N_D - ddGsim_N_D)), axis=0)

    np.savetxt(savedir+"/mut/ddG.dat",ddG_all)
    
    ## Then compute each term of MC2004 equation (9) in a vectorized fashion.
    termA_TS_D = np.ones(Fij.shape,float)*(sum(qij[states[1],:])/float(len(qij[states[1],:])) - sum(qij[states[0],:])/float(len(qij[states[0],:])))
    termA_N_D = np.ones(Fij.shape,float)*(sum(qij[states[2],:])/float(len(qij[states[2],:])) - sum(qij[states[0],:])/float(len(qij[states[0],:])))
    termA_all = np.concatenate((termA_TS_D,termA_N_D),axis=0)
    
    termB_all = 1. - np.concatenate((Fij,Fij),axis=0) 

    termC_TS_D = ((np.dot(np.exp(beta*dH[:,states[1]]),qij[states[1],:])).T/sum(np.exp(beta*dH[:,states[1]]).T)).T \
              - ((np.dot(np.exp(beta*dH[:,states[0]]),qij[states[0],:])).T/sum(np.exp(beta*dH[:,states[0]]).T)).T
    termC_N_D = ((np.dot(np.exp(beta*dH[:,states[2]]),qij[states[2],:])).T/sum(np.exp(beta*dH[:,states[2]]).T)).T \
              - ((np.dot(np.exp(beta*dH[:,states[0]]),qij[states[0],:])).T/sum(np.exp(beta*dH[:,states[0]]).T)).T
    termC_all = np.concatenate((termC_TS_D,termC_N_D),axis=0)

    ## Then contruct the matrix M. 
    M = -beta*( termA_all -  termB_all*termC_all)

    ## DEBUGGING
    #print "Term A", termA_all.shape
    #print "Term B", termB_all.shape
    #print "Term C", termC_all.shape
    #print "M", M.shape
    #print "ddG", ddG_all.shape

    np.savetxt(savedir+"/mut/termA.dat",termA_all)
    np.savetxt(savedir+"/mut/termB.dat",termB_all)
    np.savetxt(savedir+"/mut/termC.dat",termC_all)
    np.savetxt(savedir+"/mut/M.dat",M)
    
    return ddG_all,epsij,M


def calculate_dH_for_mutants(Model,System,append_log):
    ''' First task is to calculate the perturbations for each mutation for
        each frame in the trajectory.   May be generalized in the future or 
        moved inside Model to deal with Models with multiple parameters per
        interaction (e.g. desolvation barrier, etc.)
    '''
    
    append_log(System.subdir,"Starting: Calculating_dH")

    cwd = os.getcwd()
    sub = cwd+"/"+System.subdir+"/"+System.mutation_active_directory
    T = get_Tf_choice(sub)
    savedir = sub+"/"+T+"_agg"

    os.chdir(System.subdir)

    mutants = [ x.split()[1]+x.split()[0]+x.split()[2] for x in open("mutants/mutations.txt","r").readlines()[1:] ]

    sigij,epsij,deltaij,interaction_nums,keep_interactions,pairs,traj,traj_dist = load_eps_delta_sig_traj(savedir)
    Fij = get_mutant_fij(mutants,keep_interactions)
    qij = get_Qij(Model,traj_dist,sigij,deltaij,interaction_nums)

    for j in range(len(Fij)):
        mut = mutants[j]
        if not os.path.exists(savedir+"/dH_"+mut+".dat"):
            fij = Fij[j]
            print "    Computing dH vectorized for ", mut
            dH_k = -1.*np.array([ sum(x) for x in fij*qij ])
            print "    Saving dH for ",mut
            np.savetxt(savedir+"/dH_"+mut+".dat",dH_k)
    os.chdir(cwd)
    append_log(System.subdir,"Finished: Calculating_dH")

def calculate_phi_values(Model,System,append_log,coord):
    ''' Calculate the phi values for a trajectory. Requires only state 
        definitions and dH (energetic perturbation for each mutation).
    '''
    
    append_log(System.subdir,"Starting: Calculating_phi_values")
    cwd = os.getcwd()
    sub = cwd+"/"+System.subdir+"/"+System.mutation_active_directory
    T = get_Tf_choice(sub)
    savedir = sub+"/"+T+"_agg"
    beta = 1./(GAS_CONSTANT_KJ_MOL*float(T))
    os.chdir(System.subdir)
    if not os.path.exists(savedir+"/phi"):
        os.mkdir(savedir+"/phi")

    print "  Loading mutations.txt"
    mutants = [ x.split()[1]+x.split()[0]+x.split()[2] for x in open("mutants/mutations.txt","r").readlines()[1:] ]
    print "  Getting state bounds for coordinate:",coord
    bounds, states = get_state_bounds(savedir,coord)
    num_states = len(states)
    print "  Loading dH for mutants"
    dH = get_mutant_dH(savedir,mutants)

    ## Compute deltaG for each state. Then DeltaDelta G with 
    ## respect to the first state (assumed to be the denatured state).
    ## Units of kT.
    print "  Computing ddG and phi values..."
    dG = [ -np.log(sum(np.exp(-beta*dH[:,states[X]]).T)/float(len(dH[:,states[X]]))) for X in range(num_states) ]
    ddG = [ dG[X]-dG[0] for X in range(1,num_states) ]

    ## Compute the phi value for each mutation. Phi is the ratio
    ## of DeltaDeltaG of the transition state(s) to DeltaDelta G
    ## of the native state (assumed to be the last state). unitless.
    phi = [ ddG[X]/ddG[-1] for X in range(len(ddG)-1) ]
    
    save_phi_values(savedir,mutants,coord,bounds,dG,ddG,phi)
    os.chdir(cwd)
    append_log(System.subdir,"Finished: Calculating_phi_values")

def get_mutant_dH(path,mutants):
    ''' Load the mutant energy perturbations dH_<mut>.dat'''

    i = 0
    for mut in mutants:
        temp = np.loadtxt(path+"/dH_"+mut+".dat")
        print "    Loading:",mut
        if i == 0:
            dH = np.zeros((len(mutants),len(temp)),float)
            dH[i,:] = temp
        else:
            dH[i,:] = temp
        i += 1
    
    return dH

def get_exp_ddG():
    ''' Get experimental ddG data from mutants/mutations.dat'''

    ddG_exp_all = np.loadtxt("mutants/mutations.dat",skiprows=1,usecols=(3,4))
     
    ddG_exp_TS_D = ddG_exp_all[:,0]
    ddG_exp_N_D = ddG_exp_all[:,1]

    return ddG_exp_TS_D, ddG_exp_N_D

def get_sim_ddG(savedir,coord,bounds):
    ''' Get the saved ddG from simulation that should have been computed already.'''

    index_sim = len(bounds)+1
    num = len(bounds)-1
    
    ddG_sim_all = np.loadtxt(savedir+"/phi/"+coord+"_phi.dat",skiprows=1,usecols=(4,5))
     
    ddG_sim_TS_D = ddG_sim_all[:,0]
    ddG_sim_N_D = ddG_sim_all[:,1]

    return ddG_sim_TS_D, ddG_sim_N_D


def get_mutant_fij(mutants,keep_interactions):
    ''' Load in the fraction of contact loss for each mutation.
        The matrix needs to be filtered to include only pair 
        interactions that correspond to parameters that are going
        to be mutated.'''
    k = 0
    for mut in mutants:
        fij_temp = np.loadtxt("mutants/fij_"+mut+".dat")
        fij_all = []
        for i in range(len(fij_temp)-4):
            fij_all.extend(fij_temp[i,i+4:])
        fij = np.array(fij_all)[keep_interactions != 0]
        if k == 0:
            Fij = np.zeros((len(mutants),len(fij)),float)
            Fij[0,:] = fij
        else:
            Fij[k,:] = fij
        k += 1
        
    return Fij

def get_Qij(Model,r,sig,delta,interaction_nums):
    ''' Calculates the normalized interaction betwen nonbonded pairs.'''
    print "  Calculating Qij..."
    qij = Model.nonbond_interaction(r,sig,delta)
    return qij

def get_state_bounds(path,coord):
    ''' Get bounds for each state for specified coordinate. Return a list of boolean
        arrays that specifies if each frame is in the given state or not.'''
    #print path+"/"+coord+"_states.txt" ## DEBUGGING
    #print open(path+"/"+coord+"_states.txt","r").read() ## DEBUGGING

    statefile = open(path+"/"+coord+"_states.txt","r").readlines()[1:]
    bounds = []
    for line in statefile:
        bounds.append([line.split()[0],float(line.split()[1]),float(line.split()[2]),float(line.split()[3])])

    if coord in ["Q","Qnh","Qh"]:
        data = np.loadtxt(path+"/"+coord+".dat")
        data /= max(data)
    elif coord == "Nh":
        data = np.loadtxt(path+"/Nh.dat")
    elif coord == "Rg":
        dummy,data = np.loadtxt(path+"/radius_cropped.xvg",unpack=True)
    elif coord == "rmsd":
        dummy,data = np.loadtxt(path+"/rmsd.xvg",unpack=True)
    else:
        print "ERROR!"
        print "  No option for coordinate: ", coord
        print "  Exiting"
        raise SystemExit

    states = []
    for i in range(len(bounds)):
        print "  State: ", bounds[i][0], " is defined as between: ",bounds[i][2], bounds[i][3]
        states.append((bounds[i][2] <= data)*(data <= bounds[i][3]))

    return bounds,states

def get_Tf_choice(sub):
    if not os.path.exists(sub+"/Tf_choice.txt"):
        print "ERROR!"
        print "  Please create ",sub+"/Tf_choice.txt with your choice to do mutations at."
        print "  Exiting"
        raise SystemExit
    else:
        Tf_choice = open(sub+"/Tf_choice.txt").read().split()[0]
        print "  Calculating dH for temp ",Tf_choice
    return Tf_choice

        
def load_eps_delta_sig_traj(subdir):
    ''' Load in the info from the BeadBead.dat file. Sig_ij, eps_ij, delta_ij and
        index pairs. This information is constant for a trajectory. Filter all fields
        to keep only interactions with nonzero interaction type.

        In calculating the mutations only modify parameters that have interaction_type
        in the BeadBead.dat =/= [0,ds,ss]. 
    '''
    print "  Loading BeadBead.dat"
    beadbead = np.loadtxt(subdir+"/BeadBead.dat",dtype=str) 
    sigij = beadbead[:,5].astype(float)
    epsij = beadbead[:,6].astype(float)
    deltaij = beadbead[:,7].astype(float)
    interaction_numbers = beadbead[:,4].astype(str)
    pairs = beadbead[:,:2].astype(int) 
    pairs -= np.ones(pairs.shape,int)

    keep_interactions = np.zeros(len(interaction_numbers),int)
    for i in range(len(interaction_numbers)):
        if interaction_numbers[i] in ["ds","ss"]:
            pass
        else:
            keep_interactions[i] = int(interaction_numbers[i])

    #print keep_interactions != 0       ## DEBUGGING
    #print sum((keep_interactions != 0).astype(int))      ## DEBUGGING
    sigij = sigij[keep_interactions != 0]
    epsij = epsij[keep_interactions != 0]
    deltaij = deltaij[keep_interactions != 0]
    pairs = pairs[keep_interactions != 0]

    print "  Only modifying ",sum((keep_interactions != 0).astype(int)), " parameters out of ", len(keep_interactions)
    ## Use mdtraj to compute the distances between pairs.
    print "  Loading traj.xtc with mdtraj..."
    traj = md.load(subdir+"/traj.xtc",top=subdir+"/Native.pdb")
    print "  Computing distances with mdtraj..."
    traj_dist = md.compute_distances(traj,pairs)

    return sigij,epsij,deltaij,interaction_numbers,keep_interactions,pairs,traj,traj_dist

def save_phi_values(savedir,mutants,coord,bounds,dG,ddG,phi):
    ''' Save the calculated dG, ddG, and phi values for states'''

    header_string = "# mut" 
    for i in range(len(bounds)):
        header_string += "     dG_"+bounds[i][0]+"(kT)"
    header_string += "    "
    for i in range(1,len(bounds)):
        header_string += " ddG_"+bounds[i][0]+"-"+bounds[0][0]+"(kT)"
    for i in range(1,len(bounds)-1):
        header_string += "   Phi_"+bounds[i][0]+"/"+bounds[-1][0]

    data_string = ''
    for j in range(len(mutants)):
        line = "%6s"%mutants[j]
        for i in range(len(dG)):
            line += "  %10.5f " % dG[i][j]
        for i in range(len(ddG)):
            line += "  %10.5f " % ddG[i][j]
        for i in range(len(phi)):
            line += "  %10.5f  " % phi[i][j]
        data_string += line+"\n"
    print "ddG and Phi values:"
    print header_string
    print data_string

    outputfile = open(savedir+"/phi/"+coord+"_phi.dat","w")
    outputfile.write(header_string+"\n"+data_string)
    outputfile.close()


if __name__ == '__main__':
    ## TESTING calculating ddG, phi values.
    def dummy_func(sub,string):
        pass 
    
    subdirs = ["r15"]
    Models = models.load_models(subdirs,dryrun=True)
    Systems = systems.load_systems(subdirs)
    Model = Models[0]
    System = Systems[0]
    path = System.subdir+"/"+System.mutation_active_directory+"/131.17_agg"

    '''
    #bounds, states = get_state_bounds(path,"Q") ## DEBUGGING
    dH, states = calculate_phi_values(Model,System,dummy_func)
    '''

    delta_eps = calculate_thermodynamic_perturbation(Model,System,dummy_func)
