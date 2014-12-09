""" Solve for new parameters

Description:

    Submodule to solve for the new parameter set using the thermodynamic
perturbation technique outlined in Matysiak Clementi 2004. Solutions are
found with the Truncated Singular Value Decomposition (TSVD) algorithm

References:

(1) Matysiak, S.; Clementi, C. Optimal Combination of Theory and Experiment for
the Characterization of the Protein Folding Landscape of S6: How Far Can a
Minimalist Model Go? J. Mol. Biol. 2004, 343, 235-248.
"""

import matplotlib.pyplot as plt
import numpy as np

def find_solutions(model,method):

    target_feature = np.loadtxt("target_feature.dat")
    target_feature_err = np.loadtxt("target_feature_err.dat")
    sim_feature = np.loadtxt("sim_feature.dat")
    sim_feature_err = np.loadtxt("sim_feature_err.dat")
    Jacobian = np.loadtxt("Jacobian.dat")
    Jacobian_err = np.loadtxt("Jacobian_err.dat")
    J = Jacobian

    epsilons = model.contact_epsilons
    norm_eps = np.linalg.norm(epsilons)

    ## Normalize the target step. 
    df = target_feature - sim_feature
    #df /= np.linalg.norm(df)   ## Arbitrary scaling

    u,s,v = np.linalg.svd(J)
    temp  = list(0.5*(s[:-1] + s[1:])) + [0.0]
    temp.reverse()
    cutoffs = np.array(temp)

    nrm_soln = []
    nrm_resd = []
    condition_number = []
    solutions = []
    Taus = []
    for i in range(len(cutoffs)):
        S = np.zeros(J.shape) 
        tau = cutoffs[i] 
        s_use = s[s > tau]
        S[np.arange(len(s_use)),np.arange(len(s_use))] = 1./s_use
        J_pinv = np.dot(v.T,np.dot(S.T,u.T))
        x_soln = np.dot(J_pinv,df)

        #x_soln = np.dot(v.T,np.dot(S.T,np.dot(u.T,df)))

        residual = np.dot(J,x_soln) - df

        nrm_soln.append(np.linalg.norm(x_soln))
        nrm_resd.append(np.linalg.norm(residual))
        solutions.append(x_soln)
        Taus.append(tau)

        J_use = np.dot(v.T,np.dot(S.T,u.T))     ## This isn't right
        cond_num = np.linalg.norm(J_use)*np.linalg.norm(J_pinv)
        condition_number.append(cond_num)

    save_solution_data(solutions,Taus,nrm_soln,nrm_resd,norm_eps,condition_number,s)


def save_solution_data(solutions,Taus,nrm_soln,nrm_resd,norm_eps,condition_number,s):
    """ Save and plot data for choosing solution and damping parameter 

    
    Description:
        
        A solution should be chosen that 
    """

    ## Save data
    for i in range(len(solutions)):
        np.savetxt("xp_%d.dat" % i,solutions[i])

    np.savetxt("taus.dat",Taus)
    np.savetxt("solution_norms.dat",nrm_soln)
    np.savetxt("residual_norms.dat",nrm_resd)
    np.savetxt("perturbation_norms.dat",np.array(nrm_soln)/norm_eps)
    np.savetxt("singular_values.dat",s)
    np.savetxt("condition_num.dat",condition_number)

    plot_Lcurve(nrm_resd,nrm_soln,Taus)
    plt.savefig("Lcurve.png")
    plt.savefig("Lcurve.pdf")

    plot_Lcurve_curvature(nrm_resd,nrm_soln,Taus)
    plt.savefig("Lcurve_curvature.png")
    plt.savefig("Lcurve_curvature.pdf")

    plot_condition_number(Taus,condition_number)
    plt.savefig("condition_vs_tau.png")
    plt.savefig("condition_vs_tau.pdf")

    plot_solutions(Taus,solutions)
    plt.savefig("solutions.png")
    plt.savefig("solutions.pdf")
    
def plot_solutions(Taus,solutions):

    solutions = np.array(solutions)

    plt.figure()
    plt.plot(Taus,solutions) 
    plt.title("solution values")
    plt.xlabel("$\\tau$",fontsize=16)
    plt.ylabel("$\\delta\\epsilon_{\\tau}$",fontsize=16)
    plt.ylim(-1,5)

def plot_Lcurve(nrm_resd,nrm_soln,Taus,skip=30):
    x = np.log10(nrm_resd)
    y = np.log10(nrm_soln)
    y2 = np.log10(Taus)

    ## Plot the L-curve
    fig, ax1 =  plt.subplots()
    ax2 = ax1.twinx()
    line1, = ax1.plot(x,y,'b')
    line2, = ax2.plot(x,y2,'r')
    ax1.set_xlabel("$\\log{||J\\delta\\epsilon_{\\tau} - \\delta f||}$",fontsize=16)
    ax1.set_ylabel("$\\log{||\\delta\\epsilon_{\\tau}}||}$",fontsize=16)
    ax2.set_ylabel("$\\log{\\tau}$",fontsize=16)
    ticks = ax1.get_xticks()
    ax2.set_xticks(ticks)
    ax2.grid(True)
    plt.title("L-curve for choosing damping parameter $\\tau$")

    ax1.yaxis.label.set_color(line1.get_color())
    ax2.yaxis.label.set_color(line2.get_color())

    ax1.tick_params(axis='y', colors=line1.get_color())
    ax2.tick_params(axis='y', colors=line2.get_color())

def plot_Lcurve_curvature(nrm_resd,nrm_soln,Taus):
    x = np.log10(nrm_resd)
    y = np.log10(nrm_soln)
    y2 = np.log10(Taus)
    ## Calculate the second derivative of the L-curve.
    x2 = np.array([ 0.5*(x[i] + x[i+1]) for i in range(len(x)-1)])
    x3 = np.array([ 0.5*(x2[i] + x2[i+1]) for i in range(len(x2)-1)])

    dy = np.diff(y)
    dx = np.diff(x)
    dx2 = np.diff(x2)
    dydx = dy/dx
    ddy = np.diff(dydx)
    ddyddx = ddy/dx2

    ## Plot the second derivative of the L-curve to help the user select the
    ## optimum damping parameter, tau_choice. tau_choice should be 
    ## the tau where the L-curve curvature is the most positive.
    fig, ax1 =  plt.subplots()
    ax2 = ax1.twinx()
    line1, = ax1.plot(x3,ddyddx,'b')
    ax1.xaxis.grid(True)
    ax1.yaxis.grid(True)
    line2, = ax2.plot(x,y2,'r')
    ax1.set_xlabel("$\\log{||J\\delta\\epsilon_{\\tau} - \\delta f||}$",fontsize=16)
    ax1.set_ylabel("$\\frac{d^2}{dx^2}\\log{||\\delta\\epsilon_{\\tau}}||}$",fontsize=16)
    ax2.set_ylabel("$\\log{\\tau}$",fontsize=16)
    plt.title("$\\tau_{choice}$ should be $\\tau$ of max positive curvature")

    ax1.yaxis.label.set_color(line1.get_color())
    ax2.yaxis.label.set_color(line2.get_color())
    ax1.tick_params(axis='y', colors=line1.get_color())
    ax2.tick_params(axis='y', colors=line2.get_color())

def plot_condition_number(Taus,condition_number):
    y2 = np.log10(Taus)
    ## Plot condition number versus damping/regularization parameter tau.
    fig, ax1 =  plt.subplots()
    ax2 = ax1.twinx()
    line1, = ax1.plot(y2,condition_number,'b')
    ax1.xaxis.grid(True)
    ax1.yaxis.grid(True)
    line2, = ax2.plot(y2,np.log(condition_number),'r')
    ax1.set_xlabel("$\\log{\\tau}$",fontsize=16)
    ax1.set_ylabel("condition number",fontsize=14)
    ax2.set_ylabel("$\\kappa = \\log{cond}$",fontsize=16)
    plt.title("Condition number versus damping parameter")

    ax1.yaxis.label.set_color(line1.get_color())
    ax2.yaxis.label.set_color(line2.get_color())
    ax1.tick_params(axis='y', colors=line1.get_color())
    ax2.tick_params(axis='y', colors=line2.get_color())


if __name__ == '__main__':
    import os
    import argparse
    import model_builder as mdb

    parser = argparse.ArgumentParser(description='.')
    parser.add_argument('--name', type=str, required=True, help='Name of protein to plot.')
    parser.add_argument('--iteration', type=int, required=True, help='Iteration to plot.')
    args = parser.parse_args()

    name = args.name
    iteration = args.iteration
    
    model = mdb.check_inputs.load_model("%s" % name,dry_run=True)
    model.iteration = iteration
    model.contact_epsilons = np.ones(model.n_contacts,float)
    method = "ddG_MC2004"

    cwd = os.getcwd()
    if not os.path.exists("%s/iteration_%d/test" % (name,iteration)):
        os.mkdir("%s/iteration_%d/test" % (name,iteration))
    os.chdir("%s/iteration_%d/test" % (name,iteration))

    find_solutions(model,method)

    os.chdir(cwd)



