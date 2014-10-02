""" Solve for new parameters

Description:

    Submodule to solve for the new parameter set using the thermodynamic
perturbation technique outlined in Matysiak Clementi 2004. Solutions are
found with the Levenberg-Marquardt algorithm, an algorithm for damped 
least-squares.


References:

(1) Matysiak, S.; Clementi, C. Optimal Combination of Theory and Experiment for
the Characterization of the Protein Folding Landscape of S6: How Far Can a
Minimalist Model Go? J. Mol. Biol. 2004, 343, 235-248.
"""

import matplotlib.pyplot as plt
import numpy as np

global GAS_CONSTANT_KJ_MOL
GAS_CONSTANT_KJ_MOL = 0.0083144621

def Levenberg_Marquardt_solution(model,method,scaling=False):
    """ Solve for new parameters using the Levenberg-Marquardt algorithm 

    Description:
    
        Solves the ill-posed (rectangular) system for the new parameters
    using the damped least-squares, Levenberg-Marquardt algorithm. Ill-posed
    systems usually have some small singular values that would blow the 
    solution up if the matrix inverse was calculated. The idea is to reduce 
    the influence of very small singular values, which have a lower signal 
    to noise ratio, on the solution. 

        The usual least-squares solution seeks to minimize the norm of the
    residual, 
                    Jx = df   =>  x s.t. min||Jx - df||

    damped least-squares modifies the problem to,

                A = J^T*J

    (A + lambda*diag(A))x = J^T*df  =>  x s.t. min{||Jx - b|| + lambda*||x||}

    Now a solution is sought that simultaneously minimizes the residual 
    and the norm of the solution. The choice of the damping parameter lambda
    can be determined

    """

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
    df /= np.linalg.norm(df)
    JTdf = np.dot(J.T,df)
    JTJ = np.dot(J.T,J)

    ## The damping parameter, lambda, is scanned from below the smallest
    ## singular value to above the largest singular value.
    ## The contributions from singular values less than the damping 
    ## parameter are 'damped' out by a filter function.
    n_rows = JTJ.shape[0]
    Levenberg = np.identity(n_rows)

    if scaling:
        ## Scale the diagonal by the curvature. This makes it the Levenberg-Marquardt method
        Levenberg[(np.arange(n_rows),np.arange(n_rows))] = np.diag(JTJ)

    u,s,v = np.linalg.svd(JTJ)
    if np.log10(min(s)) < -10:
        smin = -10
    else:
        smin = np.log10(min(s))
    Lambdas = np.logspace(smin,np.log10(max(s)),num=200)

    nrm_soln = []
    nrm_resd = []
    condition_number = []
    solutions = []
    for i in range(len(Lambdas)):
        Lambda = Lambdas[i]
        lhs = JTJ + Lambda*Levenberg
        x_soln = np.linalg.solve(lhs,JTdf) 
        residual = np.dot(J,x_soln) - df

        nrm_soln.append(np.linalg.norm(x_soln))
        nrm_resd.append(np.linalg.norm(residual))
        solutions.append(x_soln)

        condition_number.append(np.linalg.norm(lhs)*np.linalg.norm(np.linalg.inv(lhs)))

    save_solution_data(solutions,Lambdas,nrm_soln,nrm_resd,norm_eps,condition_number,s)

def save_solution_data(solutions,Lambdas,nrm_soln,nrm_resd,norm_eps,condition_number,s):
    """ Save and plot data for choosing solution and damping parameter 

    
    Description:
        
        A solution should be chosen that 
    """

    ## Save data
    for i in range(len(solutions)):
        np.savetxt("xp_%d.dat" % i,solutions[i])

    np.savetxt("lambdas.dat",Lambdas)
    np.savetxt("solution_norms.dat",nrm_soln)
    np.savetxt("residual_norms.dat",nrm_resd)
    np.savetxt("perturbation_norms.dat",np.array(nrm_soln)/norm_eps)
    np.savetxt("singular_values.dat",s)
    np.savetxt("condition_num.dat",condition_number)

    plot_Lcurve(nrm_resd,nrm_soln,Lambdas)
    plt.savefig("Lcurve.png")
    plt.savefig("Lcurve.pdf")

    plot_Lcurve_curvature(nrm_resd,nrm_soln,Lambdas)
    plt.savefig("Lcurve_curvature.png")
    plt.savefig("Lcurve_curvature.pdf")

    plot_condition_number(Lambdas,condition_number)
    plt.savefig("condition_vs_lambda.png")
    plt.savefig("condition_vs_lambda.pdf")

    plot_solutions(Lambdas,solutions)
    plt.savefig("solutions.png")
    plt.savefig("solutions.pdf")
    
def plot_solutions(Lambdas,solutions):

    solutions = np.array(solutions)

    plt.figure()
    plt.plot(Lambdas,solutions) 
    plt.title("solution values")
    plt.xlabel("$\\lambda$",fontsize=16)
    plt.ylabel("$\\delta\\epsilon_{\\lambda}$",fontsize=16)
    plt.ylim(-1,5)

def plot_Lcurve(nrm_resd,nrm_soln,Lambdas,skip=30):
    x = np.log10(nrm_resd)
    y = np.log10(nrm_soln)
    y2 = np.log10(Lambdas)

    ## Plot the L-curve
    fig, ax1 =  plt.subplots()
    ax2 = ax1.twinx()
    line1, = ax1.plot(x,y,'b')
    line2, = ax2.plot(x,y2,'r')
    ax1.set_xlabel("$\\log{||J\\delta\\epsilon_{\\lambda} - \\delta f||}$",fontsize=16)
    ax1.set_ylabel("$\\log{||\\delta\\epsilon_{\\lambda}}||}$",fontsize=16)
    ax2.set_ylabel("$\\log{\\lambda}$",fontsize=16)
    ticks = ax1.get_xticks()
    ax2.set_xticks(ticks)
    ax2.grid(True)
    plt.title("L-curve for choosing damping parameter $\\lambda$")

    ax1.yaxis.label.set_color(line1.get_color())
    ax2.yaxis.label.set_color(line2.get_color())

    ax1.tick_params(axis='y', colors=line1.get_color())
    ax2.tick_params(axis='y', colors=line2.get_color())

def plot_Lcurve_curvature(nrm_resd,nrm_soln,Lambdas):
    x = np.log10(nrm_resd)
    y = np.log10(nrm_soln)
    y2 = np.log10(Lambdas)
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
    ## optimum damping parameter, lambda_choice. lambda_choice should be 
    ## the lambda where the L-curve curvature is the most positive.
    fig, ax1 =  plt.subplots()
    ax2 = ax1.twinx()
    line1, = ax1.plot(x3,ddyddx,'b')
    ax1.xaxis.grid(True)
    ax1.yaxis.grid(True)
    line2, = ax2.plot(x,y2,'r')
    ax1.set_xlabel("$\\log{||J\\delta\\epsilon_{\\lambda} - \\delta f||}$",fontsize=16)
    ax1.set_ylabel("$\\frac{d^2}{dx^2}\\log{||\\delta\\epsilon_{\\lambda}}||}$",fontsize=16)
    ax2.set_ylabel("$\\log{\\lambda}$",fontsize=16)
    plt.title("$\\lambda_{choice}$ should be $\\lambda$ of max positive curvature")

    ax1.yaxis.label.set_color(line1.get_color())
    ax2.yaxis.label.set_color(line2.get_color())
    ax1.tick_params(axis='y', colors=line1.get_color())
    ax2.tick_params(axis='y', colors=line2.get_color())

def plot_condition_number(Lambdas,condition_number):
    y2 = np.log10(Lambdas)
    ## Plot condition number versus damping/regularization parameter lambda.
    fig, ax1 =  plt.subplots()
    ax2 = ax1.twinx()
    line1, = ax1.plot(y2,condition_number,'b')
    ax1.xaxis.grid(True)
    ax1.yaxis.grid(True)
    line2, = ax2.plot(y2,np.log(condition_number),'r')
    ax1.set_xlabel("$\\log{\\lambda}$",fontsize=16)
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
    model.Mut_iteration = iteration
    model.contact_epsilons = np.ones(model.n_contacts,float)
    method = "ddG_MC2004"

    cwd = os.getcwd()
    if not os.path.exists("%s/Mut_%d/test" % (name,iteration)):
        os.mkdir("%s/Mut_%d/test" % (name,iteration))
    os.chdir("%s/Mut_%d/test" % (name,iteration))

    Levenberg_Marquardt_solution(model,method)

    os.chdir(cwd)



