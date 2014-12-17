import matplotlib.pyplot as plt
import numpy as np


def save_solution_data(solutions,Lambdas,nrm_soln,nrm_resd,norm_eps,condition_number,s):
    ''' Save and plot data for choosing solution and damping parameter 

    
    Description:
        
        A solution should be chosen that 
    '''

    ## Save data
    for i in range(len(solutions)):
        np.savetxt("xp_%d.dat" % i,solutions[i])

    np.savetxt("lambdas.dat",Lambdas)
    np.savetxt("solution_norms.dat",nrm_soln)
    np.savetxt("residual_norms.dat",nrm_resd)
    np.savetxt("perturbation_norms.dat",np.array(nrm_soln)/norm_eps)
    np.savetxt("singular_values.dat",s)
    #np.savetxt("condition_num.dat",condition_number)

    plot_singular_values(s)
    plt.savefig("singular_values.png")
    plt.savefig("singular_values.pdf")

    plot_log_singular_values(s)
    plt.savefig("log_singular_values.png")
    plt.savefig("log_singular_values.pdf")

    plot_Lcurve(nrm_resd,nrm_soln,Lambdas)
    plt.savefig("Lcurve.png")
    plt.savefig("Lcurve.pdf")

    plot_Lcurve_curvature(nrm_resd,nrm_soln,Lambdas)
    plt.savefig("Lcurve_curvature.png")
    plt.savefig("Lcurve_curvature.pdf")

    #plot_condition_number(Lambdas,condition_number)
    #plt.savefig("condition_vs_lambda.png")
    #plt.savefig("condition_vs_lambda.pdf")

    plot_solutions(Lambdas,solutions)
    plt.savefig("solutions.png")
    plt.savefig("solutions.pdf")

def plot_singular_values(s):

    plt.figure()
    plt.plot(s/max(s),'ro') 
    plt.title("singular values")
    plt.xlabel("i",fontsize=16)
    plt.ylabel("singular value",fontsize=16)

def plot_log_singular_values(s):

    plt.figure()
    plt.plot(np.log10(s/max(s)),'ro') 
    plt.title("singular values")
    plt.xlabel("i",fontsize=16)
    plt.ylabel("log(singular value)",fontsize=16)
    
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
    ax2.tick_params(axis='y', colors=line2.get_color())

