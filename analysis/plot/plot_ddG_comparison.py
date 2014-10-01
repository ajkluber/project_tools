import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import argparse
import os
import os.path

'''                                                                                                                               
Author: Fernando Yrazu.
                                                                                                                                 
Created: May 2014              
                                                                                                                                 
Description: 
 This program plots the ddG0, ddGdag and Phi-value comparison from simulation and experiment 
                                                                                                                                  
Procedure:                                                                                                                        
   Execute the program with the following options:
   --prot (protein_name)
   --iter (iteration_number)
   The iteration number should start from zero (meaning, raw simulation without experimental adjustment)
                                                    
Changelog:                                                                                                                       
May 2014 Created                                                                                                                 June 2014 Added phi-value comparison                                                                                             July 2014 Updated to new directory format, eliminated references to excel docs that were unnecessary                
'''
def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--prot', type=str, required=True, help='Input protein name')
    parser.add_argument('--iter', type=str, required=True, help='Input iteration number')
    args = parser.parse_args()

    return args


def plot_ddGs(protein, current_dir, iteration, select_temp):
    select_path_1 = current_dir+'/'+protein+'/Mut_'+iteration+'/'+select_temp+'_1/phi/'
    select_path_2 = current_dir+'/'+protein+'/Mut_'+iteration+'/'+select_temp+'_2/phi/'
    select_path_3 = current_dir+'/'+protein+'/Mut_'+iteration+'/'+select_temp+'_3/phi/'

    ddGdag_sim_1 = np.loadtxt(select_path_1+'Q_phi.dat', usecols=(4,))
    ddG0_sim_1 = np.loadtxt(select_path_1+'Q_phi.dat', usecols=(5,))
    ddGdag_sim_2 = np.loadtxt(select_path_2+'Q_phi.dat', usecols=(4,))
    ddG0_sim_2 = np.loadtxt(select_path_2+'Q_phi.dat', usecols=(5,))
    ddGdag_sim_3 = np.loadtxt(select_path_3+'Q_phi.dat', usecols=(4,))
    ddG0_sim_3 = np.loadtxt(select_path_3+'Q_phi.dat', usecols=(5,))

    ddGdag_sim = np.mean(np.column_stack((ddGdag_sim_1,ddGdag_sim_2,ddGdag_sim_3)), axis=1)
    ddG0_sim = np.mean(np.column_stack((ddG0_sim_1,ddG0_sim_2,ddG0_sim_3)), axis=1)
    err_ddGdag_sim = np.abs(np.amax(np.column_stack((ddGdag_sim_1,ddGdag_sim_2,ddGdag_sim_3)), axis=1)-
                            np.amin(np.column_stack((ddGdag_sim_1,ddGdag_sim_2,ddGdag_sim_3)), axis=1))
    err_ddG0_sim = np.abs(np.amax(np.column_stack((ddG0_sim_1,ddG0_sim_2,ddG0_sim_3)), axis=1)-
                          np.amin(np.column_stack((ddG0_sim_1,ddG0_sim_2,ddG0_sim_3)), axis=1))

    ddGdag_exp_raw = np.loadtxt(current_dir+'/'+protein+'_calculated_ddG.dat', usecols=(6,))
    err_ddGdag_exp_raw = np.loadtxt(current_dir+'/'+protein+'_calculated_ddG.dat', usecols=(7,))
    ddG0_exp_raw = np.loadtxt(current_dir+'/'+protein+'_calculated_ddG.dat', usecols=(4,))
    err_ddG0_exp_raw = np.loadtxt(current_dir+'/'+protein+'_calculated_ddG.dat', usecols=(5,))

    index_raw = np.loadtxt(current_dir+'/'+protein+'_calculated_ddG.dat', usecols=(1,))
    usable = np.loadtxt(current_dir+'/'+protein+'_calculated_ddG.dat', usecols=(8,), dtype=str)
    location = np.loadtxt(current_dir+'/'+protein+'_calculated_ddG.dat', usecols=(0,), dtype=str)
    exclude = np.loadtxt(current_dir+'/'+protein+'_calculated_ddG.dat', usecols=(11,), dtype=int)

    ddGdag_exp = []
    err_ddGdag_exp = []
    ddG0_exp = []
    err_ddG0_exp = []
    index =[]
    for i in range(len(usable)):
        if all([usable[i]=='True',exclude[i]==0]):
            ddGdag_exp.append(ddGdag_exp_raw[i])
            err_ddGdag_exp.append(err_ddGdag_exp_raw[i])
            ddG0_exp.append(ddG0_exp_raw[i])
            err_ddG0_exp.append(err_ddG0_exp_raw[i])
            index.append(index_raw[i])
  
    #Calculate Phi values        
    phi_exp = []
    phi_sim = []
    phi_exp_rel = []
    phi_sim_rel = []
    phi_exp_raw = []
    phi_sim_raw = []
    # Approximate value in kT that makes a ddG0 usable for the Phi-value calculation  
    ddG0_usable_minimum = 1.30
 
    for i in range(len(ddG0_exp)):
        exp = ddGdag_exp[i]/ddG0_exp[i]
        sim = ddGdag_sim[i]/ddG0_sim[i]
        # Approximate value in kT that 
        if ddG0_exp[i]<ddG0_usable_minimum:
            phi_exp_rel.append(exp)
            phi_sim_rel.append(sim)
        else:
            phi_exp.append(exp)
            phi_sim.append(sim)
        #these are for making the bar plot later on
        phi_exp_raw.append(exp)
        phi_sim_raw.append(sim)
            
    #First plot
    #Finding the maxima and minima to size plot and draw y=x line
    c = np.min(ddGdag_sim)
    d = np.min(ddGdag_exp)
    e = np.max(ddGdag_sim)
    f = np.max(ddGdag_exp)
    g = np.min(ddG0_sim)
    h = np.min(ddG0_exp)
    m = np.max(ddG0_sim)
    n = np.max(ddG0_exp)
    upper = np.max([e,f,m,n]) + 1
    lower = np.min([c,d,g,h,0]) - 1
    fig = plt.figure()
    ax1 = fig.add_subplot(111)

#    ax1.plot(ddGdag_sim, ddGdag_exp, 'ro', label="ddGdag ")                                                                    
    ax1.errorbar(ddGdag_sim, ddGdag_exp, fmt='ro', label="ddGdag ", xerr = err_ddGdag_sim, yerr= err_ddGdag_exp)
#    ax1.errorbar(ddGdag_sim, ddGdag_exp, fmt='ro', label="ddGdag ", xerr =0.08 , yerr=0.05)                                     
#    ax1.plot(ddG0_sim, ddG0_exp, 'go', label="ddG0")                                                                               
    ax1.errorbar(ddG0_sim, ddG0_exp, fmt='go', label="ddG0" , xerr = err_ddG0_sim, yerr=err_ddG0_exp)

    #Plot the tendency line of the entire sim vs exp ddGs (as in MatysiakClementi(20040                                             
    ddG_sim = np.concatenate((ddGdag_sim,ddG0_sim))
    ddG_exp = np.concatenate((ddGdag_exp,ddG0_exp))
    fit = np.polyfit(ddG_sim,ddG_exp, 1)
    fit_fn = np.poly1d(fit)
    ax1.plot([lower,upper], fit_fn([lower,upper]), '-k')
    # Plot the 45 degree line                                                                                                       
    ax1.plot([lower,upper],[lower,upper],'b-', lw=2)
    plt.title("Comparison of simulation and experimental ddG for "+protein+", iteration "+iteration)
    plt.xlim([lower,upper])
    plt.ylim([lower,upper])
    plt.xlabel('ddG from simulation')
    plt.ylabel('ddG from experiment')
    plt.legend(loc='lower right')
    plt.savefig(current_dir+'/metrics/ddG_comparison/'+protein+'_ddGs_'+iteration+'_'+select_temp+'.pdf')
    plt.clf()



    #Third plot: Phi values exp vs sim
    #Finding the maxima and minima to size plot and draw y=x line                                                                                     
    #Using latex for phi
#    plt.rc('text', usetex=True)    
                                                      
    c = np.min(phi_sim)
    d = np.min(phi_exp)
    if phi_exp_rel:
        e = np.max(phi_sim_rel)
        f = np.max(phi_exp_rel)
    else:
        e = 0.
        f = 0.

    upper = np.max([e,f]) + 1
    lower = np.min([c,d,0]) - 1

    plt.plot(phi_sim, phi_exp, 'ro')
    if phi_exp_rel:
        plt.plot(phi_sim_rel, phi_exp_rel, 'yo')

    #Plot the 45 degree line                                                         
    plt.plot([lower,upper],[lower,upper],'b-', lw=2)

    plt.title("Comparison of simulation and experimental $\Phi$ values for "+protein+", iteration "+iteration)
    plt.xlim([lower,upper])
    plt.ylim([lower,upper])
    plt.xlabel('$\Phi$ values from simulation')
    plt.ylabel('$\Phi$ values from experiment')
    plt.savefig(current_dir+'/metrics/ddG_comparison/'+protein+'_phi_values_'+iteration+'_'+select_temp+'.pdf')
    plt.clf()
    
    #Fourth plot: Phi values exp vs sim according to sequence                                 
    #Finding the maxima and minima to size plot and draw y=x line                                                                                                
    c = np.min(phi_sim)
    d = np.min(phi_exp)
    if phi_exp_rel:
        e = np.max(phi_sim+phi_sim_rel)
        f = np.max(phi_exp+phi_exp_rel)
    else:
        e = np.max(phi_sim)
        f = np.max(phi_exp)

    upper = np.max([f]) + 1
    lower = np.min([d,0]) - 1
    upper_all = np.max([e,f]) + 1
    lower_all = np.min([c,d,0]) - 1

    colors_exp = []
    colors_sim = []
    for i in range(len(ddG0_exp)):
        if ddG0_exp[i]<ddG0_usable_minimum:
            colors_exp.append('y')
            colors_sim.append('y')
        else:
            colors_exp.append('r')
            colors_sim.append('b')
  
    f, ax = plt.subplots(3,sharex=True)

    ax[0].bar(index,phi_exp_raw, color=colors_exp)
    ax[1].bar(index,phi_sim_raw, color=colors_sim)
    ax[2].bar(index,phi_sim_raw, color=colors_sim)
    ax[0].set_title("Comparison of simulation and experimental $\Phi$ values for "+protein+", iteration "+iteration+'\n according to sequence')
    ax[0].axhline(y=0, xmin=0, xmax=1)
    ax[1].axhline(y=0, xmin=0, xmax=1)
    ax[2].axhline(y=0, xmin=0, xmax=1)
    ax[0].set_ylim([lower,upper])
    ax[1].set_ylim([lower,upper])
    ax[2].set_ylim([lower_all,upper_all])
    ax[0].set_ylabel('$\Phi$ values (exp)')
    ax[1].set_ylabel('$\Phi$ values (sim)')
    ax[2].set_ylabel('$\Phi$ values (sim)')
#    f.subplots_adjust(hspace=0)
#    plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
    plt.savefig(current_dir+'/metrics/ddG_comparison/'+protein+'_phi_value_sequence_'+iteration+'_'+select_temp+'.pdf')
    plt.clf()

    #Fifth plot: Eigenvalues of M matrix for decomposition                                                                    
    select_path_mut = current_dir+'/'+protein+'/Mut_'+iteration+'/mut/singular_values_norm.dat'
    protein_dict={'r15':'ro', 'r16':'go', 'r17':'bo'}
    if protein in protein_dict:
        chosen_color = protein_dict[protein]
    else:
        chosen_color = 'ro'

    if os.path.exists(select_path_mut):
        ev = np.loadtxt(select_path_mut)
        plt.plot(ev, chosen_color)
        plt.axhline(y=0.5, xmin=0, xmax=1, color='m')
        plt.axhline(y=0.25, xmin=0, xmax=1, color='m')
        plt.axhline(y=0.1, xmin=0, xmax=1, color='m')
        plt.title("M-matrix normalized eigenvalues for "+protein)
        plt.xlabel('Ordinal #')
        plt.ylabel('Normalized eigenvalues')
        plt.savefig(current_dir+'/metrics/ddG_comparison/'+protein+'_'+select_temp+'_M_eigenvalues.pdf')

def main():
    protein = get_args().prot
    iteration = get_args().iter
    current_dir = os.getcwd()    
    
    if os.path.isdir(current_dir+'/metrics')==False:
        os.mkdir('metrics')
    os.chdir(current_dir+'/metrics/')
    if os.path.isdir(current_dir+'/metrics/ddG_comparison')==False:
        os.mkdir('ddG_comparison')
    
    temps_file = open(current_dir+'/'+protein+'/Mut_'+iteration+'/T_array_last.txt').readlines()
    select_temp = temps_file[0].split("_")[0]
    plot_ddGs(protein, current_dir, iteration,select_temp)

if __name__=="__main__":
    main()
