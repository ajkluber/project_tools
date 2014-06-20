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
May 2014 Created                                                                                                                  
June 2014 Added phi-value comparison                                                                                                             
'''
def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--prot', type=str, required=True, help='Input protein name')
    parser.add_argument('--iter', type=str, required=True, help='Input iteration number')
    args = parser.parse_args()

    return args


def plot_ddGs(protein, current_dir, iteration):
    select_temp = open(current_dir+'/'+protein+'/Mut_'+iteration+'/Tf_choice.txt').readline().split()[0]
    select_path = current_dir+'/'+protein+'/Mut_'+iteration+'/'+select_temp+'_agg/phi/'

    ddGdag_sim = np.loadtxt(select_path+'Q_phi.dat', usecols=(4,))
    ddG0_sim = np.loadtxt(select_path+'Q_phi.dat', usecols=(5,))
    ddGdag_exp_raw = np.loadtxt(current_dir+'/'+protein+'_calculated_ddG.dat', usecols=(6,))
    ddG0_exp_raw = np.loadtxt(current_dir+'/'+protein+'_calculated_ddG.dat', usecols=(4,))
    index_raw = np.loadtxt(current_dir+'/'+protein+'_calculated_ddG.dat', usecols=(1,))
    usable = np.loadtxt(current_dir+'/'+protein+'_calculated_ddG.dat', usecols=(8,), dtype=str)
    location = np.loadtxt(current_dir+'/'+protein+'_calculated_ddG.dat', usecols=(0,), dtype=str)
    exclude = np.loadtxt(current_dir+'/'+protein+'_calculated_ddG.dat', usecols=(11,), dtype=int)

    if os.path.exists(select_path+'mutations_excel.dat'):

        ddGdag_excel = np.loadtxt(select_path+'mutations_excel.dat', usecols=(3,))
        ddG0_excel = np.loadtxt(select_path+'mutations_excel.dat', usecols=(4,))
   
    ddGdag_exp = []
    ddG0_exp = []
    index =[]
    for i in range(len(usable)):
        if all([usable[i]=='True', location[i]=='core',exclude[i]==0]):
            ddGdag_exp.append(ddGdag_exp_raw[i])
            ddG0_exp.append(ddG0_exp_raw[i])
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
    upper = np.max([e,f]) + 1
    lower = np.min([c,d,0]) - 1
    
    plt.plot(ddGdag_sim, ddGdag_exp, 'ro')
    # Plot the 45 degree line
    plt.plot([lower,upper],[lower,upper],'b-', lw=2)

    plt.title("Comparison of simulation and experimental ddGdag for "+protein+", iteration "+iteration)
    plt.xlim([lower,upper])
    plt.ylim([lower,upper])
    plt.xlabel('ddGdag from simulation')
    plt.ylabel('ddGdag from experiment')
    plt.savefig(current_dir+'/metrics/ddG_comparison/'+protein+'_ddGdag_'+iteration+'.pdf')
    plt.clf()

    #Second plot
    #Finding the maxima and minima to size plot and draw y=x line  
    c = np.min(ddG0_sim)
    d = np.min(ddG0_exp)
    e = np.max(ddG0_sim)
    f = np.max(ddG0_exp)
    upper = np.max([e,f]) + 1
    lower = np.min([c,d,0]) - 1
    
    plt.plot(ddG0_sim, ddG0_exp, 'go')
    #Plot the 45 degree line
    plt.plot([lower,upper],[lower,upper],'b-', lw=2)

    plt.title("Comparison of simulation and experimental ddG0 for "+protein+", iteration "+iteration)
    plt.xlim([lower,upper])
    plt.ylim([lower,upper])
    plt.xlabel('ddG0 from simulation')
    plt.ylabel('ddG0 from experiment')
    plt.savefig(current_dir+'/metrics/ddG_comparison/'+protein+'_ddG0_'+iteration+'.pdf')
    plt.clf()

    if os.path.exists(select_path+'mutations_excel.dat'):
        c = np.min(ddG0_sim)
        d = np.min(ddG0_excel)
        e = np.max(ddG0_sim)
        f = np.max(ddG0_excel)
        upper = np.max([e,f]) + 1
        lower = np.min([c,d,0]) - 1

        plt.plot(ddG0_sim, ddG0_excel, 'go')
        #Plot the 45 degree line                                                                                                 
        plt.plot([lower,upper],[lower,upper],'b-', lw=2)
        
        plt.title("Comparison of simulation and experimental ddG0 for "+protein+", iteration "+iteration+' from EXCEL')
        plt.xlim([lower,upper])
        plt.ylim([lower,upper])
        plt.xlabel('ddG0 from simulation')
        plt.ylabel('ddG0 from experiment')
        plt.savefig(current_dir+'/metrics/ddG_comparison/'+protein+'_ddG0_'+iteration+'_excel.pdf')
        plt.clf()

        c = np.min(ddGdag_sim)
        d = np.min(ddGdag_excel)
        e = np.max(ddGdag_sim)
        f = np.max(ddGdag_excel)
        upper = np.max([e,f]) + 1
        lower = np.min([c,d,0]) - 1

        plt.plot(ddGdag_sim, ddGdag_excel, 'go')
        #Plot the 45 degree line                                                                                               
    
        plt.plot([lower,upper],[lower,upper],'b-', lw=2)

        plt.title("Comparison of simulation and experimental ddGdag for "+protein+", iteration "+iteration+' from EXCEL')
        plt.xlim([lower,upper])
        plt.ylim([lower,upper])
        plt.xlabel('ddGdag from simulation')
        plt.ylabel('ddGdag from experiment')
        plt.savefig(current_dir+'/metrics/ddG_comparison/'+protein+'_ddGdag_'+iteration+'_excel.pdf')
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
    plt.savefig(current_dir+'/metrics/ddG_comparison/'+protein+'_phi_values_'+iteration+'.pdf')
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
    plt.savefig(current_dir+'/metrics/ddG_comparison/'+protein+'_phi_value_sequence_'+iteration+'.pdf')
    plt.clf()

    #Fifth plot: Eigenvalues of M matrix for decomposition                                                                    
    select_path_mut = current_dir+'/'+protein+'/Mut_'+iteration+'/'+select_temp+'_agg/mut/singular_values_norm.dat'
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
        plt.savefig(current_dir+'/metrics/ddG_comparison/'+protein+'_M_eigenvalues.pdf')

def main():
    protein = get_args().prot
    iteration = get_args().iter
    current_dir = os.getcwd()    
    
    if os.path.isdir(current_dir+'/metrics')==False:
        os.mkdir('metrics')
    os.chdir(current_dir+'/metrics/')
    if os.path.isdir(current_dir+'/metrics/ddG_comparison')==False:
        os.mkdir('ddG_comparison')
    
    plot_ddGs(protein, current_dir, iteration)

if __name__=="__main__":
    main()
