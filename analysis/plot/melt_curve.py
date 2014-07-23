import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import os
import argparse
import model_builder.models as models
'''                                                                                                                                           
Author: Fernando Yrazu
Created: June 2014                
Description:                     
    This program plots the melting and Cv curves of the proteins and saves the corresponding .pdf files.
    It should be executed from the parent directory of the simulation. 
    It will store the files in a new directory called /metrics/cv_melting within this parent directory.                       
Procedure:                                                                                                                                              
        1. Execute the program with the following options:  
        --prot = (Required) Input the proteins that you want to evaluate  
        --iter = (Required) Input the iteration number you are working with. Defaults to 0.                                   
Changelog:                      
June 2014 Created
July 2014 Modified to account for new directory structure                  
'''

def get_args(proteins_list):
    # User input to define which metrics are to be calculated, from a list displayed in the main() function                                          
    parser = argparse.ArgumentParser(description='Select')
    parser.add_argument('--prot', type=str, required=True, choices=proteins_list, nargs='+', help='Select proteins')
    parser.add_argument('--iter', type=str, required=True, help='Select iteration number')
    args = parser.parse_args()

    return args

def cv_melting(protein, iteration_number, current_dir):
    path = current_dir+'/'+protein+'/Tf_'+iteration_number+'/whamQ/'
    model = models.load_model(protein,dryrun=True)

    Cv = np.loadtxt(path+"cv",usecols=(0,1))
    QvsT = np.loadtxt(path+"Q_vs_T",dtype=float)
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twinx()
    ax1.plot(Cv[:,0],Cv[:,1],'r')
    ax1.set_xlabel("Temperature (K)")
    ax1.set_ylabel("Heat Capacity (kJ/mol K)")
    ax1.set_title("$C_v(T)$ and $\\left< Q \\right>(T)$ for %s" % model.name)

    ax2.plot(QvsT[:,0],QvsT[:,1]/model.n_contacts,'b')
    ax2.set_ylim(0,1)
    ax2.set_ylabel("$\\left< Q \\right>(T)$")
    plt.savefig('metrics/cv_melting/cv_and_melt_'+protein+'_iter_'+iteration_number+'.pdf')


def old_cv_melting(protein, iteration_number, current_dir):
    # This is the old version of the function, now rendered obsolete by new arrangement of directories 22-Jul-2014
    path = current_dir+'/'+protein+'/Tf_'+iteration_number+'/'
    
    # First plot the heat capacity curve
    cv = np.loadtxt(path+'whamQ/Heat_rmsd_Rg.dat', usecols=([0,2])).astype(float)
    plt.plot(cv[:,0],cv[:,1])
    plt.title('Heat capacity curve for '+protein+', iteration '+iteration_number)
    plt.xlabel('T [K]')
    plt.ylabel('Cv')
    plt.savefig(current_dir+'/metrics/cv_melting/cv_'+protein+'_iter_'+iteration_number+'.pdf')
    plt.clf()
    
    # Then plot the melting curve
    t_array = open(path+'T_array.txt').readlines()
    t_array = [x.split()[0] for x in t_array]
    temps = [int(x[:-2]) for x in t_array]
    # Open every directory in t_array and extract the average Q, store it in qs
    qs = []
    for t in t_array:
        qs.append(np.average(np.loadtxt(path+t+'/Q.dat').astype(float)))
    
    plt.plot(temps, qs)
    plt.title('Melting curve for '+protein+', iteration '+iteration_number)
    plt.xlabel('T [K]')
    plt.ylabel('Average Q')
    plt.savefig(current_dir+'/metrics/cv_melting/melting_curve_'+protein+'_iter_'+iteration_number+'.pdf')
    plt.clf()


def main():
    #List of possible proteins and functions to choose from                                            
    proteins_list = ['r15', 'r16', 'r17','1SHG', '1RIS', '1TEN', '1K85','1E0G','1E41']

    proteins_choice = [x for x in get_args(proteins_list).prot]

   #Mutation number (defaults to 0)                                                                    
    iteration_number = get_args(proteins_list).iter
    if iteration_number == None:
        iteration_number = '0'

    # Determine if the metrics directory exists and create it if not                                   
    current_dir = os.getcwd()

    # We will create a separate "metrics" directory to store all the data                             \
                                                                                                       
    if os.path.isdir(current_dir+'/metrics')==False:
        os.mkdir('metrics')

    os.chdir(current_dir+'/metrics/')
    if os.path.isdir(current_dir+'/metrics/cv_melting')==False:
        os.mkdir('cv_melting')
    os.chdir(current_dir)
    
    for y in proteins_choice:
        cv_melting(y, iteration_number, current_dir)
        
if __name__=='__main__':
    main()
