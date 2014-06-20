import argparse
import shutil
import os
import model_builder.models.HeterogeneousGoModel as hetgm
import model_builder.models.CalphaBase as cab
import model_builder.models.HomogeneousGoModel as homgm
import model_builder.models.__init__ as init

'''                                                                                                                                                                                      
Author: Fernando Yrazu                                                                                                                                                                   
Created: June 2014                                                                                                                                                                       
Description:                                                                                                                                                                             
    This program checks for negative epsilons in the NewBeadBead.dat file generated after a ddG minimization run.
    It will save the previous .dat file and generate a new one, turning epsilons positive if they are negative and changing the corresponding delta to zero.
Procedure:                                                                                                                                                              
        1. Execute the program with the following options:                                                                                                                  
        --prot = (Required) Input the proteins that you want to evaluate                                                                                                       
        --iter = (Required) Input the iteration number you are working with                                                   
Changelog:                                                                                                                                                                               
June 2014 Created                                                                                                                                                                        
'''



def get_args(proteins_list):
    parser = argparse.ArgumentParser(description='Select')
    parser.add_argument('--prot', type=str, required=True, choices=proteins_list, help='Select protein')
    parser.add_argument('--iter', type=str, required=True, help='Select iteration number')
    args = parser.parse_args()
    return args

def flip_epsilons(model, subdir):
#    a = cab.CalphaBase()
#    indices,atoms,residues,coords = cab.CalphaBase.dissect_clean_pdb(a,subdir)
    #Read the old Newbeadbead file, archive it and flip the negative epsilons and their respective deltas
    directory = '/'.join(model.contact_energies.split('/')[:-1])
    old_bb_name = 'previous_'+model.contact_energies.split('/')[-1]
    old_path = directory + '/'+ old_bb_name
    shutil.move(model.contact_energies, old_path)
    new_path = directory + '/' + model.contact_energies.split('/')[-1]
    newbb = open(new_path, 'w')

    with open(old_path, 'r') as oldbb:
        for oldline in oldbb:
            line = oldline.split()
            eps = float(line[6])
            delta = float(line[7])
            if eps <0.:
                eps= -eps
                delta= 0
            if eps <0.1:
                eps = 0.1
            newbb.write('%5d%5d%8s%8s%5s%16.8E%16.8E%16.8E\n'%(int(line[0]),int(line[1]),line[2],line[3],line[4],float(line[5]), eps,delta))
                
#    model.contact_energies = new_path

       #hetgm.HeterogeneousGoModel.get_nonbonded_itp_strings(c, indices,atoms,residues,coords)            

def main():
    #List of possible proteins and functions to choose from                                                                       
    proteins_list = ['r15', 'r16', 'r17']
    protein_choice = get_args(proteins_list).prot
    iteration_number = get_args(proteins_list).iter
    current_dir = os.getcwd()
    subdir = current_dir + '/'+protein_choice
    model = init.load_model(subdir)
    flip_epsilons(model, subdir)
    

if __name__=='__main__':
    main()
