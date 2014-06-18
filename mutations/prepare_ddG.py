import numpy as np
import csv
import os
import sys
import argparse

'''                                                                                                                               
Author: Fernando Yrazu. Adapted from DMC model code from Brad Lambeth                                                                                                          
Created: May 2014                                                                                                            
                                                                                                                                  
Description:                                                                                                                      
   This module reads the published mutation data and makes the necessary calculations to yield the corresponding ddG values in units of kT. It can also start from a .csv (MS Excel) file. The only input needed is the protein name. It should be executed from the base directory and will store its results there as well.

Procedure:
 

Remaining issues:
   1) Find out how to deal with kf when kf_denaturant_conc is not zero
   2) Exception handling if wt accidentally ommited from protein_ddG.dat
   
Changelog:                                                                                                                        
May 2014 Created                                                                                                             
June 2014 Modified to include proper treatment of surface mutations                                                                                        JunJune 2014 Turned into a module 
June 2014 Added the option to incorporate ddGs directly from available data
'''

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--prot', type=str, required=True, help='Input protein name')
    args = parser.parse_args()

    return args

def translate_csv_file(protein):
    content_in_rows = []

    file_content = csv.reader(open(protein+'_ddG.csv', 'rU'))

    content_in_rows.append('# Temp:'+'\t'*8+file_content.next()[1]+'\n')
    content_in_rows.append('# Temp units (\'deg_C\' or \'K\')'+'\t'*5+file_content.next()[1]+'\n')
    content_in_rows.append('# Energy units (\'kcal\' or \'kJ\'):'+'\t'*4+file_content.next()[1]+'\n')
    content_in_rows.append('# Kinetic constants type (\'k_value\' or \'log_value\'):'+'\t'*2+file_content.next()[1]+'\n')
    content_in_rows.append('# kf denaturant concentration (M):'+'\t'*4+file_content.next()[1]+'\n')
    content_in_rows.append('# ku denaturant concentration (M):'+'\t'*4+file_content.next()[1]+'\n')
    content_in_rows.append('# Denaturant equilibrium concentration (M):'+'\t'*3+file_content.next()[1]+'\n')
    content_in_rows.append('# Type of calculation: (direct or indirect):'+'\t'*3+file_content.next()[1]+'\n')
    content_in_rows.append('# (NOTE: if the kinetic constant type is \'log_value\', please replace  \'kf, e_kf, ku and e_ku\''+
                           ' by  \'ln_kf, e_ln_kf, ln_ku and  e_ln_ku\'    in the units line.\n')
    content_in_rows.append('# (NOTE: if the energy values are in kJ/mol, please replace the corresponding units)\n')
    content_in_rows.append('# (NOTE: \'-999\' is the placeholder for missing data)\n')
    content_in_rows.append('# (NOTE: In the case of the location, the valid options are \'none\', \'core\' or \'surf\'\n')
    content_in_rows.append('# location\tindex\tresA\tresB\tkf\te_kf\tku\te_ku\tmkf\te_mkf\tmku\te_mku\tddG0\te_ddG0\tphi_f\te_phi_f\texclude\n\
')
    content_in_rows.append('# \t-\t-\t-\t-\t(1/s)\t(1/s)\t(1/s)\t(1/s)\t(1/M)\t(1/M)\t(1/M)\t(1/M)\t(kcal/mol)\t-\t-\t-\n')

    for row in file_content:
        if list(row[0])[0]=='#':
            continue
        line_to_append ='\t'+'\t'.join(row)+'\n'
        content_in_rows.append(line_to_append)

    with open(protein+'_ddG.dat', 'w') as datfile:
        for mutation in content_in_rows:
            datfile.write(mutation)

class ddG_data:
    def __init__(self, location=[],index=[], resA=[], resB=[], ddG0=[], e_ddG0=[], ddGdag=[], e_ddGdag=[], usable=[], phi=[], e_phi=[], exclude=[]):
        # units of kT
        self.location = location
        self.index = index
        self.resA = resA
        self.resB = resB
        self.ddG0 = ddG0
        self.e_ddG0 = e_ddG0
        self.ddGdag = ddGdag
        self.e_ddGdag = e_ddGdag
        self.usable = usable
        self.phi = phi
        self.e_phi = e_phi
        self.exclude = exclude

# All values in units of kT
def calc_ddGdag(wt_ln_kf, ln_kf, wt_mkf, mkf, CM, calculation_type):
    if calculation_type=='direct':
        ddGdag = wt_ln_kf - ln_kf
    else:
        ddGdag = wt_ln_kf - ln_kf - CM * (wt_mkf - mkf)
    return ddGdag

def calc_e_ddGdag(wt_e_ln_kf, e_ln_kf, wt_e_mkf, e_mkf, CM, calculation_type):
    # Using the arithmetic version of the error
    if calculation_type=='direct':
        e_ddGdag = wt_e_ln_kf + e_ln_kf
    else:
        e_ddGdag = wt_e_ln_kf + e_ln_kf + CM *(wt_e_mkf + e_mkf)
    return e_ddGdag

def calc_ddG0(ln_kf, ln_ku, mkf, mku, CM, ku_denaturant_conc, ddG0_eq_f, calculation_type, temperature, kB, placeholder_list):
    if calculation_type=='direct':
        if ddG0_eq_f in placeholder_list:
            ddG0 = 0.
        else:
            ddG0 = ddG0_eq_f/(temperature * kB)
    else:
        ddG0 = ln_ku - ln_kf + CM * (mku + mkf) - ku_denaturant_conc * mku
    return ddG0

def calc_e_ddG0(e_ln_kf, e_ln_ku, e_mkf, e_mku, CM, ku_denaturant_conc, e_ddG0_eq_f, calculation_type, temperature, kB, placeholder_list):
    # Using the arithmetic version of the error  
    if calculation_type=='direct':
        e_ddG0 = e_ddG0_eq_f/(temperature * kB)
        if e_ddG0_eq_f in placeholder_list:
            e_ddG0 = 0.
    else:
        e_ddG0 = e_ln_ku + e_ln_kf + CM * (e_mku + e_mkf) + ku_denaturant_conc * e_mku
    return e_ddG0

def calc_phi(item_ddGdag, item_ddG0, phi_f, calculation_type, placeholder_list):
    if calculation_type=='direct' and phi_f not in placeholder_list:
        phi = phi_f
    else:
        # Just for informative purposes, we will include all phi values (even if ddG0 is very low)
        if item_ddG0!=0:
            phi = item_ddGdag / item_ddG0
        else:
            phi = 0.
    
    return phi
    
def calc_e_phi(item_ddGdag, item_e_ddGdag, item_ddG0, item_e_ddG0, e_phi_f, calculation_type, placeholder_list):
    # Using the arithmetic version of the error  
    if calculation_type=='direct' and e_phi_f not in placeholder_list:
        e_phi = e_phi_f
    else:
        if item_e_ddG0 != 0.:
            e_phi = abs((item_e_ddGdag * item_ddG0 + item_e_ddG0 * item_ddGdag)/ np.square(item_ddG0))
        else:
            e_phi = 0.

    return e_phi

def import_ddG_file(protein):
    file_name = protein + '_ddG.dat'

    with open(file_name, 'r') as ddG:
        ddG_file = ddG.readlines()
    ddG.close()

    # Some exception handling provisions
    units_dict_kcal = {'index_number':'', 'resA':'', 'resB':'', 'kf': '1/s', 'e_kf':'1/s', 'ku':'1/s', 'e_ku':'1/s',
                       'mkf':'1/M','e_mkf':'1/M', 'mku':'1/M', 'e_mku':'1/M', 'ddG0_eq':'kcal/mol', 'e_ddG0_eq':'kcal/mol'}
    units_dict_kJ = {'index_number':'', 'resA':'', 'resB':'', 'kf': '1/s', 'e_kf':'1/s', 'ku':'1/s', 'e_ku':'1/s',
                       'mkf':'1/M','e_mkf':'1/M', 'mku':'1/M', 'e_mku':'1/M', 'ddG0_eq':'kJ/mol', 'e_ddG0_eq':'kJ/mol'}
    temperature_list = ['deg_C', 'K']
    placeholder_list = ['-999', -999, -999.000000]
    k_type_options = ['k_value', 'log_value']
    location_options = ['none','core','surf']
    calculation_type_options = ['direct', 'indirect']
    
    
    # Read experimental temperature
    if ddG_file[1].split()[-1] =='K':
        temperature = float(ddG_file[0].split()[-1])
    elif ddG_file[1].split()[-1] == 'deg_C':
        temperature = 273.15 + float(ddG_file[0].split()[-1])
    else:
        sys.exit('Please give a valid option for the temperature units')
        
    # Read energy units to determine kB
    if ddG_file[2].split()[-1] =='kcal':
        kB = 0.0019872  # Boltzmann constant in kcal/(mol.K)
        units_dict = units_dict_kcal
    elif ddG_file[2].split()[-1] =='kJ':
        kB = 0.0083145  # in kJ/(mol.K)
        units_dict = units_dict_kJ
    else: 
        sys.exit('Please give a valid option for the energy units')

    # Determine if kinetic constants are given in terms of their log or not ('k_value' or 'log_value'):
    k_type =ddG_file[3].split()[-1]
   
    if k_type not in k_type_options:
        sys.exit('Please give a valid option for the kinetic constants type')

    
    # Read kf and ku denaturant concentration:
    kf_denaturant_conc = float(ddG_file[4].split()[-1])
    ku_denaturant_conc = float(ddG_file[5].split()[-1])

    # Read denaturant concentration at vertex of chevron plot
    CM = float(ddG_file[6].split()[-1])
    calculation_type = ddG_file[7].split()[-1]
    if calculation_type not in calculation_type_options:
        sys.exit('Please give a valid option for the calculation type')

    location = []
    index =[]
    resA = []
    resB = []
    ddG0 = []
    e_ddG0 = []
    ddGdag = []
    e_ddGdag = []
    usable = []
    phi = []
    e_phi = []
    exclude = []
    
    for line in ddG_file:
        # Obviate the first few rows that do not contain mutations data
        if line.split()[0]=='#':
            continue
    
        values = line.split()
        location_value = values[0]
        if location_value not in location_options:
            sys.exit('Please verify that all mutations have a valid location option')
            
        index_number = values[1]
        if index_number == 'wt':
            # Store wild-type data first, with missing data check for uncertainties
            wt_kf = float(values[4])
            wt_e_kf = float(values[5])
            if wt_e_kf in placeholder_list:
                wt_e_kf = 0.
            wt_ku = float(values[6])
            wt_e_ku = float(values[7])
            if wt_e_ku in placeholder_list:
                wt_e_ku = 0.
            wt_mkf = float(values[8])
            wt_e_mkf = float(values[9])
            if wt_e_mkf in placeholder_list:
                wt_e_mkf = 0.
            wt_mku = float(values[10])
            wt_e_mku = float(values[11])
            if wt_e_mku in placeholder_list:
                wt_e_mku = 0.
            

            if any(x in placeholder_list for x in [wt_kf,wt_mkf]) and calculation_type=='indirect':
                sys.exit('The program cannot run without kf or ln(kf) and mkf for the wild type')
            if wt_kf in placeholder_list and calculation_type=='direct':
                sys.exit('The program cannot run without kf for the wild type')
               
            if k_type =='log_value':
                wt_ln_kf = wt_kf
                wt_e_ln_kf = wt_e_kf
            else: 
                wt_ln_kf = np.log(wt_kf)
                wt_e_ln_kf = wt_e_kf/wt_kf
            continue
        
        else:
            #All mutations
            residue_A =  values[2]
            residue_B = values[3]
            kf = float(values[4])
            e_kf = float(values[5])
            ku = float(values[6])
            e_ku = float(values[7])
            mkf = float(values[8])
            e_mkf = float(values[9])
            mku = float(values[10])
            e_mku = float(values[11])
            #The following values are only used for the direct calculation_type
            ddG0_eq_f = float(values[12])
            e_ddG0_eq_f = float(values[13])
            phi_f = float(values[14])
            e_phi_f = float(values[15])
            item_exclude = int(values[16])
        # We will append all mutations that are not the wild type 
        location.append(location_value)
        index.append(index_number)
        resA.append(residue_A)
        resB.append(residue_B)

        # Exception handling: kf, mkf, ku and mku (or their logarithmic equivalent) must be given for mutation to be usable
        if any(x in placeholder_list for x in [kf, mkf, ku, mku]) and calculation_type=='indirect':
            print 'Mutation '+ residue_A + index_number + residue_B + ' cannot be used because either kf, ku, mkf or mku are missing'
            ddG0.append(0)
            e_ddG0.append(0)
            ddGdag.append(0)
            e_ddGdag.append(0)
            usable.append(False)
            phi.append(0)
            e_phi.append(0)
            exclude.append(item_exclude)             
            continue
        
        if any(x in placeholder_list for x in [kf]) and calculation_type=='direct':
            print 'Mutation '+ residue_A + index_number + residue_B + ' cannot be used because kf is missing'
            ddG0.append(0)
            e_ddG0.append(0)
            ddGdag.append(0)
            e_ddGdag.append(0)
            usable.append(False)
            phi.append(0)
            e_phi.append(0)
            exclude.append(item_exclude)
            continue
        
        if all([any(x in placeholder_list for x in [ddG0_eq_f]), calculation_type=='direct', location_value not in ['surf']]):
            print 'Mutation '+ residue_A + index_number + residue_B + ' cannot be used because kf is missing'
            ddG0.append(0)
            e_ddG0.append(0)
            ddGdag.append(0)
            e_ddGdag.append(0)
            usable.append(False)
            phi.append(0)
            e_phi.append(0)
            exclude.append(item_exclude)
            continue


        # Mutations that pass the previous filter are usable
        usable.append(True)

        # Calculate the logarithm of the kinetic constants if given as k_value
        if k_type =='log_value':
            ln_kf = kf
            e_ln_kf = e_kf
            ln_ku = ku
            e_ln_ku = e_ku
        else: 
            ln_kf = np.log(kf)
            e_ln_kf = e_kf/kf
            try:    
                ln_ku = np.log(ku)
            except: 
                ln_ku = 0.
            try:
                e_ln_ku = e_ku/ku
            except:
                e_ln_ku = 0.

        # Error calculation is not indispensable, but then we must reset the respective uncertainties to zero if not available
        if e_ln_kf in placeholder_list:
            e_ln_kf = 0.

        if e_ln_ku in placeholder_list:
            e_ln_ku = 0.

        if e_mkf in placeholder_list:
            e_mkf = 0.

        if e_mku in placeholder_list:
            e_mku = 0.
        
        if e_ddG0_eq_f in placeholder_list:
            e_ddG_eq_f = 0.
        
        if e_phi_f in placeholder_list:
            e_phi_f = 0.

        # Calculate parameters to yield in final .dat file
        # Variables within the scope of this module are preceded by "item_" so that they are not confused with the corresponding lists
        item_ddG0 = calc_ddG0(ln_kf, ln_ku, mkf, mku, CM, ku_denaturant_conc, ddG0_eq_f, calculation_type, temperature, kB, placeholder_list)
        item_e_ddG0 = calc_e_ddG0(e_ln_kf, e_ln_ku, e_mkf, e_mku, CM, ku_denaturant_conc, e_ddG0_eq_f, calculation_type, temperature, kB, placeholder_list)
        item_ddGdag = calc_ddGdag(wt_ln_kf, ln_kf, wt_mkf, mkf, CM, calculation_type)
        item_e_ddGdag = calc_e_ddGdag(wt_e_ln_kf, e_ln_kf, wt_e_mkf, e_mkf, CM, calculation_type)
        item_phi = calc_phi(item_ddGdag, item_ddG0, phi_f, calculation_type, placeholder_list)
        item_e_phi = calc_e_phi(item_ddGdag, item_e_ddGdag, item_ddG0, item_e_ddG0, e_phi_f, calculation_type, placeholder_list)
        
        ddG0.append(item_ddG0)
        e_ddG0.append(item_e_ddG0)
        ddGdag.append(item_ddGdag)
        e_ddGdag.append(item_e_ddGdag)
        phi.append(item_phi)
        e_phi.append(item_e_phi)
        exclude.append(item_exclude)
        
    # Now, do the trick for the surface mutations:
    i = 0
    n_muts = len(location)-1
    while i < n_muts:
        if location[i]=='surf' and exclude[i]==0:
            # Mutations from A or G are already good to use, otherwise we need to calculate the difference
            if resA[i]!='A' and resA[i]!='G':
                # We need to see if the surf mutation came in the needed pairs
                if index[i]==index[i+1]:
                    # Append the A->G surface mutation at the end of the list
                    location.append('surf')
                    index.append(index[i])
                    resA.append('A')
                    resB.append('G')
                    # Check that they are in the correct order
                    if resB[i]=='A' and resB[i+1]=='G':
                        first = i
                        second = i+1
                    elif resB[i]=='G' and resB[i+1]=='A':
                        first = i+1
                        second = i
                    # Calculate the ddG of the A->G mutation
                    ddG0.append(ddG0[second]-ddG0[first])
                    e_ddG0.append(e_ddG0[second]+e_ddG0[first])
                    ddGdag.append(ddGdag[second]-ddGdag[first])
                    e_ddGdag.append(e_ddGdag[second]+e_ddGdag[first])
                    usable.append(True)
                    phi.append(calc_phi(ddGdag[len(location)-1], ddG0[len(location)-1], phi_f, 'indirect', placeholder_list))
                    e_phi.append(calc_e_phi(ddGdag[len(location)-1], e_ddGdag[len(location)-1],
                                            ddG0[len(location)-1], e_ddG0[len(location)-1], e_phi_f, 'indirect', placeholder_list))
                    exclude.append(0)
                    # Turn the original pair into usable=False
                    usable[i]=False
                    usable[i+1]=False
                    i+=2
                else:
                    usable[i]=False
                    i+=1
            else:
                i+=1
        else:
            i+=1
                    
    Calculated_ddG_list = ddG_data(location, index, resA, resB, ddG0, e_ddG0, ddGdag, e_ddGdag, usable, phi, e_phi, exclude)
    
    # Save ddG_data class object to a pickle file to be accessed later, and also to a .dat file to be human readable
    
    dat_file = open(protein+'_calculated_ddG.dat', 'w') 
    a = '  '
    b = '   '
    c = '    '
    dat_file.write('#loc\tindex\tresA\tresB'+2*c+' ddG0'+2*c+b+'e_ddG0'+2*c+' ddGdag'+2*c+a+'e_ddGdag'+3*c+'usable'+2*c+'phi'+3*c+'e_phi'+c+a+'exclude\n')
    for i in range(len(Calculated_ddG_list.location)):
        dat_file.write('%s\t%s\t%s\t%s\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t\t%s\t%10.4f\t%10.4f\t%d\n'% (Calculated_ddG_list.location[i], Calculated_ddG_list.index[i],
                                                                                                  Calculated_ddG_list.resA[i], Calculated_ddG_list.resB[i], 
                                                                                                  Calculated_ddG_list.ddG0[i], Calculated_ddG_list.e_ddG0[i], 
                                                                                                  Calculated_ddG_list.ddGdag[i], Calculated_ddG_list.e_ddGdag[i],
                                                                                                  Calculated_ddG_list.usable[i],
                                                                                                  Calculated_ddG_list.phi[i], Calculated_ddG_list.e_phi[i],
                                                                                                  Calculated_ddG_list.exclude[i]))
    dat_file.close()
    
def main():
    #find out protein from input... should be replaced by reading system.info                                                     
    protein = get_args().prot
    #See if the input .csv file exists
    try:
        translate_csv_file(protein)
    finally:
        import_ddG_file(protein)

if __name__ == '__main__':
    main()
