''' A recipe to apply the Matysiak Clementi 2004 algorithm


Description:

    This recipes performs the algorithm in reference (1) to add energetic
heterogeneity to a Go-model.


Reference:

(1) Matysiak, S.; Clementi, C. Optimal Combination of Theory and Experiment for
the Characterization of the Protein Folding Landscape of S6: How Far Can a
Minimalist model Go? J. Mol. Biol. 2004, 343, 235-248.
'''

import os
import argparse
import numpy as np

from project_tools import simulation, analysis, parameter_fitting
from recipe_manager import ProjectManager
import model_builder as mdb


class MatysiakClementi2004(ProjectManager):
    
    ''' A project manager to reproduce Matysiak Clementi 2004 algorithm. 


    Description:

        This project manager is intended to automate the generation of a
    heterogeneous Go-model via the thermodynamic perturbation algorithm
    put forth by Matysiak & Clementi (1).
        

    Recipe:
    1. Run constant temperature runs to calculate heat capacity Cv(T). Identify
        the folding temperature Tf as the peak in Cv(T).
    2. Run longer simulations around Tf to witness folding/unfolding transitions. 
    3. Use equilibrium data and wham to get F(Q) at different temperatures. Find
        a more accurate folding temperature as temperature at which the folded
        and unfolded state have same free energy.
    3. Prepare the experimental data to get <name>_calculated_ddG.dat input file.
    3. Mutate pdbs with MODELLER and calculate fraction of heavy atom contact
        loss for each contact and each mutation.

    Reference:

    (1) Matysiak, S.; Clementi, C. Optimal Combination of Theory and Experiment for
    the Characterization of the Protein Folding Landscape of S6: How Far Can a
    Minimalist model Go? J. Mol. Biol. 2004, 343, 235-248.
    '''


    def logical_flowchart_starting(self,model,task):
        sub = model.subdir
        if task == "Tf_loop_iteration":
            print "Checking if Tf_loop_iteration completed..."
            simulation.constant_temp.check_completion(model,self.append_log)
            lasttime2,action2,task2 = self.check_modelbuilder_log(sub)
            if action2 == "Finished:":
                print "Finished Tf_loop_iteration..."
                print "Starting Tf_loop_analysis..."
                analysis.constant_temp.analyze_temperature_array(model,self.append_log)
        elif task == "Tf_loop_analysis":
            print "Checking if Tf_loop_analysis completed..."
            analysis.constant_temp.check_completion(model,self.append_log)
        elif task == "Equil_Tf":
            print "Starting to check if Equil_Tf completed..."
            simulation.constant_temp.check_completion(model,self.append_log,long=True)
            lasttime2,action2,task2 = self.check_modelbuilder_log(sub)
            if action2 == "Finished:":
                print "Finished Equil_Tf_iteration..."
                print "Starting Equil_Tf_analysis..."
                analysis.constant_temp.analyze_temperature_array(model,self.append_log,long=True)
        elif task == "Equil_Tf_analysis":
            print "Starting to check if Equil_Tf_analysis completed..."
            analysis.constant_temp.check_completion(model,self.append_log,long=True)
        else:
            print "ERROR!"
            print "  Couldn't find next option for task:",task
            print "  Please check that things are ok."
            print "  Exiting."
            raise SystemExit

    def logical_flowchart_finished(self,model,task):
        sub = model.subdir
        if task == "Tf_loop_iteration":
            print "Finished Tf_loop_iteration..."
            print "Starting Tf_loop_analysis..."
            analysis.constant_temp.analyze_temperature_array(model,self.append_log)
        elif task == "Tf_loop_analysis":
            print "Finished Tf_loop_analysis..."
            flag = analysis.constant_temp.run_wham_heat_capacity(model,self.append_log)
            if flag == 1:
                pass 
            else:
                print "Starting Tf_loop_iteration..."
                simulation.constant_temp.folding_temperature_loop(model,self.append_log)
        elif task == "Tf_wham":
            print "Starting equilibrium simulations at Tf..."
            simulation.constant_temp.run_equilibrium_simulations(model,self.append_log)
        elif task == "Equil_Tf":
            print "Starting Equil_Tf_analysis..."
            analysis.constant_temp.analyze_temperature_array(model,self.append_log,long=True)
        elif task == "Equil_Tf_analysis":
            ## Use the following sub module to plot PMFS of coordinates:
            ## analysis.plot.pmfs
            ## Run heat capacity for equilibrium runs. Cv(T), F(Q)
            analysis.constant_temp.run_wham_heat_capacity(model,self.append_log,long=True)
        elif task == "Equil_Tf_wham":
            print "Starting calculating feature vector and Jacobian..."
            parameter_fitting.prepare_newtons_method(model,"ddG_MC2004",self.append_log)
        elif task == "Calculating_Jacobian":
            print "Solving for solutions with Levenberg-Marquardt method..."
            parameter_fitting.solve_newtons_method(model,"ddG_MC2004",self.append_log)
        elif task == "Solving_Newtons_Method":
            ## Write new parameter file
            ## Start the next round of simulations with new parameters.
            parameter_fitting.save_new_parameters(model,"ddG_MC2004",self.append_log)
            simulation.constant_temp.start_next_Tf_loop_iteration(model,self.append_log)
        else:
            print "ERROR!"
            print "  Couldn't find next option for task:",task
            print "  Please check that things are ok."
            print "  Exiting."
            raise SystemExit

    def new_project(self,args,modeloptions):
        ''' Start a new simulation project'''

        subdirs = [ x[:-4] for x in args.pdbs ]
        for sub in subdirs:
            if os.path.exists(sub) == False:
                os.mkdir(sub)
            else:
                print "Subdirectory: %s already exists! just fyi" % sub

        print "Starting a new simulation project..."
        Models = mdb.check_inputs.new_models(subdirs,modeloptions)

        self.save_model_info(Models)
        if args.temparray != None:
            for n in range(len(subdirs)):
                Models[n].initial_T_array = args.temparray

        for k in range(len(Models)):
            model = Models[k]
            open("%s/Native.pdb" % model.subdir,"w").write(model.cleanpdb)
            open("%s/clean.pdb" % model.subdir,"w").write(model.cleanpdb_full)
            open("%s/clean_noH.pdb" % model.subdir,"w").write(model.cleanpdb_full_noH)
            open("%s/%s.pdb" % (model.subdir,model.subdir),"w").write(model.cleanpdb_full_noH)
            np.savetxt("%s/contact_map.dat" % (model.subdir),model.Qref,delimiter=" ",fmt="%1d")
            np.savetxt("%s/contacts.dat" % (model.subdir),model.contacts,delimiter=" ",fmt="%4d")

        for k in range(len(Models)):
            model = Models[k]
            print "Starting Tf_loop_iteration for subdirectory: ", model.subdir
            simulation.constant_temp.folding_temperature_loop(model,self.append_log,new=True)

        self.save_model_info(Models)
        print "Success"


def get_args():
    ''' Get command line arguments '''

    parser = argparse.ArgumentParser(description='Options for MatysiakClementi2004 recipe.')
    sp = parser.add_subparsers(dest='action')

    ## Options for initializing a new simulation project.
    new_parser = sp.add_parser('new')
    new_parser.add_argument('--pdbs', type=str, required=True, nargs='+',help='PDBs to start simulations.')
    new_parser.add_argument('--pairwise_params_file', type=str, default="None", help='Optional, specify pairwise interactions.')
    new_parser.add_argument('--model_params_file', type=str, default="None", help='Optional, specify model parameters.')
    new_parser.add_argument('--epsilon_bar', type=float,help='Optional, average strength of contacts. epsilon bar.')
    new_parser.add_argument('--fitting_solver', type=str, default="Levenberg", help='Optional, specify solution algorithm for fitting.')
    new_parser.add_argument('--fitting_includes', type=str, nargs='+', default="None", help='Optional, specify directories included in fitting.')
    new_parser.add_argument('--fitting_allowswitch', type=str, default="False", help='Optional, allow contacts to switch between attractive/repulsive in fitting.')
    new_parser.add_argument('--disulfides', type=int, nargs='+', help='Optional pairs of disulfide linked residues.')
    new_parser.add_argument('--temparray', type=int, nargs='+',help='Optional initial temp array: T_min T_max deltaT. Default: predicts now')
    new_parser.add_argument('--dry_run', action='store_true', help='Add this option for dry run. No simulations started.')

    ## Options for continuing from a previously saved simulation project.
    run_parser = sp.add_parser('continue')
    run_parser.add_argument('--subdirs', type=str, nargs='+', help='Subdirectories to continue',required=True)
    run_parser.add_argument('--dry_run', action='store_true', help='Dry run. No simulations started.')

    ## Options for manually adding a temperature array.
    add_parser = sp.add_parser('add')
    add_parser.add_argument('--subdirs', type=str, nargs='+', help='Subdirectories to add temp array',required=True)
    add_parser.add_argument('--short_temps', type=int, nargs='+', help='T_initial T_final dT for new temps')
    add_parser.add_argument('--long_temps', type=float, nargs='+', help='List of temps for new equilibrium simulations')
    add_parser.add_argument('--dry_run', action='store_true', help='Dry run. No simulations started.')

    ## Options for manually extending some temperatures.
    ext_parser = sp.add_parser('extend')
    ext_parser.add_argument('--subdirs', type=str, nargs='+', help='Subdirectories to add temp array',required=True)
    ext_parser.add_argument('--factor', type=float, help='Factor by which you want to extend simulations. e.g. --factor 2 doubles length',required=True)
    ext_parser.add_argument('--short_temps', type=float, nargs='+', help='Temperatures that you want extended')
    ext_parser.add_argument('--long_temps', type=float, nargs='+', help='T_initial T_final dT for new mutational sims array')
    ext_parser.add_argument('--dry_run', action='store_true', help='Dry run. No simulations started.')

    args = parser.parse_args()

    args.fitting_data = "ddG_MC2004"
    args.model_code = "HetGo"
    args.bead_model = "CA"

    if args.action == "new":
        modeloptions = mdb.check_inputs.new_args(args)
    else:
        modeloptions = []
        
    return args, modeloptions

if __name__ == "__main__":
    
    args, options = get_args()
    MatysiakClementi2004(args,options)
