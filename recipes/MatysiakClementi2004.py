""" A recipe to apply the Matysiak Clementi 2004 algorithm


Description:

    This recipes performs the algorithm in reference (1) to add energetic
heterogeneity to a Go-model.


Reference:

(1) Matysiak, S.; Clementi, C. Optimal Combination of Theory and Experiment for
the Characterization of the Protein Folding Landscape of S6: How Far Can a
Minimalist Model Go? J. Mol. Biol. 2004, 343, 235-248.
"""

import os
import argparse

from project_tools.manager import ProjectManager
from project_tools import simulation, analysis, mutations
from model_builder import models
from model_builder import systems


class MatysiakClementi2004(ProjectManager):
    
    """ A project manager to reproduce Matysiak Clementi 2004 algorithm. 


    Description:

        This project manager is intended to automate the generation of a
    heterogeneous Go-model via the thermodynamic perturbation algorithm
    put forth by Matysiak & Clementi (1).
        

    Reference:

    (1) Matysiak, S.; Clementi, C. Optimal Combination of Theory and Experiment for
    the Characterization of the Protein Folding Landscape of S6: How Far Can a
    Minimalist Model Go? J. Mol. Biol. 2004, 343, 235-248.
    """


    def logical_flowchart_starting(self,System,Model,sub,task):
        if task == "Tf_loop_iteration":
            print "Checking if Tf_loop_iteration completed..."
            simulation.Tf_loop.check_completion(System,self.append_log)
            lasttime2,action2,task2 = self.check_modelbuilder_log(sub)
            if action2 == "Finished:":
                print "Finished Tf_loop_iteration..."
                print "Starting Tf_loop_analysis..."
                analysis.Tf_loop.analyze_temperature_array(System,self.append_log)
        elif task == "Tf_loop_analysis":
            print "Checking if Tf_loop_analysis completed..."
            analysis.Tf_loop.check_completion(System,self.append_log)
        elif task == "wham_Cv":
            print "Starting to check if wham_Cv completed..."
            analysis.Tf_loop.continue_wham(System,self.append_log)
        elif task == "wham_FreeEnergy":
            print "Starting Equil_Tf..."
            simulation.Tf_loop.run_equilibrium_simulations(Model,System,self.append_log)
        elif task == "Equil_Tf":
            print "Starting to check if Equil_Tf completed..."
            simulation.Tf_loop.check_completion(System,self.append_log,equil=True)
            lasttime2,action2,task2 = self.check_modelbuilder_log(sub)
            if action2 == "Finished:":
                print "Finished Equil_Tf_iteration..."
                print "Starting Equil_Tf_analysis..."
                analysis.Tf_loop.analyze_temperature_array(System,self.append_log,equil=True)
        elif task == "Equil_Tf_analysis":
            print "Starting to check if Equil_Tf_analysis completed..."
            analysis.Tf_loop.check_completion(System,self.append_log,equil=True)
        elif task == "Calculating_MC2004":
            print "ERROR!"
            print "  ",task, " should have finished!"
            print "  Please check that things are ok."
            print "  Exiting."
            raise SystemExit
        else:
            print "ERROR!"
            print "  Couldn't find next option for task:",task
            print "  Please check that things are ok."
            print "  Exiting."
            raise SystemExit

    def logical_flowchart_finished(self,System,Model,sub,task):
        if task == "Tf_loop_iteration":
            print "Finished Tf_loop_iteration..."
            print "Starting Tf_loop_analysis..."
            analysis.Tf_loop.analyze_temperature_array(System,self.append_log)
        elif task == "Tf_loop_analysis":
            print "Finished Tf_loop_analysis..."
            flag = analysis.Tf_loop.check_if_wham_is_next(System,self.append_log)
            if flag == 1:
                pass 
            else:
                print "Starting Tf_loop_iteration..."
                simulation.Tf_loop.folding_temperature_loop(Model,System,self.append_log)
        elif task == "wham_Cv":
            print "Finished wham_Cv..."
            print "Stating wham_FreeEnergy..."
            analysis.Tf_loop.continue_wham(System,self.append_log)
        elif task == "Equil_Tf":
            print "Starting Equil_Tf_analysis..."
            analysis.Tf_loop.analyze_temperature_array(System,self.append_log,equil=True)
        elif task == "Equil_Tf_analysis":
            ## Aggregrate equil_Tf data for each temperature and plot PMFs
            print "Starting aggregate data..."
            analysis.Tf_loop.aggregate_equilibrium_runs(System,self.append_log)
            print "Plotting aggregated data PMFS..."
            analysis.plot.pmfs.plot_aggregated_data(System,self.append_log)
        elif task == "Aggregating_Equil_Runs":
            ## If plotting diddn't work before
            print "Plotting aggregated data PMFS..."
            analysis.plot.pmfs.plot_aggregated_data(System,self.append_log)
        elif task == "Plotting_Agg_Data":
            print "Starting prepping mutant pdbs..."
            mutations.mutatepdbs.prepare_mutants(System,self.append_log)
        elif task == "Preparing_Mutants":
            print "Starting calculating dH for mutants..."
            mutations.phi_values.calculate_dH_for_mutants(Model,System,self.append_log)
        elif task == "Calculating_dH":
            mutations.phi_values.calculate_phi_values(Model,System,self.append_log,"Q")
        elif task == "Calculating_phi_values":
            mutations.perturbation.calculate_MC2004_perturbation(Model,System,self.append_log)
        elif task == "Calculating_MC2004":
            ## Start the next round of simulations with new parameters.
            simulation.Tf_loop.start_next_Tf_loop_iteration(Model,System,self.append_log)
        else:
            print "ERROR!"
            print "  Couldn't find next option for task:",task
            print "  Please check that things are ok."
            print "  Exiting."
            raise SystemExit

    def new_project(self,args,modeloptions):
        """ Start a new simulation project"""

        subdirs = [ x[:-4] for x in args.pdbs ]
        for sub in subdirs:
            if os.path.exists(sub) == False:
                os.mkdir(sub)
            else:
                print "Subdirectory: ", sub, " already exists! just fyi"

        print "Starting a new simulation project..."
        Models = models.new_models(subdirs,modeloptions)
        Systems = systems.new_systems(subdirs)

        self.prepare_systems(Models,Systems)
        self.save_model_system_info(Models,Systems,subdirs)

        if args.temparray != None:
            for n in range(len(subdirs)):
                Systems[n].initial_T_array = args.temparray

        ## The first step depends on the type of model.
        for k in range(len(subdirs)):
            print "Starting Tf_loop_iteration for subdirectory: ", subdirs[k]
            Model = Models[k]
            System = Systems[k]
            simulation.Tf_loop.folding_temperature_loop(Model,System,self.append_log,new=True)

        self.save_model_system_info(Models,Systems,subdirs)
        print "Success"


def get_args():
    """ Get command line arguments """

    parser = argparse.ArgumentParser(description='Run .')
    sp = parser.add_subparsers(dest='action')

    ## Options for initializing a new simulation project.
    new_parser = sp.add_parser('new')
    new_parser.add_argument('--pdbs', type=str, required=True, nargs='+',help='PDBs to start simulations.')
    new_parser.add_argument('--dryrun', action='store_true', help='Add this option for dry run. No simulations started.')
    new_parser.add_argument('--temparray', type=int, nargs='+',help='Optional initial temp array: Ti Tf dT. Default: 50 350 50')

    ## Options for continuing from a previously saved simulation project.
    run_parser = sp.add_parser('continue')
    run_parser.add_argument('--subdirs', type=str, nargs='+', help='Subdirectories to continue',required=True)
    run_parser.add_argument('--dryrun', action='store_true', help='Dry run. No simulations started.')

    ## Options for manually adding a temperature array.
    add_parser = sp.add_parser('add')
    add_parser.add_argument('--subdirs', type=str, nargs='+', help='Subdirectories to add temp array',required=True)
    add_parser.add_argument('--temparray', type=int, nargs='+', help='T_initial T_final dT for new temp array',required=True)
    add_parser.add_argument('--mutarray', type=int, nargs='+', help='T_initial T_final dT for new mutational sims array')
    add_parser.add_argument('--dryrun', action='store_true', help='Dry run. No simulations started.')

    args = parser.parse_args()

    if args.dryrun != False:
        options = {"Dry_Run":True}
    else:
        options = {"Dry_Run":False}

    options["Model_Code"] = "HetGo"
    options["Bead_Model"] = "CA"
    options["Solvent"] = None
    options["R_CD"] = None
    options["Epsilon_Bar"] = 1.0
    options["Disulfides"] = None
    options["Contact_Energies"] = "MC2004"

    modeloptions = models.check_options(options)

    return args, modeloptions

if __name__ == "__main__":
    
    args, options = get_args()
    MatysiakClementi2004(args,options)
