""" A recipe to apply the FRETFit algorithm


Description:

    This recipes performs an algorithm to add energetic heterogeneity to a
Go-model in order to reproduce a FRET pairwise distance distribution.


Reference:
"""

import os
import argparse

from project_tools import simulation, analysis, parameter_fitting
from recipe_manager import ProjectManager
import model_builder as mdb


class FRETFit(ProjectManager):
    
    """ A recipe to apply the FRETFit algorithm


    Description:

        This recipes performs an algorithm to add energetic heterogeneity to a
    Go-model in order to reproduce a FRET pairwise distance distribution. The
    algorithm was inspired by Matysiak


    Reference:
    """


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
            simulation.constant_temp.check_completion(model,self.append_log,equil=True)
            lasttime2,action2,task2 = self.check_modelbuilder_log(sub)
            if action2 == "Finished:":
                print "Finished Equil_Tf_iteration..."
                print "Starting Equil_Tf_analysis..."
                analysis.constant_temp.analyze_temperature_array(model,self.append_log,equil=True)
        elif task == "Equil_Tf_analysis":
            print "Starting to check if Equil_Tf_analysis completed..."
            analysis.constant_temp.check_completion(model,self.append_log,equil=True)
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
        #elif task == "Tf_loop_analysis":
        #    print "Finished Tf_loop_analysis..."
        #    flag = analysis.constant_temp.run_wham_heat_capacity(model,self.append_log)
        #    if flag == 1:
        #        pass 
        #    else:
        #        print "Starting Tf_loop_iteration..."
        #        simulation.constant_temp.folding_temperature_loop(model,self.append_log)
        #elif task == "Tf_wham":
        #    print "Starting equilibrium simulations at Tf..."
        #    simulation.constant_temp.run_equilibrium_simulations(model,self.append_log)
        #elif task == "Equil_Tf":
        #    print "Starting Equil_Tf_analysis..."
        #    analysis.constant_temp.analyze_temperature_array(model,self.append_log,equil=True)
        #elif task == "Equil_Tf_analysis":
        ### Use the following sub module to plot PMFS of coordinates:
        ### analysis.plot.pmfs
        #    ## Run heat capacity for equilibrium runs. Cv(T), F(Q)
        #    analysis.constant_temp.run_wham_heat_capacity(model,self.append_log,Mut=True)
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
        Models = mdb.models.new_models(subdirs,modeloptions)

        self.save_model_info(Models)
        if args.temparray != None:
            for n in range(len(subdirs)):
                Models[n].initial_T_array = args.temparray

        for k in range(len(Models)):
            model = Models[k]
            print "Starting Tf_loop_iteration for subdirectory: ", model.subdir
            simulation.constant_temp.folding_temperature_loop(model,self.append_log,new=True)

        self.save_model_info(Models)
        print "Success"


def get_args():
    """ Get command line arguments """

    parser = argparse.ArgumentParser(description='Run .')
    sp = parser.add_subparsers(dest='action')

    ## Options for initializing a new simulation project.
    new_parser = sp.add_parser('new')
    new_parser.add_argument('--pdbs', type=str, required=True, nargs='+',help='PDBs to start simulations.')
    new_parser.add_argument('--epsilon_bar', type=float, help='Optional, average strength of contacts. epsilon bar.')
    new_parser.add_argument('--disulfides', type=int, nargs='+', help='Optional pairs of disulfide linked residues.')
    new_parser.add_argument('--temparray', type=int, nargs='+',help='Optional initial temp array: Ti Tf dT. Default: 50 350 50')
    new_parser.add_argument('--dryrun', action='store_true', help='Add this option for dry run. No simulations started.')

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

    ## Options for manually extending some temperatures.
    ext_parser = sp.add_parser('extend')
    ext_parser.add_argument('--subdirs', type=str, nargs='+', help='Subdirectories to add temp array',required=True)
    ext_parser.add_argument('--factor', type=float, help='Factor by which you want to extend simulations. e.g. --factor 2 doubles length',required=True)
    ext_parser.add_argument('--Tf_temps', type=float, nargs='+', help='Temperatures that you want extended')
    ext_parser.add_argument('--Mut_temps', type=float, nargs='+', help='T_initial T_final dT for new mutational sims array')
    ext_parser.add_argument('--dryrun', action='store_true', help='Dry run. No simulations started.')


    args = parser.parse_args()

    if args.dryrun != False:
        options = {"Dry_Run":True}
    else:
        options = {"Dry_Run":False}

    if args.action == "new":
        if args.epsilon_bar != False:
            options["Epsilon_Bar"] = args.epsilon_bar
        else:
            options["Epsilon_Bar"] = None
        if args.disulfides != False:
            options["Disulfides"] = args.disulfides
        else:
            options["Disulfides"] = None
    else:
        options["Epsilon_Bar"] = None
        options["Disulfides"] = None

    options["Fitting_Data"] = "FRET"
    options["Model_Code"] = "HetGo"
    options["Bead_Model"] = "CA"
    options["Contact_Energies"] = "FRETFit"

    modeloptions = mdb.models.check_options(options)

    return args, modeloptions

if __name__ == "__main__":
    
    args, options = get_args()
    FRETFit(args,options)
