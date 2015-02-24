"""A recipe to reproduce Matysiak Clementi 2004 algorithm

Description:
------------

    This project manager is intended to automate the generation of a
heterogeneous Go-model via the thermodynamic perturbation algorithm
put forth by Matysiak & Clementi (1). The recipe is as follows:

1. Run constant temperature runs to calculate heat capacity Cv(T). Identify the
    folding temperature Tf as the peak in Cv(T).
2. Run longer simulations around Tf to witness folding/unfolding transitions. 
3. Use equilibrium data and wham to get F(Q) at different temperatures. Find a
    more accurate folding temperature as temperature at which the folded and
    unfolded state have same free energy.
3. Prepare the experimental data to get <name>_calculated_ddG.dat input file.
3. Mutate pdbs with MODELLER and calculate fraction of heavy atom contact loss
    for each contact and each mutation.

Reference:
----------
(1) Matysiak, S.; Clementi, C. Optimal Combination of Theory and Experiment for
the Characterization of the Protein Folding Landscape of S6: How Far Can a
Minimalist model Go? J. Mol. Biol. 2004, 343, 235-248.
"""

import os
import argparse
import logging
from numpy import savetxt

from project_tools import simulation, analysis, parameter_fitting
import model_builder as mdb

#############################################################################
# Start a new project
#############################################################################
def new_project(args):
    """ Start a new simulation project"""

    names = args.names
    for sub in names:
        if not os.path.exists(sub):
            os.mkdir(sub)
        else:
            raise IOError("Subdirectory: %s already exists!" % sub)

    print "Starting a new simulation project"
    Models,Fittingopts = mdb.inputs.new_models(names)

    if args.temparray != None:
        for n in range(len(names)):
            Models[n].initial_T_array = args.temparray

    for k in range(len(Models)):
        model = Models[k]
        model.dry_run = args.dry_run
        os.chdir("%s" % model.name)
        open("Native.pdb","w").write(model.cleanpdb)
        open("clean.pdb","w").write(model.cleanpdb_full)
        open("clean_noH.pdb","w").write(model.cleanpdb_full_noH)
        open("%s.pdb" % model.name,"w").write(model.cleanpdb_full_noH)
        savetxt("contact_map.dat",model.Qref,delimiter=" ",fmt="%1d")
        savetxt("contacts.dat",model.pairs,delimiter=" ",fmt="%5d")
        savetxt("pairs.dat",model.pairs,delimiter=" ",fmt="%5d")
        os.chdir("..")

    for k in range(len(Models)):
        model = Models[k]
        fitopts = Fittingopts[k]
        #set default fitopts
        if fitopts["iteration"] == None:
            fitopts["iteration"] = 0
        iteration = fitopts["iteration"]
        print "Starting Tf_loop_iteration for %s: " % model.name
        simulation.constant_temp.folding_temperature_loop(model,iteration,new=True)
        fitopts["last_completed_task"] = "Starting: Tf_loop_iteration"
        mdb.inputs.save_model(model,fitopts)

#############################################################################
# Continue an existing project
#############################################################################
def continue_project(args):
    """ Checks where something left off and continues it."""
    names = args.names
    Models,Fittingopts = mdb.inputs.load_models(names,dry_run=args.dry_run)

    for i in range(len(names)):
        model = Models[i]
        fitopts = Fittingopts[i]
        print "Checking progress for directory:  %s" % model.name
        print "Last task was %s" % fitopts["last_completed_task"]
        action = fitopts["last_completed_task"].split()[0]
        task = fitopts["last_completed_task"].split()[1]
        if action == "Starting:":
            logical_flowchart_starting(model,fitopts,task)
        elif action == "Finished:":
            logical_flowchart_finished(model,fitopts,task)
        else:
            raise IOError("Not valid action: %s" % action)
        mdb.inputs.save_model(model,fitopts)


def logical_flowchart_starting(model,fitopts,task):
    sub = model.name
    iteration = fitopts["iteration"]
    if task == "Tf_loop_iteration":
        print "Checking if Tf_loop_iteration completed"
        simulation.constant_temp.check_completion(model,iteration)
        fitopts["last_completed_task"] = "Starting: Tf_loop_iteration"
        analysis.constant_temp.analyze_temperature_array(model,iteration)
        fitopts["last_completed_task"] = "Starting: Tf_loop_analysis"
        print "Starting Tf_loop_analysis"
    elif task == "Tf_loop_analysis":
        print "Checking if Tf_loop_analysis completed"
        analysis.constant_temp.check_completion(model,iteration)
        fitopts["last_completed_task"] = "Finished: Tf_loop_analysis"
    elif task == "Equil_Tf":
        print "Starting to check if Equil_Tf completed"
        simulation.constant_temp.check_completion(model,iteration,long=True)
        fitopts["last_completed_task"] = "Finished: Equil_Tf"
        analysis.constant_temp.analyze_temperature_array(model,iteration,long=True)
        fitopts["last_completed_task"] = "Starting: Equil_Tf_analysis"
    elif task == "Equil_Tf_analysis":
        print "Starting to check if Equil_Tf_analysis completed"
        analysis.constant_temp.check_completion(model,iteration,long=True)
        fitopts["last_completed_task"] = "Finished: Equil_Tf_analysis"
    else:
        raise ValueError("  Couldn't find next option for task: %s" % task)

def logical_flowchart_finished(model,fitopts,task):
    sub = model.name
    iteration = fitopts["iteration"]
    if task == "Tf_loop_iteration":
        print "Finished: Tf_loop_iteration"
        print "Starting: Tf_loop_analysis"
        analysis.constant_temp.analyze_temperature_array(model,iteration)
        fitopts["last_completed_task"] = "Starting: Tf_loop_analysis"
    elif task == "Tf_loop_analysis":
        print "Finished: Tf_loop_analysis"
        flag = analysis.constant_temp.run_wham_heat_capacity(model,iteration)
        fitopts["last_completed_task"] = "Finished: Tf_wham"
        if flag == 1:
            pass 
        else:
            simulation.constant_temp.folding_temperature_loop(model,iteration)
            print "Starting Tf_loop_iteration"
            fitopts["last_completed_task"] = "Starting: Tf_loop_iteration"
    elif task == "Tf_wham":
        print "Starting equilibrium simulations at Tf"
        simulation.constant_temp.run_equilibrium_simulations(model,iteration)
        fitopts["last_completed_task"] = "Starting: Equil_Tf"
    elif task == "Equil_Tf":
        print "Starting Equil_Tf_analysis"
        analysis.constant_temp.analyze_temperature_array(model,iteration,long=True)
        fitopts["last_completed_task"] = "Starting: Equil_Tf_analysis"
    elif task == "Equil_Tf_analysis":
        # Run heat capacity for equilibrium runs. Cv(T), F(Q)
        analysis.constant_temp.run_wham_heat_capacity(model,iteration,long=True)
        fitopts["last_completed_task"] = "Finished: Equil_Tf_wham"
    elif task == "Equil_Tf_wham":
        print "Starting calculating feature vector and Jacobian"
        parameter_fitting.prepare_newtons_method(model,fitopts)
        fitopts["last_completed_task"] = "Finished: Solving_Newtons_Method"
    elif task == "Calculating_Jacobian":
        print "Solving for solutions with Levenberg-Marquardt method"
        fitopts["last_completed_task"] = "Starting: Solving_Newtons_Method"
        parameter_fitting.solve_newtons_method(model,fitopts)
        fitopts["last_completed_task"] = "Finished: Solving_Newtons_Method"
    elif task == "Solving_Newtons_Method":
        # Write new parameter file
        # Start the next round of simulations with new parameters.
        parameter_fitting.save_new_parameters(model,fitopts)
        fitopts["iteration"] += 1
        fitopts["last_completed_task"] = "Finished: Saving_New_Params"
    elif task == "Saving_New_Params":
        simulation.constant_temp.start_next_Tf_loop_iteration(model,iteration)
        fitopts["last_completed_task"] = "Starting: Tf_loop_iteration"
    else:
        raise ValueError("  Couldn't find next option for task: %s" % task)

#############################################################################
# Add some simulations to existing project 
#############################################################################
def add_temperature_array(args):
    """ Adds manually adds a temperature array."""

    names = args.names
    Models,Fittingopts = mdb.inputs.load_models(names,dry_run=args.dry_run)

    if args.short_temps != None:
        T_min = args.short_temps[0] 
        T_max = args.short_temps[1] 
        deltaT = args.short_temps[2] 
        long = False
    elif args.long_temps != None:
        temps = args.long_temps
        long = True
    else:
        print "ERROR! Must use --short_temps or --long_temps with this option"
        print " Exiting."
        raise SystemExit

    for i in range(len(Models)):
        model = Models[i]
        iteration = Fittingopts[i]["iteration"]
        if not os.path.exists("%s/iteration_%d" % (model.name,iteration)):
            os.mkdir("%s/iteration_%d" % (model.name,iteration))
        if long == False:
            print "Manually adding temperature array Ti=%d Tf=%d dT=%d" % (T_min,T_max,deltaT)
            print "Starting constant_temp_iteration"
            simulation.constant_temp.manually_add_temperature_array(model,iteration,T_min,T_max,deltaT)
            Fittingopts[i]["last_completed_task"] = "Starting: Tf_loop_iteration"
        elif long == True:
            print "Manually adding equilibrium sims ", temps
            simulation.constant_temp.manually_add_equilibrium_runs(model,iteration,temps)
            Fittingopts[i]["last_completed_task"] = "Starting: Equil_Tf"
        mdb.inputs.save_model(model,Fittingopts[i])

#############################################################################
# Extend some simulations in existing project 
#############################################################################
def extend_temperatures(args):
    """ Manually extends."""
    factor = args.factor

    names = args.names
    Models,Fittingopts = mdb.inputs.load_models(names,dry_run=args.dry_run)

    if (args.short_temps != None) and (args.long_temps != None):
        print "ERROR!"
        print "  Specify either --short_temps or --long_temps. Not both!"
        print "  Exiting"
        raise SystemExit
    if (args.short_temps == None) and (args.long_temps == None):
        print "ERROR!"
        print "  Specify either --short_temps or --long_temps."
        print "  Exiting"
        raise SystemExit
    if (args.short_temps != None):
        temps = args.short_temps
        method = "short"
        print "Extending short temps",temps
    if (args.long_temps != None):
        temps = args.long_temps
        method = "long"
        print "Extending longational temps",temps

    for i in range(len(names)):
        model = Models[i]
        iteration = Fittingopts[i]["iteration"]
        simulation.constant_temp.manually_extend_temperatures(model,iteration,method,temps,factor)
        mdb.inputs.save_model(model,Fittingopts[i])

#############################################################################
# Get command line arguments 
#############################################################################
def get_args():
    """Get command line arguments """
    common = argparse.ArgumentParser(description='common arguments for all parsers', add_help=False)
    
    common.add_argument('--dry_run',
        action='store_true',
        help='Dry run. No simulations submitted.')
        
    parser = argparse.ArgumentParser(description='Options for MatysiakClementi2004 recipe.')
    sp = parser.add_subparsers(dest='action')

    # Options for initializing a new project.
    new_parser = sp.add_parser('new', parents=[common])
    new_parser.add_argument('--names',
        type=str,
        required=True,
        nargs='+',
        help='Names corresponding to pdbs to start simulations.')

    new_parser.add_argument('--temparray',
        type=int,
        nargs='+',
        help='Optional initial temp array: T_min T_max deltaT. Default: predicts now')


    # Options for continuing from a previously saved simulation project.
    run_parser = sp.add_parser('continue', parents=[common])
    run_parser.add_argument('--names',
        type=str,
        required=True,
        nargs='+',
        help='Subdirectories to continue')
    
    # Options for manually adding a temperature array.
    add_parser = sp.add_parser('add', parents=[common])
    add_parser.add_argument('--names',
        type=str,
        required=True,
        nargs='+',
        help='Subdirectories to add temp array')

    add_parser.add_argument('--short_temps',
        type=int,
        nargs='+',
        help='T_initial T_final dT for new temps')

    add_parser.add_argument('--long_temps',
        type=float,
        nargs='+',
        help='List of temps for new equilibrium simulations')

    # Options for manually extending some temperatures.
    ext_parser = sp.add_parser('extend', parents=[common])
    ext_parser.add_argument('--names',
        type=str,
        required=True,
        nargs='+',
        help='Subdirectories to add temp array')

    ext_parser.add_argument('--factor',
        type=float,
        required=True,
        help='Factor by which you want to extend simulations. e.g. --factor 2 doubles length')

    ext_parser.add_argument('--short_temps',
        type=float,
        nargs='+',
        help='Temperatures that you want extended')

    ext_parser.add_argument('--long_temps',
        type=float,
        nargs='+',
        help='Temperatures that you want extended')
        

    args = parser.parse_args()

    return args

if __name__ == "__main__":
    args = get_args()

    path = os.getcwd()
    if args.action == 'new':
        new_project(args)
    elif args.action == 'continue':
        continue_project(args)
    elif args.action == 'add':
        add_temperature_array(args)
    elif args.action == 'extend':
        extend_temperatures(args)
    print "Success"
