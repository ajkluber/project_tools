""" Basic project manager class automates projects at top level.

Description:

ProjectManager
    The ProjectManager class is top-level helper class to track and execute
varied procedures for structure-based modeling. At this level the tasks that
need to be done are abstract (such as 'determining folding temperature' or
'refining parameter'). 
    ProjectManager relies on model_builder to generate the proper Gromacs input
files for the desired coarse-grain model, but doesn't need to know of the
details.
    The goal is to conceal as much of the details as possible away from the
user, so that the user can focus on top-level questions. For this reason any
function that requires manipulating the data is best moved to a different 
module. 

    
See Also: 

    development_notes.txt for chronological list of changes and development
plan.

"""

import argparse
import os
import shutil
import subprocess as sb
import time
import numpy as np

import simulation
import analysis
import mutations

from model_builder import models
from model_builder import systems

global GAS_CONSTANT_KJ_MOL
GAS_CONSTANT_KJ_MOL = 0.0083144621

def print_header():

    print "---------------------------- Model Builder v0.1 ------------------------------"
    print " Your using model_builder!  A helper module for prepping CG simulations for"
    print " Gromacs. More coming soon!"
    print " Version 0.1 "
    print " Words to live by:\n"
    print "             'If you can calculate it, you should calculate it' - PGW \n"
    #print "               'One never notices what has been done; "
    #print "                one can only see what remains to be done' - Marie Curie \n"
    #print "  'The best way to have a good idea is to have a lot of ideas' - Linus Pauling\n"
    #print "'A ship in port is safe, but that's not what ships are built for' - Grace Hopper'\n"
    #print "                 'Science and everyday life cannot and"
    #print "                  should not be separated' - Rosalind Franklin\n"
    print "-------------------------------- Good Luck! ----------------------------------"

def get_args():
    ''' Get command line arguments and check that they are consistent/allowed.
    '''
    print_header()

    parser = argparse.ArgumentParser(description='Build a model of a system.')
    sp = parser.add_subparsers(dest='action')

    ## Options for initializing a new simulation project.
    new_parser = sp.add_parser('new')
    new_parser.add_argument('--model', type=str, required=True, help='Choose model type: HetGo, HomGo, DMC')
    new_parser.add_argument('--beads', type=str, required=True, help='Choose model beads: CA, CACB.')
    new_parser.add_argument('--pdbs', type=str, required=True, nargs='+',help='PDBs to start simulations.')
    new_parser.add_argument('--contact_energies', type=str, help='HetGo Contact energies: MJ, Bach, MC2004, from file.')
    new_parser.add_argument('--temparray', type=int, nargs='+',help='Optional initial temp array: Ti Tf dT. Default: 50 350 50')
    new_parser.add_argument('--solvent', action='store_true', help='Add this option for solvent.')
    new_parser.add_argument('--dryrun', action='store_true', help='Add this option for dry run. No simulations started.')
    new_parser.add_argument('--R_CD', type=float, help='Optional specific ratio of contact to dihedral energy.')
    new_parser.add_argument('--epsilon_bar', type=float, help='Optional, average strength of contacts. epsilon bar.')
    new_parser.add_argument('--cutoff', type=float, help='Optional cutoff for heavy atom determination of native contacts.')
    new_parser.add_argument('--disulfides', type=int, nargs='+', help='Optional pairs of disulfide linked residues.')
    new_parser.add_argument('--email', type=str, help='Optional email address for PBS to send sim details.')

    ## Options for continuing from a previously saved simulation project.
    run_parser = sp.add_parser('continue')
    run_parser.add_argument('--subdirs', type=str, nargs='+', help='Subdirectories to continue',required=True)
    run_parser.add_argument('--dryrun', action='store_true', help='Dry run. No simulations started.')
    run_parser.add_argument('--email', type=str, help='Optional email address for PBS to send sim details.')

    ## Options for manually adding a temperature array.
    add_parser = sp.add_parser('add')
    add_parser.add_argument('--subdirs', type=str, nargs='+', help='Subdirectories to add temp array',required=True)
    add_parser.add_argument('--temparray', type=int, nargs='+', help='T_initial T_final dT for new temp array',required=True)
    add_parser.add_argument('--mutarray', type=int, nargs='+', help='T_initial T_final dT for new mutational sims array')
    add_parser.add_argument('--dryrun', action='store_true', help='Dry run. No simulations started.')
    add_parser.add_argument('--email', type=str, help='Optional email address for PBS to send sim details.')

    args = parser.parse_args()

    if args.dryrun != False:
        options = {"Dry_Run":True}
    else:
        options = {"Dry_Run":False}

    if args.action == "new":
        ## Check if options for model make sense.
        options["Model_Code"] = args.model
        options["Bead_Model"] = args.beads
        options["Solvent"] = args.solvent
        options["R_CD"] = args.R_CD
        options["Epsilon_Bar"] = args.epsilon_bar
        options["Disulfides"] = args.disulfides
        options["Contact_Energies"] = args.contact_energies
        modeloptions = models.check_options(options)
    elif args.action == "continue":
        modeloptions = options
    else:
        modeloptions = options
    return args, modeloptions


class ProjectManager(object):
    """ A shell class to 


    """

    def __init__(self,args,modeloptions):
        self.path = os.getcwd()
        
        if args.action == 'new':
            self.new_project(args,modeloptions)
        elif args.action == 'continue':
            self.continue_project(args)
        elif args.action == 'add':
            self.add_temperature_array(args)

    def append_log(self,sub,string):
        open(self.path+'/'+sub+'/modelbuilder.log','a').write(self.append_time_label()+' '+string+'\n')

    def append_time_label(self):
        now = time.localtime()
        now_string = "%s:%s:%s:%s:%s" % (now.tm_year,now.tm_mon,now.tm_mday,now.tm_hour,now.tm_min)
        return now_string

    def create_subdirs(self,System):
        ''' Create the subdirectories for the system. Skip if they already exist.'''
        for i in range(len(System.subdirs)): 
            sub = System.subdirs[i]
            if os.path.exists(sub) != True:
                os.mkdir(sub)
                os.mkdir(sub+"/"+System.Tf_active_directory[i])
                #open(sub+"/Tf_active_directory.txt","w").write(System.Tf_active_directory[i])
                self.append_log(sub,"Creating new subdirectory: %s" % sub)

    def add_temperature_array(self,args):
        ''' Adds manually adds a temperature array.'''
        args.pdbs = [ name+'.pdb' for name in  args.subdirs ]

        subdirs = args.subdirs
        Models = models.load_models(subdirs,dryrun=args.dryrun)
        Systems = systems.load_systems(subdirs)
        self.prepare_systems(Models,Systems)

        Ti = args.temparray[0] 
        Tf = args.temparray[1] 
        dT = args.temparray[2] 

        for i in range(len(subdirs)):
            sub = subdirs[i]
            Model = Models[i]
            System = Systems[i]
            lasttime,action,task = self.check_modelbuilder_log(sub)
            print "Checking progress for directory:  ", sub
            print "Last task was %s %s at %s" % (action,task,lasttime) ## DEBUGGING
            print "Manually adding temperature array Ti=%d Tf=%d dT=%d" % (Ti,Tf,dT)
            print "Starting Tf_loop_iteration..."
            simulation.Tf_loop.manually_add_temperature_array(Model,System,self.append_log,Ti,Tf,dT)
            
        self.save_model_system_info(Models,Systems,subdirs)
        print "Success"

    def check_modelbuilder_log(self,sub):
        ''' Gets last line of sub/modelbuilder.log to determine where the 
            program last left off. Could probably be more robust.'''
        modelbuilder = open(sub+'/modelbuilder.log','r').readlines()
        lasttime, action,task = modelbuilder[-1].split()
        return lasttime,action,task

    def continue_project(self,args):
        ''' Checks where something left off and continues it.'''
        subdirs = args.subdirs
        Models = models.load_models(subdirs,dryrun=args.dryrun)
        Systems = systems.load_systems(subdirs)
        self.prepare_systems(Models,Systems)

        if args.dryrun == True:
            print "  Dry run complete. Saving and Exiting."
        else:
            for i in range(len(subdirs)):
                Model = Models[i]        
                System = Systems[i]      
                subdir = subdirs[i]
                
                lasttime,action,task = self.check_modelbuilder_log(subdir)
                print "Checking progress for directory:  ", subdir
                print "Last task was %s %s at %s" % (action,task,lasttime)
                if action == "Starting:":
                    self.logical_flowchart_starting(System,Model,subdir,task)
                elif action == "Finished:":
                    self.logical_flowchart_finished(System,Model,subdir,task)
                elif action == "Error":
                    pass

        self.save_model_system_info(Models,Systems,subdirs)
        print "Success"

    def save_model_system_info(self,Models,Systems,subdirs):
        ''' Save the model and system info strings.'''
        print "Saving system.info progress..."
        for i in range(len(subdirs)):
            open(subdirs[i]+"/system.info","w").write(Systems[i].__repr__())
            open(subdirs[i]+"/model.info","w").write(Models[i].__repr__())

    def prepare_systems(self,Models,Systems):
        ''' New style of preparing files: on subdirectory basis.'''
        for i in range(len(Models)):
            if not os.path.exists(Systems[i].path+"/"+Systems[i].subdir+"/"+Systems[i].subdir+".pdb"):
                shutil.copy(Systems[i].subdir+".pdb",Systems[i].subdir)
            if not os.path.exists(Systems[i].path+"/"+Systems[i].subdir+"/Qref_shadow"):
                os.mkdir(Systems[i].path+"/"+Systems[i].subdir+"/Qref_shadow")
            Models[i].prepare_system(Systems[i])
            print "Done preparing systems."

def main():
    args, modeloptions = get_args()
    ProjectManager(args,modeloptions)

if __name__ == '__main__':
    main()
