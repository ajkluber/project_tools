import argparse
import os
import shutil
import subprocess as sb
import time
import numpy as np

import systems
import analysis
import simulation
import models

'''
ModelBuilder Class

Purpose:
    The ModelBuilder class is top-level helper class to track and execute
varied procedures for structure-based modeling. At this level the tasks that
need to be done are abstract (such as 'determining folding temperature' or
'refining parameter'), so ModelBuilder relies on the internal machinery of
the Model and System classes to handle the details.
    The goal is to conceal as much of the details as possible away from the
user, so that the user can focus on top-level questions. For this reason any
function that requires manipulating the data is best moved to a different 
module. 
    
Description:
    This module contains the logical flow for running simulations for 
several coarse-grain models.
    
See development_notes.txt for chronological list of changes and development
plan.

'''

def print_header():

    print "---------------------------- Model Builder v0.1 ------------------------------"
    print " Your using model_builder!  A helper module for prepping CG simulations for"
    print " Gromacs. More coming soon!"
    print " Version 0.1 "
    print " Words to live by:\n"
    #print "             'If you can calculate it, you should calculate it' - PGW \n"
    #print "               'One never notices what has been done; "
    #print "                one can only see what remains to be done' - Marie Curie \n"
    print "  'The best way to have a good idea is to have a lot of ideas' - Linus Pauling\n"
    print "-------------------------------- Good Luck! ----------------------------------"

def get_args():
    ''' Get command line arguments and check that they are consistent/allowed.
        Maybe also do a cursory check of stored data for the 'continue' option
        to ensure integrity of data. 
    '''
    print_header()

    parser = argparse.ArgumentParser(description='Build a model of a system.')
    sp = parser.add_subparsers(dest='action')

    ## Options for initializing a new simulation project.
    new_parser = sp.add_parser('new')
    #new_parser.add_argument('--name', type=str, required=True, help='Name of system.')
    new_parser.add_argument('--type', type=str, required=True, help='Choose model type: HetGo, HomGo, DMC')
    new_parser.add_argument('--beads', type=str, required=True, help='Choose model beads: CA, CACB.')
    new_parser.add_argument('--pdbs', type=str, required=True, nargs='+',help='PDBs to start simulations.')
    new_parser.add_argument('--contact_energies', type=str, help='HetGo Contact energies: MJ, Bach, MC2004, from file.')
    new_parser.add_argument('--temparray', type=int, nargs='+',help='Optional initial temp array: Ti Tf dT. Default: 50 350 50')
    new_parser.add_argument('--solvent', action='store_true', help='Add this option for solvent.')
    new_parser.add_argument('--dryrun', action='store_true', help='Add this option for dry run. No simulations started.')
    new_parser.add_argument('--R_CD', type=float, help='Optional specific ratio of contact to dihedral energy.')
    new_parser.add_argument('--cutoff', type=float, help='Optional cutoff for heavy atom determination of native contacts.')
    new_parser.add_argument('--disulfides', type=int, nargs='+', help='Optional pairs of disulfide linked residues.')

    ## Options for continuing from a previously saved simulation project.
    run_parser = sp.add_parser('continue')
    run_parser.add_argument('--subdirs', type=str, nargs='+', help='Subdirectories to continue',required=True)
    run_parser.add_argument('--dryrun', action='store_true', help='Dry run. No simulations started.')

    ## Options for manually adding a temperature array.
    add_parser = sp.add_parser('add')
    add_parser.add_argument('--subdirs', type=str, nargs='+', help='Subdirectories to add temp array',required=True)
    add_parser.add_argument('--temparray', type=int, nargs='+', help='T_initial T_final dT for new temp array',required=True)
    add_parser.add_argument('--dryrun', action='store_true', help='Add this option for solvent.')

    args = parser.parse_args()

    if args.action == "new":
        ## Check if options for model make sense.
        options = {}
        options["Model_Code"] = args.type
        options["Bead_Model"] = args.beads
        options["Solvent"] = args.solvent
        options["R_CD"] = args.R_CD
        options["Disulfides"] = args.disulfides
        options["Contact_Energies"] = args.contact_energies
        modeloptions = models.check_options(options)

    elif args.action == "continue":
        modeloptions = {}
        pass
        ##
    else:
        modeloptions = {}

    if args.dryrun != False:
        modeloptions["Dry_Run"] = True
    else:
        modeloptions["Dry_Run"] = False

    return args, modeloptions


class ModelBuilder(object):
    def __init__(self,args,modeloptions):
        ''' Model Initialization.''' 
        self.path = os.getcwd()
        
        ## Not used 
        procedure = {'HomGo':["Prepping system","Submitting T_array",
                              "Analyzing T_array"],
                     'HetGo':["Prepping system","Submitting T_array",
                              "Analyzing T_array","Mutations"],
                     'DMC':["Prepping system","Submitting T_array",
                              "Analyzing T_array","Mutations"]}

        go_model_procedure = ["Prepping system","Submitting T_array","Analyzing T_array"]

        if args.action == 'new':
            ## This option creates a new project with the specifications given
            ## on the command line. Then saves all that info in a format that 
            ## can be automatically loaded the next time around.
            self.new_project(args,modeloptions)
        elif args.action == 'continue':
            ## project where you don't remember where it left off. Outputs a nice
            ## summary of the state of the simulation.
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

        System = system.System(args)
        self.load_model_system_info(System)
        #Models = models.load_models(System.subdirs)
        modelinfo = open(args.subdirs[0]+'/model.info','r').readlines()
        modeltype = modelinfo[3].split()[0]
        Model = models.get_model(modeltype)
        self.prepare_system(Model,System)

        Ti = args.temparray[0] 
        Tf = args.temparray[1] 
        dT = args.temparray[2] 

        ## DEBUGGING 
        #print Ti, Tf, dT       
        #print type(Ti), type(Tf), type(dT)
        #raise SystemExit

        for i in range(len(System.subdirs)):
            sub = System.subdirs[i]
            lasttime,action,task = self.check_modelbuilder_log(sub)
            print "Checking progress for directory:  ", sub
            print "Last task was %s %s at %s" % (action,task,lasttime) ## DEBUGGING
            print "Manually adding temperature array Ti=%d Tf=%d dT=%d" % (Ti,Tf,dT)
            print "Starting Tf_loop_iteration..."
            simulation.Tf_loop.manually_add_temperature_array(Model,System,i,self.append_log,Ti,Tf,dT)
            
        self.save_model_system_info(Model,System,subdirs)
        print "Success"

    def check_modelbuilder_log(self,sub):
        ''' Gets last line of sub/modelbuilder.log to determine where the 
            program last left off. Could probably be more robust.'''
        modelbuilder = open(sub+'/modelbuilder.log','r').readlines()
        lasttime, action,task = modelbuilder[-1].split()
        return lasttime,action,task

    def continue_project(self,args):
        ''' Checks where something left off and continues it.'''
        args.pdbs = [ name+'.pdb' for name in  args.subdirs ]
        subdirs = args.subdirs

        ## Read in options for each directory.
        Models = models.load_models(subdirs)
        Systems = systems.load_systems(subdirs)
        Model = Models[0]       ## Temporary for backwards compatibility.
        System = Systems[0]     ## Temporary for backwards compatibility.

        self.prepare_systems(Models,Systems)
        self.save_model_system_info(Model,System,subdirs)

        raise SystemExit

        #####   OLD CODE vvvvv
        #System = systems.system.System(args)
        #self.load_model_system_info(System)
        #modelinfo = open(args.subdirs[0]+'/model.info','r').readlines()
        #modeltype = modelinfo[3].split()[0]
        #Model = models.get_model(modeltype)
        #if len(System.R_CD) != 0:
        #    if System.R_CD[0] != None:
        #        self.prepare_system(Model,System,R_CD=System.R_CD[0])
        #    else:
        #        self.prepare_system(Model,System)
        #else:
        #    self.prepare_system(Model,System)
        #####   OLD CODE ^^^^^

        if args.dryrun == True:
            print "Dry run complete. Exiting."
        else:
            for i in range(len(System.subdirs)):
                ## Desired format:
                ## Model = Models[i]        # Access Model from list
                ## System = Systems[i]      # Access System from list
                sub = System.subdirs[i]
                lasttime,action,task = self.check_modelbuilder_log(sub)
                print "Checking progress for directory:  ", sub
                print "Last task was %s %s at %s" % (action,task,lasttime) ## DEBUGGING
                if action == "Starting:":
                    logical_flowchart_starting(System,Model,i,sub,task)
                elif action == "Finished:":
                    logical_flowchart_finished(System,Model,i,sub,task)
                elif action == "Error":
                    pass

        self.save_model_system_info(Model,System,subdirs)
        print "Success"

    def logical_flowchart_finished(System,Model,i,sub,task):
        if task == "Tf_loop_iteration":
            print "Finished Tf_loop_iteration..."
            print "Starting Tf_loop_analysis..."
            analysis.Tf_loop.analyze_temperature_array(System,i,self.append_log)
        elif task == "Tf_loop_analysis":
            print "Finished Tf_loop_analysis..."
            flag = analysis.Tf_loop.check_if_wham_is_next(System,i,self.append_log)
            if flag == 1:
                pass 
            else:
                print "Starting Tf_loop_iteration..."
                simulation.Tf_loop.folding_temperature_loop(Model,System,i,self.append_log)
        elif task == "wham_Cv":
            print "Finished wham_Cv..."
            print "Stating wham_FreeEnergy..."
            analysis.Tf_loop.continue_wham(System,i,self.append_log)
        elif task == "Equil_Tf":
            print "Starting Equil_Tf_analysis..."
            analysis.Tf_loop.analyze_temperature_array(System,i,self.append_log,equil=True)

    def logical_flowchart_starting(System,Model,i,sub,task):
        if task == "Tf_loop_iteration":
            print "Starting to check if Tf_loop_iteration completed..."
            simulation.Tf_loop.check_completion(System,i,self.append_log)
            lasttime2,action2,task2 = self.check_modelbuilder_log(sub)
            if action2 == "Finished:":
                print "Finished Tf_loop_iteration..."
                print "Starting Tf_loop_analysis..."
                analysis.Tf_loop.analyze_temperature_array(System,i,self.append_log)
        elif task == "Tf_loop_analysis":
            print "Starting to check if Tf_loop_analysis completed..."
            analysis.Tf_loop.check_completion(System,i,self.append_log)
        elif task == "wham_Cv":
            print "Starting to check if wham_Cv completed..."
            analysis.Tf_loop.continue_wham(System,i,self.append_log)
        elif task == "wham_FreeEnergy":
            ## Start equilibrium runs.
            print "Starting Equil_Tf..."
            simulation.Tf_loop.run_equilibrium_simulations(Model,System,i,self.append_log)
        elif task == "Equil_Tf":
            print "Starting to check if Equil_Tf completed..."
            simulation.Tf_loop.check_completion(System,i,self.append_log,equil=True)
            lasttime2,action2,task2 = self.check_modelbuilder_log(sub)
            if action2 == "Finished:":
                print "Finished Equil_Tf_iteration..."
                print "Starting Equil_Tf_analysis..."
                analysis.Tf_loop.analyze_temperature_array(System,i,self.append_log,equil=True)
        elif task == "Equil_Tf_analysis":
            print "Starting to check if Equil_Tf_analysis completed..."
            analysis.Tf_loop.check_completion(System,i,self.append_log,equil=True)

    def new_project(self,args,modeloptions):
        ''' Starting a new simulation project.'''
        subdirs = [ x[:-4] for x in args.pdbs ]

        for sub in subdirs:
            if os.path.exists(sub) == False:
                os.mkdir(sub)
            else:
                print "Subdirectory: ", sub, " already exists! just fyi"

        print "Starting a new simulation project..."
        ## Transitioning to using list of Model objects instead of one singluar
        ## Model object. 3-10-14 AK
        Models = models.new_models(subdirs,modeloptions)
        Model = Models[0]
        ## Transistioning to using a list of System objects
        Systems = systems.new_systems(subdirs)
        System = Systems[0]

        self.prepare_systems(Models,Systems)

        self.save_model_system_info(Model,System,subdirs)
        #print dir(System) ## DEBUGGING
        #print System.topology_files.keys() ## DEBUGGING
        raise SystemExit

        System = systems.system.System(args)

        self.create_subdirs(System)
        if args.cutoff != None:
            print "Using cutoff", args.cutoff
            cutoff = args.cutoff
        else:
            cutoff = 5.5
        if args.R_CD != None:
            print "Using R_CD = ",args.R_CD
            self.prepare_system(Model,System,R_CD=args.R_CD,cutoff=cutoff)
        else:
            self.prepare_system(Model,System,cutoff=cutoff)
        self.save_model_system_info(Model,System,subdirs)
        self.load_model_system_info(System)
        if args.temparray != None:
            System.initial_T_array = args.temparray

        ## The first step depends on the type of model.
        if args.type in ["HomGo","HetGo"]:
            for k in range(len(System.subdirs)):
                ## To Do: Prepare each Model System pair. 
                print "Starting the first Tf_loop_iteration..."
                if args.dryrun == True:
                    print "Dry run complete. Exiting."
                else:
                    simulation.Tf_loop.folding_temperature_loop(Model,System,k,self.append_log)
        elif args.type == "DMC":
            pass

        #self.append_log("Analysis T_array")
        #analysis.analyze_temperature_array(System)
        #self.append_log("Finished: T_array")
        print "Success"

    def load_model_system_info(self,System):
        ''' Save the model and system info strings.'''
        print "Loading system.info and model.info ..."
        for i in range(len(System.subdirs)):
            #print System.subdirs[i]
            System.load_info_file(i)
            #modelname = ''
        #print System.__repr__(i)

    def save_model_system_info(self,Model,System,subdirs):
        ''' Save the model and system info strings.'''
        print "Saving system.info progress..."
        for i in range(len(subdirs)):
            open(subdirs[i] + "/system.info","w").write(System.__repr__())
            open(subdirs[i] + "/model.info","w").write(Model.__repr__())

    def prepare_systems(self,Models,Systems):
        ''' New style of preparing files: on subdirectory basis.'''
        for i in range(len(Models)):
            if os.path.exists(Systems[i].path + "/" + Systems[i].subdir + "/" + Systems[i].subdir + ".pdb") == False:
                shutil.copy(Systems[i].subdir + ".pdb", Systems[i].subdir)
            #print os.path.exists(Systems[i].path+"/"+Systems[i].subdir+"/Qref_shadow")     ## DEBUGGING

            if os.path.exists(Systems[i].path+"/"+Systems[i].subdir+"/Qref_shadow") == False:
                os.mkdir(Systems[i].path+"/"+Systems[i].subdir+"/Qref_shadow")
            Models[i].new_prepare_system(Systems[i])
            print "Done preparing systems."

    def prepare_system(self,Model,System,R_CD=None,cutoff=5.5):
        ''' Extract all the topology files from Model. 
            SOON TO BE MOVED INTO THE MODEL CLASS.
        '''
        print "Preparing files..."
        prots_Qref = System.shadow_contacts()
        System.write_Native_pdb_CA()
        if R_CD != None:
            for i in range(len(System.subdirs)):
                Nc = float(sum(sum(prots_Qref[i])))
                Nd = float(len(prots_Qref[i])-4)
                System.nonbond_params.append((R_CD*Nd/Nc)*Model.backbone_param_vals["Kd"])
                System.R_CD.append(R_CD)
        else:
            for i in range(len(System.subdirs)):
                N = len(prots_Qref[i])
                Nc = float(sum(sum(prots_Qref[i])))
                Nd = float(len(prots_Qref[i])-4)
                print "Num contacts per residue: ",Nc/N
                System.nonbond_params.append(Model.nonbond_param)
                System.R_CD.append(None)
        prots_indices, prots_residues, prots_coords = System.get_atom_indices(Model.beadmodel)
        prots_ndxs = Model.get_index_string(prots_indices)
        topology_files = Model.get_itp_strings(prots_indices, prots_residues, prots_coords,prots_ndxs,prots_Qref,R_CD=R_CD)
        System.topology_files = topology_files

    
def main():
    args, modeloptions = get_args()
    ModelBuilder(args,modeloptions)

if __name__ == '__main__':
    main()
