import argparse
import os
import subprocess as sb
import time
import numpy as np
import cPickle 

import system
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
    
To Do:
- Decide on a logging format that can also be read to determine the
  the last successfully finshed step.
  a. Implement command line option to check last step finished.
  b. Decide on format to store the meta-procedure for each type of 
    simulation. Probably store this in Model class.
  c. Implement logging of every important task.  

- Ideas for new command line functionality. 

- Code procedure for running many temperatures.

- Think about analysis tools. Probably should be their own module.

February 15 2014
Alexander Kluber



'''


class ModelBuilder(object):
    def __init__(self,args):
        ''' Model Initialization.''' 
        self.path = os.getcwd()
        
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
            self.new_project(args)
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
        ''' Checks where something left off and continues it.'''
        args.pdbs = [ name+'.pdb' for name in  args.subdirs ]

        System = system.System(args)
        self.load_model_system_info(System)
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
            
        self.save_model_system_info(Model,System)
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

        System = system.System(args)
        self.load_model_system_info(System)
        modelinfo = open(args.subdirs[0]+'/model.info','r').readlines()
        modeltype = modelinfo[3].split()[0]
        Model = models.get_model(modeltype)
        if len(System.R_CD) != 0:
            if System.R_CD[0] != None:
                self.prepare_system(Model,System,R_CD=System.R_CD[0])
            else:
                self.prepare_system(Model,System)
        else:
            self.prepare_system(Model,System)

        #print System.__repr__(0) ## DEBUGGING
        #print Model.__repr__()   ## DEBUGGING
        #raise SystemExit

        if args.dryrun == True:
            print "Dry run complete. Exiting."
        else:
            for i in range(len(System.subdirs)):
                sub = System.subdirs[i]
                lasttime,action,task = self.check_modelbuilder_log(sub)
                print "Checking progress for directory:  ", sub
                print "Last task was %s %s at %s" % (action,task,lasttime) ## DEBUGGING
                if action == "Starting:":
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
                elif action == "Finished:":
                    if task == "Tf_loop_iteration":
                        print "Finished Tf_loop_iteration..."
                        print "Starting Tf_loop_analysis..."
                        analysis.Tf_loop.analyze_temperature_array(System,i,self.append_log)
                    elif task == "Tf_loop_analysis":
                        print "Finished Tf_loop_analysis..."
                        flag = analysis.Tf_loop.check_if_wham_is_next(System,i,self.append_log)
                        if flag == 1:
                            continue
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
                elif action == "Error":
                    pass

        self.save_model_system_info(Model,System)
        print "Success"
        

    def new_project(self,args):
        ''' Starting a new simulation project.'''
        #self.append_log("Project %s started" % args.name)
        print "Starting a new project..."
        Model = models.get_model(args.type)
        System = system.System(args)
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
        self.save_model_system_info(Model,System)
        self.load_model_system_info(System)
        if args.temparray != None:
            System.initial_T_array = args.temparray

        ## The first step depends on the type of model.
        if args.type in ["HomGo","HetGo"]:
            for k in range(len(System.subdirs)):
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

    def save_model_system_info(self,Model,System):
        ''' Save the model and system info strings.'''
        print "Saving system.info progress..."
        for i in range(len(System.subdirs)):
            System.write_info_file(i)
            Model.nonbond_param = System.nonbond_params[i]
            Model.write_info_file(System.subdirs[i])

    def prepare_system(self,Model,System,R_CD=None,cutoff=5.5):
        ''' Extract all the topology files from Model.'''
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
    parser = argparse.ArgumentParser(description='Build a model of a system.')
    sp = parser.add_subparsers(dest='action')

    ## Initializing a new simulation project.
    new_parser = sp.add_parser('new')
    #new_parser.add_argument('--name', type=str, required=True, help='Name of system.')
    new_parser.add_argument('--type', type=str, required=True, help='Choose model type: HetGo, HomGo, DMC')
    new_parser.add_argument('--beads', type=str, required=True, help='Choose model beads: CA, CACB.')
    new_parser.add_argument('--pdbs', type=str, required=True, nargs='+',help='PDBs to start simulations.')
    new_parser.add_argument('--temparray', type=int, nargs='+',help='Optional initial temp array: Ti Tf dT. Default: 50 350 50')
    new_parser.add_argument('--solvent', action='store_true', help='Add this option for solvent.')
    new_parser.add_argument('--dryrun', action='store_true', help='Add this option for dry run. No simulations started.')
    new_parser.add_argument('--R_CD', type=float, help='Optional specific ratio of contact to dihedral energy.')
    new_parser.add_argument('--cutoff', type=float, help='Optional cutoff for heavy atom determination of native contacts.')

    ## Checking on a simulation project.
    run_parser = sp.add_parser('continue')
    run_parser.add_argument('--subdirs', type=str, nargs='+', help='Subdirectories to continue',required=True)
    run_parser.add_argument('--dryrun', action='store_true', help='Dry run. No simulations started.')

    ## Manually adding a temperature array.
    add_parser = sp.add_parser('add')
    add_parser.add_argument('--subdirs', type=str, nargs='+', help='Subdirectories to add temp array',required=True)
    add_parser.add_argument('--temparray', type=int, nargs='+', help='T_initial T_final dT for new temp array',required=True)
    add_parser.add_argument('--dryrun', action='store_true', help='Add this option for solvent.')

    args = parser.parse_args()
    
    
    #print args
    #raise SystemExit
    ModelBuilder(args)

if __name__ == '__main__':
    main()
