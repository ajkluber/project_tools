import argparse
import os
import subprocess as sb
import time
import numpy as np
import cPickle 

import system
import analysis.analysis as analysis
import simulation
from models.HomogeneousGoModel import HomogeneousGoModel
from models.HeterogeneousGoModel import HeterogeneousGoModel
from models.DMCModel import DMCModel

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
            self.new_project(args)
        elif args.action == 'check':
            self.check_project(args)
            ## I would like this action to check the state of the current 
            ## simulation. Each simulation style has a defined procedure (more 
            ## or less) so this action will check how many steps the procedure has 
            ## successfully executed. This will be useful when returning to a 
            ## project where you don't remember where it left off. Outputs a nice
            ## summary of the state of the simulation.
            pass
        elif args.action == 'clean':
            ## Not sure what to put 
            pass
        elif args.action == 'continue':
            
            ## I would like this action to continue with the next step in the
            ## procedure, picking up right after the last succesfful step.
            #self.load_project()
            #System = system.System("system.info")
            #Model = model.Model("model.info")
            pass

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
                open(sub+"/Tf_active_directory.txt","w").write(System.Tf_active_directory[i])
                self.append_log(sub,"Creating new subdirectory: %s" % sub)

    def new_project(self,args):
        ''' Starting a new simulation project.'''
        newmodel = {'HomGo':HomogeneousGoModel, 'HetGo':HeterogeneousGoModel, 
                      'DMC':DMCModel}[args.type]
        #self.append_log("Project %s started" % args.name)
        Model = newmodel(self.path)
        System = system.System(args)
        self.create_subdirs(System)
        self.prepare_system(Model,System)
        self.save_model_system_info(Model,System)
        self.load_model_system_info(Model,System)

        ## The first step depends on the type of model.
        if args.type in ["HomGo","HetGo"]:
            simulation.Tf_loop.folding_temperature_loop(Model,System,self.append_log)
        elif args.type == "DMC":
            pass

        #self.append_log("Analysis T_array")
        #analysis.analyze_temperature_array(System)
        #self.append_log("Finished: T_array")
        print "Success"

    def load_model_system_info(self,Model,System):
        ''' Save the model and system info strings.'''
        for i in range(len(System.subdirs)):
            System.load_info_file(i)
            modelname = 

    def save_model_system_info(self,Model,System):
        ''' Save the model and system info strings.'''
        for i in range(len(System.subdirs)):
            System.write_info_file(i)
            Model.write_info_file(System.subdirs[i])

    def prepare_system(self,Model,System):
        ''' Extract all the topology files from Model.'''
        System.clean_pdbs()
        System.write_Native_pdb_CA()
        prots_indices, prots_residues, prots_coords = System.get_atom_indices(Model.beadmodel)
        prots_ndxs = Model.get_index_string(prots_indices)
        topology_files = Model.get_itp_strings(prots_indices, prots_residues, prots_coords,prots_ndxs)
        System.topology_files = topology_files

def main():
    parser = argparse.ArgumentParser(description='Build a model of a system.')
    sp = parser.add_subparsers(dest='action')

    ## Initializing a new simulation project.
    new_parser = sp.add_parser('new')
    new_parser.add_argument('--name', type=str, required=True, help='Name of system.')
    new_parser.add_argument('--type', type=str, required=True, help='Choose model type: HetGo, HomGo, DMC')
    new_parser.add_argument('--beads', type=str, required=True, help='Choose model beads: CA, CACB.')
    new_parser.add_argument('--pdbs', type=str, required=True, nargs='+',help='PDBs to start simulations.')
    new_parser.add_argument('--solvent', action='store_true', help='Add this option for solvent.')

    ## Checking on a simulation project.
    run_parser = sp.add_parser('check')
    run_parser.add_argument('--subdir', nargs='+', help='Subdirectories to check',required=True)
    args = parser.parse_args()
    
    
    #print args
    #raise SystemExit
    ModelBuilder(args)

if __name__ == '__main__':
    main()
