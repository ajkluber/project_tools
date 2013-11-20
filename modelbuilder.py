import argparse
import os
import subprocess as sb
import time
import numpy as np
import cPickle 

import system
import analysis
import sim
from HomogeneousGoModel import HomogeneousGoModel
from HeterogeneousGoModel import HeterogeneousGoModel
from DMCModel import DMCModel

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
        newmodel = {'HomGo':HomogeneousGoModel, 'HetGo':HeterogeneousGoModel, 
                      'DMC':DMCModel}[args.type]
        
        procedure = {'HomGo':["Prepping system","Submitting T_array",
                              "Analyzing T_array"],
                     'HetGo':["Prepping system","Submitting T_array",
                              "Analyzing T_array","Mutations"],
                     'DMC':["Prepping system","Submitting T_array",
                              "Analyzing T_array","Mutations"]}

        go_model_rocedure = ["Prepping system","Submitting T_array","Analyzing T_array"]

        if args.action == 'new':
            self.append_log("Project %s started" % args.name)
            Model = newmodel(self.path)
            System = system.System(args)
            self.create_subdirs(System.subdirs)
            self.append_log("Prepping system")
            self.prep_system(System,Model)
            self.append_log("Finished: Prepping system")
            #print self.Model  ## DEBUGGING
            #print self.System  ## DEBUGGING

            self.append_log("Submitting T_array")
            sim.run_temperature_array(Model,System)
            self.append_log("Finished: Submitting T_array")
            #self.append_log("Analys T_array")
            #analysis.analyze_temperature_array(System)
            #self.append_log("Finished: T_array")
            #self.save_project(Model,System)
            print "Success"
            raise SystemExit
            
            #self.write_info_file()
        elif args.action == 'check':
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

    def append_log(self,string):
        logfile = open(self.path + '/modelbuilder.log','a').write(self.append_time_label()+' '+string+'\n')

    def append_time_label(self):
        now = time.localtime()
        now_string = "%s:%s:%s:%s" % (now.tm_mon,now.tm_mday,now.tm_hour,now.tm_min)
        return now_string

    def create_subdirs(self,subdirs):
        ''' Create the subdirectories for the system.'''
        pdbs = ''
        for i in range(len(subdirs)): 
            sub = subdirs[i]
            try:
                os.mkdir(sub)
            except:
                pass
            pdbs += sub + "/ "
        self.append_log("Creating new subdirectories: %s" % pdbs)

    def write_info_file(self):
        ''' Writes model.info file for new simulation project. 
            DOESN"T WORK. NEEDS TO BE UPDATED.'''

        template_info = open('/home/ajk8/projects/dmc_model/gmx/%s.params' % self.modeltype, 'r').read()
        pdbs = ''
        for pdb in self.pdbs: pdbs += pdb + " "
        info = open('model.info','w')
        info.write("Simulation project initialized on: %s\n" % time.asctime())
        info.write("[ System ]\n")
        info.write("SystemName:    %s\n"  % self.systemname )
        info.write("PDBs:    %s\n"  % pdbs )
        info.write("\n")
        info.write("[ Model ]\n")
        info.write(template_info)
        info.write("Solvent: %s" % self.solvent)
        info.close()

    def prep_system(self,System,Model):
        ''' Extract all the topology files from Model.'''
        System.clean_pdbs()
        System.write_Native_pdb_CA()
        prots_indices, prots_residues, prots_coords = System.get_atom_indices(Model.beadmodel)
        prots_ndxs = Model.get_index_string(System.subdirs,prots_indices)
        topology_files = Model.get_itp_strings(prots_indices, prots_residues, prots_coords,prots_ndxs)
        System.topology_files = topology_files

    def load_project(self):
        ''' Load Model and System objects for later use. Doesn't work.'''
        mdlpkl = open("model.pkl","rb")
        Model = cPickle.load(mdlpkl)
        mdlpkl.close()
        syspkl = open("system.pkl","rb")
        System = cPickle.load(syspkl)
        syspkl.close()
        return Model, System
        
    def save_project(self,Model,System):
        ''' Save the Model and System objects for later use. Currently doesn't work.'''
        print Model.__dict__
        raise SystemExit
        syspkl = open(System.path+"/system.pkl","wb")
        cPickle.dump(System,syspkl,protocol=cPickle.HIGHEST_PROTOCOL)
        mdlpkl = open(System.path+"/model.pkl","wb")
        cPickle.dump(Model,mdlpkl,protocol=cPickle.HIGHEST_PROTOCOL)
        mdlpkl.close()
        syspkl.close()
            
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
    #run_parser = sp.add_parser('check')
    #run_parser.add_argument('--output', type=str, help='Method of calculation',required=True)
    #run_parser.add_argument('--Ti', type=int, help='Initial Temperature')
    args = parser.parse_args()
    
    
    #print args
    #raise SystemExit
    ModelBuilder(args)

if __name__ == '__main__':
    main()
