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
        
        procedure = {'HomGo':["Prepping system","Submitting T_array",
                              "Analyzing T_array"],
                     'HetGo':["Prepping system","Submitting T_array",
                              "Analyzing T_array","Mutations"],
                     'DMC':["Prepping system","Submitting T_array",
                              "Analyzing T_array","Mutations"]}

        go_model_rocedure = ["Prepping system","Submitting T_array","Analyzing T_array"]

        if args.action == 'new':
            self.new_project(args)
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

    def append_log(self,sub,string):
        logfile = open(self.path+'/'+sub+'/modelbuilder.log','a').write(self.append_time_label()+' '+string+'\n')

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

    def folding_temperature_loop(self,Model,System):
        ''' The "folding temperature loop" is one of the several large-scale 
            logical structures in modelbuilder. It is entered anytime we want
            to determine the folding temperature. This could be when we have
            started a new project, refined the paramters, or returned to a 
            project in progress. The folding temperature loop successively 
            narrows in on the folding temperature.'''

        for k in range(len(System.subdirs)):
            sub = System.path +"/"+ System.subdirs[k] +"/"+ System.Tf_active_directory[k]
            print sub  ## DEBUGGING
            if (not os.path.exists(sub)):
                os.mkdir(sub)
            ## Check to see if the folding temperature has been found. If yes, then continue.
            if (not os.path.exists(sub+"/Tf.txt")):
                ## Check to see if a previous temperature range was used.
                os.chdir(sub)
                self.folding_temperature_loop_extension(Model,System,k)
            else:
                ## Folding temperature has been found. Continuing.
                continue

    def folding_temperature_loop_extension(self,Model,System,k):
        if (not os.path.exists("Ti_Tf_dT.txt")):
            ## For initial exploration use very broad temperature increments.
            self.append_log(System.subdirs[k],"Starting: Submitting T_array iteration %d ; refinement %d" % \
                            (System.Tf_iteration[k],System.Tf_refinements[k][System.Tf_iteration[k]]))
            self.append_log(System.subdirs[k],"  Ti = %d , Tf = %d , dT = %d" % (100, 300, 50))
            sim.run_temperature_array(Model,System,k,Ti=100,Tf=300,dT=50)
            self.append_log(System.subdirs[k],"Finished: Submitting T_array iteration %d ; refinement %d" % \
                            (System.Tf_iteration[k],System.Tf_refinements[k][System.Tf_iteration[k]]))
            #System.Tf_refinements[k][System.Tf_iteration[k]] += 1
        else:
            ## Use previous range to determine new range. 
            ## DOESN'T WORK YET. USE ANALYSIS TO ESTIMATE THE BRACKETING TEMPS.
            Ti,Tf,dT = open("Ti_Tf_dT.txt","r").read().split()
            Ti = int(Ti); Tf = int(Tf); dT = int(dT)
            lowerT, upperT = open("T_brackets.txt","r").read().split()
            lowerT = int(lowerT); upperT = int(upperT)
            newdT = float(dT)/5.
            ## If new dT is less than 1 then don't do it.
            if newdT <= 3.:
                newdT = 1
                midT = int(0.5*(float(lowerT)+upperT))
                newTi = midT - 5
                newTf = midT + 5
            else:
                newTi = lowerT + newdT
                newTf = upperT - newdT
            self.append_log(System.subdirs[k],"Starting: Submitting T_array iteration %d ; refinement %d" % \
                            (System.Tf_iteration[k],System.Tf_refinements[k][System.Tf_iteration[k]]))
            self.append_log(System.subdirs[k],"  Ti = %d , Tf = %d , dT = %d" % (newTi, newTf, dT))
            sim.run_temperature_array(Model,System,Ti=newTi,Tf=newTf,dT=newdT)
            self.append_log(System.subdirs[k],"Finished: Submitting T_array iteration %d ; refinement %d" % \
                            (System.Tf_iteration[k],System.Tf_refinements[k][System.Tf_iteration[k]]))
            System.Tf_refinements[k][System.Tf_iteration[k]] += 1
            open("Ti_Tf_dT.txt","w").write("%d %d %d" % (newTi,newTf,newdT))

    def new_project(self,args):
        ''' Starting a new simulation project.'''
        newmodel = {'HomGo':HomogeneousGoModel, 'HetGo':HeterogeneousGoModel, 
                      'DMC':DMCModel}[args.type]
        #self.append_log("Project %s started" % args.name)
        Model = newmodel(self.path)
        System = system.System(args)
        self.create_subdirs(System)
        self.prepare_system(System,Model)

        ## The first step depends on the type of model.
        if args.type in ["HomGo","HetGo"]:
            self.folding_temperature_loop(Model,System)
        elif args.type == "DMC":
            pass
        #self.append_log("Analys T_array")
        #analysis.analyze_temperature_array(System)
        #self.append_log("Finished: T_array")
        #self.save_project(Model,System)
        print "Success"
        #raise SystemExit

    def save_project_status(self,Model,System):

        open("system.info","w").write(System.__repr__())
        open("model.info","w").write(Model.__repr__())

    def prepare_system(self,System,Model):
        ''' Extract all the topology files from Model.'''
        System.clean_pdbs()
        System.write_Native_pdb_CA()
        prots_indices, prots_residues, prots_coords = System.get_atom_indices(Model.beadmodel)
        prots_ndxs = Model.get_index_string(System.subdirs,prots_indices)
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
    #run_parser = sp.add_parser('check')
    #run_parser.add_argument('--output', type=str, help='Method of calculation',required=True)
    #run_parser.add_argument('--Ti', type=int, help='Initial Temperature')
    args = parser.parse_args()
    
    
    #print args
    #raise SystemExit
    ModelBuilder(args)

if __name__ == '__main__':
    main()
