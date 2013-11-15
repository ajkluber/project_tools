import argparse
import os
import subprocess as sb
import time
import numpy as np

import model
import system
import mdp

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
        
        if args.action == 'new':
            Model = model.get_model(args.type,self.path)
            System = system.System(args)
            self.create_subdirs(System.subdirs)
            #print self.Model  ## DEBUGGING
            #print self.System  ## DEBUGGING
            Model.prep_system(System)
            self.run_temperature_array(Model,System)
            print "Success"
            raise SystemExit
            
            #self.write_info_file()
            #self.append_log("Project %s started on %s" % (args.name,time.asctime()))
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
            pass

    def append_log(self,string):
        logfile = open(self.path + '/modelbuilder.log','a').write(self.time_label()+' '+string+'\n')

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
        self.append_log("Creating new subdirectories: %s\n" % pdbs)

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
        
    def run_temperature_array(self,Model,System):
        ''' Run many constant temperature runs over a range of temperatures to
            find the folding temperature. '''
    
        Temperatures = range(100,200,10)
        for i in range(len(System.subdirs)):
            self.append_log("Starting Temperature array for protein: %s" % System.subdirs[i])
            T_string = ''
            for T in Temperatures:
                T_string += "%d_0\n" % T
                self.append_log("  running T=%d" % T)
                self.run_constant_temp(Model,System,T,i)
                ## Need to append logfile.
            open(self.path+"/"+System.subdirs[i]+"/T_array.txt","w").write(T_string)

    def run_constant_temp(self,Model,System,T,prot_num):
        ''' Start a constant temperature simulation with Gromacs. First it has
            to write the gromacs files stored in the System object, then it
            calls a function to submit the job.'''
        grompp_mdp = mdp.get_constant_temperature_mdp(Model,T)
        os.chdir(self.path)
        simpath = System.subdirs[prot_num]+"/"+str(T)+"_0"
        try:
            os.mkdir(simpath)
        except:
            pass
        for filename in System.topology_files[prot_num].iterkeys():
            #print "Writing: ", filename    ## DEBUGGING
            open(simpath+"/"+filename,"w").write(System.topology_files[prot_num][filename])
        open(simpath+"/grompp.mdp","w").write(grompp_mdp)
        open(simpath+"/Native.pdb","w").write(System.native_pdbs[prot_num])
        for m in range(len(Model.interaction_groups)):
            tablefile = "table_%s.xvg" % Model.interaction_groups[m]
            np.savetxt(simpath+"/"+tablefile,Model.tables[m],fmt="%16.15e",delimiter=" ")
        np.savetxt(simpath+"/table.xvg",Model.other_table,fmt="%16.15e",delimiter=" ")
        ## Start simulation
        os.chdir(simpath)
        self.submit_run(System.subdirs[prot_num]+"_"+str(T))
        
    def submit_run(self,jobname,walltime="24:00:00",queue="serial"):
        ''' Executes the constant temperature runs.'''

        prep_step1 = 'trjconv -f Native.pdb -o Native.xtc'
        prep_step2 = 'grompp -n index.ndx -c Native.pdb -maxwarn 1'
        prep_step3 = 'echo -e "System\\n" | trjconv -f Native.pdb -o conf.gro'
        #prep_step4 = 'grompp -n index.ndx -maxwarn 1'

        sb.call(prep_step1.split())
        sb.call(prep_step2.split())
        sb.call(prep_step3,shell=True)
        #sb.call(prep_step4.split())

        pbs_string = "#!/bin/bash \n"
        pbs_string +="### Number of nodes and procs/node \n"
        pbs_string +="#PBS -l nodes=1:ppn=1,walltime=%s \n" % walltime
        pbs_string +="###PBS -W group_list=pbc \n"
        pbs_string +="#PBS -q %s \n" % queue
        pbs_string +="#PBS -V \n"
        pbs_string +="### output files \n"
        pbs_string +="#PBS -o out \n"
        pbs_string +="#PBS -e err \n"
        pbs_string +="### Job Name (max 15 chars.) \n"
        pbs_string +="#PBS -N %s \n\n" % jobname
        pbs_string +="cd $PBS_O_WORKDIR\n"
        pbs_string +="mdrun -nt 1"

        open("run.pbs","w").write(pbs_string)
        qsub = "qsub run.pbs"
        sb.call(qsub.split())

    def time_label(self):
        now = time.localtime()
        now_string = "%s:%s:%s:%s" % (now.tm_mon,now.tm_mday,now.tm_hour,now.tm_min)
        return now_string
            
def main():
    ''' Two possible branches: 1. Calculate reference matrix, 2. Calculate Q '''
    parser = argparse.ArgumentParser(description='Calculate the (Non)Native contact matrix')
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
