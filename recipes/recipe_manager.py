''' Top-level class that automates the execution of a recipe

Description:

ProjectManager
    The ProjectManager class is a task manager that executes a preset procedure
specified in a recipe. Recipes hold the logical flow to for example perform
parameter fitting of a structure-based model.


    
'''

import os
import time

from project_tools import simulation, analysis, parameter_fitting
import model_builder as mdb

global GAS_CONSTANT_KJ_MOL
GAS_CONSTANT_KJ_MOL = 0.0083144621  ## Gas constant in kJ/(mole K)

def print_header():

    print "------------------------------ project_tools ---------------------------------"
    print " Your using project_tools, a multiscale toolbox from the Clementi lab"
    print " Version 0.0 \n"
    #print " Words to live by:\n"
    print "             'If you can calculate it, you should calculate it' - PGW \n"
    #print "               'One never notices what has been done; "
    #print "                one can only see what remains to be done' - Marie Curie \n"
    #print "  'The best way to have a good idea is to have a lot of ideas' - Linus Pauling\n"
    #print "'A ship in port is safe, but that's not what ships are built for' - Grace Hopper'\n"
    #print "                 'Science and everyday life cannot and"
    #print "                  should not be separated' - Rosalind Franklin\n"
    print "------------------------------------------------------------------------------"

class ProjectManager(object):
    ''' A shell class to handle the simulations for a project.

    Description:

        A class that can handle simulation projects.
    '''

    def __init__(self,args,modeloptions):
        self.path = os.getcwd()
        
        if args.action == 'new':
            self.new_project(args,modeloptions)
        elif args.action == 'continue':
            self.continue_project(args)
        elif args.action == 'add':
            self.add_temperature_array(args)
        elif args.action == 'extend':
            self.extend_temperatures(args)

    def append_log(self,sub,string,subdir=False):
        if subdir:
            open('%s/%s/%s.log' % (self.path,sub,sub),'a').write("%s %s\n" % (self.append_time_label(),string))
        else:
            open('%s/%s/modelbuilder.log' % (self.path,sub),'a').write("%s %s\n" % (self.append_time_label(),string))

    def append_time_label(self):
        now = time.localtime()
        now_string = "%s:%s:%s:%s:%s" % (now.tm_year,now.tm_mon,now.tm_mday,now.tm_hour,now.tm_min)
        return now_string

    def create_subdirs(self,Models):
        ''' Create the subdirectories for the system. Skip if they already exist.'''
        for i in range(len(Models)): 
            sub = Models[i].subdirs
            if os.path.exists(sub) != True:
                os.mkdir(sub)
                os.mkdir("%s/iteration_%d" % (sub,Models[i].iteration))
                self.append_log(sub,"Creating new subdirectory: %s" % sub)

    def add_temperature_array(self,args):
        ''' Adds manually adds a temperature array.'''

        subdirs = args.subdirs
        Models = mdb.check_inputs.load_models(subdirs,dry_run=args.dry_run)
    
        if args.temparray != None:
            T_min = args.temparray[0] 
            T_max = args.temparray[1] 
            deltaT = args.temparray[2] 
            long = False
        elif args.mutarray != None:
            temps = args.mutarray
            long = True
        else:
            print "ERROR! Must use --temparray or --mutarray with this option"
            print " Exiting."
            raise SystemExit

        for i in range(len(Models)):
            model = Models[i]
            sub = model.subdir
            if long == False:
                print "Manually adding temperature array Ti=%d Tf=%d dT=%d" % (T_min,T_max,deltaT)
                print "Starting constant_temp_iteration..."
                simulation.constant_temp.manually_add_temperature_array(model,self.append_log,T_min,T_max,deltaT)
            elif long == True:
                print "Manually adding equilibrium sims ", temps
                simulation.constant_temp.manually_add_equilibrium_runs(model,self.append_log,temps)
            
            
        self.save_model_info(Models)
        print "Success"

    def extend_temperatures(self,args):
        ''' Manually extends.'''
        factor = args.factor

        subdirs = args.subdirs
        Models = mdb.check_inputs.load_models(subdirs,dry_run=args.dry_run)

        if (args.Tf_temps != None) and (args.Mut_temps != None):
            print "ERROR!"
            print "  Specify either --Tf_temps or --Mut_temps. Not both!"
            print "  Exiting"
            raise SystemExit
        if (args.Tf_temps == None) and (args.Mut_temps == None):
            print "ERROR!"
            print "  Specify either --Tf_temps or --Mut_temps."
            print "  Exiting"
            raise SystemExit
        if (args.Tf_temps != None):
            temps = args.Tf_temps
            method = "short"
            print "Extending Tf temps",temps
        if (args.Mut_temps != None):
            temps = args.Mut_temps
            method = "long"
            print "Extending Mutational temps",temps

        for i in range(len(subdirs)):
            model = Models[i]
            sub = model.subdir
            simulation.constant_temp.manually_extend_temperatures(model,self.append_log,method,temps,factor)
            
        self.save_model_info(Models)
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
        Models = mdb.check_inputs.load_models(subdirs,dry_run=args.dry_run)

        if args.dry_run == True:
            print "  Dry run complete. Exiting."
        else:
            for i in range(len(subdirs)):
                model = Models[i]        
                subdir = model.subdir
                
                lasttime,action,task = self.check_modelbuilder_log(subdir)
                print "Checking progress for directory:  ", subdir
                print "Last task was %s %s at %s" % (action,task,lasttime)
                if action == "Starting:":
                    self.logical_flowchart_starting(model,task)
                elif action == "Finished:":
                    self.logical_flowchart_finished(model,task)
                elif action == "Error":
                    pass

        self.save_model_info(Models)
        print "Success"

    def save_model_info(self,Models):
        ''' Save the model.info strings.'''
        print "Saving model.info ..."
        for i in range(len(Models)):
            open(Models[i].subdir+"/model.info","w").write(Models[i].get_model_info_string())
