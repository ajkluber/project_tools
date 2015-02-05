""" Top-level class that automates the execution of a recipe

Description:

ProjectManager
    The ProjectManager class is a task manager that executes a preset procedure
specified in a recipe. Recipes hold the logical flow to for example perform
parameter fitting of a structure-based model.


    
"""

import os
import time

from project_tools import simulation, analysis, parameter_fitting
import model_builder as mdb

global GAS_CONSTANT_KJ_MOL
GAS_CONSTANT_KJ_MOL = 0.0083144621  # Gas constant in kJ/(mole K)

class ProjectManager(object):
    """ A shell class to handle the simulations for a project.

    Description:

        A class that can handle simulation projects.
    """

    def __init__(self,args):
        self.path = os.getcwd()
        
        if args.action == 'new':
            self.new_project(args)
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

    def add_temperature_array(self,args):
        """ Adds manually adds a temperature array."""

        subdirs = args.subdirs
        Models = mdb.inputs.load_models(subdirs,dry_run=args.dry_run)
    
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
        """ Manually extends."""
        factor = args.factor

        subdirs = args.subdirs
        Models = mdb.inputs.load_models(subdirs,dry_run=args.dry_run)

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

        for i in range(len(subdirs)):
            model = Models[i]
            sub = model.subdir
            simulation.constant_temp.manually_extend_temperatures(model,self.append_log,method,temps,factor)
            
        self.save_model_info(Models)
        print "Success"


    def continue_project(self,args):
        """ Checks where something left off and continues it."""
        names = args.names
        Models,Fittingopts = mdb.inputs.load_models(names,dry_run=args.dry_run)

        for i in range(len(names)):
            model = Models[i]
            fitopts = Fittingopts[i]
            print "Checking progress for directory:  %s" % model.name
            print "Last task was %s" % fitopts["last_completed_task"]
            action = fitopts["last_completed_task"].split()[0]
            if action == "Starting:":
                self.logical_flowchart_starting(model,fitopts)
            elif action == "Finished:":
                self.logical_flowchart_finished(model,fitopts)
            else:
                raise IOError("Not valid action: %s" % action)
            mdb.inputs.save_model(model,fitopts)
        print "Success"
