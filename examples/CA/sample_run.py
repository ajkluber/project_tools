"""
Author: Alexander Kluber ajkluber@rice.edu
Date: June 2014

               Running a C-alpha structure-based model

    This script shows how to create and run a C-alpha structure-based model
starting from just a pdb structure and the list of contacts.

"""

import os
import numpy as np

import model_builder as mdb
import project_tools as pjt

# Creates coarse-grain model from the SH3.ini file in this directory. All
# files needed for simulation are automatically generated. 
name = "SH3"
model = mdb.inputs.new_model_from_config(name)

# Uncomment to Start simulations for a range of temperatures.
#nsteps = "400000"
#T_min = 50
#T_max = 150
#deltaT = 5
#cwd = os.getcwd()
#for T in range(T_min,T_max+deltaT,deltaT):
#    if not os.path.exists(str(T)):
#        os.mkdir(str(T))
#        print "  Running ",T
#    # Write simulation files to subdirectory and submit job
#    os.chdir(str(T))
#    model.save_simulation_files()
#    pjt.simulation.constant_temp.run_constant_temp(model,T,nsteps=nsteps)
#    os.chdir(cwd)

