"""
Author: Alexander Kluber ajkluber@rice.edu
Date: June 2014

               Running a C-alpha structure-based model

    This script shows how to create and run a C-alpha structure-based model
starting from just a pdb structure.

Directions:
    1. Run this script on the command line within a directory that contains 
1SHG.pdb.

"""

import os
import numpy as np
import model_builder as mdb
import project_tools as pjt

name = "1FMK"
pdb = "%s.pdb" % name
nsteps = "400000"
T_min = 50
T_max = 150
deltaT = 5


contacts = np.loadtxt("%s.contacts" % name)
LJtype = np.ones(len(contacts),int)

## Initialize a C-alpha Go-model with a pdb. All topology files needed for
## simulation are automatically generated. 
model = mdb.models.SmogCalpha.SmogCalpha(pdb,contacts=contacts,LJtype=LJtype)


### Start simulations for a range of temperatures.
#cwd = os.getcwd()
#for T in range(T_min,T_max+deltaT,deltaT):
#    if not os.path.exists(str(T)):
#        os.mkdir(str(T))
#        print "  Running ",T
#    ## Write simulation files to subdirectory and submit job
#    os.chdir(str(T))
#    pjt.simulation.constant_temp.run_constant_temp(model,T,nsteps=nsteps)
#    os.chdir(cwd)

