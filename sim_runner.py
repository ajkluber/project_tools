import numpy as np


import modelbuilder


'''

Purpose:

    Build simulation files for coarse-grain simulations.
Purpose:
    sim_runner is a top-level helper class to track and execute varied
procedures for coarse-grain modeling. At this level the tasks that need to be
done are abstract (such as 'determining folding temperature' or 'refining
parameter'). sim_runner relies on the modelbuilder to prepare the input
files. 
    The goal is to conceal as much of the details as possible away from the
user, so that the user can focus on top-level questions. For this reason any
function that requires manipulating the data is best moved to a different 
module. 
    

'''
