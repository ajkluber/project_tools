import numpy as np
import os
import subprocess as sb


'''
    Feb 5 2014
    Alexander Kluber

        Starting to implement the mutations to the heterogeneous Go model.
'''





if __name__ == '__main__':

    ## Start by reading the mutations file. Should be an attribute of System
    ## object.

    muts = open("mutations.txt","r").readlines()[1:]

    

    ## If they don't exist, use MODELLER to create the mutated pdbs.


    ## Use shadow map to create all-atom contact map. For each mutated pdb
    ## determine the 
