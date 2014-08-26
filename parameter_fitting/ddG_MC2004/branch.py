""" Copy only necessary sub tree to proceed with mutation 

Use if you want to continue simulation

To Do:
- Write helper function that converts old naming convention to 
  new naming convention (mut -> newton) for backwards 
  compatibility.

"""

import argparse
import os
import shutil
from glob import glob

def get_args():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--name', type=str, required=True, help='')
    parser.add_argument('--iteration', type=int, required=True, help='')
    parser.add_argument('--dest', type=str, required=True, help='')
    args = parser.parse_args()
    return args

if __name__ == "__main__":

    args = get_args()
    name = args.name
    iteration = args.iteration
    destination = args.dest

    cwd = os.getcwd()
    path = "%s/Mut_%d/" % (name,iteration)
    sub = "%s/%s/Mut_%d" % (cwd,name,iteration)

    paths = ["/Mut_%d/mut" % iteration,"/mutants","/Qref_shadow"]

    for P in paths:
        if not os.path.exists("%s/%s%s" % (destination,name,P)):
            print " Making directory %s" % P
            os.makedirs("%s/%s%s" % (destination,name,P))

    files = [ \
    "/mutants/calculated_ddG.dat",
    "/Qref_cryst.dat",
    "/contacts.dat",
    "/model.info",
    "/modelbuilder.log",
    "/Mut_%d/T_array_last.txt" % iteration]

    shutil.copy("%s.pdb" % name, destination+"/")
    shutil.copy("%s_calculated_ddG.dat" % name, destination+"/")
    for file in files:
        if not os.path.exists("%s/%s%s" % (destination,name,file)):
            print " copying ",file
            shutil.copy("%s%s" % (name,file),"%s/%s%s" % (destination,name,file))

    
    mut_data = glob("%s/Mut_%d/mut/*" % (name,iteration))
    print " copying %s/Mut_%d/mut/*" % (name,iteration)
    for data in mut_data:
        shutil.copy(data,"%s/%s/Mut_%d/newton" % (destination,name,iteration))

    mutants = glob("%s/mutants/fij*.dat" % name)
    print " copying mutant fij_*dat"
    for mut in mutants:
        shutil.copy(mut,"%s/%s/mutants" % (destination,name))
