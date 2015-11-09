''' Copy only necessary sub tree to proceed with mutation 

Use if you want to continue simulation

To Do:
- Write helper function that converts old naming convention to 
  new naming convention (mut -> newton) for backwards 
  compatibility.

- Update to new directory organization scheme.
  Instead of name/iteration_# -> name/iteration_#

'''

import argparse
import os
import shutil
from glob import glob

def branch(args):
    name = args.name
    iteration = args.iteration
    destination = args.dest

    cwd = os.getcwd()
    path = "%s/iteration_%d/" % (name,iteration)
    sub = "%s/%s/iteration_%d" % (cwd,name,iteration)

    paths = ["/iteration_%d/mut" % iteration,"/mutants","/Qref_shadow"]

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
    "/iteration_%d/T_array_last.txt" % iteration]

    shutil.copy("%s.pdb" % name, destination+"/")
    shutil.copy("%s_calculated_ddG.dat" % name, destination+"/")
    for file in files:
        if not os.path.exists("%s/%s%s" % (destination,name,file)):
            print " copying ",file
            shutil.copy("%s%s" % (name,file),"%s/%s%s" % (destination,name,file))

    
    mut_data = glob("%s/iteration_%d/mut/*" % (name,iteration))
    print " copying %s/iteration_%d/mut/*" % (name,iteration)
    for data in mut_data:
        shutil.copy(data,"%s/%s/iteration_%d/newton" % (destination,name,iteration))

    mutants = glob("%s/mutants/fij*.dat" % name)
    print " copying mutant fij_*dat"
    for mut in mutants:
        shutil.copy(mut,"%s/%s/mutants" % (destination,name))

def update_names(args):
    name = args.name
    iteration = args.iteration

    cwd = os.getcwd()
    sub = "%s/%s/iteration_%d" % (cwd,name,iteration)

    if not os.path.exists("%s/newton" % sub):
        os.mkdir("%s/newton" % sub) 

    shutil.copy("%s/mut/ddGexp.dat" % sub, "%s/newton/target_feature.dat" % sub)
    shutil.copy("%s/mut/ddGexp_err.dat" % sub, "%s/newton/target_feature_err.dat" % sub)
    shutil.copy("%s/mut/ddGsim.dat" % sub, "%s/newton/sim_feature.dat" % sub)
    shutil.copy("%s/mut/ddGsim_err.dat" % sub, "%s/newton/sim_feature_err.dat" % sub)
    shutil.copy("%s/mut/M.dat" % sub, "%s/newton/Jacobian.dat" % sub)
    shutil.copy("%s/mut/Merr.dat" % sub, "%s/newton/Jacobian_err.dat" % sub)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='')
    sp = parser.add_subparsers(dest='action')

    branch_parser = sp.add_parser("branch")
    branch_parser.add_argument('--name', type=str, required=True, help='')
    branch_parser.add_argument('--iteration', type=int, required=True, help='')
    branch_parser.add_argument('--dest', type=str, required=True, help='')

    update_parser = sp.add_parser("update")
    update_parser.add_argument('--name', type=str, required=True, help='')
    update_parser.add_argument('--iteration', type=int, required=True, help='')

    args = parser.parse_args()

    if args.action == "branch":
        branch(args)
    elif args.action == "update":
        update_names(args)
    else:
        pass
