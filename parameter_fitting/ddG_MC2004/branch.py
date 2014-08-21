""" Copy only necessary sub tree to proceed with mutation


"""

import argparse
import os
import shutil
from glob import glob

import model_builder as mdb

import phi_values as phi
import mutatepdbs as mut


def get_args():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--subdir', type=str, required=True, help='')
    parser.add_argument('--dest', type=str, required=True, help='')
    args = parser.parse_args()
    return args

if __name__ == "__main__":

    args = get_args()
    subdir = args.subdir
    destination = args.dest

    Models = mdb.models.load_models([subdir],dryrun=True)
    model = Models[0]
    path = model.subdir+"/Mut_0/"

    cwd = os.getcwd()
    sub = cwd+"/"+model.subdir+"/Mut_0"

    paths = ["/Mut_0/mut","/mutants","/Qref_shadow"]

    for P in paths:
        if not os.path.exists(destination+"/"+subdir+P):
            print " Making directory ",P
            os.makedirs(destination+"/"+subdir+P)

    files = [ \
    "/mutants/calculated_ddG.dat",
    "/Qref_cryst.dat",
    "/contacts.dat",
    "/model.info",
    "/modelbuilder.log",
    "/Mut_0/T_array_last.txt"]

    shutil.copy(subdir+".pdb",destination+"/")
    shutil.copy(subdir+"_calculated_ddG.dat",destination+"/")
    for file in files:
        if not os.path.exists(destination+"/"+subdir+file):
            print " copying ",file
            shutil.copy(subdir+file,destination+"/"+subdir+file)

    
    mut_data = glob(subdir+"/Mut_0/mut/*")
    print " copying Mut_0/mut/*"
    for data in mut_data:
        shutil.copy(data,destination+"/"+subdir+"/Mut_0/mut/")

    mutants = glob(subdir+"/mutants/fij*.dat")
    print " copying mutant fij_*dat"
    for mut in mutants:
        shutil.copy(mut,destination+"/"+subdir+"/mutants/")
