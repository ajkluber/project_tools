""" Copy only necessary sub tree to proceed with mutation


"""

import argparse
import os
import shutil
from glob import glob

from model_builder import models
from model_builder import systems


def get_args():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--subdir', type=str, required=True, help='')
    parser.add_argument('--dest', type=str, required=True, help='')
    args = parser.parse_args()
    return args

if __name__ == "__main__":

    def dummy_func(sub,string):
        pass 

    args = get_args()
    subdir = args.subdir
    destination = args.dest

    Tf_choice = open(subdir+"/Mut_0/Tf_choice.txt","r").read()[:-1]
    Models = models.load_models([subdir],dryrun=True)
    Systems = systems.load_systems([subdir])
    Model = Models[0]
    System = Systems[0]
    path = System.subdir+"/"+System.mutation_active_directory+"/"+Tf_choice+"_agg"

    cwd = os.getcwd()
    sub = cwd+"/"+System.subdir+"/"+System.mutation_active_directory
    T = phi.get_Tf_choice(sub)
    savedir = sub+"/"+T+"_agg"

    paths = [ \
    "/Mut_0/"+Tf_choice+"_agg/mut",
    "/Mut_0/"+Tf_choice+"_agg/phi",
    "/mutants",
    "/Qref_shadow"]

    for P in paths:
        if not os.path.exists(destination+"/"+subdir+P):
            print " Making directory ",P
            os.makedirs(destination+"/"+subdir+P)

    files = [ \
    "/mutants/calculated_ddG.dat",
    "/Qref_cryst.dat",
    "/model.info",
    "/system.info",
    "/modelbuilder.log",
    "/Mut_0/Tf_choice.txt",
    "/Mut_0/"+Tf_choice+"_agg/mut/M.dat",
    "/Mut_0/"+Tf_choice+"_agg/mut/eps.dat",
    "/Mut_0/"+Tf_choice+"_agg/mut/ddG.dat",
    "/Mut_0/"+Tf_choice+"_agg/phi/Q_phi.dat" ]

    shutil.copy(subdir+".pdb",destination+"/")
    for file in files:
        if not os.path.exists(destination+"/"+subdir+file):
            print " copying ",file
            shutil.copy(subdir+file,destination+"/"+subdir+file)

    mutants = glob(subdir+"/mutants/fij*.dat")
    print " copying mutant fij_*dat"
    for mut in mutants:
        shutil.copy(mut,destination+"/"+subdir+"/mutants/")
