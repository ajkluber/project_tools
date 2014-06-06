""" Copy only necessary sub tree to proceed with mutation


"""

import argparse
import os
import shutil
from glob import glob

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
    Models = models.load_models(subdirs,dryrun=True)
    Systems = systems.load_systems(subdirs)
    Model = Models[0]
    System = Systems[0]
    path = System.subdir+"/"+System.mutation_active_directory+"/"+Tf_choice+"_agg"

    cwd = os.getcwd()
    sub = cwd+"/"+System.subdir+"/"+System.mutation_active_directory
    T = phi.get_Tf_choice(sub)
    savedir = sub+"/"+T+"_agg"

    os.makedirs(destination+"/"+subdir+"/Mut_0/"+Tf_choice+"_agg/mut")
    os.makedirs(destination+"/"+subdir+"/Mut_0/"+Tf_choice+"_agg/phi")
    os.makedirs(destination+"/"+subdir+"/mutants")
    os.makedirs(destination+"/"+subdir+"/Qref_shadow")


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
    "/Mut_0/"+Tf_choice+"_agg/phi/Q_phi.dat",
    ]

    shutil.copy(subdir+".pdb",destination+"/")
    for file in files:
        shutil.copy(subdir+file,destination+"/"+subdir+file)

    mutants = glob(subdir+"/mutants/fij*.dat")
    for mut in mutants:
        shutil.copy(mut,destination+"/"+subdir+"/mutants/")
