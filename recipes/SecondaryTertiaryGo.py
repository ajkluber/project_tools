""" A recipe to run the Secondary Tertiary Go-model


CURRENTLY NOT RIGHT.

UPDATE TO REMOVE "model_builder.systems" REFERENCE


Description:

    This recipes performs an algorithm to add energetic heterogeneity to a
Go-model considering secondary interactions (backbone H-bonds) and tertiary
interactions (sidechain contacts) separately.

    The hope is that by treating backbone H-bonding by a set of weaker,
triangulating contacts between adjacent residues the secondary structure will
become more 'cooperative'. Cooperative in this case is in the sense of
helix-coil theory where forming one helical residue somewhat preorganizes the
next H-bond to add on.

Reference:
"""

import os
import argparse

from project_tools import simulation, analysis, parameter_fitting
from recipe_manager import ProjectManager
from model_builder import models


class SecondaryTertiaryGo(ProjectManager):
    
    """ A recipe to run the Secondary Tertiary Go-model


    Description:

        This recipes performs an algorithm to add energetic heterogeneity to a
    Go-model considering secondary interactions (backbone H-bonds) and tertiary
    interactions (sidechain contacts) separately.

        The hope is that by treating backbone H-bonding by a set of weaker,
    triangulating contacts between adjacent residues the secondary structure will
    become more 'cooperative'. Cooperative in this case is in the sense of
    helix-coil theory where forming one helical residue somewhat preorganizes the
    next H-bond to add on.

    Reference:
    """

def get_args():
    """ Get command line arguments """

    parser = argparse.ArgumentParser(description='Run .')
    sp = parser.add_subparsers(dest='action')

    ## Options for initializing a new simulation project.
    new_parser = sp.add_parser('new')
    new_parser.add_argument('--pdbs', type=str, required=True, nargs='+',help='PDBs to start simulations.')
    new_parser.add_argument('--dryrun', action='store_true', help='Add this option for dry run. No simulations started.')
    new_parser.add_argument('--temparray', type=int, nargs='+',help='Optional initial temp array: Ti Tf dT. Default: 50 350 50')

    ## Options for continuing from a previously saved simulation project.
    run_parser = sp.add_parser('continue')
    run_parser.add_argument('--subdirs', type=str, nargs='+', help='Subdirectories to continue',required=True)
    run_parser.add_argument('--dryrun', action='store_true', help='Dry run. No simulations started.')

    ## Options for manually adding a temperature array.
    add_parser = sp.add_parser('add')
    add_parser.add_argument('--subdirs', type=str, nargs='+', help='Subdirectories to add temp array',required=True)
    add_parser.add_argument('--temparray', type=int, nargs='+', help='T_initial T_final dT for new temp array',required=True)
    add_parser.add_argument('--mutarray', type=int, nargs='+', help='T_initial T_final dT for new mutational sims array')
    add_parser.add_argument('--dryrun', action='store_true', help='Dry run. No simulations started.')

    args = parser.parse_args()

    if args.dryrun != False:
        options = {"Dry_Run":True}
    else:
        options = {"Dry_Run":False}

    options["Model_Code"] = "HetGo"
    options["Bead_Model"] = "CA"
    options["Solvent"] = None
    options["R_CD"] = None
    options["Epsilon_Bar"] = 1.0
    options["Disulfides"] = None
    options["Contact_Energies"] = "SecTer"

    modeloptions = models.check_options(options)

    return args, modeloptions

if __name__ == "__main__":
    
    args, options = get_args()
    SecondaryTertiaryGo(args,options)
