
'''
November 11 2013
Purpose:
    Easy way to generate grompp.mdp files for Gromacs simulations. When called
this module returns a grompp.mdp as a string. The requested variables will be 
already set (e.g. temperature, timestep, friction constant).

To Do:
- Add different functions for different styles of simulations. E.g. simulated 
  annealing simulations. (potentially other types in the future, like pulling
  simulations)
  a. Add many options to each style of parameter file. 
  b. There are probably many options that we haven't thought of modifying 
    before.

'''

def get_constant_temperature_mdp(model,T,nsteps=400000000):
    ''' Generate grompp.mdp file string. Default is a ~24hr run.'''
    mdp_string = '; RUN CONTROL PARAMETERS \n'
    mdp_string += 'integrator               = sd  \n'
    mdp_string += 'tinit                    = 0 \n'
    mdp_string += 'dt                       = 0.005 \n'
    mdp_string += 'nsteps                   = %d \n\n' % nsteps
    mdp_string += '; For exact run continuation or redoing part of a run \n'
    mdp_string += 'init_step                = 0 \n'
    mdp_string += 'comm_mode                = angular \n'
    mdp_string += 'nstcalcenergy            = 10 \n'
    mdp_string += 'nstcomm                  = 10 \n'
    mdp_string += 'comm_grps                = System \n\n'
    mdp_string += '; OUTPUT CONTROL OPTIONS \n'
    mdp_string += 'nstxout                  = 0 \n'
    mdp_string += 'nstvout                  = 0 \n'
    mdp_string += 'nstfout                  = 0 \n'
    mdp_string += 'nstlog                   = 50000 \n'
    mdp_string += 'nstenergy                = 10000 \n'
    mdp_string += 'nstxtcout                = 10000 \n'
    mdp_string += 'xtc_grps                 = System \n'
    mdp_string += 'energygrps               = %s \n\n' % model.energygrps
    mdp_string += '; NEIGHBORSEARCHING PARAMETERS \n'
    mdp_string += 'nstlist                  = 10 \n'
    mdp_string += 'ns-type                  = grid \n'
    mdp_string += 'pbc                      = no \n'
    mdp_string += 'periodic_molecules       = no \n'
    mdp_string += 'rlist                    = 2.0 \n\n'
    mdp_string += '; OPTIONS FOR ELECTROSTATICS AND VDW \n'
    mdp_string += '; Method for doing electrostatics \n'
    mdp_string += 'coulombtype              = User \n'
    mdp_string += 'rcoulomb-switch          = 0   \n'
    mdp_string += 'rcoulomb                 = 2.0 \n'
    mdp_string += '; Method for doing Van der Waals \n'
    mdp_string += 'vdw-type                 = User \n'
    mdp_string += 'rvdw-switch              = 0  \n'
    mdp_string += 'rvdw                     = 2.0  \n'
    mdp_string += '; Apply long range dispersion corrections for Energy and Pressure \n'
    mdp_string += 'DispCorr                 = no \n'
    mdp_string += '; Extension of the potential lookup tables beyond the cut-off \n'
    mdp_string += 'table-extension          = 1.0 \n\n'
    mdp_string += '; Seperate tables between energy group pairs \n'
    mdp_string += 'energygrp_table          = %s \n\n' % model.energygrp_table
    mdp_string += '; IMPLICIT SOLVENT ALGORITHM \n'
    mdp_string += 'implicit_solvent         = No \n\n'
    mdp_string += '; OPTIONS FOR WEAK COUPLING ALGORITHMS \n'
    mdp_string += 'Tcoupl                   = no \n'
    mdp_string += 'ld_seed                  = -1 \n'
    mdp_string += 'tc-grps                  = System \n\n'
    mdp_string += 'tau_t                    = 0.2 \n'
    mdp_string += 'ref_t                    = %d \n\n' % T
    mdp_string += 'Pcoupl                   = no \n\n'
    mdp_string += '; OPTIONS FOR BONDS     \n'
    mdp_string += 'constraints              = none \n\n'
    mdp_string += '; GENERATE VELOCITIES FOR STARTUP RUN \n'
    mdp_string += 'gen-vel                  = yes \n'
    mdp_string += 'gen_temp                 = %f \n' % T
    mdp_string += 'gen_seed                 = -1 \n'
    return mdp_string
