
"""
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

"""

def constant_temperature(T,nsteps):
    """ Generate grompp.mdp file string."""
    mdp_string = '; RUN CONTROL PARAMETERS \n'
    mdp_string += 'integrator               = sd  \n'
    mdp_string += 'dt                       = 0.0005 \n'
    mdp_string += 'nsteps                   = %s \n\n' % nsteps
    mdp_string += '; OUTPUT CONTROL OPTIONS \n'
    mdp_string += 'nstxout                  = 0 \n'
    mdp_string += 'nstvout                  = 0 \n'
    mdp_string += 'nstfout                  = 0 \n'
    mdp_string += 'nstlog                   = 5000 \n'
    mdp_string += 'nstenergy                = 1000 \n'
    mdp_string += 'nstxtcout                = 1000 \n'
    mdp_string += 'xtc_grps                 = system \n'
    mdp_string += 'energygrps               = system \n\n' 
    mdp_string += '; NEIGHBORSEARCHING PARAMETERS \n'
    mdp_string += 'nstlist                  = 20 \n'
    mdp_string += 'ns-type                  = grid \n'
    mdp_string += 'pbc                      = no \n'
    mdp_string += 'periodic_molecules       = no \n'
    mdp_string += 'rlist                    = 2.0 \n'
    mdp_string += 'rcoulomb                 = 2.0 \n'
    mdp_string += 'rvdw                     = 2.0  \n\n'
    mdp_string += '; OPTIONS FOR ELECTROSTATICS AND VDW \n'
    mdp_string += 'coulombtype              = User \n'
    mdp_string += 'vdw-type                 = User \n'
    mdp_string += 'table-extension          = 1.0 \n\n'
    mdp_string += '; OPTIONS FOR TEMP COUPLING \n'
    mdp_string += 'Tcoupl                   = no \n'
    mdp_string += 'ld_seed                  = -1 \n'
    mdp_string += 'tc-grps                  = system \n'
    mdp_string += 'tau_t                    = 1 \n'
    mdp_string += 'ref_t                    = %s \n' % T
    mdp_string += 'Pcoupl                   = no \n\n'
    mdp_string += '; GENERATE VELOCITIES FOR STARTUP RUN \n'
    mdp_string += 'gen_vel                  = yes \n'
    mdp_string += 'gen_temp                 = %s \n' % T
    mdp_string += 'gen_seed                 = -1 \n\n'
    mdp_string += '; REMOVE CENTER OF MASS\n'
    mdp_string += 'comm_mode                = angular \n'
    mdp_string += 'comm_grps                = System \n'
    return mdp_string
