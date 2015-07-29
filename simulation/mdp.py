
"""
November 11 2013
Purpose:
    Easy way to generate grompp.mdp files for Gromacs simulations. When called
this module returns a grompp.mdp as a string. The requested variables will be 
already set (e.g. temperature, timestep, friction constant).

#TODO(alex):
- mdp for simulated annealing

"""

def constant_temperature(T,nsteps,nstout="1000"):
    """ Generate grompp.mdp file string. Gromacs 4.5 """
    mdp_string = "; Run control parameters \n"
    mdp_string += "integrator               = sd  \n"
    mdp_string += "dt                       = 0.0005 \n"
    mdp_string += "nsteps                   = %s \n\n" % nsteps
    mdp_string += "; output control options \n"
    mdp_string += "nstxout                  = 0 \n"
    mdp_string += "nstvout                  = 0 \n"
    mdp_string += "nstfout                  = 0 \n"
    mdp_string += "nstlog                   = 5000 \n"
    mdp_string += "nstenergy                = %s \n" % nstout
    mdp_string += "nstxtcout                = %s \n" % nstout
    mdp_string += "xtc_grps                 = system \n"
    mdp_string += "energygrps               = system \n\n" 
    mdp_string += "; neighborsearching parameters \n"
    mdp_string += "nstlist                  = 20 \n"
    mdp_string += "ns-type                  = grid \n"
    mdp_string += "pbc                      = no \n"
    mdp_string += "periodic_molecules       = no \n"
    mdp_string += "rlist                    = 2.0 \n"
    mdp_string += "rcoulomb                 = 2.0 \n"
    mdp_string += "rvdw                     = 2.0  \n\n"
    mdp_string += "; options for electrostatics and vdw \n"
    mdp_string += "coulombtype              = User \n"
    mdp_string += "vdw-type                 = User \n"
    mdp_string += "table-extension          = 1.0 \n\n"
    mdp_string += "; options for temp coupling \n"
    mdp_string += "Tcoupl                   = no \n"
    mdp_string += "ld_seed                  = -1 \n"
    mdp_string += "tc-grps                  = system \n"
    mdp_string += "tau_t                    = 1 \n"
    mdp_string += "ref_t                    = %s \n" % T
    mdp_string += "Pcoupl                   = no \n\n"
    mdp_string += "; generate velocities for startup run \n"
    mdp_string += "gen_vel                  = yes \n"
    mdp_string += "gen_temp                 = %s \n" % T
    mdp_string += "gen_seed                 = -1 \n\n"
    mdp_string += "; remove center of mass\n"
    mdp_string += "comm_mode                = angular \n"
    mdp_string += "comm_grps                = System \n"
    return mdp_string

def simulated_annealing(Tlist,pslist,nsteps,nstout="1000"):
    """ Generate grompp.mdp file string. Gromacs 4.6 """
    # Should we assert that the total number of steps is longer than the
    # annealing schedule? 

    annealtimes = ""
    annealtemps = ""
    for i in range(len(Tlist)):
        annealtimes += "%.2f " % (pslist[i] + sum(pslist[:i]))
        annealtemps += "%.2f " % Tlist[i]

    mdp_string = "; Run control parameters \n"
    mdp_string += "integrator               = sd  \n"
    mdp_string += "dt                       = 0.0005 \n"
    mdp_string += "nsteps                   = %s \n\n" % nsteps
    mdp_string += "; output control options \n"
    mdp_string += "nstxout                  = 0 \n"
    mdp_string += "nstvout                  = 0 \n"
    mdp_string += "nstfout                  = 0 \n"
    mdp_string += "nstlog                   = 5000 \n"
    mdp_string += "nstenergy                = %s \n" % nstout
    mdp_string += "nstxtcout                = %s \n" % nstout
    mdp_string += "xtc_grps                 = system \n"
    mdp_string += "energygrps               = system \n\n" 
    mdp_string += "; neighborsearching parameters \n"
    mdp_string += "nstlist                  = 20 \n"
    mdp_string += "ns-type                  = grid \n"
    mdp_string += "pbc                      = no \n"
    mdp_string += "periodic_molecules       = no \n"
    mdp_string += "rlist                    = 2.0 \n"
    mdp_string += "rcoulomb                 = 2.0 \n"
    mdp_string += "rvdw                     = 2.0  \n\n"
    mdp_string += "; options for electrostatics and vdw \n"
    mdp_string += "coulombtype              = User \n"
    mdp_string += "vdw-type                 = User \n"
    mdp_string += "table-extension          = 1.0 \n\n"
    mdp_string += "; options for temp coupling \n"
    mdp_string += "Tcoupl                   = no \n"
    mdp_string += "ld_seed                  = -1 \n"
    mdp_string += "tc-grps                  = system \n"
    mdp_string += "tau_t                    = 1 \n"
    mdp_string += "ref_t                    = %s \n" % Tlist[0]
    mdp_string += "Pcoupl                   = no \n\n"
    mdp_string += "; generate velocities for startup run \n"
    mdp_string += "gen_vel                  = yes \n"
    mdp_string += "gen_temp                 = %s \n" % Tlist[0]
    mdp_string += "gen_seed                 = -1 \n\n"
    mdp_string += "; remove center of mass\n"
    mdp_string += "comm_mode                = angular \n"
    mdp_string += "comm_grps                = System \n"
    mdp_string += "; simulated annealing schedule\n"
    mdp_string += "annealing                = single\n"
    mdp_string += "annealing-npoints        = %d\n" % len(Tlist)
    mdp_string += "annealing-time           = %s\n" % annealtimes
    mdp_string += "annealing-temp           = %s\n" % annealtemps
    return mdp_string
