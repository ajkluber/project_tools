""" Parameter fitting using Newton's method

Description:

    This submodule performs parameter fitting to match simulation features to
experimental (or designed) data. This parameter fitting is an extension of the
Matysiak, Clementi 2004 (see reference (1)). It uses the multivariate
Newton-Raphson method to correct simulation parameters in order to fit some
inputted target features. See references (2),(3) for more on Netwon's method in
general.


Submodules:

  newton 
      Solves for new parameters by Newton's method: by inverting the 
    Jacobian.

  ddG_MC2004
      Calculates Jacobian and feature vector for Delta Delta G's to match
    experimental data from phi-value analysis.

  contact_Qi
      Calculates Jacobian and feature vector for contact probabilities.


References:

1. Matysiak, S.; Clementi, C. Optimal Combination of Theory and Experiment for
the Characterization of the Protein Folding Landscape of S6: How Far Can a
Minimalist model Go? J. Mol. Biol. 2004, 343, 235-248.

2. J.E. Dennis; Robert B. Schnabel. "Numerical Methods for Unconstrained
Optimization and Nonlinear Equations". SIAM. 1996.

3. Jorge Nocedal; Stephen Wright. "Numerical Optimization". Springer. 2000

"""

import os 
import numpy as np
import logging

import newton_solver
import ddG_MC2004
import FRET
import RMSF
import contact_Qi

def prepare_newtons_method(model,fitopts):
    """ Prepare the files to do newtons method """

    name = model.name
    iteration = fitopts["iteration"]
    method = fitopts["data_type"]

    logging.basicConfig(filename="%s.log" % name,level=logging.INFO)
    logger = logging.getLogger("parameter_fitting")

    available_methods = ["ddG_MC2004","FRET","RMSF","contact_Qi"]
    modules = {"ddG_MC2004":ddG_MC2004,"FRET":FRET,"RMSF":RMSF,"contact_Qi":contact_Qi}

    if method not in available_methods:
        raise ValueError("Method %s not in available methods %s" % (method, available_methods.__repr__()) )
    else:
        submodule = modules[method]

    if not os.path.exists("%s/iteration_%d/newton/Jacobian.dat" % (name,iteration)):
        logger.info(" Starting: Calculating_Jacobian")

        sim_feature_avg, sim_feature_err, Jacobian_avg, Jacobian_err = submodule.compute_Jacobian.calculate_average_Jacobian(model,fitopts)
        target_feature, target_feature_err = submodule.compute_Jacobian.get_target_feature(model)

        if not os.path.exists("%s/iteration_%d/newton" % (name,iteration)):
            os.mkdir("%s/iteration_%d/newton" % (name,iteration))
        if not os.path.exists("%s/iteration_%d/%s" % (name,iteration,method)):
            os.mkdir("%s/iteration_%d/%s" % (name,iteration,method))

        print "  Saving feature vector and Jacobian in %s/iteration_%d/%s" % (name,iteration,method)
        np.savetxt("%s/iteration_%d/%s/target_feature.dat" % (name,iteration,method), target_feature)
        np.savetxt("%s/iteration_%d/%s/target_feature_err.dat" % (name,iteration,method), target_feature_err)
        np.savetxt("%s/iteration_%d/%s/sim_feature.dat" % (name,iteration,method), sim_feature_avg)
        np.savetxt("%s/iteration_%d/%s/sim_feature_err.dat" % (name,iteration,method), sim_feature_err)
        np.savetxt("%s/iteration_%d/%s/Jacobian.dat" % (name,iteration,method), Jacobian_avg)
        np.savetxt("%s/iteration_%d/%s/Jacobian_err.dat" % (name,iteration,method) ,Jacobian_err)

        ## To Do:
        ##  - Code Fitting_Includes option
        ##      - Collect Jacobian rows from all fitting_includes directories.
        ##      - Map columns (parameters) to match those of the first directory. Stack the rows.
        ##      - Save in the first fitting directory.

        print "  Saving feature vector and Jacobian in %s/iteration_%d/newton" % (name,iteration)
        np.savetxt("%s/iteration_%d/newton/target_feature.dat" % (name,iteration), target_feature)
        np.savetxt("%s/iteration_%d/newton/target_feature_err.dat" % (name,iteration), target_feature_err)
        np.savetxt("%s/iteration_%d/newton/sim_feature.dat" % (name,iteration), sim_feature_avg)
        np.savetxt("%s/iteration_%d/newton/sim_feature_err.dat" % (name,iteration), sim_feature_err)
        np.savetxt("%s/iteration_%d/newton/Jacobian.dat" % (name,iteration), Jacobian_avg)
        np.savetxt("%s/iteration_%d/newton/Jacobian_err.dat" % (name,iteration) ,Jacobian_err)

        logger.info(" Finished: Calculating_Jacobian")

    solve_newtons_method(model,fitopts)

def solve_newtons_method(model,fitopts):
    """ Solve the newton problem """
    name = model.name
    iteration = model.iteration
    method = fitopts["data_type"]

    solver_opts = { "Levenberg":newton_solver.Levenberg_Marquardt,
                    "TSVD":newton_solver.Truncated_SVD,
                    "TSVD_Cplex":newton_solver.Truncated_SVD_cplex}

    if fitopts["solver"] not in solver_opts.keys():
        print "ERROR! requested solver algorithm %s not found!" % fitopts["solver"]
        print " Choose from available solvers: ", solver_opts.keys()
        print " Exiting."
        raise SystemExit
    else:
        solver = solver_opts[fitopts["solver"]]

    logging.basicConfig(filename="%s.log" % name,level=logging.INFO)
    logger = logging.getLogger("parameter_fitting")
    logger.info(" Starting: Solving_Newtons_Method")
    print "  Solving for solutions with %s method" % fitopts["solver"]
    cwd = os.getcwd()
    os.chdir("%s/iteration_%d/newton" % (name,iteration))
    solver.find_solutions(model)
    os.chdir(cwd)
    logger.info(" Finished: Solving_Newtons_Method")

def save_new_parameters(model,fitopts):
    """ Save new parameters """
    available_methods = ["ddG_MC2004","FRET","RMSF","contact_Qi"]
    modules = {"ddG_MC2004":ddG_MC2004,"FRET":FRET,"RMSF":RMSF,"contact_Qi":contact_Qi}

    method = fitopts["data_type"]
    if method not in available_methods:
        raise ValueError("Requested fitting data %s not in %s" % (method,available_methods.__repr__()))

    submodule = modules[method]

    name = model.name
    iteration = model.iteration
    if not os.path.exists("%s/iteration_%d/newton/Lambda_index.txt" % (name,iteration)):
        raise IOError("%s/iteration_%d/newton/Lambda_index.txt  must exist!" % (name,iteration))
    else:
        soln_index = int(open("%s/iteration_%d/newton/Lambda_index.txt" % (name,iteration),"r").read().rstrip("\n"))

    cwd = os.getcwd()
    os.chdir("%s/iteration_%d/newton" % (name,iteration))
    submodule.save_new_parameters.save(model,soln_index)
    os.chdir(cwd)

if __name__ == "__main__":
    pass
