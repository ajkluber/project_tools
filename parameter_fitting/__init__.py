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

import newton
import ddG_MC2004
import FRET
import RMSF
import contact_Qi

def prepare_newtons_method(model,method,append_log):
    """ Prepare the files to do newtons method """


    available_methods = ["ddG_MC2004","FRET","RMSF","contact_Qi"]
    modules = {"ddG_MC2004":ddG_MC2004,"FRET":FRET,"RMSF":RMSF,"contact_Qi":contact_Qi}

    if method not in available_methods:
        print "ERROR! requested method not found %s" % method
        print " Choose from available_methods: ",available_methods
        print " Exiting."
        raise SystemExit

    submodule = modules[method]

    name = model.subdir
    iteration = model.Mut_iteration

    append_log(name,"Starting: Calculating_Jacobian")

    target_feature, target_feature_err = submodule.compute_Jacobian.get_target_feature(model)
    sim_feature_avg, sim_feature_err, Jacobian_avg, Jacobian_err = submodule.compute_Jacobian.calculate_average_Jacobian(model)

    if not os.path.exists("%s/Mut_%d/newton" % (name,iteration)):
        os.mkdir("%s/Mut_%d/newton" % (name,iteration))
    if not os.path.exists("%s/Mut_%d/%s" % (name,iteration,method)):
        os.mkdir("%s/Mut_%d/%s" % (name,iteration,method))

    print "  Saving feature vector and Jacobian in %s/Mut_%d/%s" % (name,iteration,method)
    np.savetxt("%s/Mut_%d/%s/target_feature.dat" % (name,iteration,method), target_feature)
    np.savetxt("%s/Mut_%d/%s/target_feature_err.dat" % (name,iteration,method), target_feature_err)
    np.savetxt("%s/Mut_%d/%s/sim_feature.dat" % (name,iteration,method), sim_feature_avg)
    np.savetxt("%s/Mut_%d/%s/sim_feature_err.dat" % (name,iteration,method), sim_feature_err)
    np.savetxt("%s/Mut_%d/%s/Jacobian.dat" % (name,iteration,method), Jacobian_avg)
    np.savetxt("%s/Mut_%d/%s/Jacobian_err.dat" % (name,iteration,method) ,Jacobian_err)


    ## To Do:
    ##  - Collect Jacobian rows from all fitting_includes directories.
    ##  - Map columns (parameters) to match those of the first directory. Stack the rows.
    ##  - Save in the first fitting directory.

    print "  Saving feature vector and Jacobian in %s/Mut_%d/newton" % (name,iteration)
    np.savetxt("%s/Mut_%d/newton/target_feature.dat" % (name,iteration), target_feature)
    np.savetxt("%s/Mut_%d/newton/target_feature_err.dat" % (name,iteration), target_feature_err)
    np.savetxt("%s/Mut_%d/newton/sim_feature.dat" % (name,iteration), sim_feature_avg)
    np.savetxt("%s/Mut_%d/newton/sim_feature_err.dat" % (name,iteration), sim_feature_err)
    np.savetxt("%s/Mut_%d/newton/Jacobian.dat" % (name,iteration), Jacobian_avg)
    np.savetxt("%s/Mut_%d/newton/Jacobian_err.dat" % (name,iteration) ,Jacobian_err)

    append_log(name,"Finished: Calculating_Jacobian")
    append_log(name,"Starting: Solving_Newtons_Method")


    ## Find solutions with Levenbeg_Marquardt algorithm
    print "  Solving for solutions with Levenberg-Marquardt method"
    cwd = os.getcwd()
    os.chdir("%s/Mut_%d/newton" % (name,iteration))
    newton.solver.Levenberg_Marquardt_solution(model,method)
    os.chdir(cwd)



if __name__ == "__main__":
    ## Module tested. Works!
    import model_builder as mdb
    def dummy(this,that):
        pass

    name = "1RIS"
    method = "ddG_MC2004"

    model = mdb.models.load_model(name) 
    model.Mut_iteration = 0

    prepare_newtons_method(model,method,dummy)
