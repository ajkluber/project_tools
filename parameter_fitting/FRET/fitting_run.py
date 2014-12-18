"""Script for running a specialized fitting sequence"""

import compute_Jacobian as compJ
import save_new_paramter as snp
import truncated_SVD_FRET as fitting
import numpy as np

def prepare_newtons_method(model,method,append_log)
    name = model.subdir
    iteration = model.iteration
    if not os.path.exists("%s/iteration_%d/newton/Jacobian.dat" % (name,iteration)):
        append_log(name,"Starting: Calculating_Jacobian")

        sim_feature_avg, sim_feature_err, Jacobian_avg, Jacobian_err = compJ.calculate_average_Jacobian(model)
        target_feature, target_feature_err = compJ.get_target_feature(model)
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
        
        print "  Saving feature vector and Jacobian in %s/iteration_%d/newton" % (name,iteration)
        np.savetxt("%s/iteration_%d/newton/target_feature.dat" % (name,iteration), target_feature)
        np.savetxt("%s/iteration_%d/newton/target_feature_err.dat" % (name,iteration), target_feature_err)
        np.savetxt("%s/iteration_%d/newton/sim_feature.dat" % (name,iteration), sim_feature_avg)
        np.savetxt("%s/iteration_%d/newton/sim_feature_err.dat" % (name,iteration), sim_feature_err)
        np.savetxt("%s/iteration_%d/newton/Jacobian.dat" % (name,iteration), Jacobian_avg)
        np.savetxt("%s/iteration_%d/newton/Jacobian_err.dat" % (name,iteration) ,Jacobian_err)
        
        append_log(name,"Finished: Calculating_Jacobian")
        solve_newtons_method(model,method,append_log)
    
def solve_newtons_method(model,method,append_log):
    name = model.subdir
    iteration = model.iteration
    append_log(name,"Starting: Solving_Newtons_Method")
    cwd = os.getcwd()
    os.chdir("%s/iteration_%d/newton" % (name,iteration))
    fitting.find_solutions(model,method)
    os.chdir(cwd)
    append_log(name,"Finished: Solving_Newtons_Method")

def save_new_parameters(model,method,append_log):
    name = model.subdir
    iteration = model.iteration
    if not os.path.exists("%s/iteration_%d/newton/solution.dat" % (name,iteration)):
        print "ERRROR, MISSING SOLUTION"
    cwd = os.getcwd()
    os.chdir("%s/iteration_%d/newton" % (name,iteration))
    snp.savej(model)

    os.chdir(cwd)
    
    
    
if __name__ == "__main__":
    asdfasdf
