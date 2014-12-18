"""Module for finding the solution to a Jacobian fitting.

Made especially for the FRET fitting package"""

import numpy as np
import matplotlib.pyplot as plt
import os

def find_solutions(model,method):
    jump = 1000 #factor by which a difference is considered a jump down. truncate here
    
    target_feature = np.loadtxt("target_feature.dat")
    target_feature_err = np.loadtxt("target_feature_err.dat")
    sim_feature = np.loadtxt("sim_feature.dat")
    sim_feature_err = np.loadtxt("sim_feature_err.dat")
    Jacobian = np.loadtxt("Jacobian.dat")
    Jacobian_err = np.loadtxt("Jacobian_err.dat")
    J = Jacobian
    df = target_feature - sim_feature
    u,s,v = np.linalg.svd(J)
    umat = np.matrix(u)
    vmat = np.matrix(v)
    df = np.matrix(df)
    eps = model.model_param_values
    
    truncate = False
    nums = np.shape(s)[0]
    count = 1
    trun_indx = nums
    print "Current singular values are", s
    while (not truncate) and count < nums:
        if s[count-1]/s[count] > jump:
            truncate = True
            trun_indx = count 
        count += 1
    
    smat = np.zeros((np.shape(umat)[1],np.shape(vmat)[0]))
    smat = np.matrix(smat)
    
    for i in range(nums):
        if i < trun_indx:
            smat[i,i] = s[i]**(-1)
    umat = np.matrix(umat)
    vmat = np.matrix(vmat)
    df = np.matrix(df)
    
    #Transpose to solve
    
    umat = np.transpose(umat)
    smat = np.transpose(smat)
    vmat = np.transpose(vmat)
    df = np.transpose(df)
    
    print np.shape(df)
    print np.shape(umat)
    print np.shape(smat)
    print np.shape(vmat)
    print np.shape(umat*df)
    
    deps = vmat*(smat*(umat*df))
    os.mkdir("mytest")
    cwd = os.getcwd()
    os.chdir("mytest")
    save_solutions(deps,s,trun_indx)
    os.chdir(cwd)
    
def save_solutions(deps, s, trun_indx):
    np.savetxt("solution.dat", deps)
    np.savetxt("singular_values.dat", s)
    np.savetxt("truncation_index.dat", trun_indx)
    plt.plot(np.arange(np.shape(s)[0]), s, markerstyle="o", linewidth=2)
    plt.title("Singular Values", fontsize=20)
    plt.set_yscale("log")
    plt.xlabel("index", fontsize=20)
    plt.ylabel("Value", fontsize=20)
    
