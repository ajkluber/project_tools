import numpy as np
import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import os

import model_builder.analysis.contacts as ct


'''
Created 2-12-2014
Alexander Kluber

    A module under development that plots the 1D, 2D pmfs of equilibrium data
as well as the contact maps as a function of several coordinates. Calculating
contact maps over equilibrium trajectories is aided by the
model_builder.analysis.contacts submodule, but takes really long so it should
be done from within a PBS script.

To Do:
- Create command line functionality that can be called from a PBS script.
- Write a function that submits PBS scripts to plot quantities.

'''

def get_equil_data(name,T,mutnum=0,numsub=4):
    ''' Accumulate data from equilibrium temperature directories.'''
    path = name +"/Mut_"+str(mutnum)
    for i in range(1,numsub):
        q = np.loadtxt(path+"/"+T+"_"+str(i)+"/Qprob.dat",dtype=float)
        qh = np.loadtxt(path+"/"+T+"_"+str(i)+"/Qhprob.dat",dtype=float)
        qnh = np.loadtxt(path+"/"+T+"_"+str(i)+"/Qnhprob.dat",dtype=float)
        nh = np.loadtxt(path+"/"+T+"_"+str(i)+"/Nh.dat",dtype=float)
        if i == 1:
            Q = q
            Qh = qh
            Qnh = qnh
            Nh = nh
        else:
            Q = np.hstack([Q,q])  
            Qh = np.hstack([Qh,qh])  
            Qnh = np.hstack([Qnh,qnh])  
            Nh = np.hstack([Nh,nh])
    Q /= float(max(Q))
    Qh /= float(max(Qh))
    Qnh /= float(max(Qnh))
    Nh /= float(max(Nh))
    return Q,Qh,Qnh,Nh

def plot_2D_equil_pmfs(name,T):
    ''' Plot 2D pmfs.'''

    path = name +"/Mut_0/"+T+"_pmfs"
    if os.path.exists(path) == False:
        os.mkdir(path)
    #Tf = open(name+"/Tf_0/Tf.txt","r").read()[:-1]
    #Tf = "%.2f" % float(Tf)
    Q,Qh,Qnh,Nh = get_equil_data(name,T)

    qvals = np.unique(list(Q)+[1.0,])
    qhvals = np.unique(list(Qh)+[1.0,])
    qnhvals = np.unique(list(Qnh)+[1.0,])
    nhvals = np.unique(list(Nh)+[1.0,])

    q_qh_hist,xedges,yedges = np.histogram2d(Qh,Q,bins=[qhvals,qvals],normed=True)
    q_qh_hist[q_qh_hist == 0.0] == 1.0
    pmf = -np.log(q_qh_hist)
    pmf -= min(pmf.ravel())

    levels = np.arange(0,6,0.5)
    X,Y = np.meshgrid(yedges[:-1],xedges[:-1])
    plt.figure()
    plt.subplot(1,1,1,aspect=1)
    plt.contourf(X,Y,pmf,levels=levels)
    #plt.contour(X,Y,pmf)
    plt.xlim(0,1)
    plt.ylim(0,1)
    plt.xlabel("$Q$")
    plt.ylabel("$Q_h$")
    plt.title("$F(Q,Q_h)$ for "+name+" T="+T)
    cbar = plt.colorbar()
    cbar.set_label("$F / kT$")
    #plt.savefig(path+"/Q_Qh_"+T+"_pmf.pdf")
    plt.savefig(path+"/Q_Qh_pmf.pdf")

    qnh_qh_hist,xedges,yedges = np.histogram2d(Qh,Qnh,bins=[qhvals,qnhvals],normed=True)
    qnh_qh_hist[qnh_qh_hist == 0.0] == 1.0
    pmf = -np.log(qnh_qh_hist)
    pmf -= min(pmf.ravel())

    levels = np.arange(0,6,0.5)
    X,Y = np.meshgrid(yedges[:-1],xedges[:-1])
    plt.figure()
    plt.subplot(1,1,1,aspect=1)
    plt.contourf(X,Y,pmf,levels=levels)
    #plt.contour(X,Y,pmf)
    plt.xlim(0,1)
    plt.ylim(0,1)
    plt.xlabel("$Q_{nh}$")
    plt.ylabel("$Q_h$")
    plt.title("$F(Q_{nh},Q_h)$ for "+name+" T="+T)
    cbar = plt.colorbar()
    cbar.set_label("$F / kT$")
    #plt.savefig(path+"/Qnh_Qh_"+T+"_pmf.pdf")
    plt.savefig(path+"/Qnh_Qh_pmf.pdf")

    q_qnh_hist,xedges,yedges = np.histogram2d(Qnh,Q,bins=[qhvals,qnhvals],normed=True)
    q_qnh_hist[q_qnh_hist == 0.0] == 1.0
    pmf = -np.log(q_qnh_hist)
    pmf -= min(pmf.ravel())

    levels = np.arange(0,6,0.5)
    X,Y = np.meshgrid(yedges[:-1],xedges[:-1])
    plt.figure()
    plt.subplot(1,1,1,aspect=1)
    plt.contourf(X,Y,pmf,levels=levels)
    #plt.contour(X,Y,pmf)
    plt.xlim(0,1)
    plt.ylim(0,1)
    plt.xlabel("$Q$")
    plt.ylabel("$Q_{nh}$")
    plt.title("$F(Q,Q_{nh})$ for "+name)
    cbar = plt.colorbar()
    cbar.set_label("$F / kT$")
    #plt.savefig(path+"/Q_Qnh_"+T+"_pmf.pdf")
    plt.savefig(path+"/Q_Qnh_pmf.pdf")

    qnh_nh_hist,xedges,yedges = np.histogram2d(Nh,Qnh,bins=[nhvals,qnhvals],normed=True)
    qnh_nh_hist[qnh_nh_hist == 0.0] == 1.0
    pmf = -np.log(qnh_nh_hist)
    pmf -= min(pmf.ravel())

    levels = np.arange(0,6,0.5)
    X,Y = np.meshgrid(yedges[:-1],xedges[:-1])
    plt.figure()
    plt.subplot(1,1,1,aspect=1)
    plt.contourf(X,Y,pmf,levels=levels)
    #plt.contour(X,Y,pmf)
    plt.xlim(0,1)
    plt.ylim(0,1)
    plt.xlabel("$Q_{nh}$")
    plt.ylabel("$N_{h}$")
    plt.title("$F(Q_{nh},N_{h})$ for "+name)
    cbar = plt.colorbar()
    cbar.set_label("$F / kT$")
    #plt.savefig(path+"/Qnh_Nh_"+T+"_pmf.pdf")
    plt.savefig(path+"/Qnh_Nh_pmf.pdf")

def plot_1D_equil_pmfs(name,T):
    ''' Plot 1D pmfs.'''

    path = name+"/Mut_0/"+T+"_pmfs"
    if os.path.exists(path) == False:
        os.mkdir(path)
    #Tf = open(name+"/Tf_0/Tf.txt","r").read()[:-1]
    #Tf = "%.2f" % float(Tf)
    Q,Qh,Qnh,Nh = get_equil_data(name,T)

    qvals = np.unique(list(Q)+[1.0,])
    qhvals = np.unique(list(Qh)+[1.0,])
    qnhvals = np.unique(list(Qnh)+[1.0,])
    nhvals = np.unique(list(Nh)+[1.0,])

    nq,binsq = np.histogram(Q,bins=qvals,density=True)
    nqh,binsqh = np.histogram(Qh,bins=qhvals,density=True)
    nqnh,binsqnh = np.histogram(Qnh,bins=qnhvals,density=True)
    nnh,binsnh = np.histogram(Nh,bins=nhvals,density=True)
    #nq,binsq = np.histogram(Q,bins=50,density=True)
    #nqh,binsqh = np.histogram(Qh,bins=10,density=True)
    #nqnh,binsqnh = np.histogram(Qnh,bins=35,density=True)

    pmfq = -np.log(nq)
    pmfqh = -np.log(nqh)
    pmfqnh = -np.log(nqnh)
    pmfnh = -np.log(nnh)

    pmfq -= min(pmfq)
    pmfqh -= min(pmfqh)
    pmfqnh -= min(pmfqnh)
    pmfnh -= min(pmfnh)

    plt.figure()
    plt.plot(binsq[1:]/max(binsq),pmfq)
    plt.xlabel("$Q$")
    plt.ylabel("$F(Q) / kT$")
    plt.title("$F(Q)$ for "+name+" $T_f=%s$" % T)
    plt.ylim(0,6)
    plt.xlim(0,1)
    #plt.savefig(path+"/Q_pmf_"+name+".pdf")
    #plt.savefig(path+"/Q_pmf_"+T+".pdf")
    plt.savefig(path+"/Q_pmf.pdf")

    plt.figure()
    plt.plot(binsqh[1:]/max(binsqh),pmfqh)
    plt.xlabel("$Q_h$")
    plt.ylabel("$F(Q_h) / kT$")
    plt.title("$F(Q_h)$ for "+name+" $T_f=%s$" % T)
    plt.ylim(0,6)
    plt.xlim(0,1)
    #plt.savefig(path+"/Qh_pmf_"+name+".pdf")
    #plt.savefig(path+"/Qh_pmf_"+T+".pdf")
    plt.savefig(path+"/Qh_pmf.pdf")

    plt.figure()
    plt.plot(binsqnh[1:]/max(binsqnh),pmfqnh)
    plt.xlabel("$Q_{nh}$")
    plt.ylabel("$F(Q_{nh}) / kT$")
    plt.title("$F(Q_{nh})$ for "+name+" $T_f=%s$" % T)
    plt.ylim(0,6)
    plt.xlim(0,1)
    #plt.savefig(path+"/Qnh_pmf_"+name+".pdf")
    #plt.savefig(path+"/Qnh_pmf_"+T+".pdf")
    plt.savefig(path+"/Qnh_pmf.pdf")

    plt.figure()
    plt.plot(binsnh[1:]/max(binsnh),pmfnh)
    plt.xlabel("$N_{h}$")
    plt.ylabel("$F(N_{h}) / kT$")
    plt.title("$F(N_{h})$ for "+name+" $T_f=%s$" % T)
    plt.ylim(0,6)
    plt.xlim(0,1)
    #plt.savefig(path+"/Qnh_pmf_"+name+".pdf")
    #plt.savefig(path+"/Nh_pmf_"+T+".pdf")
    plt.savefig(path+"/Nh_pmf.pdf")

def plot_contact_maps(name,T,numtemps):
    ''' Plot contact maps as a function of several reaction
        coordinates. '''

    coords = ["Q","Qh","Qnh","Nh"]
    Q,Qh,Qnh,Nh = get_equil_data(name,T)
    data = [Q,Qh,Qnh,Nh]

    cwd = os.getcwd()
    os.chdir(name+"/Mut_0")
    path = T+"_maps"
    if os.path.exists(path) == False:
        print "Creating subdirectory ",path
        os.mkdir(path)
    for i in range(len(data)):
        coord = coords[i]
        datum = data[i]
        n, bins = np.histogram(datum,bins=10)
        print "Temperature:   ",T
        print "Bin occupancy: ",n
        print "Bin edges:     ",bins

        framestate = np.zeros(len(datum),int)
        tempstates = [ -1*np.ones(len(datum),int) for j in range(len(bins)-1) ]
        numstates = len(bins)-1  
        ## Assign each frame to a bin. Works.
        for n in range(len(datum)):
            for i in range(numstates):
                if (bins[i] <= datum[n] < datum[i+1]):
                    framestate[n] = i
                    tempstates[i][n] = i

        print framestate
        if list(framestate).count(-1) != 0:
            print "Framestate is not right. Exiting"
            os.chdir(cwd)
            raise SystemExit
        savedir = path+"/"+coord+"_maps"
        
        if os.path.exists(savedir) == False:
            os.mkdir(savedir)
        np.savetxt(savedir+"/counts.dat",bins)
        np.savetxt(savedir+"/bins.dat",bins)

        ct.equil_contacts_for_states(framestate,numstates,savedir,T,numtemps)
    os.chdir(cwd)

def plot_equil_data(T,name,walltime="00:20:00"):
    ''' PBS script to call plot_equil.py to plot 1D and 2D pmfs as well as 
        contact maps as a function of multiple coordinates.'''
    ##   **** NOT DONE YET ***
    maps_pbs = "#!/bin/bash\n"
    maps_pbs +="#PBS -N maps_"+name+"\n"
    maps_pbs +="#PBS -q serial\n"
    maps_pbs +="#PBS -l nodes=1:ppn=1,walltime=%s\n" % walltime
    maps_pbs +="#PBS -j oe\n"
    maps_pbs +="#PBS -e maps_%s.err\n" % T
    maps_pbs +="#PBS -o maps_%s.out\n" % T
    maps_pbs +="#PBS -V\n\n"
    maps_pbs +="cd $PBS_O_WORKDIR\n"
    maps_pbs +='python -m model_builder.plot equil --T %s --subdir %s\n' % (T,name)

    open("maps.pbs","w").write(energy_pbs)
    qsub = "qsub maps.pbs"
    sb.call(qsub.split(),stdout=open("maps.out","w"),stderr=open("pbs.err","w"))




#name = "sh3/1FMK"
#name = "s6/1RIS"
name = "r17"
Tf = open(name+"/Tf_0/Tf.txt","r").read()[:-1]
Tf = "%.2f" % float(Tf)
for i in range(5):
    T = float(Tf) + float(Tf)*(0.003*i)
    T = "%.2f" % T
    print T
    plot_contact_maps(name,T,5)
    #plot_1D_equil_pmfs(name,T)
    #plot_2D_equil_pmfs(name,T)

