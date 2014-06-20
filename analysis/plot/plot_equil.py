""" Plots data from 'equilibrium' simulations (i.e. really long).

**** Probably soon to be DEPRECATED 

Description:

    Collects and plots data from really long (hopefully equilibrium) simulations 
around the folding temperature. Plots 1D and 2D pmfs (1).


References:
(1) http://en.wikipedia.org/wiki/Potential_of_mean_force
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import os
import argparse
import subprocess as sb

import project_tools.analysis.contacts as ct

def main():
    ''' Two possible branches: 1. Calculate reference matrix, 2. Calculate Q '''
    parser = argparse.ArgumentParser(description='Calculate the (Non)Native contact matrix')
    sp = parser.add_subparsers(dest='action')

    map_parser = sp.add_parser('maps')
    map_parser.add_argument('--T', type=str, required=True, help='Temperature')
    map_parser.add_argument('--num', type=int, required=True, help='Number of temperature')
    map_parser.add_argument('--subdir', type=str, required=True, help='Directory holding Mut_0')

    TS_parser = sp.add_parser('TS')
    TS_parser.add_argument('--T', type=str, required=True, help='Temperature')
    TS_parser.add_argument('--num', type=int, required=True, help='Number of temperature')
    TS_parser.add_argument('--subdir', type=str, required=True, help='Directory holding Mut_0')
    TS_parser.add_argument('--Qbrackets', type=float, nargs="+", required=True, help='Values of Q that bracket the TS region.')
    args = parser.parse_args()
    
    if args.action == 'maps':
        plot_contact_maps(args.subdir,args.T,args.num)
    elif args.action == "TS":
        plot_TS_contact_map(args.subdir,args.T,args.num,args.Qbrackets)

def get_equil_data(T,path=".",numsub=3):
    ''' Accumulate data from equilibrium temperature directories.'''
    #path = name +"/Mut_"+str(mutnum)
    #path = "."
    for i in range(1,numsub+1):
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

def get_data(coord):
    if coord in ["Q","Qh","Qnh"]:
        x = np.loadtxt(coord+"prob.dat")
        x /= max(x)
    elif coord == "Nh":
        x = np.loadtxt(coord+".dat") 
        x /= max(x)
    elif coord == "rmsd":
        dummy, x = np.loadtxt(coord+".xvg",unpack=True)
    elif coord == "Rg":
        dummy, x = np.loadtxt("radius_cropped.xvg",unpack=True)
    else:
        print "ERROR!"
        print "  Coordinate: ",coord," not found"
    return x

def plot_2D_pmf(coord1,coord2,T,):

    path = os.getcwd()
    savedir = path+"/pmfs"
    if os.path.exists(savedir) == False:
        os.mkdir(savedir)

    x = get_data(coord1)
    y = get_data(coord2)

    xvals = np.unique(list(x))[::4]
    yvals = np.unique(list(y))[::4]

    hist,xedges,yedges = np.histogram2d(x,y,bins=[xvals,yvals],normed=True)
    hist[hist == 0.0] == 1.0
    pmf = -np.log(hist)
    pmf -= min(pmf.ravel())
    X,Y = np.meshgrid(yedges[:-1],xedges[:-1])

    levels = np.arange(0,6,0.5)

    plt.figure()
    plt.subplot(1,1,1,aspect=1)
    plt.contourf(X,Y,pmf,levels=levels)
    if (coord1.startswith("Q")) or (coord2.startswith("Q")) or (coord1 == "Nh") or (coord2 == "Nh"):
        plt.xlim(0,1)
        plt.ylim(0,1)
    else: 
        plt.xlim(0,max(x))
        plt.ylim(0,max(y))
    plt.xlabel(coord1,fontsize="xx-large")
    plt.ylabel(coord2,fontsize="xx-large")
    plt.title("F("+coord1+","+coord2+")",fontsize="xx-large")
    cbar = plt.colorbar()
    cbar.set_label("F / kT",fontsize="xx-large")
    #plt.savefig(path+"/Q_Qh_"+T+"_pmf.pdf")
    plt.savefig(savedir+"/"+coord1+"_"+coord2+"_pmf.pdf")

def plot_2D_equil_pmfs(name,T,num):
    ''' Plot 2D pmfs.'''

    path = name +"/Mut_0/"+T+"_pmfs"
    if os.path.exists(path) == False:
        os.mkdir(path)
    #Tf = open(name+"/Tf_0/Tf.txt","r").read()[:-1]
    #Tf = "%.2f" % float(Tf)
    Q,Qh,Qnh,Nh = get_equil_data(T,path=name+"/Mut_0",numsub=num)

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
    plt.xlabel("$Q$",fontsize="xx-large")
    plt.ylabel("$Q_h$",fontsize="xx-large")
    plt.title("$F(Q,Q_h)$ for "+name+" T="+T,fontsize="xx-large")
    cbar = plt.colorbar()
    cbar.set_label("$F / kT$",fontsize="xx-large")
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
    plt.xlabel("$Q_{nh}$",fontsize="xx-large")
    plt.ylabel("$Q_h$",fontsize="xx-large")
    plt.title("$F(Q_{nh},Q_h)$ for "+name+" T="+T,fontsize="xx-large")
    cbar = plt.colorbar()
    cbar.set_label("$F / kT$",fontsize="xx-large")
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
    plt.xlabel("$Q$",fontsize="xx-large")
    plt.ylabel("$Q_{nh}$",fontsize="xx-large")
    plt.title("$F(Q,Q_{nh})$ for "+name,fontsize="xx-large")
    cbar = plt.colorbar()
    cbar.set_label("$F / kT$",fontsize="xx-large")
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
    plt.xlabel("$Q_{nh}$",fontsize="xx-large")
    plt.ylabel("$N_{h}$",fontsize="xx-large")
    plt.title("$F(Q_{nh},N_{h})$ for "+name,fontsize="xx-large")
    cbar = plt.colorbar()
    cbar.set_label("$F / kT$",fontsize="xx-large")
    #plt.savefig(path+"/Qnh_Nh_"+T+"_pmf.pdf")
    plt.savefig(path+"/Qnh_Nh_pmf.pdf")

def plot_1D_equil_pmfs(name,T,num):
    ''' Plot 1D pmfs.'''

    path = name+"/Mut_0/"+T+"_pmfs"
    if os.path.exists(path) == False:
        os.mkdir(path)
    #Tf = open(name+"/Tf_0/Tf.txt","r").read()[:-1]
    #Tf = "%.2f" % float(Tf)
    Q,Qh,Qnh,Nh = get_equil_data(T,path=name+"/Mut_0",numsub=num)

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

    np.savetxt(path+"/Q_n.dat",nq,delimiter=" ",fmt="%.4f")
    np.savetxt(path+"/Q_bins.dat",binsq,delimiter=" ",fmt="%.4f")
    np.savetxt(path+"/Qh_n.dat",nqh,delimiter=" ",fmt="%.4f")
    np.savetxt(path+"/Qh_bins.dat",binsqh,delimiter=" ",fmt="%.4f")
    np.savetxt(path+"/Qnh_n.dat",nqnh,delimiter=" ",fmt="%.4f")
    np.savetxt(path+"/Qnh_bins.dat",binsqnh,delimiter=" ",fmt="%.4f")
    np.savetxt(path+"/Nh_n.dat",nnh,delimiter=" ",fmt="%.4f")
    np.savetxt(path+"/Nh_bins.dat",binsnh,delimiter=" ",fmt="%.4f")

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
    plt.xlabel("$Q$",fontsize="xx-large")
    plt.ylabel("$F(Q) / kT$",fontsize="xx-large")
    plt.title("$F(Q)$ for "+name+" $T_f=%s$" % T,fontsize="xx-large")
    plt.ylim(0,6)
    plt.xlim(0,1)
    #plt.savefig(path+"/Q_pmf_"+name+".pdf")
    #plt.savefig(path+"/Q_pmf_"+T+".pdf")
    plt.savefig(path+"/Q_pmf.pdf")

    plt.figure()
    plt.plot(binsqh[1:]/max(binsqh),pmfqh)
    plt.xlabel("$Q_h$",fontsize="xx-large")
    plt.ylabel("$F(Q_h) / kT$",fontsize="xx-large")
    plt.title("$F(Q_h)$ for "+name+" $T_f=%s$" % T,fontsize="xx-large")
    plt.ylim(0,6)
    plt.xlim(0,1)
    #plt.savefig(path+"/Qh_pmf_"+name+".pdf")
    #plt.savefig(path+"/Qh_pmf_"+T+".pdf")
    plt.savefig(path+"/Qh_pmf.pdf")

    plt.figure()
    plt.plot(binsqnh[1:]/max(binsqnh),pmfqnh)
    plt.xlabel("$Q_{nh}$",fontsize="xx-large")
    plt.ylabel("$F(Q_{nh}) / kT$",fontsize="xx-large")
    plt.title("$F(Q_{nh})$ for "+name+" $T_f=%s$" % T,fontsize="xx-large")
    plt.ylim(0,6)
    plt.xlim(0,1)
    #plt.savefig(path+"/Qnh_pmf_"+name+".pdf")
    #plt.savefig(path+"/Qnh_pmf_"+T+".pdf")
    plt.savefig(path+"/Qnh_pmf.pdf")

    plt.figure()
    plt.plot(binsnh[1:]/max(binsnh),pmfnh)
    plt.xlabel("$N_{h}$",fontsize="xx-large")
    plt.ylabel("$F(N_{h}) / kT$",fontsize="xx-large")
    plt.title("$F(N_{h})$ for "+name+" $T_f=%s$" % T,fontsize="xx-large")
    plt.ylim(0,6)
    plt.xlim(0,1)
    #plt.savefig(path+"/Qnh_pmf_"+name+".pdf")
    #plt.savefig(path+"/Nh_pmf_"+T+".pdf")
    plt.savefig(path+"/Nh_pmf.pdf")

def plot_TS_contact_map(name,T,numtemps,brackets):
    ''' Plot contact maps as a function of several reaction coordinates. '''

    coords = ["Q","Qh","Qnh","Nh"]
    Q,Qh,Qnh,Nh = get_equil_data(T,numsub=numtemps)
    data = [Q,Qh,Qnh,Nh]
    bins = np.array(brackets)

    cwd = os.getcwd()
    #os.chdir(name+"/Mut_0")
    path = T+"_maps"
    if os.path.exists(path) == False:
        print "Creating subdirectory ",path
        os.mkdir(path)
    ## Testing first on just Q.
    #for i in range(len(data)):
    for i in [0]:
        coord = coords[i]
        datum = data[i]
        #n, bins = np.histogram(datum,bins=10)
        print "Calculating transition state contact map"
        print "Temperature:   ",T
        print "Bin edges:     ",bins

        framestate = -1*np.ones(len(datum),int)
        tempstates = [ -1*np.ones(len(datum),int) for j in range(len(bins)-1) ]
        numstates = len(bins)-1  
        ## Assign each frame to a bin. Works.
        for n in range(len(datum)):
            if (bins[0] <= datum[n] <= bins[1]):
                framestate[n] = 0
                tempstates[i][n] = 0

        print framestate

        savedir = path+"/"+coord+"_TS_map"
        if os.path.exists(savedir) == False:
            os.mkdir(savedir)
        np.savetxt(savedir+"/counts.dat",np.array([list(framestate).count(1)]))
        np.savetxt(savedir+"/bins.dat",bins)

        ct.equil_contacts_for_states(framestate,numstates,savedir,T,numtemps)
    os.chdir(cwd)

def plot_contact_maps(name,T,numtemps):
    ''' Plot contact maps as a function of several reaction coordinates. '''

    coords = ["Q","Qh","Qnh","Nh"]
    Q,Qh,Qnh,Nh = get_equil_data(T,numsub=numtemps)
    data = [Q,Qh,Qnh,Nh]

    cwd = os.getcwd()
    #os.chdir(name+"/Mut_0")
    path = T+"_maps"
    if os.path.exists(path) == False:
        print "Creating subdirectory ",path
        os.mkdir(path)
    ## Testing first on just Q.
    #for i in range(len(data)):
    for i in [0]:
        coord = coords[i]
        datum = data[i]
        counts, bins = np.histogram(datum,bins=10)
        print "Temperature:   ",T
        print "Bin occupancy: ",counts
        print "Bin edges:     ",bins

        framestate = np.zeros(len(datum),int)
        tempstates = [ -1*np.ones(len(datum),int) for j in range(len(bins)-1) ]
        numstates = len(bins)-1  
        ## Assign each frame to a bin. Works.
        for n in range(len(datum)):
            for k in range(numstates):
                if (bins[k] <= datum[n] < datum[k+1]):
                    framestate[n] = k
                    tempstates[k][n] = k

        print framestate
        if list(framestate).count(-1) != 0:
            print "Framestate is not right. Exiting"
            os.chdir(cwd)
            raise SystemExit
        savedir = path+"/"+coord+"_maps"
        
        if os.path.exists(savedir) == False:
            os.mkdir(savedir)
        np.savetxt(savedir+"/counts.dat",counts)
        np.savetxt(savedir+"/bins.dat",bins)

        ct.equil_contacts_for_states(framestate,numstates,savedir,T,numtemps)
    os.chdir(cwd)

def submit_contact_map_calculator(subdir,T,num,mutnum=0,walltime="00:40:00"):
    ''' PBS script to call plot_equil.py to plot contact maps as a function
        of multiple coordinates.'''
    ##   **** NOT DONE YET ***
    maps_pbs = "#!/bin/bash\n"
    maps_pbs +="#PBS -N maps_"+subdir+"\n"
    maps_pbs +="#PBS -q serial\n"
    maps_pbs +="#PBS -l nodes=1:ppn=1,walltime=%s\n" % walltime
    maps_pbs +="#PBS -j oe\n"
    maps_pbs +="#PBS -e maps_%s.err\n" % T
    maps_pbs +="#PBS -o maps_%s.out\n" % T
    maps_pbs +="#PBS -V\n\n"
    maps_pbs +="cd $PBS_O_WORKDIR\n"
    maps_pbs +='python -m model_builder.plot.plot_equil maps --subdir %s --T %s --num %d\n' % (subdir,T,num)

    open(subdir+"/Mut_0/maps.pbs","w").write(maps_pbs)
    cwd = os.getcwd()
    os.chdir(subdir+"/Mut_0")
    qsub = "qsub maps.pbs"
    sb.call(qsub.split(),stdout=open("maps.out","w"),stderr=open("maps.err","w"))
    os.chdir(cwd)

def submit_TS_contact_map_calculator(subdir,T,num,brackets,mutnum=0,walltime="00:40:00"):
    ''' PBS script to call plot_equil.py to plot contact maps as a function
        of multiple coordinates.'''
    ##   **** NOT DONE YET ***
    TS_pbs = "#!/bin/bash\n"
    TS_pbs +="#PBS -N TS_"+subdir+"\n"
    TS_pbs +="#PBS -q serial\n"
    TS_pbs +="#PBS -l nodes=1:ppn=1,walltime=%s\n" % walltime
    TS_pbs +="#PBS -j oe\n"
    TS_pbs +="#PBS -e TS_%s.err\n" % T
    TS_pbs +="#PBS -o TS_%s.out\n" % T
    TS_pbs +="#PBS -V\n\n"
    TS_pbs +="cd $PBS_O_WORKDIR\n"
    TS_pbs +='python -m model_builder.plot.plot_equil TS --subdir %s --T %s --num %d --Qbrackets %.4f %.4f\n' % (subdir,T,num,brackets[0],brackets[1])

    open(subdir+"/Mut_0/TS.pbs","w").write(TS_pbs)
    cwd = os.getcwd()
    os.chdir(subdir+"/Mut_0")
    qsub = "qsub TS.pbs"
    sb.call(qsub.split(),stdout=open("TS.out","w"),stderr=open("TS.err","w"))
    os.chdir(cwd)

def plot_equil_data(subdir,T,num,inc=0.003):
    #name = "sh3/1FMK"
    #name = "s6/1RIS"
    #name = "r17"
    #Tf = open(subdir+"/Tf_0/Tf.txt","r").read()[:-1]
    #Tf = "%.2f" % float(Tf)
    coord = "Nh"

    ## Plot equilbirium pmfs.
    print "Plotting equilibrium data for %s at temperature %s" % (subdir,T)
    print "  1D pmfs for Q,Qh,Qnh,Nh..."
    plot_1D_equil_pmfs(subdir,T,num)
    print "  2D pmfs: Qh vs Q; Qh vs Qnh; Nh vs Qnh..."
    plot_2D_equil_pmfs(subdir,T,num)

    ## Plot equilibrium contact maps as a function of reactin coords.
    #print "  contact maps as a function of Q,Qh,Qnh,Nh..."
    #submit_contact_map_calculator(subdir,T,num)

    #cwd = os.getcwd()
    #for i in range(1,num+1):
        #T = float(Tf) + float(Tf)*(inc*i)
        #T = "%.2f" % T


        #print "  contact maps as a function of Q,Qh,Qnh,Nh..."
        #os.chdir(subdir+"/Mut_0")
        #submit_contact_map_calculator(subdir,T,num,name)
        #os.chdir(cwd)

def plot_transition_states(subdir,T,num,brackets):
    ''' Just plots the states for defined range of Q.'''
    submit_TS_contact_map_calculator(subdir,T,num,brackets)

if __name__ == '__main__':
    main()
