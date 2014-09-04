import matplotlib as mpl
mpl.use('Agg')
#import matplotlib.pylab as py
import matplotlib.pyplot as plt
import numpy as np
import os 
import argparse

#from project_tools.analysis.whamdata import WhamData

"""

    Author: Alexander Kluber
    Date: August 24 2013

    Code to replace Contour.m, for making contour plots of the free energy 
    data that comes from WHAM.

    Problem is that some matplotlib functionality doesn't work quite right.
    This could be because I compiled numpy with the Intel mkl libraries on
    installation. The error occurs when plt.colorbar or plt.contourf is used.
 
"""

global R
R = 0.00831415

def smooth_xyz(X,Y,Z):
    ''' WORKS '''
    winsize = 3
    shp = (Z.shape[0]-winsize,Z.shape[1]-winsize)
    newX = np.zeros(shp,float)
    newY = np.zeros(shp,float)
    newZ = np.zeros(shp,float)

    for i in range(shp[0]):
        for j in range(shp[1]):
            newX[i,j] = np.mean(X[i,j:j+winsize])
            newY[i,j] = np.mean(Y[i:i+winsize,j])
            #newZ[i,j] = np.mean(Z[i:i+winsize,j:j+winsize])
            #newZ[i,j] = np.mean(Z[i,j:j+winsize])
            newZ[i,j] = np.mean(Z[i:i+winsize,j])
            if j >= winsize:
                newZ[i,j-winsize] = np.mean(newZ[i,j-winsize:j])
        
    return newX,newY,newZ

def individual_FreeEnergy(path,name):
    ''' For use in a directory containing different scales.'''
    fileendings = ["rmsd_Rg","Q_Rg","Qh_Qnh","Q_Qh"]
    for ending in fileendings:
        plt.figure()
        data = WhamData(path=path,inputfilename="input_%s_F.wham" % ending)
        F = np.loadtxt(path+"/wham/FreeEnergy_%s.dat" % ending,usecols=(0,1,2))

        X = F[:,0].reshape((data.nbins_d1,data.nbins_d2)) 
        Y = F[:,1].reshape((data.nbins_d1,data.nbins_d2)) 
        Z = F[:,2].reshape((data.nbins_d1,data.nbins_d2)) 

        maxF = max(Z[Z != 999.])
        minF = min(Z[Z != 999.])
        Z[Z == 999.] = maxF + 0.1
        Z -= minF
        Z /= R*data.T_out_min
        #X,Y,Z = smooth_xyz(X,Y,Z)

        if ending in ["Qh_Qnh","Q_Qh","Q_Rg"]:
            if ending == "Q_Rg":
                X /= X.max()
            else:
                X /= X.max()
                Y /= Y.max()

        levels = np.arange(0,6,0.5)

        print "Plotting FreeEnergy_%s" % ending
        if ending == "Qh_Qnh":
            CS = plt.contourf(Y,X,Z,levels=levels)
            plt.xlabel("%s" % ending.split("_")[1])
            plt.ylabel("%s" % ending.split("_")[0])
            plt.xlim(0,1)
            plt.ylim(0,1)
        else:
            if ending in ["Q_Qh"]:
                CS = plt.contourf(X,Y,Z,levels=levels)
                plt.xlabel("%s" % ending.split("_")[0])
                plt.ylabel("%s" % ending.split("_")[1])
                plt.xlim(0,1)
                plt.ylim(0,1)
            else: 
                CS = plt.contourf(X,Y,Z,50)
                plt.xlabel("%s" % ending.split("_")[0])
                plt.ylabel("%s" % ending.split("_")[1])
        plt.colorbar()
        plt.title("Free Energy %s %s" % (ending,name))
        plt.savefig(path+"/wham/FreeEnergy_%s.pdf" % ending)

def individual_1DFreeEnergy(path,name):
    ''' For use in a directory containing different scales.'''
    
    fileendings = ["Q","Qh","Qnh","Rg","rmsd"]
    for ending in fileendings:
        plt.figure()
        data = WhamData(path=path,inputfilename="input_%s_1D_F.wham" % ending)
        F = np.loadtxt(path+"/wham/FreeEnergy_%s_1D.dat" % ending,usecols=(0,1,2))

        X = F[:,0].reshape((data.nbins_d1,data.nbins_d2)) 
        Y = F[:,1].reshape((data.nbins_d1,data.nbins_d2)) 
        Z = F[:,2].reshape((data.nbins_d1,data.nbins_d2)) 

        maxF = max(Z[Z != 999.])
        minF = min(Z[Z != 999.])
        #print "Max F ", maxF
        Z[Z == 999.] = maxF + 0.1
        Z -= minF
        Z /= R*data.T_out_min
        #Z /= R*data.T_out_min

        if ending in ["Q","Qh","Qnh"]:
            X /= X.max()
        
        #print X.shape
        #print Z.shape
        #raise SystemExit

        #X,Y,Z = smooth_xyz(X,Y,Z)
        plt.plot(X[:,0],Z)
        #plt.set_ylim(0,1.25)
        if ending in ["Q","Qh","Qnh"]:
            plt.xlim(0,1)
            plt.ylim(0,6)

        print "Plotting FreeEnergy_%s" % ending
        plt.xlabel("%s" % ending)
        plt.ylabel("Free Energy / $k_{B}T$")
        plt.title("Free Energy %s %s" % (ending,name))
        plt.savefig(path+"/wham/FreeEnergy_%s.pdf" % ending)

def individual_HeatCap(path,name):
    ''' For use in a directory containing different scales.'''
    
    plt.figure()
    T,Cv = np.loadtxt(path+"/wham/Heat_rmsd_Rg.dat",usecols=(0,2),unpack=True)

    print "Plotting Heat Capacity"
    plt.plot(T,Cv)
    plt.xlabel("Temperature (K)")
    plt.ylabel("Heat Capacity kJ/mol")
    plt.title("Heat Capcity %s" % name)
    plt.savefig(path+"/wham/Cv.pdf")

def subplot_free_energy():
    ''' For use in a directory containing different scales.'''
    #Tf = [[371,354],[332,305],[290,266]]
    #fileendings = ["rmsd_Rg","Q_Rg","Qh_Qnh","Q_Qh"]
    fileendings = ["Q_Rg","Qh_Qnh","Q_Qh"]
    #ending=fileendings[3]
    scales = [[1,0.9],[0.8,0.7],[0.6,0.5]]

    R = 0.008314461

    for m in range(len(fileendings)):
        ending = fileendings[m]
        print "Plotting subplots for %s " % ending
        plt.figure()
        fig, axes = plt.subplots(3,2,sharex=True,sharey=True)
        for i in range(3):
            for j in range(2):

                scale = str(scales[i][j])
                Tf = int(open(scale+"/Tf.txt").read())
                data = WhamData(path=str(scale)+"/",inputfilename="input_%s_F.wham" % ending)
                F = np.loadtxt(scale+"/wham/FreeEnergy_%s.dat" % ending,usecols=(0,1,2))

                X = F[:,0].reshape((data.nbins_d1,data.nbins_d2)) 
                Y = F[:,1].reshape((data.nbins_d1,data.nbins_d2)) 
                Z = F[:,2].reshape((data.nbins_d1,data.nbins_d2)) 

                maxF = max(Z[Z != 999.])
                minF = min(Z[Z != 999.])
                print "Max F ", maxF
                Z[Z == 999.] = maxF + 0.1
                Z -= minF
                
                # Normalize free energy to units of k_b*T_f
                Z /= R*data.T_out_min
                
                X,Y,Z = smooth_xyz(X,Y,Z)
                CS = axes[i,j].contour(X,Y,Z,50)
                #axes[i,j].plot(X[0,0],Y[0][0],'k.',label=scale+" $T_f=$"+str(Tf[i][j]))
                #axes[i,j].legend(loc=1)
                #axes[i,j].title("$\\epsilon = %s T_f = %d $" % (scale,Tf[i][j]), size="xx-small")
                #axes[i,j].text(0.5,0.75,"$\\epsilon=$"+scale+" $T_f=$"+str(Tf[i][j]),bbox=dict(facecolor='white',alpha=0.5),transform=axes[i,j].transAxes)
                axes[i,j].text(0.1,0.85,"$\\epsilon=$"+scale+" $T_f=$"+str(Tf),bbox=dict(facecolor='white',alpha=0.5),transform=axes[i,j].transAxes)
                #if i == 0:
                #    axes[i,j].set_title("$\\epsilon = $"+scale+" $T_f = $"+str(Tf[i][j]), size="small")
                #else:
                #    axes[i,j].set_xlabel("$\\epsilon = $"+scale+" $T_f = $"+str(Tf[i][j]), size="small")


        plt.subplots_adjust(wspace=0,hspace=0)
        fig.text(0.5, 0.04, '%s' % ending.split("_")[0], ha='center', va='center')
        fig.text(0.06, 0.5, '%s' % ending.split("_")[1], ha='center', va='center', rotation='vertical')
        #plt.ylabel("Radius of Gyration (Rg)")
        #plt.xlabel("Native Contacts (Q)")
        #plt.title("Free Energy for scale="+str(scale))
        plt.suptitle("Free Energy at $T_f$ for non-bonded scaling",fontsize=12)
        plt.savefig("FreeEnergy_subplot_%s.pdf" % ending)

def subplot_free_energy_1D_from_2D_data():
    ''' For use in a directory containing different scales.'''
    #Tf = [[371,354],[332,305],[290,266]]
    #fileendings = ["rmsd_Rg","Q_Rg","Qh_Qnh","Q_Qh"]
    fileendings = ["Q_Rg","Qh_Qnh","Q_Qh"]
    scales = [[1,0.9],[0.8,0.7],[0.6,0.5]]
    R = 0.008314461

    #for m in range(len(fileendings)):
    ending = fileendings[0]
    print "Plotting subplots for %s " % ending
    plt.figure()
    fig, axes = plt.subplots(3,2,sharex=True,sharey=True)
    for i in range(3):
        for j in range(2):

            scale = str(scales[i][j])
            print "Scale = ",scale
            Tf = int(open(scale+"/Tf.txt").read())
            data = WhamData(path=str(scale)+"/",inputfilename="input_%s_F.wham" % ending)
            F = np.loadtxt(scale+"/wham/FreeEnergy_%s.dat" % ending,usecols=(0,1,2))

            X = F[:,0].reshape((data.nbins_d1,data.nbins_d2)) 
            Y = F[:,1].reshape((data.nbins_d1,data.nbins_d2)) 
            Z = F[:,2].reshape((data.nbins_d1,data.nbins_d2)) 

            maxF = max(Z[Z != 999.])
            minF = min(Z[Z != 999.])
            print "Max F ", maxF
            Z[Z == 999.] = maxF + 0.1
            Z -= minF
            
            # Normalize free energy to units of k_b*T_f
            Z /= R*data.T_out_min
            
            temp = np.exp(-Z)
            ZQ = -np.log(np.array([ sum(temp[n,:]) for n in range(len(Z)) ]))
            
            ZQ -= min(ZQ)
            
            #X,Y,Z = smooth_xyz(X,Y,Z)
            axes[i,j].plot(X[:,0],ZQ)
            axes[i,j].set_ylim(0,1.25)
            #CS = axes[i,j].contour(X,Y,Z,50)
            #axes[i,j].plot(X[0,0],Y[0][0],'k.',label=scale+" $T_f=$"+str(Tf[i][j]))
            #axes[i,j].legend(loc=1)
            #axes[i,j].title("$\\epsilon = %s T_f = %d $" % (scale,Tf[i][j]), size="xx-small")
            #axes[i,j].text(0.5,0.75,"$\\epsilon=$"+scale+" $T_f=$"+str(Tf[i][j]),bbox=dict(facecolor='white',alpha=0.5),transform=axes[i,j].transAxes)
            axes[i,j].text(0.1,0.85,"$\\epsilon=$"+scale+" $T_f=$"+str(Tf),bbox=dict(facecolor='white',alpha=0.5),transform=axes[i,j].transAxes)
            #if i == 0:
            #    axes[i,j].set_title("$\\epsilon = $"+scale+" $T_f = $"+str(Tf[i][j]), size="small")
            #else:
            #    axes[i,j].set_xlabel("$\\epsilon = $"+scale+" $T_f = $"+str(Tf[i][j]), size="small")


    print "saving"
    plt.subplots_adjust(wspace=0,hspace=0)
    fig.text(0.5, 0.04, '%s' % ending.split("_")[0], ha='center', va='center')
    fig.text(0.06, 0.5, 'Free Energy (kT)', ha='center', va='center', rotation='vertical')
    #fig.text(0.06, 0.5, '%s' % ending.split("_")[1], ha='center', va='center', rotation='vertical')
    #plt.ylabel("Radius of Gyration (Rg)")
    #plt.xlabel("Native Contacts (Q)")
    #plt.title("Free Energy for scale="+str(scale))
    plt.suptitle("Free Energy at $T_f$",fontsize=12)
    plt.savefig("FreeEnergy_1D_subplot_%s.pdf" % ending.split("_")[0])

def subplot_free_energy_1D():
    ''' For use in a directory containing different scales.'''
    #Tf = [[371,354],[332,305],[290,266]]
    #fileendings = ["rmsd_Rg","Q_Rg","Qh_Qnh","Q_Qh"]
    fileendings = ["Q","Qh","Qnh","Rg","rmsd"]
    scales = [[1,0.9],[0.8,0.7],[0.6,0.5]]
    R = 0.008314461

    for m in range(len(fileendings)):
        ending = fileendings[m]
        print "Plotting subplots for %s " % ending
        plt.figure()
        fig, axes = plt.subplots(3,2,sharex=True,sharey=True)
        for i in range(3):
            for j in range(2):

                scale = str(scales[i][j])
                print "     Scale = ",scale
                Tf = int(open(scale+"/Tf.txt").read())
                data = WhamData(path=str(scale)+"/",inputfilename="input_%s_1D_F.wham" % ending)
                F = np.loadtxt(scale+"/wham/FreeEnergy_%s_1D.dat" % ending,usecols=(0,1,2))

                X = F[:,0].reshape((data.nbins_d1,data.nbins_d2)) 
                Y = F[:,1].reshape((data.nbins_d1,data.nbins_d2)) 
                Z = F[:,2].reshape((data.nbins_d1,data.nbins_d2)) 

                maxF = max(Z[Z != 999.])
                minF = min(Z[Z != 999.])
                #print "Max F ", maxF
                Z[Z == 999.] = maxF + 0.1
                Z -= minF
                Z /= R*data.T_out_min
                
                #print X.shape
                #print Z.shape
                #raise SystemExit

                #X,Y,Z = smooth_xyz(X,Y,Z)
                axes[i,j].plot(X[:,0],Z)
                axes[i,j].set_ylim(0,1.25)
                #CS = axes[i,j].contour(X,Y,Z,50)
                #axes[i,j].plot(X[0,0],Y[0][0],'k.',label=scale+" $T_f=$"+str(Tf[i][j]))
                #axes[i,j].legend(loc=1)
                #axes[i,j].title("$\\epsilon = %s T_f = %d $" % (scale,Tf[i][j]), size="xx-small")
                #axes[i,j].text(0.5,0.75,"$\\epsilon=$"+scale+" $T_f=$"+str(Tf[i][j]),bbox=dict(facecolor='white',alpha=0.5),transform=axes[i,j].transAxes)
                axes[i,j].text(0.1,0.85,"$\\epsilon=$"+scale+" $T_f=$"+str(Tf),bbox=dict(facecolor='white',alpha=0.5),transform=axes[i,j].transAxes)
                #if i == 0:
                #    axes[i,j].set_title("$\\epsilon = $"+scale+" $T_f = $"+str(Tf[i][j]), size="small")
                #else:
                #    axes[i,j].set_xlabel("$\\epsilon = $"+scale+" $T_f = $"+str(Tf[i][j]), size="small")


        print "saving as: FreeEnergy_1D_subplot_%s.pdf" % ending
        plt.subplots_adjust(wspace=0,hspace=0)
        fig.text(0.5, 0.04, '%s' % ending, ha='center', va='center')
        fig.text(0.06, 0.5, 'Free Energy (kT)', ha='center', va='center', rotation='vertical')
        #fig.text(0.06, 0.5, '%s' % ending.split("_")[1], ha='center', va='center', rotation='vertical')
        #plt.ylabel("Radius of Gyration (Rg)")
        #plt.xlabel("Native Contacts (Q)")
        #plt.title("Free Energy for scale="+str(scale))
        plt.suptitle("Free Energy at $T_f$",fontsize=12)
        plt.savefig("FreeEnergy_1D_subplot_%s.pdf" % ending)

def subplot_free_energy_R151617():
    ''' For use in a directory containing different scales.'''
    #Tf = [[371,354],[332,305],[290,266]]
    fig, axes = plt.subplots(3,2,sharex=True,sharey=True)
    fileendings = ["rmsd_Rg","Q_Rg","Qh_Qnh","Q_Qh"]
    ending=fileendings[2]
    scales = [[1,0.9],[0.8,0.7],[0.6,0.5]]
    proteins = ["r15","r16","r17"]
    cwd = os.getcwd()
    for prot in proteins:
        print "Plotting all FreeEnergies for: ", prot
        os.chdir(cwd + "/" + prot)
        subplot_free_energy()

def main():
    parser = argparse.ArgumentParser(description='Plotting utility.')
    sp = parser.add_subparsers(dest='action')

    single_parser = sp.add_parser('single')
    single_parser.add_argument('--name', type=str, required=True, help='Name of system.')
    single_parser.add_argument('--path', type=str, required=True, help='Path that holds the wham directory.')
    single_parser.add_argument('--FreeEnergy1D', action='store_true', help='.')
    single_parser.add_argument('--FreeEnergy', action='store_true', help='.')
    single_parser.add_argument('--HeatCap', action='store_true', help='.')
    

    subplot_parser = sp.add_parser('subplot')
    subplot_parser.add_argument('--name', type=str, required=True, help='Name of system.')
    subplot_parser.add_argument('--path', type=str, required=True, help='Path that holds the wham directory.')

    args = parser.parse_args()
    
    #print args
    #raise SystemExit
    #if args.action == 'single':
    #    if args.FreeEnergy == True:
    #        individual_FreeEnergy(args.path,args.name)
    #    elif args.FreeEnergy1D == True:
    #        individual_1DFreeEnergy(args.path,args.name)
    #    elif args.HeatCap == True:
    #        individual_HeatCap(args.path,args.name)
    individual_FreeEnergy(args.path,args.name)
    individual_1DFreeEnergy(args.path,args.name)
    individual_HeatCap(args.path,args.name)


if __name__ == '__main__':
   
    main()

    #individual_plots()
    #subplot_free_energy()
    #subplot_free_energy_1D()
    #subplot_free_energy()
