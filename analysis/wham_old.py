""" Python wrapper for WHAM code
    Author: Alexander Kluber
    Created: Aug 20, 2013

    Purpose: Have easy and flexible Python code for running WHAM.
    Description: 
        Python wrapper for Cecilia's fortran90 wham code. Preps input files and
    submits PBS jobs.
        I started writing this because it was getting difficult to figure out
    an error I was getting using the fortran WHAM code. Although I found that
    my error was a trivial input error, this code could be useful to continue
    to develop.
        This script helps prep the input files for the fortran WHAM code that
    we use. Most of this code is just formating strings so that the WHAM code
    will read
    it correctly.
        Note that this module requires WhamData.py to function. WhamData is a
    class that contains all the WHAM variables. There is also a BASH script I
    use to call this module. Also note that the core calculations are done 
    with Cecilia's fortran WHAM code. This program doesn't do any calculations
    itself, without the WHAM fortran code ths program.

    Procedure:
        The procedure requires several executions to acheive a final result.
    There is first a preparation step followed by a job submission step. The
    procedure, outlined below, will first calculate the heat capacity (Cv) as a
    function of temperature. The folding temperature can be taken as the peak
    in he heat capacity, which is subsequently used to generate 1D and 2D free
    energy plots as a function of different reaction coordinates.

    1. First, use options: prep --output HeatCap.
    2. Then, use options:  calc --output HeatCap
        Wait for the heat capacity runs to finish. There will be a Heat.dat
        file written when it finishes. This needs to exist before the free
        energy is calculated.
    3. After Heat.dat is created:
        i. For 2D free energy data 
            ia. Use option:      prep --output FreeEnergy
            ib. Then use option: calc --output FreeEnergy
        ii. For 1D free energy data 
            ia. Use option:      prep --output 1DFreeEnergy
            ib. Then use option: calc --output 1DFreeEnergy

    Changelog:
    August 20, 2013     Created
    September 13, 2013  Updated and expanded
    - Put the class definition WhamData into a separate file.
    November 16, 2013   Refactoring
    - Breaking large code blocks into functions to increase clarity and
    readability.
    - Changed the arguments of the two primary functions prep_input_files and
    run_wham so that they can be imported and used by other programs. Added
    prep_input_files_command_line and run_wham_command_line for backwards 
    compatibility. 

"""

import numpy as np
import subprocess as sb
import argparse
import os

from whamdata import WhamData


def main():
    """ Two possible branches: 1. Prep WHAM inputs; 2. Submit WHAM PBS job."""
    parser = argparse.ArgumentParser(description='Perform WHAM.')
    sp = parser.add_subparsers(dest='action')

    prep_parser = sp.add_parser('prep')
    prep_parser.add_argument('--output', type=str, help='Method',required=True)
    prep_parser.add_argument('--Ti', type=int, help='Initial Temperature',required=True)
    prep_parser.add_argument('--dT', type=int, help='Temperature step',required=True)
    prep_parser.add_argument('--Tf', type=int, help='Final Temperature',required=True)
    prep_parser.add_argument('--Tguess', type=int, help='Guess of Folding Temp')
    prep_parser.add_argument('--path', type=str, default=".", help='Path that holds temp dirs.')

    run_parser = sp.add_parser('run')
    run_parser.add_argument('--output', type=str, help='Method of calculation',required=True)
    run_parser.add_argument('--Ti', type=int, help='Initial Temperature')
    run_parser.add_argument('--dT', type=int, help='Temperature step')
    run_parser.add_argument('--Tf', type=int, help='Final Temperature')
    run_parser.add_argument('--Tguess', type=int, help='Folding Temperature')
    run_parser.add_argument('--path', type=str, default=".", help='Path that holds the temp dirs')

    args = parser.parse_args()
    
    F = {'prep':prep_input_files_command_line, 'run':run_wham_command_line}[args.action]
    F(args)

def prep_input_files_command_line(args):
    """ This function collates the proper data files into the Input_for_WHAM
        files. Then prepares the input.wham file for the desired calculation.
        Requires     
    """
    temps = range(args.Ti,args.Tf+args.dT,args.dT)
    N = len(temps)
    inputs = ("rmsd","Rg","Q","Qh","Qnh","Epot")
    numinputs = len(inputs)
    if os.path.exists(path +"/maxs_mins_rmsd_Rg_Q_Qh_Qnh_Et.dat"):
        maxesmins = np.loadtxt(args.path +"/maxs_mins_rmsd_Rg_Q_Qh_Qnh_Et.dat")
        maxs = maxesmins[0,:] 
        mins = maxesmins[1,:] 
        calcflag = 0
    else:
        calcflag = 1  
        maxs = -1000000*np.ones((numinputs),float)
        mins =  1000000*np.ones((numinputs),float)

    firstblock = str(N) + "\n"
    for t in temps:
        ## Prep the first block of the WHAM input file. As well as the
        ## Input_for_WHAM_###.dat files
        if calcflag == 1:
            maxs, mins = prepare_Input_for_WHAM_for_temperature(t,args.path,maxs,mins,numinputs)
        firstblock += str(t) + "\n"
        firstblock += "Input_for_WHAM_"+str(t)+".dat"+"\n"

    if calcflag == 1:
        np.savetxt(args.path +"/maxs_mins_rmsd_Rg_Q_Qh_Qnh_Et.dat",np.array([maxs,mins]))
    ## The second block of the input.wham file depends on what calculation is 
    ## desired. The three choices are: HeatCap, FreeEnergy, 1DFreeEnergy.
    if args.output == "HeatCap":
        ## Write the input file for heat capactiy calculation.
        secondblock = second_block_string_Cv(inputs,maxs,mins,0,1,len(inputs)-1,temps,args.Tguess)
        open(path+"/wham/input_rmsd_Rg_Cv.wham","w").write(firstblock + secondblock)
    elif args.output == "FreeEnergy" or args.output == "1DFreeEnergy":
        ## Write the input files for free energy calculations.
        prepare_free_energy_inputs(args.path,args.output,inputs,maxs,mins,firstblock)
    else:
        print "Select a proper output option: --output = [ HeatCap | FreeEnergy ]"

def prep_input_files(Ti,Tf,dT,path,outputstyle):
    """ This function collates the proper data files into the Input_for_WHAM
        files. Then prepares the input.wham file for the desired calculation.
        Requires reaction coordinate data (with specific filenames) to exist 
        or will crash. 
    """
    temps = range(Ti,Tf+dT,dT)
    N = len(temps)
    inputs = ("rmsd","Rg","Q","Qh","Qnh","Epot")
    numinputs = len(inputs)
    if os.path.exists(path +"/maxs_mins_rmsd_Rg_Q_Qh_Qnh_Et.dat"):
        maxesmins = np.loadtxt(path +"/maxs_mins_rmsd_Rg_Q_Qh_Qnh_Et.dat")
        maxs = maxesmins[0,:] 
        mins = maxesmins[1,:] 
        calcflag = 0
    else:
        calcflag = 1  
        maxs = -1000000*np.ones((numinputs),float)
        mins =  1000000*np.ones((numinputs),float)

    firstblock = str(N) + "\n"
    for t in temps:
        ## Prep the first block of the WHAM input file. As well as the
        ## Input_for_WHAM_###.dat files
        print "  Prepping Input_for_WHAM_"+str(t)+".dat "
        if calcflag == 1:
            maxs, mins = prepare_Input_for_WHAM_for_temperature(t,path,maxs,mins,numinputs)
        firstblock += str(t) + "\n"
        firstblock += "Input_for_WHAM_"+str(t)+".dat"+"\n"

    if calcflag == 1:
        np.savetxt(path +"/maxs_mins_rmsd_Rg_Q_Qh_Qnh_Et.dat",np.array([maxs,mins]))
    ## The second block of the input.wham file depends on what calculation is 
    ## desired. The three choices are: HeatCap, FreeEnergy, 1DFreeEnergy.
    if outputstyle == "HeatCap":
        ## Write the input file for heat capactiy calculation.
        secondblock = second_block_string_Cv(inputs,maxs,mins,0,1,len(inputs)-1,temps,300)
        open(path+"/wham/input_rmsd_Rg_Cv.wham","w").write(firstblock + secondblock)
    elif outputstyle == "FreeEnergy" or outputstyle == "1DFreeEnergy":
        ## Write the input files for free energy calculations.
        prepare_free_energy_inputs(path,outputstyle,inputs,maxs,mins,firstblock)
    else:
        print "Select a proper output option: --output = [ HeatCap | FreeEnergy ]"

def prepare_Input_for_WHAM_for_temperature(t,path,maxs,mins,numinputs):
    """ Get reaction coordinates for temperature "t". Note that all these files
        must exist or the program will crash. Format is as follows:
        Input_for_WHAM_###.dat
        rmsd  Rg  Q  Qh  Qnh  Epot

        Where,
        rmsd = Root-mean-squared-deviation from native structure.
        Rg = Radius of gyration. 
        Q = Number of native contacts.
        Qh = Number of helical (i, i+4) contacts.
        Qnh = Number of non-helical contacts. Just Qnh = Q - Qh
        Epot = Total potential energy. 

    """ 
    time,Eb,Eang,Ed,Enonbond,Etot = np.loadtxt(path+"/"+str(t)+"_0/energyterms.xvg",unpack=True)
    time, rmsd = np.loadtxt(path+"/"+str(t)+"_0/rmsd.xvg",unpack=True)
    frames, rg = np.loadtxt(path+"/"+str(t)+"_0/radius_cropped.xvg",unpack=True)
    Q = np.loadtxt(path+"/"+str(t)+"_0/Qprob.dat")
    Qh = np.loadtxt(path+"/"+str(t)+"_0/Qhprob.dat")
    Qnh = Q - Qh

    numframes = len(Q)
    data = np.zeros((numframes,numinputs),float)
    data[:,0] = rmsd
    data[:,1] = rg
    data[:,2] = Q
    data[:,3] = Qh
    data[:,4] = Qnh
    data[:,5] = Etot
    np.savetxt(path+"/wham/Input_for_WHAM_"+str(t)+".dat",data,delimiter=" ",fmt="%10.5f")

    for j in range(numinputs):
        if max(data[1:,j]) > maxs[j]:
            maxs[j] = max(data[1:,j])
        if min(data[1:,0]) < mins[j]:
            mins[j] = min(data[1:,j])
    return maxs,mins


def prepare_free_energy_inputs(path,outputstyle,inputs,maxs,mins,firstblock):
    """ Note that this requires that the Heat.dat has already finished 
        running. It requires the the Heat capacity curve to find the 
        folding temperature. """
    T, E, Cv, F = np.loadtxt(path+"/wham/Heat_rmsd_Rg.dat",unpack=True)
    maxCv = max(Cv)
    Tf = int(T[list(Cv).index(maxCv)])
    Tf_rounded = int(Tf) 
    open(path+"/Tf.txt","w").write("%.2f" % Tf)
    open(path+"/Tf_rounded.txt","w").write("%d" % Tf_rounded)

    if outputstyle == "FreeEnergy":
        ## FreeEnergy option is for 2D free energy calculations.
        print "Writing input_*_F.wham files..."
        reaction_coord_pairs = ((0,1),(2,1),(3,4),(2,3))
        for idx1,idx2 in reaction_coord_pairs:
            print "Writing input file for c1 = %s  c2 = %s ..." % (inputs[idx1],inputs[idx2])
            secondblock = second_block_string_F(inputs,maxs,mins,idx1,idx2,len(inputs)-1,Tf)
            open(path+"/wham/input_%s_%s_F.wham" % 
                (inputs[idx1],inputs[idx2]) ,"w").write(firstblock + secondblock)
    else:
        ## 1DFreeEnergy option is for 1D free energy calculations.
        print "Writing input_*_1D_F.wham files..."
        reaction_coords = (0,1,2,3,4)
        for idx1 in reaction_coords:
            print "Writing input file for c1 = %s ..." % inputs[idx1]
            secondblock = second_block_string_F(inputs,maxs,mins,idx1,1,len(inputs)-1,Tf,nbins2=1,oneD=True)
            open(path+"/wham/input_%s_1D_F.wham" % inputs[idx1] ,"w").write(firstblock + secondblock)

def second_block_string_Cv(inputs,maxs,mins,idx1,idx2,Eidx,temps,Tguess,nbins1=150,nbins2=150,Enbins=150):
    """ Second block string
    Ex:
    2 3 4
    20 1 210
    1.3 0.042 100
    -435 8 156
    FreeEnergy.2d.dat
    Heat_371.dat
    1 371 1
    0
    23
    0
    """
    stp1 = (maxs[idx1] - mins[idx1])/nbins1
    stp2 = (maxs[idx2] - mins[idx2])/nbins2
    Estp = (maxs[Eidx] - mins[Eidx])/Enbins
    block = "%d %d %d\n" % (idx1+1,idx2+1,Eidx+1)
    block += "%f %f %d\n" % (mins[idx1],stp1,nbins1)
    block += "%f %f %d\n" % (mins[idx2],stp2,nbins2)
    block += "%f %f %d\n" % (mins[Eidx],Estp,Enbins)
    block += "FreeEnergy_%s_%s.dat\n" % (inputs[idx1],inputs[idx2])
    block += "Heat_%s_%s.dat\n" % (inputs[idx1],inputs[idx2])
    block += "%d %d %.1f\n" % (2*(temps[-1]-temps[0]),temps[0],0.5)
    block += "0\n"
    block += "%d\n" % Tguess ## Need guess for folding temperature. FINISH LATER
    block += "0\n"
    return block

def second_block_string_F(inputs,maxs,mins,idx1,idx2,Eidx,Tf,nbins1=150,nbins2=150,Enbins=150,oneD=False):
    """ Second block string 
    Ex:
    20 1 210
    1.3 0.042 100
    -435 8 156
    FreeEnergy.2d.dat
    Heat_371.dat
    1 371 1
    0
    23
    0
    """
    Qlist = ["Q","Qh","Qnh"] 
    #if (inputs[idx1] in Qlist) or (inputs[idx2] in Qlist):
    if idx1 == 2 or idx1 == 3:
        print "Q step size = 1"
        stp1 = 1
        nbins1 = maxs[idx1] - mins[idx1] + 5
        if idx2 == 4 or idx2 == 3:
            print "%s step size = 1" % inputs[idx2]
            stp2 = 1
            nbins2 = maxs[idx2] - mins[idx2] + 5
        else:
            stp2 = (maxs[idx2] - mins[idx2])/nbins2
    else:
        stp1 = (maxs[idx1] - mins[idx1])/nbins1
        stp2 = (maxs[idx2] - mins[idx2])/nbins2
    Estp = (maxs[Eidx] - mins[Eidx])/Enbins
    block = "%d %d %d\n" % (idx1+1,idx2+1,Eidx+1)
    block += "%f %f %d\n" % (mins[idx1],stp1,nbins1)
    block += "%f %f %d\n" % (mins[idx2],stp2,nbins2)
    block += "%f %f %d\n" % (mins[Eidx],Estp,Enbins)
    if oneD == True:
        block += "FreeEnergy_%s_1D.dat\n" % inputs[idx1]
    else:
        block += "FreeEnergy_%s_%s.dat\n" % (inputs[idx1],inputs[idx2])
    block += "Heat_%s_%s_dummy.dat\n" % (inputs[idx1],inputs[idx2])
    block += "1 %d 1\n" % Tf
    block += "0\n"
    block += "%d\n" % int(round(Tf/10.)*10) 
    block += "0\n"
    return block

def run_wham_command_line(args):
    """ Creates and submits PBS jobs to execute WHAM. This should be in the 
        directory that contains the "wham/" directory. All input and output
        files are written to the "wham/" directory. For the free energy 
        option, a default list of coordinates is used, called "inputs". """
    inputs = ("rmsd","Rg","Q","Qh","Qnh","Epot")
    #reaction_coord_pairs = (("rmsd","Rg"),("Q","Rg"),("Qh","Qnh"))
    #reaction_coords = ("rmsd","Rg","Q","Qh","Qnh")
    reaction_coords = (0,1,2,3,4)
    reaction_coord_pairs = ((0,1),(2,1),(3,4),(2,3))

    cwd = os.getcwd()
    #copy_WHAM_executable = "cp /projects/cecilia/ajk8/dmc_model/free_energy_analysis/WHAM " + cwd + "/wham/"
    #copy_WHAM_executable = "cp /projects/cecilia/ajk8/model_builder/analysis/WHAM " + cwd + "/wham/"
    projects = os.environ["PROJECTS"]
    copy_WHAM_executable = "cp "+projects+"/project_tools/analysis/WHAM " + cwd + "/wham/"
    sb.call(copy_WHAM_executable.split())

    if args.output == "HeatCap":
        submit_heat_capacity_job()
    elif args.output == "FreeEnergy":
        submit_free_energy_job(reaction_coord_pairs,inputs)
    elif args.output == "1DFreeEnergy":
        submit_free_energy_1D_job(reaction_coords,inputs)
    print "Finished..."

def run_wham(outputstyle):
    """ Creates and submits PBS jobs to execute WHAM. This should be in the 
        directory that contains the "wham/" directory. All input and output
        files are written to the "wham/" directory. For the free energy 
        option, a default list of coordinates is used, called "inputs". """
    inputs = ("rmsd","Rg","Q","Qh","Qnh","Epot")
    reaction_coords = (0,1,2,3,4)
    reaction_coord_pairs = ((0,1),(2,1),(3,4),(2,3))
    cwd = os.getcwd()
    #copy_WHAM_executable = "cp /projects/cecilia/ajk8/dmc_model/free_energy_analysis/WHAM " + cwd + "/wham/"
    #copy_WHAM_executable = "cp /projects/cecilia/ajk8/model_builder/analysis/WHAM " + cwd + "/wham/"
    ## Trying to remove references to ajk8 directory.
    projects = os.environ["PROJECTS"]
    copy_WHAM_executable = "cp "+projects+"/project_tools/analysis/WHAM " + cwd + "/wham/"
    sb.call(copy_WHAM_executable.split())
    if outputstyle == "HeatCap":
        submit_heat_capacity_job()
    elif outputstyle == "FreeEnergy":
        submit_free_energy_job(reaction_coord_pairs,inputs)
    elif outputstyle == "1DFreeEnergy":
        submit_free_energy_1D_job(reaction_coords,inputs)
    #print "Finished..."
    
def submit_heat_capacity_job():
    """ Writes and submits the PBS job for calculating the heat capacity as
        a function of temperature. The peak of the heat capacity can be 
        taken as the folding temperature. """
    cwd = os.getcwd()
    os.chdir(cwd + "/wham")
    print "Writing Cv PBS script..."
    wham_pbs = "#!/bin/bash \n"
    wham_pbs += "#PBS -N WHAM_Cv\n"
    wham_pbs += "#PBS -q serial\n"
    wham_pbs += "#PBS -l nodes=1:ppn=1,walltime=00:10:00\n"
    wham_pbs += "#PBS -j oe\n"
    wham_pbs += "#PBS -V\n"
    wham_pbs += "cd $PBS_O_WORKDIR\n"
    wham_pbs += "./WHAM<input_rmsd_Rg_Cv.wham &> wham_Cv.log\n"
    open("wham_Cv.pbs","w").write(wham_pbs)
    print "Submitting job..."
    qsub = "qsub wham_Cv.pbs"
    sb.call(qsub.split())
    print "WHAM job submitted... "
    os.chdir(cwd)

def submit_free_energy_job(reaction_coord_pairs,inputs):
    ''' Writes and submits multiple PBS jobs for calculating the free energy 
        as a function of different pairs of coordinates. '''
    cwd = os.getcwd()
    os.chdir(cwd + "/wham")
    for idx1,idx2 in reaction_coord_pairs:    
        wham_pbs = "#!/bin/bash \n"
        wham_pbs += "#PBS -N WHAM_%s_%s_F\n" % (inputs[idx1],inputs[idx2])
        wham_pbs += "#PBS -q serial\n"
        wham_pbs += "#PBS -l nodes=1:ppn=1,walltime=00:10:00\n"
        wham_pbs += "#PBS -j oe\n"
        wham_pbs += "#PBS -V\n"
        wham_pbs += "cd $PBS_O_WORKDIR\n"
        wham_pbs += "./WHAM<input_%s_%s_F.wham &> wham_%s_%s_F.log\n" % (inputs[idx1],inputs[idx2],inputs[idx1],inputs[idx2])
        open("wham_%s_%s_F.pbs" % (inputs[idx1],inputs[idx2]) ,"w").write(wham_pbs)
        print "Submitting job wham_%s_%s_F.pbs" % (inputs[idx1],inputs[idx2])
        qsub = "qsub wham_%s_%s_F.pbs" % (inputs[idx1],inputs[idx2])
        sb.call(qsub.split())
    os.chdir(cwd)

def submit_free_energy_1D_job(reaction_coords,inputs):
    ''' Writes and submits multiple PBS jobs for calculating 1D free energy 
        curves a function of different coordinates. '''
    cwd = os.getcwd()
    os.chdir(cwd + "/wham")
    for idx1 in reaction_coords:    
        wham_pbs = "#!/bin/bash \n"
        wham_pbs += "#PBS -N WHAM_%s_F\n" % inputs[idx1]
        wham_pbs += "#PBS -q serial\n"
        wham_pbs += "#PBS -l nodes=1:ppn=1,walltime=00:10:00\n"
        wham_pbs += "#PBS -j oe\n"
        wham_pbs += "#PBS -V\n"
        wham_pbs += "cd $PBS_O_WORKDIR\n"
        wham_pbs += "./WHAM<input_%s_1D_F.wham &> wham_%s_1D_F.log\n" % (inputs[idx1],inputs[idx1])
        open("wham_%s_1D_F.pbs" % inputs[idx1] ,"w").write(wham_pbs)
        print "Submitting job wham_%s_1D_F.pbs" % inputs[idx1]
        cmd3 = "qsub wham_%s_1D_F.pbs" % inputs[idx1]
        sb.call(cmd3.split())
    os.chdir(cwd)

if __name__ == '__main__':
    main()
