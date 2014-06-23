import numpy as np
import subprocess as sb
import os


import model_builder as mdb
import project_tools as pjt


pdb = "1SHG.pdb"
name = "1SHG"

base = mdb.models.SmogCalphaBase.SmogCalphaBase()
base.clean_pdb(pdb)
if not os.path.exists("Qref_shadow"):
    os.mkdir("Qref_shadow")

open("Qref_shadow/clean.pdb","w").write(base.cleanpdb_full)
open("Native.pdb","w").write(base.cleanpdb)

base.shadow_contacts(".")

base.dissect_native_pdb(".")
base.generate_topology()
base.generate_grofile()

base.get_interaction_table()

nsteps = "4000000"


cwd = os.getcwd()
for T in range(50,150,5):
    print "  Running ",T
    if not os.path.exists(str(T)):
        os.mkdir(str(T))

    os.chdir(str(T))

    mdp = pjt.simulation.mdp.get_constant_temperature_mdp_smog(str(T),nsteps)
    pbs_string = pjt.simulation.Tf_loop.get_pbs_string(name+str(T),"serial","1","08:00:00")

    open("Native.pdb","w").write(base.cleanpdb)
    open("clean.pdb","w").write(base.cleanpdb_full)
    open("clean_noH.pdb","w").write(base.cleanpdb_full_noH)
    open("1SHG.top","w").write(base.topology)
    open("1SHG.gro","w").write(base.grofile)
    open("nvt.mdp","w").write(mdp)
    open("run.pbs","w").write(pbs_string)
    np.savetxt("table.xvg",base.table,fmt="%20.12e",delimiter=" ")

    cmd1 = "grompp -f nvt.mdp -c %s.gro -p %s.top -o run.tpr" % (name,name)
    sb.call(cmd1,shell=True,stdout=open("prep.out","w"),stderr=open("prep.err","w"))

    cmd2 = "qsub run.pbs"
    sb.call(cmd2,shell=True)

    os.chdir(cwd)
