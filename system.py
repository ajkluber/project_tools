import os
import numpy as np
import time

'''
Purpose:
    The purpose of the System class is to hold all the data relevant to the 
particular protein(s) of interest. The Model class extracts this data to 
construct the input files which are returned to the System object.

To Do:
- Create the __repr__ method to print in a pretty way.


'''


'''
PDB File Format:
http://www.wwpdb.org/documentation/format33/sect9.html#ATOM

COLUMNS        DATA  TYPE    FIELD        DEFINITION
-------------------------------------------------------------------------------------
 1 -  6        Record name   "ATOM  "
 7 - 11        Integer       serial       Atom  serial number.
13 - 16        Atom          name         Atom name.
17             Character     altLoc       Alternate location indicator.
18 - 20        Residue name  resName      Residue name.
22             Character     chainID      Chain identifier.
23 - 26        Integer       resSeq       Residue sequence number.
27             AChar         iCode        Code for insertion of residues.
31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
55 - 60        Real(6.2)     occupancy    Occupancy.
61 - 66        Real(6.2)     tempFactor   Temperature  factor.
77 - 78        LString(2)    element      Element symbol, right-justified.
79 - 80        LString(2)    charge       Charge  on the atom.
'''

class System(object):
    def __init__(self,args):
        self.path = os.getcwd()
        self.systemname = args.name
        self.itp_files = 0.
        self.pdbs = args.pdbs
        self.subdirs = [ pdb.split('.pdb')[0] for pdb in self.pdbs ] 
        self.numprots = len(self.subdirs)
        self.Tf_iteration = [ 0 for i in range(self.numprots) ]
        self.Tf_active_directory = [ 'Tf_0' for i in range(self.numprots) ]
        self.Tf_refinements = [ {0:0} for i in range(self.numprots) ]
        self.mutation_iteration = [ 0 for i in range(self.numprots) ]
        self.mutation_active_directory = [ '' for i in range(self.numprots) ]
        
    def __repr__(self):
        reprstring = "[ System_Name ]\n"
        reprstring += "%s\n" % self.systemname
        reprstring += "[ Main_Path ]\n"
        reprstring += "%s\n" % self.path
        reprstring += "[ PDBs ]\n"
        reprstring += "%s\n" % self.subdirs.__repr__()
        reprstring += "[ Tf_iteration ]\n"
        reprstring += "%s\n" % self.Tf_iteration.__repr__()
        reprstring += "[ Tf_refinements ]\n"
        reprstring += "%s\n" % self.Tf_refinements.__repr__()
        reprstring += "[ Tf_active_directories ]\n"
        reprstring += "%s\n" % self.Tf_active_directories.__repr__()
        reprstring += "[ mutation_iteration ]\n"
        reprstring += "%s\n" % self.mutation_iteration.__repr__()
        reprstring += "[ mutation_active_directories ]\n"
        reprstring += "%s\n" % self.mutation_active_directories.__repr__()
        return reprstring

    def append_log(self,sub,string):
        logfile = open(self.path+'/'+sub,'a').write(self.append_time_label()+' '+string+'\n')

    def append_time_label(self):
        now = time.localtime()
        now_string = "%s:%s:%s:%s:%s" % (now.tm_year,now.tm_mon,now.tm_mday,now.tm_hour,now.tm_min)
        return now_string

    def clean_pdbs(self):
        ''' Copy.'''
        
        for i in range(len(self.pdbs)):
            sub = self.subdirs[i]
            cleanpdb = self.clean_pdb(self.pdbs[i])
            open(sub+'/clean.pdb','w').write(cleanpdb)


    def clean_pdb(self,pdb):
        ''' Create Native.pdb in subdirectory. Currently only for CA. Expand later for
            CACB. PDB fixed-width column format is given by:
        ATOM     44  C   ALA A  11      12.266  21.667  20.517  1.00 28.80           C  
        '''
        atomid = 0
        resid = 1
        cleanpdb = ''
        for line in open(pdb,'r'):
            if line[:3] in ['TER','END']:
                break
            else:
                #if line[:4] == 'ATOM' and line[13:15] == 'CA':
                if line[:4] == 'ATOM':
                    newline = 'ATOM%7s  %-4s%3s A%4d%s' % \
                            (atomid,line[13:17],line[17:20],resid,line[26:])
                    atomid += 1
                    resid += 1
                    cleanpdb += newline
        cleanpdb += 'TER'
        return cleanpdb

    def get_atom_indices(self,beadmodel):
        ''' Extract the atom indices for desired atoms. Also get the residue names
            and the native state coordinates.'''
        prots_indices = []
        prots_residues = []
        prots_coords = []
        for sub in self.subdirs:
            atomindices = dict((atomtype,[]) for atomtype in beadmodel )
            residues = []
            coords = []
            for line in open(sub+"/Native.pdb","r"):
                if line[13:15] == "CA":
                    residues.append(line[17:20])
                for atomtype in beadmodel:
                    if line[13:15] == atomtype:
                        atomindices[atomtype].append(int(line[6:13]))
                        #print "**%s**%s**%s" % (line[31:39], line[39:47], line[47:55]) ## DEBUGGING
                        coords.append([float(line[31:39]),float(line[39:47]),float(line[47:55])]) 
            prots_indices.append(atomindices)
            prots_residues.append(residues)
            prots_coords.append(np.array(coords))
        return prots_indices, prots_residues, prots_coords

    def write_Native_pdb(self,beadmodel):
        ''' Depending on the beadmodel that is input (from Model)
            write the corresponding Native.pdb '''
        if len(beadmodel) == 1 and beadmodel[0] == "CA":
            self.write_Native_pdb_CA()
        elif len(beadmodel) == 2:
            self.write_Native_pdb_CACB()

    def write_Native_pdb_CA(self):
        ''' Write the Native.pdb for CA topology.'''
        self.native_pdbs = []
        for sub in self.subdirs:
            resid = 1
            nativepdb = ''
            for line in open(sub+"/clean.pdb","r"):
                if line[:3] == "TER":
                    nativepdb += line
                    break
                else:
                    #print "**%s**" % line[13:15] ## DEBUGGING
                    if line[13:15] == "CA":
                        newline = "%s%5d%s%4d%s" % (line[:6],resid,line[11:22],resid,line[26:])
                        nativepdb += newline
                        resid += 1
            open(sub+"/Native.pdb","w").write(nativepdb)
            self.native_pdbs.append(nativepdb)

    def write_Native_pdb_CACB(self):
        ''' Write the Native.pdb for CACB topology. '''
        for sub in self.subdirs:
            atmid = 1
            resid = 1
            nativepdb = ''
            for line in open(sub+"/clean.pdb","r"):
                if line[:3] != "TER":
                    nativepdb += line
                    break
                else:
                    if line[13:15] == "CA":
                        #newline = "%6s%5d%11s%4d%s" % (line[:6],resid,line[11:23],resid,line[])
                        newline = "%s%5d%s%4d%s\n" % (line[:6],atmid,line[11:23],resid,line[26:])
                        atmid += 1
                        resid += 1
                    elif line[13:15] == "CB":
                        newline = "%s%5d%s%4d%s\n" % (line[:6],atmid,line[11:23],resid,line[26:])

            open(sub+"/Native.pdb","w").write(nativepdb)
        return nativepdb
