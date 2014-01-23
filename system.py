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
        self.itp_files = 0.
        self.pdbs = args.pdbs
        self.subdirs = [ pdb.split('.pdb')[0] for pdb in self.pdbs ] 
        self.numprots = len(self.subdirs)
        self.error = [ 0 for i in range(self.numprots) ]
        self.Tf_iteration = [ 0 for i in range(self.numprots) ]
        self.Tf_active_directory = [ 'Tf_0' for i in range(self.numprots) ]
        self.Tf_refinements = [ {0:0} for i in range(self.numprots) ]
        self.mutation_iteration = [ 0 for i in range(self.numprots) ]
        self.mutation_active_directory = [ '' for i in range(self.numprots) ]
        
    def __repr__(self,i):
        reprstring = "[ Main_Path ]\n"
        reprstring += "%s\n" % self.path
        reprstring += "[ PDB ]\n"
        reprstring += "%s\n" % self.pdbs[i]
        reprstring += "[ subdir ]\n"
        reprstring += "%s\n" % self.subdirs[i]
        reprstring += "[ Tf_iteration ]\n"
        reprstring += "%s\n" % self.Tf_iteration[i]
        reprstring += "[ Tf_refinements ]\n"
        for key in self.Tf_refinements[i].iterkeys():
            reprstring += "%d %d\n" % (key, self.Tf_refinements[i][key])
        reprstring += "[ Tf_active_directory ]\n"
        reprstring += "%s\n" % self.Tf_active_directory[i]
        reprstring += "[ mutation_iteration ]\n"
        reprstring += "%s\n" % self.mutation_iteration[i]
        reprstring += "[ mutation_active_directory ]\n"
        reprstring += "%s\n" % self.mutation_active_directory[i]
        return reprstring

    def load_info_file(self,i):
        ''' Load in the dynamic system info. File must have exact
            format. '''
        info_file = open(self.subdirs[i]+'/system.info','r')
        line = info_file.readline()
        while line != '':
            #print line[:-1] ## DEBUGGING
            #print line.split() ## DEBUGGING
            value = line.split()[1]
            if value == "Tf_iteration":
                self.Tf_iteration[i] = int(info_file.readline())
            elif value == "Tf_refinements":
                line = info_file.readline()
                while line[0] != "[":
                    key,val = line.split()
                    self.Tf_refinements[i][int(key)] = int(val)
                    line = info_file.readline()
                continue
            elif value == "Tf_active_directory":
                self.Tf_active_directory[i] = info_file.readline()[:-1]
            elif value == "mutation_iteration":
                self.mutation_iteration[i] = int(info_file.readline())
            elif value == "mutation_active_directory":
                self.mutation_active_directory[i] = info_file.readline()[:-1]
            else:
                if line[0] == "[":
                    line = info_file.readline()
            line = info_file.readline()

    def append_log(self,sub,string):
        logfile = open(self.path+'/'+sub+'/'+sub+'.log','a').write(self.append_time_label()+' '+string+'\n')

    def append_time_label(self):
        now = time.localtime()
        now_string = "%s:%s:%s:%s:%s" % (now.tm_year,now.tm_mon,now.tm_mday,now.tm_hour,now.tm_min)
        return now_string

    def clean_pdbs(self):
        ''' Clean the PDB files by writing only the ATOM lines up until the first
            TER or END. Doesn't remove any atoms.'''
        
        for i in range(len(self.pdbs)):
            sub = self.subdirs[i]
            cleanpdb = self.clean_pdb(self.pdbs[i])
            open(sub+'/clean.pdb','w').write(cleanpdb)
            open(sub+'/'+sub+'.pdb','w').write(cleanpdb)

    def clean_pdb(self,pdb):
        ''' Create Native.pdb in subdirectory. Currently only for CA. Expand later for
            CACB. PDB fixed-width column format is given by:
        ATOM     44  C   ALA A  11      12.266  21.667  20.517  1.00 28.80           C  
        '''
        atomid = 1
        first_flag = 0
        cleanpdb = ''
        for line in open(pdb,'r'):
            if line[:3] in ['TER','END']:
                break
            else:
                if line[:4] == 'ATOM':
                    if first_flag == 0:
                        if line[16] in ["A"," "]:
                            newline = 'ATOM%7s  %-4s%3s A%4d%s' % \
                                    (atomid,line[13:17],line[17:20],1,line[26:])
                            atomid += 1
                            first_flag = 1
                            first_index = int(line[22:26]) - 1
                            cleanpdb += newline
                    else:
                        if line[16] in ["A"," "]:
                            newline = 'ATOM%7s  %-4s%3s A%4d%s' % \
                                    (atomid,line[13:17],line[17:20],int(line[22:26])-first_index,line[26:])
                            atomid += 1
                            cleanpdb += newline
        
        #cleanpdb += 'TER'
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

    def get_native_contacts(self,heavy_atoms,atoms_per_res):
        ''' Calculate contact map based on comparing inter-residue heavy atom 
            distances. WORKS!'''

        N = len(heavy_atoms)
        D = np.zeros((N,N))
        Q = np.zeros((len(atoms_per_res),len(atoms_per_res)))
        
        for i in range(1,N):
            ## DEBUGGING
            #print heavy_atoms[0:-i,:] 
            #print heavy_atoms[i:,:]
            #print heavy_atoms[0:-i,:].shape
            #print heavy_atoms[i:,:].shape
    
            diff = heavy_atoms[0:-i,:] - heavy_atoms[i:,:]
            d = np.sqrt(diff[:,0]**2 + diff[:,1]**2 + diff[:,2]**2 )
            diag = (np.arange(0,N-i),np.arange(i,N))
            D[diag] = d
        C = (D < 5.5).astype(int)

        for i in range(len(atoms_per_res)):
            N = sum(atoms_per_res[:i+1])
            n = atoms_per_res[i]
            for j in range(i+4,len(atoms_per_res)):
                M = sum(atoms_per_res[:j+1])
                m = atoms_per_res[j]
                
                res_contact = C[N:N+n,M:M+m].any()
                Q[i,j] = int(res_contact)
                
        #np.savetxt("heavy_atom_dist.dat",D,delimiter=" ",fmt="%5.3f")       ## DEBUGGING
        #np.savetxt("heavy_atom_contacts.dat",C,delimiter=" ",fmt="%1d")       ## DEBUGGING
        #np.savetxt("Native_contact.dat",Q,delimiter=" ",fmt="%1d")       ## DEBUGGING

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
        prots_heavy_atoms = []
        prots_heavy_atoms_per_res = []
        for sub in self.subdirs:
            atomid = 1
            prev_resid = 1
            nativepdb = ''
            heavy_atoms = []
            num_heavy_atoms = []
            temp_num_atoms = 0
            for line in open(sub+"/clean.pdb","r"):
                if line[:3] == "TER":
                    nativepdb += line
                    break
                else:
                    #print "**%s**" % line[13:15] ## DEBUGGING

                    if line[13:15] == "CA":
                        #newline = "%s%5d%s%s%s%4d%s" % (line[:6],resid,line[11:16]," ",line[17:22],resid,line[26:])
                        newline = "%s%5d%s%s%s%4d%s" % (line[:6],atomid,line[11:16]," ",line[17:22],atomid,line[26:])
                        nativepdb += newline
                        atomid += 1

                    if line[77] != "H":
                        ## Collect heavy atoms for calculation of native contacts.
                        resid = int(line[22:26])
                        if resid != prev_resid:
                            num_heavy_atoms.append(temp_num_atoms)
                            temp_num_atoms = 1
                            prev_resid = resid
                        else:
                            temp_num_atoms += 1
                        
                        heavy_atoms.append(np.array([float(line[30:38]),float(line[38:46]),float(line[46:54])]))

            #print np.array(heavy_atoms[:15])        ## DEBUGGING
            #print num_heavy_atoms[:10]              ## DEBUGGING
            self.get_native_contacts(np.array(heavy_atoms),num_heavy_atoms)

            prots_heavy_atoms.append(heavy_atoms)
            prots_heavy_atoms_per_res.append(num_heavy_atoms)
            open(sub+"/Native.pdb","w").write(nativepdb)
            self.native_pdbs.append(nativepdb)


        
        #raise SystemExit        ## DEBUGGING
            

        

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

    def write_info_file(self,i):
        ''' Writes model.info file in subdirectory. The data of the Model object    
            is static, so this is straightforward documentation.'''
        open(self.subdirs[i]+"/system.info","w").write(self.__repr__(i))
