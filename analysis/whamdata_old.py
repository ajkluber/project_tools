import numpy as np

'''
    Author: Alexander Kluber
    Created: August 20, 2013
    Purpose: Have easy and flexible Python code for running WHAM.

    Description:
        This class just serves as a (not very elegant) container for the
    data needed for WHAM. Its primary function is parsing the input.wham
    files that are used by Cecilia's fortran WHAM code. If the WHAM
    computation was done natively in Python, then the process would 
    probably be a lot simpler. For now this will do.
        This isn't that readable, but thats because its mainly just
    formating.

'''

class WhamData(object):
    def __init__(self,path=".",inputfilename="input.wham"):
        self.path = path
        self.inputfilename = inputfilename
        self.itermax = 1500
        self.tol = 0.001
        self.tolz = 0.
        self.error_tol = 0.
        self.Tscale = 0.0083144621
        self.get_inputs()

    def get_inputs(self):
        inputfile = open(self.path+"/wham/"+self.inputfilename,"r").readlines()
        n_input_temp = int(inputfile[0])

        self.temp = [ int(x) for x in inputfile[1:2*n_input_temp:2] ]
        self.name_input_file = \
            [ x[:-1] for x in inputfile[2:2*n_input_temp:2] ]
        self.ncol_d1, self.ncol_d2, self.ncol_ene =  \
            [ int(x) for x in inputfile[2*n_input_temp + 1].split() ]

        self.d1_min, self.d1_step, self.nbins_d1 =   \
            [ float(x) for x in inputfile[2*n_input_temp + 2].split() ]

        self.d2_min, self.d2_step, self.nbins_d2 =  \
            [ float(x) for x in inputfile[2*n_input_temp + 3].split() ]

        self.e_min, self.e_step, self.nbins_e = \
            [ float(x) for x in inputfile[2*n_input_temp + 4].split() ]

        self.nbins_d1 = int(self.nbins_d1)
        self.nbins_d2 = int(self.nbins_d2)
        self.nbins_e = int(self.nbins_e)

        self.e_max = self.e_min + self.e_step*self.nbins_e

        self.name_output_file1 = inputfile[2*n_input_temp + 5][:-1]
        self.name_output_file2 = inputfile[2*n_input_temp + 6][:-1]

        self.n_temp_out, self.T_out_min, self.T_step = \
            [ float(x) for x in inputfile[2*n_input_temp + 7].split() ]

        self.T_out_max = self.T_out_min + self.T_step*self.n_temp_out

        flag = int(inputfile[2*n_input_temp + 8])
        if flag == 0:
            zold = np.zeros(n_input_temp,float)
            znew = np.zeros(n_input_temp,float)
        else:
            print "Using guess for F(T). NOT USED. EXITING"
            raise SystemExit

        self.nfix = int(inputfile[2*n_input_temp + 9])
        self.e_shift = int(inputfile[2*n_input_temp + 10])
        
        self.dens = np.zeros((n_input_temp,self.nbins_d1,self.nbins_d2,self.nbins_e),int)
        self.n_input_temp = n_input_temp

        ############### DEBUGGING #####################
        #print self.temp, self.name_input_file
        #print self.ncol_d1, self.ncol_d2, self.ncol_ene
        #print self.d2_min, self.d2_step, self.nbins_d2
        #print self.e_min, self.e_step, self.nbins_e
        #print self.name_output_file1
        #print self.name_output_file2
        #raise SystemExit

    
        self.temp = np.array(self.temp)*self.Tscale

    def read_temp_data(self):


        self.nconf = np.zeros(self.n_input_temp,int)
        self.ncount = np.zeros((self.nbins_d1,self.nbins_d2),int)
        for idx in range(self.n_input_temp):

            xdum = np.loadtxt("Input_for_WHAM_"+str(self.temp[idx])+".dat")
            d1 = xdum[:,self.ncol_d1 - 1]       
            d2 = xdum[:,self.ncol_d2 - 1]       
            ene = xdum[:,self.ncol_ene - 1]/4.184

            ########### DEBUGGING #####################
            #print xdum.shape
            #print self.ncol_d1 - 1
            #print self.ncol_d2 - 1
            #print self.ncol_ene - 1
            #raise SystemExit
    
            print "reading in temperature ", self.temp[idx]
            for n in range(len(d1)):
                n_d1 = int((d1[n] - self.d1_min)/self.d1_step)
                n_d2 = int((d2[n] - self.d2_min)/self.d2_step)
                n_e = int((ene[n] - self.e_min)/self.e_step)

                try:
                    self.dens[idx,n_d1,n_d2,n_e] += 1
                    self.nconf[idx] += 1
                    self.ncount[n_d1,n_d2] += 1
                except IndexError:
                    print "warning: value out of range "
        
        print self.nconf
        print self.ncount
        #print self.dens  ## DEBUGGING
