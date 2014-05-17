import os

from project_tools.manager import ProjectManager
from project_tools import simulation, analysis, mutations


class MatysiakClementi2004(ProjectManager):
    
    ''' A project manager to reproduce Matysiak Clementi 2004 algorithm. 

    Description:
        This project manager is intended to automate the generation of a
    heterogeneous Go-model via the thermodynamic perturbation algorithm
    put forth by Matysiak & Clementi (1).
        
    Reference:
    (1) Matysiak, S.; Clementi, C. Optimal Combination of Theory and Experiment for
    the Characterization of the Protein Folding Landscape of S6: How Far Can a
    Minimalist Model Go? J. Mol. Biol. 2004, 343, 235–248.
    '''


    def logical_flowchart_starting(self,System,Model,sub,task):
        if task == "Tf_loop_iteration":
            print "Checking if Tf_loop_iteration completed..."
            simulation.Tf_loop.check_completion(System,self.append_log)
            lasttime2,action2,task2 = self.check_modelbuilder_log(sub)
            if action2 == "Finished:":
                print "Finished Tf_loop_iteration..."
                print "Starting Tf_loop_analysis..."
                analysis.Tf_loop.analyze_temperature_array(System,self.append_log)
        elif task == "Tf_loop_analysis":
            print "Checking if Tf_loop_analysis completed..."
            analysis.Tf_loop.check_completion(System,self.append_log)
        elif task == "wham_Cv":
            print "Starting to check if wham_Cv completed..."
            analysis.Tf_loop.continue_wham(System,self.append_log)
        elif task == "wham_FreeEnergy":
            print "Starting Equil_Tf..."
            simulation.Tf_loop.run_equilibrium_simulations(Model,System,self.append_log)
        elif task == "Equil_Tf":
            print "Starting to check if Equil_Tf completed..."
            simulation.Tf_loop.check_completion(System,self.append_log,equil=True)
            lasttime2,action2,task2 = self.check_modelbuilder_log(sub)
            if action2 == "Finished:":
                print "Finished Equil_Tf_iteration..."
                print "Starting Equil_Tf_analysis..."
                analysis.Tf_loop.analyze_temperature_array(System,self.append_log,equil=True)
        elif task == "Equil_Tf_analysis":
            print "Starting to check if Equil_Tf_analysis completed..."
            analysis.Tf_loop.check_completion(System,self.append_log,equil=True)
        else:
            print "ERROR!"
            print "  Couldn't find next option for task:",task
            print "  Please check that things are ok."
            print "  Exiting."
            raise SystemExit

    def logical_flowchart_finished(self,System,Model,sub,task):
        if task == "Tf_loop_iteration":
            print "Finished Tf_loop_iteration..."
            print "Starting Tf_loop_analysis..."
            analysis.Tf_loop.analyze_temperature_array(System,self.append_log)
        elif task == "Tf_loop_analysis":
            print "Finished Tf_loop_analysis..."
            flag = analysis.Tf_loop.check_if_wham_is_next(System,self.append_log)
            if flag == 1:
                pass 
            else:
                print "Starting Tf_loop_iteration..."
                simulation.Tf_loop.folding_temperature_loop(Model,System,self.append_log)
        elif task == "wham_Cv":
            print "Finished wham_Cv..."
            print "Stating wham_FreeEnergy..."
            analysis.Tf_loop.continue_wham(System,self.append_log)
        elif task == "Equil_Tf":
            print "Starting Equil_Tf_analysis..."
            analysis.Tf_loop.analyze_temperature_array(System,self.append_log,equil=True)
        elif task == "Equil_Tf_analysis" or task == "Aggregating_Equil_Runs":
            ## Aggregrate equil_Tf data for each temperature and plot PMFs
            print "Starting aggregate data..."
            analysis.Tf_loop.aggregate_equilibrium_runs(System,self.append_log)
            print "Plotting aggregated data PMFS..."
            analysis.plot.pmfs.plot_aggregated_data(System,self.append_log)
        elif task == "Plotting_Agg_Data":
            if Model.modelnameshort in ["HomGo","HetGo","DMC"]:
                print "Starting prepping mutant pdbs..."
                #mutations.preppdbs.prep_mutants(System,self.append_log)
                print "Starting calculating dH for mutants..."
                mutations.phi_values.calculate_dH_for_mutants(Model,System,self.append_log)
        elif task == "Calculating_dH":
            mutations.phi_values.calculate_phi_values(Model,System,self.append_log,"Q")
            #mutations.phi_values.calculate_new_epsilons(Model,System,self.append_log)
        else:
            print "ERROR!"
            print "  Couldn't find next option for task:",task
            print "  Please check that things are ok."
            print "  Exiting."
            raise SystemExit

    def new_project(self,args,modeloptions):
        ''' Starting a new simulation project.'''
        subdirs = [ x[:-4] for x in args.pdbs ]
        for sub in subdirs:
            if os.path.exists(sub) == False:
                os.mkdir(sub)
            else:
                print "Subdirectory: ", sub, " already exists! just fyi"

        print "Starting a new simulation project..."
        Models = models.new_models(subdirs,modeloptions)
        Systems = systems.new_systems(subdirs)
        Model = Models[0]
        System = Systems[0]

        self.prepare_systems(Models,Systems)
        self.save_model_system_info(Models,Systems,subdirs)

        ## Not implemented yet.
        if args.temparray != None:
            System.initial_T_array = args.temparray

        ## The first step depends on the type of model.
        if args.type in ["HomGo","HetGo"]:
            for k in range(len(subdirs)):
                print "Starting Tf_loop_iteration for subdirectory: ", subdirs[k]
                ## To Do: Prepare each Model System pair. 
                Model = Models[k]
                System = Systems[k]
                simulation.Tf_loop.folding_temperature_loop(Model,System,self.append_log,new=True)
        elif args.type == "DMC":
            pass

        self.save_model_system_info(Models,Systems,subdirs)
        print "Success"

if __name__ == "__main__":
    pass
