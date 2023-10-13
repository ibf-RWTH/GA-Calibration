import os
import sys, getopt
import time
import numpy as np
import pandas as pd
import ast
pd.options.mode.chained_assignment = None  # default='warn'
import shutil
import matplotlib.pyplot as plt
from sklearn.metrics import mean_squared_error
import subprocess
import json
import math
from geneticalgorithm2 import geneticalgorithm2 as ga
from io import StringIO
import warnings
import configparser

warnings.simplefilter(action='ignore', category=FutureWarning)


class MatParams:
    MatID2MatProps = [

            ['pw_fl','brgvec','v0','rho_0','km1','km2','c3','hdrt_0','crss_0','crss_s','pw_hd','Adir','Adyn','gam_c','pw_irr','k'],
            ['pw_fl','brgvec','v0','rho_0','km1','km2','c3'],
            ['pw_fl','hdrt_0','crss_0','crss_s','pw_hd','Adir','Adyn','gam_c','pw_irr'],
            ['pw_fl','hdrt_0','shrt_0', 'crss_0','k','crss_s','pw_hd','Adir','Adyn','gam_c','pw_irr'],
            ['pw_fl','shrt_0', 'crss_0','k','crss_s','pw_hd','Adir','Adyn','gam_c','pw_irr']

    ]

    def __init__(self, root:str) -> None:
        self.root = root
        self.MatID2MatPropsBound = {}
        self.material_ids = []

    def read_varbounds(self):
        # Open and read the file
        with open(f"{self.root}/sample_files/matvarbound.txt",'r') as f:

            phase_number = -1
            for line in f:
                 # Remove leading and trailing whitespace
                line = line.strip()

                # Check if the line starts with "<:" and ends with ":>"
                if line.startswith("<:") and line.endswith(":>"):
                # Extract the section name (e.g., "Global" or "Material: 2")
                    split_line = line.split(":", 2)
                    section_name = split_line[1].strip() + split_line[2].strip()[0]

                    #get material phase numbers
                    phase_number = int(section_name[-1])
                    self.material_ids.append(phase_number)
                    self.MatID2MatPropsBound[phase_number] = {}

                else:
                    # Split the line into property and value
                    if ":" in line:
                        # check if phase number correct
                        if phase_number == -1:
                            os.system("echo Input File error: get wrong phase index, please check your matvarbound.txt input file")
                            sys.exit(-1)

                        parts = line.split(":")
                        mat_prop_name = parts[0].strip()
                        mat_prop_value = ast.literal_eval(parts[1].strip())
                        # set mat params dict value and check matvarbound.txt
                        if mat_prop_name in self.MatID2MatProps[phase_number]:
                            self.MatID2MatPropsBound[phase_number][mat_prop_name] = mat_prop_value
                        else:
                            os.system(f"echo Input File error: invalid material parameters {mat_prop_name} for material {phase_number}, please check your matvarbound.txt input file")
                            sys.exit(1)

    def get_varbounds(self)-> np.ndarray:

        varbound = []
        self.read_varbounds()
        for mat_id in self.MatID2MatPropsBound:
            mat_props_bounds = self.MatID2MatPropsBound[mat_id]
            for value in mat_props_bounds.values():
                varbound.append(value)
        varbound = np.array(varbound)

        return varbound


class Simulation:
    def __init__(self, sim_root, ex_data, job_name, sim_type, n_jobs, mat_params: MatParams):
        self.sim_root = sim_root
        self.ex_data = ex_data
        self.subroutine_dir = f'{self.sim_root}/subroutine'
        self.simulation_dir = f'{self.sim_root}/simulation'
        self.sample_files = f'{self.sim_root}/sample_files'
        self.images = f'{self.sim_root}/evaluation_images'
        self.log_dir = f'{self.sim_root}/logs'
        self.base_job_name = job_name
        self.job_name = job_name
        self.n_jobs = n_jobs
        self.sim_type = sim_type
        self.sim_type2compare_function = {
            'CP_cyclic' : self.compare_exp2sim_cyclic,
            'CP_tensile' : self.compare_exp2sim_tensile,
            'Chaboche': self.compare_exp2sim_chaboche
        }
        self.mat_params = mat_params
        self.n_phases = len(mat_params) - 1

    def create_batch_job_script(self, job_index):
        current_simulation_dir = f'{self.simulation_dir}_{job_index}'
        destination = f'{current_simulation_dir}/simulation_job_{job_index}.sh'
        if 'CP' in self.sim_type:
            source = f'{self.sample_files}/simulation_job_CP.sh'
        elif 'Chaboche' in self.sim_type:
            source = f'{self.sample_files}/simulation_job_Chaboche.sh'
        shutil.copy(source, destination)

    def blackbox(self, params):
        submitted = False
        while not submitted:
            # define path variables
            job_index = str(time.time_ns())[-4:]
            current_simulation_dir = f'{self.simulation_dir}_{job_index}'
            current_job_name = f'{self.base_job_name}_{job_index}'
            if not os.path.isdir(current_simulation_dir):
                #create directory
                self.create_job_dir(current_simulation_dir)
                if 'CP' in self.sim_type:
                    for j, mat_id in enumerate(self.mat_params.keys()):
                        if j > 0:
                            if j == 1:
                                path = self.sample_files
                            else:
                                path = current_simulation_dir
                            # create matdata.inp file according to params
                            self.manipulate_matdata(path, current_simulation_dir, mat_id, phase_index=j, optiParams=params)
                elif 'Chaboche' in self.sim_type:
                    for j, mat_id in enumerate(self.mat_params.keys()):
                        keys = self.mat_params.keys()
                        self.manipulate_matdata(self.sample_files, current_simulation_dir, mat_id, phase_index=j+1, optiParams=params)

                # create batch script
                self.create_batch_job_script(job_index=job_index)
                # submit simulation
                self.submit_batch_job(current_simulation_dir, f'simulation_job_{job_index}.sh', current_job_name)
                time.sleep(120)

                # check simulation status and wait until simulation finished
                sim_running = True
                while sim_running:
                    sim_running, complete_status = self.check_sim_status(job_name=current_job_name)
                if not complete_status:
                    os.system(f'echo {complete_status}')
                    os.system('echo in if loop')
                    # simulation crashed due to external error remove all simfiles prepare for resubmit
                    self.remove_sim_files(current_simulation_dir)
                    # delete log file only if simulations was successful
                    #if complete_status:
                    self.remove_sim_files(f'{self.log_dir}/{current_job_name}-log*')
                    submitted = False
                    self.num_props = [0]
                    continue
                # evaluate Simulation
                if 'CP' in self.sim_type:
                    sim_results = self.calcStressStrain(current_simulation_dir, current_job_name)
                elif 'Chaboche' in self.sim_type:
                    colNames = ['sim_time', 'sim_displacement', 'sim_force']
                    sim_results = pd.read_csv(current_simulation_dir+'/RF_data.txt', names=colNames)
                compare_func = self.sim_type2compare_function[self.sim_type]
                mad1, mad2, time_stamp = compare_func(sim_results)
                mad = mad1 + mad2
                # write results to files and delete simulation files
                self.write_results_file(mad_time=mad1, mad_stress=mad2, mad=mad, params=params, compare_func=compare_func.__name__, time_stamp=time_stamp)
                self.remove_sim_files(current_simulation_dir)
                # delete log file only if simulations was successful
                #if complete_status:
                self.remove_sim_files(f'{self.log_dir}/{current_job_name}-log*')

                return mad
            else:
                os.system(f'echo {current_simulation_dir} already taken')
        time.sleep(120) # if all folders are locked find again afer 120s

    def create_job_dir(self, dst_dir)  -> None:
        # create Folder for next simulation
        source_dir = f'{self.sim_root}/sample_files_{str(self.sim_type).lower()}_simulation'
        len_sample_dir = len(os.listdir(source_dir))
        if not os.path.exists(dst_dir):
                shutil.copytree(source_dir, dst_dir)
        while len(os.listdir(dst_dir)) < len_sample_dir:
            # gotta wait until all files are copied
            time.sleep(1)
        # some more waiting to make sure everything is complete
        time.sleep(10)

    def remove_sim_files(self, path_to_remove) -> None:
        os.system(f'rm -rf {path_to_remove}')
        time.sleep(10)

    def write_results_file(self, mad_time, mad_stress, mad, params, compare_func, time_stamp) -> None:
        f = open(self.log_dir + '/results.txt', 'a+')
        f.write(f'simulation type: {self.sim_type}\n')
        f.write(f'Compute Error Using: {compare_func}\n')
        f.write(f'time_Stamp: {time_stamp}\n')
        f.write(f'Error_Time: {mad_time}\n')
        f.write(f'Error_Stress: {mad_stress}\n')
        f.write(f'Error: {mad}\n')
        prop_index = 0
        phase_index = 0
        for phase_id, mat_props in self.mat_params.items():
            if phase_id == 0: #global parameters
                f.write(f'global: ')
            else:
                f.write(f'phase{phase_id}: ')
            for mat_props_name in mat_props:
                f.write(f'{mat_props_name}: {params[prop_index]}, ')
                prop_index += 1

            f.write('\n')
            phase_index += 1
        f.write('############################################################\n')
        f.close()

    def check_sim_status(self, job_name):
        # This will call a shell output in linux to ask for the job status
        stdout = None
        time.sleep(30)
        while stdout is None:
            try:
                sacct_command = [
                'sacct',
                '-o', 'JobName%100,State'
                ]

                stdout = subprocess.check_output(sacct_command, stderr=subprocess.STDOUT, text=True)
                stdout = StringIO(stdout)
            except Exception as e:
                print(f'whoops something went wrong: {e}')
                stdout = None

        # reshaping stdOutput to make it pandas readable
        try:
           lines = stdout.readlines()
           del lines[1]
           data = [line.strip().split() for line in lines]
           headers = data[0]
           df = pd.DataFrame(data[1:], columns=headers)
        except:
            # something went wrong with the sacct command retry again and return as if ev everything was still runnning
            sim_status = True
            complete_status = False
            return sim_status, complete_status

        batch_state = df.loc[df['JobName']=='batch','State'].values[-1]
        extern_state = df.loc[df['JobName']=='extern','State'].values[-1]
        try:
            job_state = df.loc[df['JobName']== f'{job_name}','State'].values[-1]
        except Exception as e:
            if os.path.exists(f'{self.sim_root}/simulation_{job_name[-4:]}/{job_name}.sta'):
                with open(f'{self.sim_root}/simulation_{job_name[-4:]}/{job_name}.sta','r') as f:
                    content = f.read()
                    if "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" in content:
                        job_state = 'COMPLETED'
                    else:
                        job_state = "PENDING"
            elif os.path.exists(f'{self.sim_root}/simulation_{job_name[-4:]}/00_Data'):
                job_state = "PENDING"
            else:
                os.system(f"echo something went wrong I had to enter the else in the except: {e}")
                sys.exit()

        # check for successful completion of simulation
        if job_state == 'COMPLETED':
            # break while loop if completed
            sim_status = False
            complete_status = True
        # check if simulation is still running
        elif job_state == 'RUNNING':
            sim_status = True
            complete_status = False
        elif job_state == 'PENDING':
            sim_status = True
            complete_status = False
        else:
            # in this case simulation is not running but also wasn't successfully
            # completed e.g. failed --> sim_status= False
            os.system(f'echo batch_state: {batch_state}')
            os.system(f'echo extern_state: {extern_state}')
            os.system(f'echo job_state: {job_state}')
            sim_status = False
            complete_status = False

        return sim_status, complete_status

    def calcStressStrain(self, current_simulation_dir, job_name):
        datadir = current_simulation_dir + '/00_Data/'
        data_grainID = datadir + "data_grain_ID_" + job_name + ".csv"
        data_time = datadir + "data_time_" + job_name + ".csv"

        grainID_df = pd.read_csv(data_grainID)
        grainID_df.sort_values(by=['element_id'], inplace=True)

        if self.n_phases > 1:
            data_matID = datadir + "data_mat_ID_" + job_name + ".csv"
            df_matID = pd.read_csv(data_matID)  # has length of amount of grains
            grainID_df.loc[:,'matID'] = 1

            for i in range(len(df_matID)):  # iterating over all grains in matID
                grainID_df.loc[grainID_df['grainID'] == i, 'matID'] = df_matID.loc[i, 'matID']

        else:
            grainID_df.loc[:,'matID'] = 1
        files = os.listdir(datadir)

        vol_files = [file for file in files if 'el_vol' in file]
        Stress_files = [file for file in files if 'Stress' in file]
        Strain_files = [file for file in files if 'Strain' in file]

        vol_files.sort()
        Stress_files.sort()
        Strain_files.sort()

        time_df = pd.read_csv(data_time)
        flowcurve_df = pd.DataFrame()

        for i in range(0, len(vol_files)):
            vol_df = pd.read_csv(datadir + vol_files[i])
            Stress_df = pd.read_csv(datadir + Stress_files[i])
            Strain_df = pd.read_csv(datadir + Strain_files[i])

            # average over intergration points in every element
            Stress_df_mean = Stress_df.groupby('element_id').mean()
            Strain_df_mean = Strain_df.groupby('element_id').mean()

            alldata = pd.merge(left=Stress_df_mean, right=grainID_df, on='element_id')
            alldata = pd.merge(alldata, Strain_df_mean, on='element_id')
            alldata = pd.merge(alldata, vol_df, on='element_id')
            grain_vol = alldata[['el_vol', 'grainID']].groupby('grainID').sum()
            grain_vol = grain_vol.rename(columns={'el_vol': 'grain_vol'})
            alldata = alldata.groupby('grainID').mean()
            alldata = pd.merge(alldata, grain_vol, on='grainID')

            alldata['vol_avg_Stress'] = alldata['Stress'].values * alldata['grain_vol'].values / sum(
                alldata['grain_vol'].values)
            alldata['vol_avg_Strain'] = alldata['Strain'].values * alldata['grain_vol'].values / sum(
                alldata['grain_vol'].values)
            if self.n_phases > 1:

                phase_1_data = alldata[alldata['matID'] == 1]
                phase_1_data.reset_index(drop=True, inplace=True)
                vol_avg_Stress_phase1 = phase_1_data['Stress'].values * \
                                        phase_1_data['grain_vol'].values / \
                                        sum(phase_1_data['grain_vol'].values)
                vol_avg_Strain_phase1 = phase_1_data['Strain'].values * \
                                        phase_1_data['grain_vol'].values / \
                                        sum(phase_1_data['grain_vol'].values)
                phase_1_data.loc[:, 'vol_avg_Stress_phase1'] = vol_avg_Stress_phase1
                phase_1_data.loc[:, 'vol_avg_Strain_phase1'] = vol_avg_Strain_phase1
                phase_1_temp_stress = sum(phase_1_data['vol_avg_Stress_phase1'].values)
                phase_1_temp_strain = sum(phase_1_data['vol_avg_Strain_phase1'].values)

                phase_2_data = alldata[alldata['matID'] == 3]

                phase_2_data.reset_index(drop=True, inplace=True)
                vol_avg_Stress_phase2 = phase_2_data['Stress'].values * \
                                        phase_2_data['grain_vol'].values / \
                                        sum(phase_2_data['grain_vol'].values)
                vol_avg_Strain_phase2 = phase_2_data['Strain'].values * \
                                        phase_2_data['grain_vol'].values / \
                                        sum(phase_2_data['grain_vol'].values)
                phase_2_data.loc[:, 'vol_avg_Stress_phase2'] = vol_avg_Stress_phase2
                phase_2_data.loc[:, 'vol_avg_Strain_phase2'] = vol_avg_Strain_phase2
                phase_2_temp_stress = sum(phase_2_data['vol_avg_Stress_phase2'].values)
                phase_2_temp_strain = sum(phase_2_data['vol_avg_Strain_phase2'].values)

                ### choose variables below ###
                temp_stress = sum(alldata['vol_avg_Stress'].values)
                temp_strain = sum(alldata['vol_avg_Strain'].values)
                ##############################
                temp_dict = {'Strain': [temp_strain], 'Stress': [temp_stress],
                             'Strain_Phase1': [phase_1_temp_strain], 'Stress_Phase1': [phase_1_temp_stress],
                             'Strain_Phase2': [phase_2_temp_strain], 'Stress_Phase2': [phase_2_temp_stress]}

            else:
                ### choose variables below ###
                temp_stress = sum(alldata['vol_avg_Stress'].values)
                temp_strain = sum(alldata['vol_avg_Strain'].values)
                ##############################
                temp_dict = {'Strain': [temp_strain], 'Stress': [temp_stress]}

            temp_df = pd.DataFrame(temp_dict)
            flowcurve_df = pd.concat([flowcurve_df, temp_df])

        if flowcurve_df.shape[0] > 0:
            results_df = time_df
            results_df['Strain'] = flowcurve_df['Strain'].values
            results_df['Stress'] = flowcurve_df['Stress'].values
            if self.n_phases > 1:
                results_df['Strain_Phase1'] = flowcurve_df['Strain_Phase1'].values
                results_df['Stress_Phase1'] = flowcurve_df['Stress_Phase1'].values
                results_df['Strain_Phase2'] = flowcurve_df['Strain_Phase2'].values
                results_df['Stress_Phase2'] = flowcurve_df['Stress_Phase2'].values
        else:
            results_df = time_df
            results_df['Strain'] = 0
            results_df['Stress'] = 0
            results_df.to_csv(self.sim_root+'results_df.csv')

        return results_df

    def compare_exp2sim_cyclic(self, simulation_df):
        assert self.n_phases == 1 or self.n_phases == 2, 'more than two phases are not yet supported'
        # in order for this to work correctly make sure the time intervals in experiment and simulation are equal!!!
        now = int(time.time())
        if simulation_df.shape[0] < 5:
            mad_time = 99999.
            mad_stress_strain = 99999
            return mad_time, mad_stress_strain, now


        experimental_df = pd.read_csv(f'{self.sample_files}/{self.ex_data}', sep='\t')

        total_time = experimental_df['sim_time'].values[-1]
        max_sim_time = simulation_df['time'].values[-1]
        mad_time = (100 * (1 - max_sim_time / total_time))**2
        comp_df = experimental_df.merge(simulation_df, left_on='sim_time',
                                        right_on='time')  # merge dfs on time for comparison

        if self.n_phases == 1:
            mad_stress = np.mean(np.abs(comp_df['stress_t'] - comp_df['Stress'])**2)
            sim_y_cols = ['Stress']
            sim_x_cols = ['Strain'] * len(sim_y_cols)
            sim_labels = ['Simulation']

            ex_y_cols = ['stress_t']
            ex_x_cols = ['strain_t'] * len(ex_y_cols)
            ex_labels = ['Experiment']

            fig_name = f'hysteresis_{now}'
            x_label = "Strain"
            y_label = "Stress(MPa)"
            self.plot_data2(fig_name=fig_name, x_label=x_label,y_label=y_label,
                            sim_data=simulation_df, ex_data=experimental_df,
                            sim_x_cols=sim_x_cols, sim_y_cols=sim_y_cols, sim_labels=sim_labels,
                            ex_x_cols=ex_x_cols, ex_y_cols=ex_y_cols, ex_labels=ex_labels)

        elif self.n_phases == 2:
            mad_stress_total = np.mean(np.abs(comp_df['stress_t'] - comp_df['Stress']))
            mad_stress_phase1 = np.mean(np.abs(comp_df['Stress_Phase1_t'] - comp_df['Stress_Phase1_y']))
            mad_stress_phase2 = np.mean(np.abs(comp_df['Stress_Phase2_t'] - comp_df['Stress_Phase2_y']))
            mad_stress = (0.8 * mad_stress_phase1 + 0.2 * mad_stress_phase2 + mad_stress_total) / 3

            self.plot_data(f'stress_vs_time_phase_1_{now}', 'time ( s )', 'Stress ( MPa )',
                           comp_df['sim_time'], comp_df['Stress_Phase1_t'], 'Experimental Data phase 1',
                           comp_df['time_y'], comp_df['Stress_Phase1_y'], 'Simulation Data phase 1')

            self.plot_data(f'stress_vs_time_phase_2_{now}', 'time ( s )', 'Stress ( MPa )',
                           comp_df['sim_time'], comp_df['Stress_Phase2_t'], 'Experimental Data phase 2',
                           comp_df['time_y'], comp_df['Stress_Phase2_y'], 'Simulation Data phase 2')

        # self.plot_data(f'hysteresis_{now}', 'Strain ( - )', 'Stress (MPa)',
        #                comp_df['strain'], comp_df['stress_t'], 'Experimental Data',
        #                comp_df['Strain'], comp_df['Stress'], 'Simulation Data')


        return mad_time, mad_stress, now

    def compare_exp2sim_tensile(self, simulation_df):
        assert self.n_phases == 2 or self.n_phases == 1
        now = int(time.time())
        if simulation_df.shape[0] < 5:
            mad_time = 99999.
            mad_stress_strain = 99999
            return mad_time, mad_stress_strain, now

        experimental_df = pd.read_csv(f'{self.sample_files}/{self.ex_data}', sep=',')
        max_exp_strain = experimental_df['strain_t'].max()
        max_sim_strain = simulation_df['Strain'].max()

        if self.n_phases == 2:

            exp_total_stress_interp = np.interp(simulation_df['Strain'], experimental_df['strain_t'], experimental_df['stress_t'])
            exp_alpha_stress_interp = np.interp(simulation_df['Strain'], experimental_df['strain_t'], experimental_df['stress_alpha'])
            exp_beta_stress_interp = np.interp(simulation_df['Strain'], experimental_df['strain_t'], experimental_df['stress_beta'])

            mad_stress_total = np.mean(np.abs(exp_total_stress_interp - simulation_df['Stress']))
            mad_stress_alpha = np.mean(np.abs(exp_alpha_stress_interp - simulation_df['Stress_Phase1']))
            mad_stress_beta = np.mean(np.abs(exp_beta_stress_interp - simulation_df['Stress_Phase2']))
            mad_strain_total = (abs(1 - max_sim_strain / max_exp_strain) * 100) **2
            mad_stress = (mad_stress_total + 0.8*mad_stress_alpha + 0.2 * mad_stress_beta)/3

            sim_y_cols = ['Stress', 'Stress_Phase1', 'Stress_Phase2']
            sim_x_cols = ['Strain'] * len(sim_y_cols)
            sim_labels = ['Simulation_Total', 'Simulation_Alpha', 'Simulation_Beta']

            ex_y_cols = ['stress_t', 'stress_alpha', 'stress_beta']
            ex_x_cols = ['strain_t'] * len(ex_y_cols)
            ex_labels = ['Experiment_Total', 'Experiment_Alpha','Experiment_Beta']

        if self.n_phases == 1:

            exp_stress_interp = np.interp(simulation_df['Strain'], experimental_df['strain_t'], experimental_df['stress_t'])
            mad_stress = np.mean(np.abs(exp_stress_interp - simulation_df['Stress']) / exp_stress_interp) * 100
            mad_strain_total = (abs(1 - max_sim_strain / max_exp_strain) * 100) **2

            sim_y_cols = ['Stress']
            sim_x_cols = ['Strain'] * len(sim_y_cols)
            sim_labels = ['Simulation']

            ex_y_cols = ['stress_t']
            ex_x_cols = ['strain_t'] * len(ex_y_cols)
            ex_labels = ['Experiment']

        fig_name = f'stress_strain_{now}'
        x_label = "Strain"
        y_label = "Stress(MPa)"
        self.plot_data2(fig_name=fig_name, x_label=x_label,y_label=y_label,
                        sim_data=simulation_df, ex_data=experimental_df,
                        sim_x_cols=sim_x_cols, sim_y_cols=sim_y_cols, sim_labels=sim_labels,
                        ex_x_cols=ex_x_cols, ex_y_cols=ex_y_cols, ex_labels=ex_labels)

        return mad_strain_total, mad_stress, now

    def compare_exp2sim_chaboche(self, simulation_df):
        # TODO: This part needs to be fixed
        # in order for this to work correctly make sure the time intervals in experiment and simulation are equal!!!
        now = int(time.time())
        if simulation_df.shape[0] < 5:
            mad_time = 99999.
            mad_stress_strain = 99999
            return mad_time, mad_stress_strain, now


        experimental_df = pd.read_csv(f'{self.sample_files}/{self.ex_data}')

        total_time = experimental_df['time'].values[-1]
        max_sim_time = simulation_df['sim_time'].values[-1]
        mad_time = (100 * (1 - max_sim_time / total_time))**2
        comp_df = experimental_df.merge(simulation_df, left_on='time',
                                        right_on='sim_time')  # merge dfs on time for comparison

        mad_force = np.mean(np.abs(comp_df['force'] - comp_df['sim_force'])*1000)
        sim_y_cols = ['sim_force']
        sim_x_cols = ['sim_time'] * len(sim_y_cols)
        sim_labels = ['Simulation']

        ex_y_cols = ['force']
        ex_x_cols = ['time'] * len(ex_y_cols)
        ex_labels = ['Experiment']

        fig_name = f'Force_{now}'
        x_label = "time (s)"
        y_label = "Force (kN)"
        self.plot_data2(fig_name=fig_name, x_label=x_label,y_label=y_label,
                        sim_data=simulation_df, ex_data=experimental_df,
                        sim_x_cols=sim_x_cols, sim_y_cols=sim_y_cols, sim_labels=sim_labels,
                        ex_x_cols=ex_x_cols, ex_y_cols=ex_y_cols, ex_labels=ex_labels)
        
        return mad_time, mad_force, now

    def plot_data(self, fig_name, x_label, y_label, x1, y1, data_label_1, x2=None, y2=None, data_label_2=None):
        plt.plot(x1, y1, label=data_label_1)

        if x2 is not None:
            assert (y2 is not None)
            assert (data_label_2 is not None)
            plt.plot(x2, y2, label='Simulation Data')
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.legend()
        plt.savefig(f'{self.images}/{fig_name}.png')
        plt.close()

    def plot_data2(self, fig_name:str, x_label:str, y_label:str, sim_data:pd.DataFrame, ex_data:pd.DataFrame,
                   sim_x_cols:list, sim_y_cols:list, sim_labels:list,
                   ex_x_cols:list, ex_y_cols:list, ex_labels:list):
        assert len(sim_x_cols) == len(sim_y_cols) and len(ex_x_cols) == len(ex_y_cols)

        sim_label_index = 0
        for x, y in zip(sim_x_cols, sim_y_cols):
            plt.plot(sim_data[x].values, sim_data[y].values, label = sim_labels[sim_label_index])
            sim_label_index += 1

        ex_label_index = 0
        for x, y in zip(ex_x_cols, ex_y_cols):
            plt.plot(ex_data[x].values, ex_data[y].values, label = ex_labels[ex_label_index])
            ex_label_index += 1

        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.legend()
        plt.savefig(f'{self.images}/{fig_name}.png')
        plt.close()

    def submit_batch_job(self, currentDir, file, jobName):
        sim_type = self.sim_type
        config = configparser.ConfigParser()
        config.read(f'{sim_root}/configs/configs.ini')
        account = config.get('JobSettings','account')
        output = config.get('MainProcessSettings','output')+jobName+'-log.%J'
        time = config.get('JobSettings','time')
        memPerCpu = config.get('JobSettings','mem-per-cpu')
        nodes = config.get('JobSettings','nodes')
        nTasks = config.get('JobSettings','ntasks')

        cmd = f'sbatch --job-name={jobName}'
        if account != 'None':
            cmd += f' --account={account}'
        cmd += f' --output={output} --time={time}' 
        cmd += f' --mem-per-cpu={memPerCpu} --nodes={nodes} --ntasks={nTasks}  {file}'

        abaqus = config.get('JobSettings','abaqus_version')
        memory = config.get('JobSettings','ABAQUS_MEM_ARG')
        threads = config.get('JobSettings','THREADS_PER_MPI')
        root = config.get('JobSettings','ROOT')
        if 'CP' in sim_type:
            pythonPath = config.get('JobSettings','PYTHON_PATH')+'/readOdb_CP.py'
        elif 'Chaboche' in sim_type:
            pythonPath = config.get('JobSettings','PYTHON_PATH')+'/readOdb_Chaboche.py'
        subroutine = config.get('JobSettings','SUBROUTINE_PATH')
        input = config.get('JobSettings', 'INPUTFILE')
        f = open(f'{currentDir}/SimJobConfig.sh', 'w+')
        f.write(f'ABAQUS={abaqus}\n')
        f.write(f'ABAQUS_MEM_ARG={memory}\n')
        f.write(f'THREADS_PER_MPI={threads}\n')
        f.write(f'ROOT={root}\n')
        f.write(f'SIM_DIR={currentDir}\n')
        f.write(f'PYTHON_PATH={pythonPath}\n')
        f.write(f'SUBROUTINE_PATH={subroutine}\n')
        f.write(f'JOBNAME={jobName}\n')
        f.write(f'INPUTFILE={input}\n')
        f.close()
        main_cwd = os.getcwd()      
        os.chdir(currentDir)
        os.system(cmd)
        os.chdir(main_cwd)

    def manipulate_matdata(self, sample_path:str, current_sim_dir:str, id:int, phase_index:int, optiParams):

        """
        sample_path: path to simulation sample files
        current_sim_dir: directory path to current simulation
        id: abaqus material id
        phase_index: internal id starting with 1 going number of phases e.g. 2
        optiParams: material parameters coming from optimizer
        """

        # Manipulate matdata.inp file according to optimizer
        # num_prev_props = np.sum(self.num_props)
        constant_params = Utils.CONSTANTPARAMS
        mat_params = Utils.EVALUATINGPARAMS
        
        if 'CP' in self.sim_type:
            mat_props_keys = mat_params[0] # Global Parameters keys/names
            mat_props_keys.extend(mat_params[id]) # adding Parameter keys/names for current phase
            # num_curr_props = len(mat_props)
            if phase_index == 1:
                mat_props_values = optiParams[:len(mat_props_keys)]
            else:
                mat_props_values = [optiParams[0]]
                mat_props_values.extend(optiParams[-(len(mat_props_keys) - 1):])

            with open(f'{sample_path}/matdata.inp', 'r') as file:
                    lines = file.readlines()

            prop_index_var = 0
            prop_index_const = 0
            for option in Utils.CONFIG.options(self.sim_type):
                if option == 'abaqus_id':
                    continue
                if option in mat_props_keys:
                    value = mat_props_values[prop_index_var]
                    prop_index_var +=1
                elif option in constant_params:
                    value = constant_params[id][prop_index_const]
                    prop_index_const +=1
                lines = [line.replace(f'%{option}%{id}%', f'{np.round(value,4)}') for line in lines]

            f = open(f'{current_sim_dir}/matdata.inp', 'w+')
            for line in lines:
                f.write(line)
            f.close()
            # self.num_props.append(num_curr_props)
        elif 'Chaboche' in self.sim_type:
            # TODO: There are errors occuring in here.....
            mat_props_values = optiParams
            mat_props_keys = mat_params

            with open(f'{sample_path}/Material.inp', 'r') as file:
                lines = file.readlines()
            prop_index_var = 0
            prop_index_const = 0
            for option in Utils.CONFIG.options(self.sim_type):
                if option == 'abaqus_id':
                    continue
                if option in mat_props_keys[id]:
                    value = mat_props_values[prop_index_var]
                    prop_index_var +=1
                elif option in constant_params[id]:
                    value = constant_params[id][option]
                    prop_index_const +=1
                lines = [line.replace(f'%{option}%', f'{np.round(value,4)}') for line in lines]

            f = open(f'{current_sim_dir}/Material.inp','w+')
            f.writelines(lines)
            f.close()

                

class Optimize:

    def __init__(self, flag, ex_data, root, name, mat_params, varbound, algorithm_param, sim_type, n_jobs):
        self.test_flag = flag
        self.ex_data = ex_data
        self.sim_root = root
        self.job_name = name
        self.mat_params = Utils.EVALUATINGPARAMS
        self.varbound = varbound
        self.algorithm_param = algorithm_param
        self.sim_type = sim_type
        self.n_jobs = n_jobs

    def init_optimizer(self):
        sim_object = Simulation(sim_root=self.sim_root, ex_data=self.ex_data, mat_params=self.mat_params,
                                job_name=self.job_name, sim_type=self.sim_type, n_jobs=self.n_jobs)
        
        if len(self.mat_params) >= 1:
            blackbox_func = sim_object.blackbox
        else:
            os.system("echo error: empty mat_params")
            sys.exit()

        model = ga(function=blackbox_func,
                   dimension=len(self.varbound),
                   variable_type='real',
                   variable_boundaries=self.varbound,
                   algorithm_parameters=self.algorithm_param,
                   function_timeout=None)
        return model, blackbox_func


class Parser:
    def __init__(self, sim_type) -> None:
        self.sim_type = sim_type

    def parse_matparams(self):
        config = Utils.CONFIG
        varbound = []
        constantParams = {}
        matparams = {}
        params_names = []
        phases = config.get('JobSettings','PHASES').split(',')
        if 'CP' in self.sim_type:
            global_params = config.options('GlobalParams')
            
            # add global params to varbound
            for global_param in global_params:
                global_param_value = ast.literal_eval(config.get('GlobalParams',global_param))

                if global_param != "abaqus_id":
                    # check params type
                    if not isinstance(global_param_value, list):
                        os.system("echo config error: global params to be calibrated must be list")
                        sys.exit()

                    params_names.append(global_param)
                    # add global params
                    varbound.append(global_param_value)
            matparams[0] = params_names
            params_names = []

            for phase in phases:
                phase_id = config.get(phase, 'abaqus_id')
                options = config.options(phase)
                constantParams[phase_id] = {}
                #add phase params to varbound
                for option in options:
                    value = config.get(phase, option)
                    if value.startswith('['):
                        value = ast.literal_eval(value)
                        params_names.append(option)
                        varbound.append(value)
                    else:
                        value = ast.literal_eval(value)
                        constantParams[phase_id][option] = value

                matparams[phase_id] = params_names
                params_names = []

            varbound = np.array(varbound)

        elif 'Chaboche' in self.sim_type:
            phase_id = config.get(phases[0], 'abaqus_id')
            constantParams[phase_id] = {}
            options = config.options(phases[0])
            #add phase params to varbound
            for option in options:
                value = config.get(phases[0], option)
                if value.startswith('['):
                    value = ast.literal_eval(value)
                    params_names.append(option)
                    varbound.append(value)
                else:
                    value = ast.literal_eval(value)
                    constantParams[phase_id][option] = value

            
            matparams[phase_id] = params_names
            varbound = np.array(varbound)

        return constantParams, matparams, varbound

class Utils:
    CONFIG = None
    CONSTANTPARAMS = None
    EVALUATINGPARAMS = None       

if __name__ == '__main__':
    #get current path
    sim_root = os.getcwd()
    os.system(f"echo sim root: {sim_root}")
    if not os.path.isdir(f'{sim_root}/evaluation_images'):
        os.mkdir(f'{sim_root}/evaluation_images')
    if not os.path.isdir(f'{sim_root}/logs'):
        os.mkdir(f'{sim_root}/logs')
    #config simulation
    config = configparser.ConfigParser()
    config.read(f'{sim_root}/configs/configs.ini')
    Utils.CONFIG = config
    sim_Type = config.get('JobSettings','sim_type')
    
    constant_params, mat_params, varbound = Parser(sim_Type).parse_matparams()
    
    Utils.CONSTANTPARAMS = constant_params
    Utils.EVALUATINGPARAMS = mat_params

    #initialize Optimizer
    algorithm_param = {'max_num_iteration': 100, \
                        'population_size': 50, \
                        'mutation_probability': 0.1, \
                        'elit_ratio': 0.1, \
                        'parents_portion': 0.3, \
                        'max_iteration_without_improv': None}

    opt = Optimize(flag=config.get('JobSettings','test_flag'),
                    ex_data=config.get('JobSettings','ex_data'),
                    root=sim_root,
                    name=config.get('JobSettings','sim_job_base_name'),
                    mat_params=mat_params,
                    varbound=varbound,
                    algorithm_param=algorithm_param,
                    sim_type=config.get('JobSettings','sim_type'),
                    n_jobs = config.get('MainProcessSettings','ntasks'))

    model, func = opt.init_optimizer()
    os.system('echo optimizer initialized starting simulations now')
    if not ast.literal_eval(config.get('JobSettings','restart_flag')):
        model.run(no_plot=True,
                    progress_bar_stream = None,
                    save_last_generation_as = f'{sim_root}/logs/lastgeneration.npz',
                    set_function=ga.set_function_multiprocess(func, n_jobs=ast.literal_eval(config.get('MainProcessSettings','ntasks'))))
    else:
        model.run(no_plot=True,
                progress_bar_stream = None,
                start_generation=f'{sim_root}/logs/lastgeneration.npz',
                set_function=ga.set_function_multiprocess(func, n_jobs=ast.literal_eval(config.get('MainProcessSettings','ntasks'))))

    f = open(sim_root + '/logs/final_results.txt', 'w')
    json.dump(model.output_dict, f, indent=4)


