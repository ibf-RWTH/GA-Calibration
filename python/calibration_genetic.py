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
import damask
import csv
from scipy.interpolate import interp1d

warnings.simplefilter(action='ignore', category=FutureWarning)

def detect_delimiter(csv_path):
    with open(csv_path, 'r') as file:
        dialect = csv.Sniffer().sniff(file.readline())
        return dialect.delimiter

class Simulation:
    def __init__(self, sim_root, ex_data, job_name, n_jobs, mat_params):
        self.sim_root = sim_root
        self.ex_data = ex_data
        self.subroutine_dir = f'{self.sim_root}/subroutine'
        self.simulation_dir = f'{self.sim_root}/simulation'
        # self.sample_files = f'{self.sim_root}/sample_files'
        self.exp_files = f'{self.sim_root}/exp_files'
        self.temp_files = f'{self.sim_root}/template_files'
        self.images = f'{self.sim_root}/evaluation_images_{job_name}'
        self.log_dir = f'{self.sim_root}/logs_{job_name}'
        self.base_job_name = job_name
        self.job_name = job_name
        self.n_jobs = n_jobs
        self.sim_type = 'tensile'
        if Utils.SOFTWARE == "abaqus":
            config = Utils.CONFIG
            sim_type = config.get('AbaqusJobSettings','sim_type')
            self.sim_type = sim_type

        self.mat_params = mat_params
        if self.sim_type == 'Chaboche':
            self.n_phases = len(mat_params)
        else:
            self.n_phases = len(mat_params) - 1

    def create_batch_job_script(self, job_index):
        current_simulation_dir = f'{self.simulation_dir}_{job_index}'
        destination = f'{current_simulation_dir}/simulation_job_{job_index}.sh'
        if 'CP' in self.sim_type:
            source = f'{self.temp_files}/simulation_job_CP.sh'
        elif 'Chaboche' in self.sim_type:
            source = f'{self.temp_files}/simulation_job_Chaboche.sh'
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
                self.create_sim_matdata(current_simulation_dir, params)
                if Utils.DEBUG:
                    if self.sim_type == 'Chaboche':
                        sim_results =pd.read_csv(f'{self.temp_files}/debug_chaboche.csv')
                    else:
                        sim_results =pd.read_csv(f'{self.temp_files}/debug_strain_stress.csv')
                else:
                    self.run_batch_job(current_simulation_dir, job_index, current_job_name)
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
                        if Utils.SOFTWARE == "abaqus":
                            continue
                        # damask simulation doesn't converge
                        elif Utils.SOFTWARE == "damask":
                            time_stamp = int(time.time_ns())
                            self.write_results_file(mad_time=99999, mad_stress=99999, mad=199998, params=params, compare_func="irrelevant", time_stamp=time_stamp)
                            return 199998

                    # evaluate Simulation
                    sim_results = self.get_sim_results(current_simulation_dir, current_job_name)

                mad1, mad2, time_stamp = self.compare_exp2sim(sim_results)
                mad = mad1 + mad2
                # write results to files and delete simulation files
                self.write_results_file(mad_time=mad1, mad_stress=mad2, mad=mad, params=params, compare_func=self.compare_exp2sim.__name__, time_stamp=time_stamp)
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
        software = Utils.CONFIG.get('MainProcessSettings','software')
        if software == 'abaqus':
            if 'CP' in self.sim_type and self.n_phases == 1:
                source_dir = f'{self.sim_root}/{software}/input_files_{str(self.sim_type).lower()}_simulation_single_phase'
            else:
                source_dir = f'{self.sim_root}/{software}/input_files_{str(self.sim_type).lower()}_simulation'
        elif software == 'damask':
            source_dir = f'{self.sim_root}/{software}/input_files_simulation'
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
        # os.system(f"echo {self.mat_params}")
        for phase_id, mat_props in self.mat_params.items():
            if phase_id == 0: #global parameters
                f.write(f'global: ')
            else:
                f.write(f'phase{phase_id}: ')
            for mat_props_name in mat_props:
                f.write(f'{mat_props_name}: {params[prop_index]}, ')
                prop_index += 1

            f.write('\n')

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
            results_df['Sim_Strain'] = flowcurve_df['Strain'].values
            results_df['Sim_Stress'] = flowcurve_df['Stress'].values
            if self.n_phases > 1:
                results_df['Sim_Strain_Phase1'] = flowcurve_df['Strain_Phase1'].values
                results_df['Sim_Stress_Phase1'] = flowcurve_df['Stress_Phase1'].values
                results_df['Sim_Strain_Phase2'] = flowcurve_df['Strain_Phase2'].values
                results_df['Sim_Stress_Phase2'] = flowcurve_df['Stress_Phase2'].values
        else:
            results_df = time_df
            results_df['Sim_Strain'] = 0
            results_df['Sim_Stress'] = 0
            results_df.to_csv(self.sim_root+'results_df.csv')

        return results_df

    def load_csv(self, csv_path):
        sep = detect_delimiter(f'{self.exp_files}/{csv_path}')
        experimental_df = pd.read_csv(f'{self.exp_files}/{csv_path}', sep=sep)
        return experimental_df

    def compare_exp2sim_cyclic(self, simulation_df):
        assert self.n_phases == 1 or self.n_phases == 2, 'more than two phases are not yet supported'
        # in order for this to work correctly make sure the time intervals in experiment and simulation are equal!!!
        now = int(time.time_ns())
        if simulation_df.shape[0] < 5:
            mad_time = 99999.
            mad_stress_strain = 99999
            return mad_time, mad_stress_strain, now


        experimental_df = pd.read_csv(f'{self.exp_files}/{self.ex_data}', sep='\t')

        total_time = experimental_df['time'].values[-1]
        max_sim_time = simulation_df['sim_time'].values[-1]
        mad_time = (100 * (1 - max_sim_time / total_time))**2

        if self.n_phases == 1:
            phase_name = config.get('JobSettings', 'phases')
            phase_id = config.get(phase_name, 'abaqus_id')
            phase_dict = {'3':'Phase1', '4':'Phase2'}
            phase = phase_dict[phase_id]
            sim_stress_interp = np.interp(simulation_df['Sim_Stress'], experimental_df['time'], experimental_df['y'])
            exp_stress_interp = experimental_df['y'].values[:len(sim_stress_interp)]
            mad_stress = np.mean(np.abs(sim_stress_interp - exp_stress_interp)**2)
            sim_y_cols = ['Sim_Stress']
            sim_x_cols = ['Sim_Strain'] * len(sim_y_cols)
            sim_labels = ['Simulation']

            ex_y_cols = [f'Stress_{phase}']
            ex_x_cols = ['strain_t'] * len(ex_y_cols)
            ex_labels = ['Experiment']

            fig_name = f'hysteresis_{now}'
            x_label = "Strain"
            y_label = "Stress(MPa)"
            self.plot_data2(fig_name=fig_name, x_label=x_label,y_label=y_label,
                            sim_data=simulation_df, ex_data=experimental_df,
                            sim_x_cols=sim_x_cols, sim_y_cols=sim_y_cols, sim_labels=sim_labels,
                            ex_x_cols=ex_x_cols, ex_y_cols=ex_y_cols, ex_labels=ex_labels)
            # Plot time sequence
            sim_y_cols = ['Sim_Stress']
            sim_x_cols = ['sim_time'] * len(sim_y_cols)
            sim_labels = ['Simulation']

            ex_y_cols = [f'Stress_{phase}']
            ex_x_cols = ['time'] * len(ex_y_cols)
            ex_labels = ['Experiment']

            fig_name = f'Force_{now}'
            x_label = "Time (s)"
            y_label = "Stress (MPa)"
            self.plot_data2(fig_name=fig_name, x_label=x_label,y_label=y_label,
                            sim_data=simulation_df, ex_data=experimental_df,
                            sim_x_cols=sim_x_cols, sim_y_cols=sim_y_cols, sim_labels=sim_labels,
                            ex_x_cols=ex_x_cols, ex_y_cols=ex_y_cols, ex_labels=ex_labels)

        elif self.n_phases == 2:
            exp_total_stress_interp = np.interp(simulation_df['sim_time'], experimental_df['time'], experimental_df['stress_t'])
            exp_phase1_stress_interp = np.interp(simulation_df['sim_time'], experimental_df['time'], experimental_df['Stress_Phase1'])
            exp_phase2_stress_interp = np.interp(simulation_df['sim_time'], experimental_df['time'], experimental_df['Stress_Phase2'])

            mad_stress_total = np.mean(np.abs(exp_total_stress_interp - simulation_df['Sim_Stress']))
            mad_stress_phase1 = np.mean(np.abs(exp_phase1_stress_interp - simulation_df['Sim_Stress_Phase1']))
            mad_stress_phase2 = np.mean(np.abs(exp_phase2_stress_interp - simulation_df['Sim_Stress_Phase2']))

            mad_stress = (mad_stress_total + 0.8*mad_stress_phase1 + 0.2 * mad_stress_phase2)/3
            # Plot Hysteresis
            sim_y_cols = ['Sim_Stress_Phase1']
            sim_x_cols = ['Sim_Strain'] * len(sim_y_cols)
            sim_labels = ['Simulation']

            ex_y_cols = [f'Stress_Phase1']
            ex_x_cols = ['strain_t'] * len(ex_y_cols)
            ex_labels = ['Experiment']

            fig_name = f'hysteresis_Phase1_{now}'
            x_label = "Strain"
            y_label = "Stress(MPa)"
            self.plot_data2(fig_name=fig_name, x_label=x_label,y_label=y_label,
                            sim_data=simulation_df, ex_data=experimental_df,
                            sim_x_cols=sim_x_cols, sim_y_cols=sim_y_cols, sim_labels=sim_labels,
                            ex_x_cols=ex_x_cols, ex_y_cols=ex_y_cols, ex_labels=ex_labels)
            # Plot Hysteresis
            sim_y_cols = ['Sim_Stress_Phase2']
            sim_x_cols = ['Sim_Strain'] * len(sim_y_cols)
            sim_labels = ['Simulation']

            ex_y_cols = [f'Stress_Phase2']
            ex_x_cols = ['strain_t'] * len(ex_y_cols)
            ex_labels = ['Experiment']

            fig_name = f'hysteresis_Phase2_{now}'
            x_label = "Strain"
            y_label = "Stress(MPa)"
            self.plot_data2(fig_name=fig_name, x_label=x_label,y_label=y_label,
                            sim_data=simulation_df, ex_data=experimental_df,
                            sim_x_cols=sim_x_cols, sim_y_cols=sim_y_cols, sim_labels=sim_labels,
                            ex_x_cols=ex_x_cols, ex_y_cols=ex_y_cols, ex_labels=ex_labels)
            # Plot time sequence
            sim_y_cols = ['Stress_Phase1']
            sim_x_cols = ['sim_time'] * len(sim_y_cols)
            sim_labels = ['Simulation']

            ex_y_cols = [f'Stress_Phase1']
            ex_x_cols = ['time'] * len(ex_y_cols)
            ex_labels = ['Experiment']

            fig_name = f'Stress_Phase1_{now}'
            x_label = "Time (s)"
            y_label = "Stress (MPa)"
            self.plot_data2(fig_name=fig_name, x_label=x_label,y_label=y_label,
                            sim_data=simulation_df, ex_data=experimental_df,
                            sim_x_cols=sim_x_cols, sim_y_cols=sim_y_cols, sim_labels=sim_labels,
                            ex_x_cols=ex_x_cols, ex_y_cols=ex_y_cols, ex_labels=ex_labels)
            # Plot time sequence
            sim_y_cols = ['Stress_Phase2']
            sim_x_cols = ['sim_time'] * len(sim_y_cols)
            sim_labels = ['Simulation']

            ex_y_cols = [f'Stress_Phase2']
            ex_x_cols = ['time'] * len(ex_y_cols)
            ex_labels = ['Experiment']

            fig_name = f'Stress_Phase2_{now}'
            x_label = "Time (s)"
            y_label = "Stress (MPa)"
            self.plot_data2(fig_name=fig_name, x_label=x_label,y_label=y_label,
                            sim_data=simulation_df, ex_data=experimental_df,
                            sim_x_cols=sim_x_cols, sim_y_cols=sim_y_cols, sim_labels=sim_labels,
                            ex_x_cols=ex_x_cols, ex_y_cols=ex_y_cols, ex_labels=ex_labels)

        return mad_time, mad_stress, now

    def compare_exp2sim_tensile(self, simulation_df):
        assert self.n_phases == 2 or self.n_phases == 1
        now = int(time.time_ns())
        if simulation_df.shape[0] < 5:
            mad_time = 99999.
            mad_stress_strain = 99999
            return mad_time, mad_stress_strain, now

        max_sim_strain = simulation_df['Sim_Strain'].max()

        if self.n_phases == 2:

            #load csv files
            ex_data_total = self.ex_data[0].replace(" ", "")
            experimental_total_df = self.load_csv(ex_data_total)
            ex_data_alpha = self.ex_data[1].replace(" ", "")
            experimental_alpha_df = self.load_csv(ex_data_alpha)
            ex_data_beta = self.ex_data[2].replace(" ", "")
            experimental_beta_df = self.load_csv(ex_data_beta)

            max_exp_strain = experimental_total_df['x'].max()

            #interp values
            exp_total_stress_interp = np.interp(simulation_df['Sim_Strain'], experimental_total_df['x'], experimental_total_df['y'])
            if Utils.SOFTWARE == 'abaqus':
                exp_alpha_stress_interp = np.interp(simulation_df['Sim_Strain'], experimental_alpha_df['x'], experimental_alpha_df['y'])
                exp_beta_stress_interp = np.interp(simulation_df['Sim_Strain'], experimental_beta_df['x'], experimental_beta_df['y'])
            elif Utils.SOFTWARE == 'damask':
                exp_alpha_stress_interp = np.interp(simulation_df['Sim_Strain_Phase1'], experimental_alpha_df['x'], experimental_alpha_df['y'])
                exp_beta_stress_interp = np.interp(simulation_df['Sim_Strain_Phase2'], experimental_beta_df['x'], experimental_beta_df['y'])

            mad_stress_total = np.mean(np.abs(exp_total_stress_interp - simulation_df['Sim_Stress']))
            mad_stress_alpha = np.mean(np.abs(exp_alpha_stress_interp - simulation_df['Sim_Stress_Phase1']))
            mad_stress_beta = np.mean(np.abs(exp_beta_stress_interp - simulation_df['Sim_Stress_Phase2']))
            mad_strain_total = (abs(1 - max_sim_strain / max_exp_strain) * 100) **2
            mad_stress = (mad_stress_total + 0.8*mad_stress_alpha + 0.2 * mad_stress_beta)/3

            fig_name = f'stress_strain_{now}'
            x_label = "Strain"
            y_label = "Stress(MPa)"

            #plot total stress strain
            self.plot_data(fig_name=fig_name, x_label=x_label,y_label=y_label,
                                             x1=simulation_df['Sim_Strain'].values, y1=simulation_df['Sim_Stress'].values, data_label_1='Sim_Total',
                                             x2=experimental_total_df['x'].values, y2=experimental_total_df['y'].values, data_label_2='Exp_Total')
            #plot alpha stress strain
            fig_name = f'stress_strain_alpha_{now}'
            self.plot_data(fig_name=fig_name, x_label=x_label,y_label=y_label,
                                             x1=simulation_df['Sim_Strain'].values, y1=simulation_df['Sim_Stress_Phase1'].values, data_label_1='Sim_Phase1',
                                             x2=experimental_alpha_df['x'].values, y2=experimental_alpha_df['y'].values, data_label_2='Exp_Phase1')
            #plot beta stress strain
            fig_name = f'stress_strain_beta_{now}'
            self.plot_data(fig_name=fig_name, x_label=x_label,y_label=y_label,
                                             x1=simulation_df['Sim_Strain'].values, y1=simulation_df['Sim_Stress_Phase2'].values, data_label_1='Sim_Phase2',
                                             x2=experimental_beta_df['x'].values, y2=experimental_beta_df['y'].values, data_label_2='Exp_Phase2')

        if self.n_phases == 1:
            ex_data = self.ex_data[0].replace(" ", "") #remove spaces
            experimental_df = self.load_csv(ex_data)
            max_exp_strain = experimental_df['x'].max()
            exp_stress_interp = np.interp(simulation_df['Sim_Strain'], experimental_df['x'], experimental_df['y'])
            mad_stress = np.mean(np.abs(exp_stress_interp - simulation_df['Sim_Stress']) / exp_stress_interp) * 100
            mad_strain_total = (abs(1 - max_sim_strain / max_exp_strain) * 100) **2

            sim_y_cols = ['Sim_Stress']
            sim_x_cols = ['Sim_Strain'] * len(sim_y_cols)
            sim_labels = ['Simulation']

            ex_y_cols = ['y']
            ex_x_cols = ['x'] * len(ex_y_cols)
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
        now = int(time.time_ns())
        if simulation_df.shape[0] < 5:
            mad_time = 99999.
            mad_stress_strain = 99999
            return mad_time, mad_stress_strain, now


        experimental_df = pd.read_csv(f'{self.exp_files}/{self.ex_data}')

        total_time = experimental_df['time'].values[-1]
        max_sim_time = simulation_df['sim_time'].values[-1]
        mad_time = (100 * (1 - max_sim_time / total_time))**2
        comp_df = experimental_df.merge(simulation_df, left_on='time',
                                        right_on='sim_time')  # merge dfs on time for comparison

        mad_force = np.mean(np.abs(comp_df['force'] - comp_df['sim_force'])*1000)

        # Plot time sequence
        sim_y_cols = ['sim_force']
        sim_x_cols = ['sim_time'] * len(sim_y_cols)
        sim_labels = ['Simulation']

        ex_y_cols = ['force']
        ex_x_cols = ['time'] * len(ex_y_cols)
        ex_labels = ['Experiment']

        fig_name = f'Force_{now}'
        x_label = "Time (s)"
        y_label = "Force (kN)"
        self.plot_data2(fig_name=fig_name, x_label=x_label,y_label=y_label,
                        sim_data=simulation_df, ex_data=experimental_df,
                        sim_x_cols=sim_x_cols, sim_y_cols=sim_y_cols, sim_labels=sim_labels,
                        ex_x_cols=ex_x_cols, ex_y_cols=ex_y_cols, ex_labels=ex_labels)

        # Plot Hysteresis
        sim_x_cols = ['sim_displacement'] * len(sim_y_cols)
        ex_x_cols = ['displacement'] * len(ex_y_cols)
        fig_name = f'Hysteresis_{now}'
        x_label = 'Displacement (mm)'

        self.plot_data2(fig_name=fig_name, x_label=x_label,y_label=y_label,
                        sim_data=simulation_df, ex_data=experimental_df,
                        sim_x_cols=sim_x_cols, sim_y_cols=sim_y_cols, sim_labels=sim_labels,
                        ex_x_cols=ex_x_cols, ex_y_cols=ex_y_cols, ex_labels=ex_labels)
        return mad_time, mad_force, now

    def compare_exp2sim_damask(self, simulation_df):
        assert self.n_phases == 2 or self.n_phases == 1
        now = int(time.time_ns())
        if simulation_df.shape[0] < 5:
            mad_time = 99999.
            mad_stress_strain = 99999
            return mad_time, mad_stress_strain, now

        experimental_df = pd.read_csv(f'{self.exp_files}/{self.ex_data}', sep=',')
        max_exp_strain = experimental_df['strain_t'].max()
        max_total_strain = simulation_df['total_strain'].max()

        exp_total_stress_interp = np.interp(simulation_df['total_strain'], experimental_df['strain_t'], experimental_df['stress_t'])
        exp_alpha_stress_interp = np.interp(simulation_df['Alpha_strain'], experimental_df['strain_t'], experimental_df['stress_alpha'])
        exp_beta_stress_interp = np.interp(simulation_df['Beta_strain'], experimental_df['strain_t'], experimental_df['stress_beta'])

        mad_stress_total = np.mean(np.abs(exp_total_stress_interp - simulation_df['total_stress']))
        mad_stress_alpha = np.mean(np.abs(exp_alpha_stress_interp - simulation_df['Alpha_stress']))
        mad_stress_beta = np.mean(np.abs(exp_beta_stress_interp - simulation_df['Beta_stress']))
        mad_strain_total = (abs(1 - max_total_strain / max_exp_strain) * 100) **2
        mad_stress = (mad_stress_total + 0.8*mad_stress_alpha + 0.2 * mad_stress_beta)/3

        sim_y_cols = ['total_stress', 'Alpha_stress', 'Beta_stress']
        sim_x_cols = ['total_strain','Alpha_strain','Beta_strain']
        sim_labels = ['Simulation_Total', 'Simulation_Alpha', 'Simulation_Beta']

        ex_y_cols = ['stress_t', 'stress_alpha', 'stress_beta']
        ex_x_cols = ['strain_t'] * len(ex_y_cols)
        ex_labels = ['Experiment_Total', 'Experiment_Alpha','Experiment_Beta']

        fig_name = f'stress_strain_{now}'
        x_label = "Strain"
        y_label = "Stress(MPa)"
        self.plot_data2(fig_name=fig_name, x_label=x_label,y_label=y_label,
                        sim_data=simulation_df, ex_data=experimental_df,
                        sim_x_cols=sim_x_cols, sim_y_cols=sim_y_cols, sim_labels=sim_labels,
                        ex_x_cols=ex_x_cols, ex_y_cols=ex_y_cols, ex_labels=ex_labels)
        return mad_strain_total, mad_stress, now

    def interp(self, x1, y1, x2, y2):
        # define interp function
        interp_func1 = interp1d(x1, y1, kind='linear', fill_value='extrapolate')
        interp_func2 = interp1d(x2, y2, kind='linear', fill_value='extrapolate')

        # define common x range
        common_x = np.linspace(min(np.min(x1), np.min(x2)), max(np.max(x1), np.max(x2)), num=max(len(x1), len(x2)))

        # interpolate
        interp_y1 = interp_func1(common_x)
        interp_y2 = interp_func2(common_x)

        return interp_y1, interp_y2

    def plot_results(self, now, simulation_df, experiment_dfs):

        if self.sim_type == "Chaboche":
            fig_name = f'Force_{now}'
            x_label = "Time (s)"
            y_label = "Force (kN)"
            self.plot_data(fig_name=fig_name, x_label=x_label,y_label=y_label,
                            x1=simulation_df['sim_time'].values, y1=simulation_df['sim_force'].values, data_label_1='Simulation',
                            x2=experiment_dfs[0]['time'].values, y2=experiment_dfs[0]['y'].values, data_label_2='Experiment')

            fig_name = f'Hysteresis_{now}'
            x_label = 'Displacement (mm)'
            self.plot_data(fig_name=fig_name, x_label=x_label,y_label=y_label,
                            x1=simulation_df['sim_displacement'].values, y1=simulation_df['sim_force'].values, data_label_1='Simulation',
                            x2=experiment_dfs[0]['x'].values, y2=experiment_dfs[0]['y'].values, data_label_2='Experiment')
        else:
            fig_names = [f'stress_strain_{now}', f'stress_strain_phase1_{now}', f'stress_strain_phase2_{now}']
            labels = ['Total', 'Phase1', 'Phase2']
            sim_xs = ['Sim_Strain', 'Sim_Strain_Phase1', 'Sim_Strain_Phase2']
            sim_ys = ['Sim_Stress', 'Sim_Stress_Phase1', 'Sim_Stress_Phase2']
            x_label = "Strain"
            y_label = "Stress(MPa)"

            n = 1 if self.n_phases == 1 else 3
            for i in range(n):
                self.plot_data(fig_name=fig_names[i], x_label=x_label,y_label=y_label,
                            x1=simulation_df[sim_xs[i]].values, y1=simulation_df[sim_ys[i]].values, data_label_1=f'Sim_{labels[i]}',
                            x2=experiment_dfs[i]['x'].values, y2=experiment_dfs[i]['y'].values, data_label_2=f'Exp_{labels[i]}')

                if self.sim_type == 'CP_Cyclic':
                    fig_names = [f'stress_{now}', f'stress_phase1_{now}', f'stress_phase2_{now}']
                    labels = ['Total', 'Phase1', 'Phase2']
                    sim_ys = ['Sim_Stress', 'Sim_Stress_Phase1', 'Sim_Stress_Phase2']
                    x_label = "Time (s)"
                    y_label = "Stress (MPa)"
                    self.plot_data(fig_name=fig_names[i], x_label=x_label,y_label=y_label,
                                x1=simulation_df[sim_xs[i]].values, y1=simulation_df[sim_ys[i]].values, data_label_1=f'Sim_{labels[i]}',
                                x2=experiment_dfs[i]['x'].values, y2=experiment_dfs[i]['y'].values, data_label_2=f'Exp_{labels[i]}')

    def compare_exp2sim(self, simulation_df):
        assert self.n_phases == 1 or self.n_phases == 2, 'more than two phases are not yet supported'
        # in order for this to work correctly make sure the time intervals in experiment and simulation are equal!!!
        now = int(time.time_ns())
        if simulation_df.shape[0] < 5:
            mad_time = 99999.
            mad_stress_strain = 99999
            return mad_time, mad_stress_strain, now

        # define independent variables
        if 'tensile' in self.sim_type:
            exp_xs = 'x'
            sim_xs = 'Sim_Strain'
        else:
            exp_xs = 'time'
            sim_xs = 'sim_time'

        if self.sim_type == 'Chaboche':
            sim_ys = 'sim_force'
        else:
            sim_ys = 'Sim_Stress'

        max_sim_xs = simulation_df[sim_xs].max()

        if self.n_phases == 1 or self.sim_type == 'Chaboche':
            ex_data = self.ex_data[0].replace(" ", "") #remove spaces
            experimental_df = self.load_csv(ex_data)
            max_exp_xs = experimental_df[exp_xs].max()
            mad_xs = (100 * (1 - max_sim_xs / max_exp_xs))**2

            exp_interp_y, sim_interp_y = self.interp(experimental_df[exp_xs].values, experimental_df['y'].values,
                                                     simulation_df[sim_xs].values,simulation_df[sim_ys].values)

            mad_ys = np.mean(np.abs((exp_interp_y - sim_interp_y) / exp_interp_y))* 100

            self.plot_results(now, simulation_df, [experimental_df])

        elif self.n_phases == 2:
             #load csv files
            ex_data_total = self.ex_data[0].replace(" ", "")
            experimental_total_df = self.load_csv(ex_data_total)
            ex_data_alpha = self.ex_data[1].replace(" ", "")
            experimental_alpha_df = self.load_csv(ex_data_alpha)
            ex_data_beta = self.ex_data[2].replace(" ", "")
            experimental_beta_df = self.load_csv(ex_data_beta)

            max_exp_xs = experimental_total_df[exp_xs].max()
            mad_xs = (100 * (1 - max_sim_xs / max_exp_xs))**2

            #interp values
            exp_interp_total_y, sim_interp_total_y = self.interp(experimental_total_df[exp_xs].values, experimental_total_df['y'].values,
                                                                  simulation_df[sim_xs].values,simulation_df[sim_ys].values)

            if Utils.SOFTWARE == 'damask' and 'tensile' in self.sim_type:
                exp_interp_alpha_y, sim_interp_alpha_y = self.interp(experimental_alpha_df[exp_xs].values, experimental_alpha_df['y'].values,
                                                                  simulation_df[f'{sim_xs}_Phase1'].values,simulation_df[f'{sim_ys}_Phase1'].values)

                exp_interp_beta_y, sim_interp_beta_y = self.interp(experimental_beta_df[exp_xs].values, experimental_beta_df['y'].values,
                                                                  simulation_df[f'{sim_xs}_Phase2'].values,simulation_df[f'{sim_ys}_Phase2'].values)
            else:
                exp_interp_alpha_y, sim_interp_alpha_y = self.interp(experimental_alpha_df[exp_xs].values, experimental_alpha_df['y'].values,
                                                                  simulation_df[sim_xs].values,simulation_df[f'{sim_ys}_Phase1'].values)

                exp_interp_beta_y, sim_interp_beta_y = self.interp(experimental_beta_df[exp_xs].values, experimental_beta_df['y'].values,
                                                                  simulation_df[sim_xs].values,simulation_df[f'{sim_ys}_Phase2'].values)

            mad_ys_total = np.mean(np.abs((exp_interp_total_y - sim_interp_total_y) / exp_interp_total_y))* 100
            mad_ys_alpha = np.mean(np.abs((exp_interp_alpha_y - sim_interp_alpha_y) / exp_interp_alpha_y))* 100
            mad_ys_beta = np.mean(np.abs((exp_interp_total_y - sim_interp_alpha_y) / exp_interp_beta_y))* 100
            mad_ys = (mad_ys_total + 0.8*mad_ys_alpha + 0.2 * mad_ys_beta) / 3

            self.plot_results(now, simulation_df, [experimental_total_df, experimental_alpha_df, experimental_beta_df])

        return mad_xs, mad_ys, now

    def plot_data(self, fig_name, x_label, y_label, x1, y1, data_label_1, x2=None, y2=None, data_label_2=None):
        plt.plot(x1, y1, label=data_label_1)

        if x2 is not None:
            assert (y2 is not None)
            assert (data_label_2 is not None)
            plt.plot(x2, y2, label=data_label_2)
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
        x_min = min(ex_data[x].values)*1.2
        x_max = max(ex_data[x].values)*1.2
        y_min = min(ex_data[y].values)*1.2
        y_max = max(ex_data[y].values)*1.2
        plt.xlim(x_min,x_max)
        plt.ylim(y_min,y_max)
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.legend()
        plt.savefig(f'{self.images}/{fig_name}.png')
        plt.close()

    def submit_batch_job(self, currentDir, file, jobName):
        sim_type = self.sim_type
        config = Utils.CONFIG
        account = config.get('AbaqusJobSettings','account')
        sim_job_base_name = config.get('AbaqusJobSettings','sim_job_base_name')
        logs_path = f'{self.sim_root}/logs_{sim_job_base_name}'
        output = logs_path +f'/{jobName}-log.%J'
        time = config.get('AbaqusJobSettings','time')
        memPerCpu = config.get('AbaqusJobSettings','mem-per-cpu')
        nodes = config.get('AbaqusJobSettings','nodes')
        nTasks = config.get('AbaqusJobSettings','ntasks')

        cmd = f'sbatch --job-name={jobName}'
        if account != 'None':
            cmd += f' --account={account}'
        cmd += f' --output={output} --time={time}'
        cmd += f' --mem-per-cpu={memPerCpu} --nodes={nodes} --ntasks={nTasks}  {file}'

        abaqus = config.get('AbaqusJobSettings','abaqus_version')
        memory = config.get('AbaqusJobSettings','ABAQUS_MEM_ARG')
        threads = config.get('AbaqusJobSettings','THREADS_PER_MPI')
        root = config.get('AbaqusJobSettings','ROOT')
        if 'CP' in sim_type:
            pythonPath = config.get('AbaqusJobSettings','PYTHON_PATH')+'/readOdb_CP.py'
        elif 'Chaboche' in sim_type:
            pythonPath = config.get('AbaqusJobSettings','PYTHON_PATH')+'/readOdb_Chaboche.py'
        subroutine = config.get('AbaqusJobSettings','SUBROUTINE_PATH')
        pythonPath = root + '/' + pythonPath
        subroutine = root + '/' + subroutine
        input = config.get('AbaqusJobSettings', 'INPUTFILE')
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
            phase = config.get('AbaqusJobSettings','phases').split(',')[phase_index - 1]
            mat_props_keys = mat_params[0].copy() # Global Parameters keys/names !!! result in mutual change here mat_props_keys = mat_params[0] then extend also results in change of mat_params[0]
            mat_props_keys.extend(mat_params[id]) # adding Parameter keys/names for current phase
            # num_curr_props = len(mat_props)
            if phase_index == 1:
                mat_props_values = optiParams[:len(mat_props_keys)]
            else:
                mat_props_values = []
                if len(mat_params[0]) > 0:
                    mat_props_values.extend(optiParams[:len(mat_params[0])])
                mat_props_values.extend(optiParams[-(len(mat_props_keys) - len(mat_params[0])):])

            with open(f'{sample_path}/matdata.inp', 'r') as file:
                    lines = file.readlines()

            prop_index_var = 0
            for gc in constant_params[0]:
                value = constant_params[0][gc]
                lines = [line.replace(f'%{gc}%{id}%', f'{np.round(value,4)}') for line in lines]

            for option in Utils.CONFIG.options(phase):
                if option == 'abaqus_id':
                    continue
                if option in mat_props_keys:
                    value = mat_props_values[prop_index_var]
                    prop_index_var +=1
                    lines = [line.replace(f'%{option}%{id}%', f'{np.round(value,4)}') for line in lines]
                elif option in constant_params[id]:
                    value = constant_params[id][option]
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
                    lines = [line.replace(f'%{option}%', f'{np.round(value,4)}') for line in lines]
                elif option in constant_params[id]:
                    value = constant_params[id][option]
                    prop_index_const +=1
                    lines = [line.replace(f'%{option}%', f'{np.round(value,4)}') for line in lines]

            f = open(f'{current_sim_dir}/Material.inp','w+')
            f.writelines(lines)
            f.close()

    def create_sim_matdata(self, current_simulation_dir, params):

        if Utils.SOFTWARE == "abaqus":
            if 'CP' in self.sim_type:
                for j, mat_id in enumerate(self.mat_params.keys()):
                    if j > 0:
                        if j == 1:
                            path = self.temp_files
                        else:
                            path = current_simulation_dir
                        # create matdata.inp file according to params
                        self.manipulate_matdata(path, current_simulation_dir, mat_id, phase_index=j, optiParams=params)
            elif 'Chaboche' in self.sim_type:
                for j, mat_id in enumerate(self.mat_params.keys()):
                    self.manipulate_matdata(self.temp_files, current_simulation_dir, mat_id, phase_index=j+1, optiParams=params)

        elif Utils.SOFTWARE == "damask":
            m = damask.ConfigMaterial()
            m = m.load(f'{self.sim_root}/damask/input_files_simulation/material.yaml')
            phases = config.get('DamaskJobSettings','phases').split(',')
            constantParams = Utils.CONSTANTPARAMS
            mat_params = Utils.EVALUATINGPARAMS

            phase_index = 1
            for phase in phases:
                prop_index = 0
                mat_props_keys = mat_params[0].copy() # Global Parameters keys/names !!! result in mutual change here mat_props_keys = mat_params[0] then extend also results in change of mat_params[0]
                mat_props_keys.extend(mat_params[phase_index])
                phase_data = m.get('phase')[phase]
                phase_id = ast.literal_eval(config.get(phase, 'phase_id'))

                # config global constant params
                for key, item in constantParams[0].items():
                    if key in phase_data['mechanical']['plastic'].keys():
                        phase_data['mechanical']['plastic'][key] = item
                    elif key in phase_data['mechanical']['elastic'].keys():
                        phase_data['mechanical']['elastic'][key] = item
                #config phase constant params
                for key, item in constantParams[phase_id].items():
                    if key in phase_data['mechanical']['plastic'].keys():
                        phase_data['mechanical']['plastic'][key] = item
                    elif key in phase_data['mechanical']['elastic'].keys():
                        phase_data['mechanical']['elastic'][key] = item
                    elif key == 'lattice':
                        phase_data[key] = item

                #config evaluation params

                if phase_index == 1:
                    mat_props_values = params[:len(mat_props_keys)]
                else:
                    mat_props_values = []
                    if len(mat_params[0]) > 0:
                        mat_props_values.extend(params[:len(mat_params[0])])
                    mat_props_values.extend(params[-(len(mat_props_keys) - len(mat_params[0])):])

                if phase_data['lattice'] == 'cF':
                    duplicate = 1
                elif phase_data['lattice'] == 'cI':
                    duplicate = 2


                for key in mat_props_keys:
                    if  key == 'xi_0_sl'  or key == 'xi_inf_sl':
                        phase_data['mechanical']['plastic'][key] = [float(mat_props_values[prop_index])] * duplicate
                    else:
                        phase_data['mechanical']['plastic'][key] = float(mat_props_values[prop_index])
                    prop_index += 1

                phase_index += 1

            m.save(f'{current_simulation_dir}/material.yaml')

    def run_batch_job(self, current_simulation_dir, job_index, current_job_name):
        software = Utils.SOFTWARE
        if software == "abaqus":
            self.create_batch_job_script(job_index=job_index)
            self.submit_batch_job(current_simulation_dir, f'simulation_job_{job_index}.sh', current_job_name)

        elif software == 'damask':
            config = Utils.CONFIG
            account = config.get('DamaskJobSettings','account')
            sim_job_base_name = config.get('DamaskJobSettings','sim_job_base_name')
            logs_path = f'{self.sim_root}/logs_{sim_job_base_name}'
            output = logs_path +f'/{current_job_name}-log.%J'
            run_time = config.get('DamaskJobSettings','time')
            memPerCpu = config.get('DamaskJobSettings','mem-per-cpu')
            nodes = config.get('DamaskJobSettings','nodes')
            cpus_per_task = config.get('DamaskJobSettings','cpus-per-task')
            main_cwd = os.getcwd()
            os.chdir(current_simulation_dir)
            os.system(f'sbatch --job-name={current_job_name} --output={output} --time={run_time} --nodes={nodes} --account={account} --mem-per-cpu={memPerCpu} --cpus-per-task={cpus_per_task} batch.sh')
            os.chdir(main_cwd)

    def damask_postproc(self, current_simulation_dir):
        result_file = f'{current_simulation_dir}/grid_load.hdf5'
        phases = Utils.CONFIG.get('DamaskJobSettings','phases').split(',')

        result = damask.Result(result_file)
        result.add_strain(F='F')
        result.add_stress_Cauchy(F='F')
        result.add_equivalent_Mises('sigma')
        result.add_equivalent_Mises('epsilon_V^0.0(F)')


        # Averaged Quantities
        strain_eq = result.place('epsilon_V^0.0(F)')
        stress_eq = result.place('sigma')

        strain_vector = list()
        stress_vector = list()

        for _, strain in strain_eq.items():
            strain_mean = np.abs(strain[:, 2, 2].mean())
            strain_vector.append(strain_mean)

        for _, stress in stress_eq.items():
            stress_mean = np.abs(stress[:, 2, 2].mean()) #Pa
            stress_vector.append(stress_mean)

        # Partitioned Flowcurves
        phase2stress = {}
        phase2strain = {}
        for phase in phases:
            phase2stress[phase] = []
            phase2strain[phase] = []

        for key, val in result.get('sigma').items():
            for subkey, subvalue in val.items():
                for phase in phases:
                    if phase in subkey:
                        phase2stress[phase].append(np.mean(subvalue[:,2,2])) #Pa

        for key, val in result.get('epsilon_V^0.0(F)').items():
            for subkey, subvalue in val.items():
                for phase in phases:
                    if phase in subkey:
                        phase2strain[phase].append(np.mean(subvalue[:,2,2]))

        data = {}
        data['Sim_Strain'] = np.array(strain_vector)
        data['Sim_Stress'] = np.array(stress_vector) /1e6 #MPa
        phase_index = 1
        for phase in phases:
            data[f'Sim_Strain_Phase{phase_index}'] = np.array(phase2strain[phase])
            data[f'Sim_Stress_Phase{phase_index}'] = np.array(phase2stress[phase])/ 1e6 #MPa
            phase_index += 1

        sim_results = pd.DataFrame(data)
        return sim_results

    def get_sim_results(self,current_simulation_dir, current_job_name):
        if Utils.SOFTWARE == 'abaqus':
            if 'CP' in self.sim_type:
                    sim_results = self.calcStressStrain(current_simulation_dir, current_job_name)
            elif 'Chaboche' in self.sim_type:
                colNames = ['sim_time', 'sim_displacement', 'sim_force']
                sim_results = pd.read_csv(current_simulation_dir+'/RF_data.txt', names=colNames)

        elif Utils.SOFTWARE == 'damask':
            sim_results = self.damask_postproc(current_simulation_dir)

        return sim_results


class Optimize:

    def __init__(self, flag, ex_data, root, name, mat_params, varbound, algorithm_param, n_jobs, software):
        self.test_flag = flag
        self.ex_data = ex_data
        self.sim_root = root
        self.job_name = name
        self.mat_params = Utils.EVALUATINGPARAMS
        self.varbound = varbound
        self.algorithm_param = algorithm_param
        self.n_jobs = n_jobs
        self.software = software

    def init_optimizer(self):
        sim_object = Simulation(sim_root=self.sim_root, ex_data=self.ex_data, mat_params=self.mat_params,
                                job_name=self.job_name, n_jobs=self.n_jobs)

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
    def __init__(self, JobSettings) -> None:
        self.JobSettings = JobSettings

    def parse_matparams(self):
        if self.JobSettings == 'AbaqusJobSettings':
          return self.parse_matparams_abaqus()
        elif self.JobSettings == 'DamaskJobSettings':
          return self.parse_matparams_damask()

    def parse_matparams_abaqus(self):
        config = Utils.CONFIG
        varbound = []
        constantParams = {}
        matparams = {}
        params_names = []
        phases = config.get(self.JobSettings,'phases').split(',')
        sim_type = config.get(JobSettings,'sim_type')
        if 'CP' in sim_type:
            global_params = config.options('GlobalParams')

            # add global params to varbound
            constantParams[0] = {}
            for global_param in global_params:
                global_param_value = config.get('GlobalParams',global_param)
                if global_param != "phase_id":
                    # add global params
                    if global_param_value.startswith('['):
                        global_param_value = ast.literal_eval(global_param_value)
                        varbound.append(global_param_value)
                        params_names.append(global_param)
                    else:
                        global_param_value = ast.literal_eval(global_param_value)
                        constantParams[0][global_param] = global_param_value
            matparams[0] = params_names
            params_names = []

            for phase in phases:
                phase_id = ast.literal_eval(config.get(phase, 'phase_id'))
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
                        if 'd' in value:
                            value = value.replace('d','e') # subroutine input file sytax

                        if option != "lattice" and option != "type":
                            value = ast.literal_eval(value)
                        constantParams[phase_id][option] = value

                matparams[phase_id] = params_names
                params_names = []

            varbound = np.array(varbound)

        elif 'Chaboche' in sim_type:
            phase_id = config.get(phases[0], 'phase_id')
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
                    if 'd' in value:
                        value = value.replace('d','e') # subroutine input file sytax
                    value = ast.literal_eval(value)
                    constantParams[phase_id][option] = value


            matparams[phase_id] = params_names
            varbound = np.array(varbound)

        return constantParams, matparams, varbound

    def parse_matparams_damask(self):
      config = Utils.CONFIG
      varbound = []
      constantParams = {}
      matparams = {}
      params_names = []
      phases = config.get(self.JobSettings,'phases').split(',')

      global_params = config.options('GlobalParams')
      # add global params to varbound
      constantParams[0] = {}
      for global_param in global_params:
          global_param_value = config.get('GlobalParams',global_param)
          if global_param != "phase_id":
              # add global params
              if global_param_value.startswith('['):
                  global_param_value = ast.literal_eval(global_param_value)
                  varbound.append(global_param_value)
                  params_names.append(global_param)
              else:
                  global_param_value = ast.literal_eval(global_param_value)
                  constantParams[0][global_param] = global_param_value
      matparams[0] = params_names
      params_names = []

      # add phase settings
      for phase in phases:
        phase_id =  ast.literal_eval(config.get(phase, 'phase_id'))
        params_names.extend(config.get(self.JobSettings, f'{phase}.variables').split(','))
        matparams[phase_id] = params_names
        constantParams[phase_id] = {}
        options = config.options(phase)
        phase_varbound = ['nan'] * len(params_names)
        for option in options:
            value = config.get(phase, option)
            if option in params_names:
                # get index of option in params_names
                opt_index = params_names.index(option)
                value = ast.literal_eval(value)
                # set value on the opt_index
                phase_varbound[opt_index] = value
            else:
                if option != "lattice" and option != "type":
                    value = ast.literal_eval(value)
                constantParams[phase_id][option] = value

        varbound.extend(phase_varbound)
        params_names = []

      os.system(f'echo varbound is: {varbound}')
      varbound = np.array(varbound)
      return constantParams, matparams, varbound

class Utils:
    CONFIG = None
    CONSTANTPARAMS = None
    EVALUATINGPARAMS = None
    SOFTWARE = None
    DEBUG = False

if __name__ == '__main__':
    #get current path
    sim_root = os.getcwd()
    os.system(f"echo sim root: {sim_root}")

    #config simulation
    config = configparser.ConfigParser(allow_no_value=True)
    config.optionxform = lambda option : option # Preserve case of the keys
    config.read(f'{sim_root}/configs/configs.ini')
    Utils.CONFIG = config
    software = config.get('MainProcessSettings','software')
    Utils.SOFTWARE = software
    JobSettings = 'AbaqusJobSettings' if software == 'abaqus' else 'DamaskJobSettings'

    ex_data=config.get(JobSettings,'ex_data').split(',')
    os.system(f'echo ex_data are: {ex_data}')
    assert len(ex_data) == 1 or len(ex_data) == 3, "wrong number of input files!"

    constant_params, mat_params, varbound = Parser(JobSettings).parse_matparams()
    Utils.CONSTANTPARAMS = constant_params
    Utils.EVALUATINGPARAMS = mat_params
    Utils.DEBUG = config.getboolean('MainProcessSettings','debug')
    os.system(f'echo debug mode: {Utils.DEBUG}')
    # read algortithm Parameters from config
    max_iters = ast.literal_eval(config.get('AlgorithmParameters','max_iters'))
    population_size = ast.literal_eval(config.get('AlgorithmParameters','population_size'))
    mutation_probability = ast.literal_eval(config.get('AlgorithmParameters','mutation_probability'))
    elit_ratio = ast.literal_eval(config.get('AlgorithmParameters','elit_ratio'))
    parents_portion = ast.literal_eval(config.get('AlgorithmParameters','parents_portion'))
    max_iteration_without_improv = ast.literal_eval(config.get('AlgorithmParameters','max_iteration_without_improv'))

    #initialize Optimizer
    algorithm_param = {'max_num_iteration': max_iters, \
                        'population_size': population_size, \
                        'mutation_probability': mutation_probability, \
                        'elit_ratio': elit_ratio, \
                        'parents_portion': parents_portion, \
                        'max_iteration_without_improv': max_iteration_without_improv}

    opt = Optimize(flag=config.get(JobSettings,'test_flag'),
                    ex_data=ex_data,
                    root=sim_root,
                    name=config.get(JobSettings,'sim_job_base_name'),
                    mat_params=mat_params,
                    varbound=varbound,
                    algorithm_param=algorithm_param,
                    n_jobs = config.get('MainProcessSettings','ntasks'),
                    software = software)

    model, func = opt.init_optimizer()
    os.system('echo optimizer initialized starting simulations now')
    name=config.get(JobSettings,'sim_job_base_name')
    if not ast.literal_eval(config.get(JobSettings,'restart_flag')):
        result = model.run(no_plot=True,
                    progress_bar_stream = None,
                    save_last_generation_as = f'{sim_root}/logs_{name}/lastgeneration.npz',
                    set_function=ga.set_function_multiprocess(func, n_jobs=ast.literal_eval(config.get('MainProcessSettings','ntasks'))))
    else:
        result = model.run(no_plot=True,
                progress_bar_stream = None,
                start_generation=f'{sim_root}/logs_{name}/lastgeneration.npz',
                set_function=ga.set_function_multiprocess(func, n_jobs=ast.literal_eval(config.get('MainProcessSettings','ntasks'))))

    f = open(sim_root + f'/logs_{name}/results.txt', 'w')
    f.write('best found solution: \n')
    prop_index = 0
    for phase_id, mat_props in mat_params.items():
        if phase_id == 0: #global parameters
            f.write(f'global: ')
        else:
            f.write(f'phase{phase_id}: ')
        for mat_props_name in mat_props:
            f.write(f'{mat_props_name}: {result.variable[prop_index]}, ')
            prop_index += 1

        f.write('\n')
    f.close()


