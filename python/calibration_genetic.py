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
    def __init__(self, sim_root, ex_data, job_name, sim_flag, n_jobs, mat_params: MatParams):
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
        self.sim_flag = sim_flag
        self.sim_flag2compare_function = {
            'cyclic' : self.compare_exp2sim_cyclic,
            'tensile' : self.compare_exp2sim_tensile,
        }
        self.mat_params = mat_params
        self.n_phases = len(mat_params)

    def create_batch_job_script(self, job_index):

        with open(f'{self.sample_files}/simulation_job.sh','r') as f:
            lines = f.readlines()

        new_lines = [line.replace('%ROOT%', f'{self.sim_root}') for line in lines]
        new_lines = [line.replace('%JOBNAME%', f'{self.job_name}_{job_index}') for line in new_lines]
        new_lines = [line.replace('simulation', f'simulation_{job_index}') for line in new_lines]
        current_simulation_dir = f'{self.simulation_dir}_{job_index}'
        with open(f'{current_simulation_dir}/simulation_job_{job_index}.sh','w+') as f:
            f.writelines(new_lines)

    def blackbox_multiphase(self, params):
        submitted = False
        while not submitted:
            # define path variables
            job_index = str(time.time_ns())[-4:]
            current_simulation_dir = f'{self.simulation_dir}_{job_index}'
            current_job_name = f'{self.base_job_name}_{job_index}'
            os.system(f'echo current directory: {current_simulation_dir}')

            if not os.path.isdir(current_simulation_dir):
                #create directory
                self.create_job_dir(current_simulation_dir)

                for j, mat_id in enumerate(self.mat_params.material_ids):
                    if j == 0:
                        path = self.sample_files
                    else:
                        path = current_simulation_dir
                    # create matdata.inp file according to params
                    self.manipulate_matdata(path, current_simulation_dir, mat_id, phase_index=j, values=params)

                # create batch script
                self.create_batch_job_script(job_index=job_index)
                # submit simulation
                self.submit_batch_job(f'{current_simulation_dir}/simulation_job_{job_index}.sh')
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
                sim_results = self.calcStressStrain(current_simulation_dir, current_job_name)
                compare_func = self.sim_flag2compare_function[self.sim_flag]
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
        source_dir = f'{self.sim_root}/sample_files_{self.sim_flag}_simulation'
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
        f.write(f'simulation type: {self.sim_flag}\n')
        f.write(f'Compute Error Using: {compare_func}\n')
        f.write(f'time_Stamp: {time_stamp}\n')
        f.write(f'Error_Time: {mad_time}\n')
        f.write(f'Error_Stress: {mad_stress}\n')
        f.write(f'Error: {mad}\n')
        prop_index = 0
        phase_index = 0
        for mat_props in mat_params.MatID2MatPropsBound.values():
            if self.mat_params.material_ids[phase_index] == 0: #global parameters
                f.write(f'global: ')
            else:
                f.write(f'phase{self.mat_params.material_ids[phase_index]}: ')
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
                os.system(f"echo {e}")
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



    @staticmethod
    def submit_batch_job(batch_file):
        os.system(f'sbatch {batch_file}')


    def manipulate_matdata(self, sample_path:str, current_sim_dir:str, id:int, phase_index:int, values):
        # Manipulate matdata.inp file according to optimizer
        num_prev_props = np.sum(self.num_props)
        mat_props = list(self.mat_params.MatID2MatPropsBound[id].keys())
        num_curr_props = len(mat_props)
        if phase_index == 0:
            mat_props_values = values[:len(mat_props)]
        else:
            mat_props_values = []
            if len(self.mat_params.MatID2MatPropsBound[0]) > 0:
                mat_props = list(self.mat_params.MatID2MatPropsBound[0].keys()) + mat_props #global values
                mat_props_values.extend(values[:self.num_props[1]]) # global values
            mat_props_values.extend(values[num_prev_props:num_prev_props+num_curr_props])


        with open(f'{sample_path}/matdata.inp', 'r') as file:
                lines = file.readlines()
        if id > 0:
            temp_lines = [line.replace(" ", "") for line in lines]
            first_line = temp_lines.index(f'<:Material:{int(id)}:>\n')
            try:
                last_line = temp_lines.index(f'<:Material:{int(id) + 1}:>\n')
            except:
                last_line = len(lines)

            prop_index = 0
            for i in range(first_line, last_line):
                if 'pw_fl' in lines[i] and 'pw_fl' in mat_props:
                    pw_fl = mat_props_values[prop_index]
                    lines[i] = f'pw_fl : {np.round(pw_fl, 4)}\n'
                    prop_index +=1

                elif 'shrt_0' in lines[i] and 'shrt_0' in mat_props:
                    shrt_0 = mat_props_values[prop_index]
                    lines[i] = f'shrt_0 : {np.round(shrt_0, 4)}\n'
                    prop_index += 1

                elif 'hdrt_0' in lines[i] and 'hdrt_0' in mat_props:
                    hdrt_0 = mat_props_values[prop_index]
                    lines[i] = f'hdrt_0 : {np.round(hdrt_0, 4)}\n'
                    prop_index += 1

                elif 'crss_0' in lines[i] and 'crss_0' in mat_props:
                    crss_0 = mat_props_values[prop_index]
                    lines[i] = f'crss_0 : {np.round(crss_0, 4)}\n'
                    prop_index += 1

                elif 'k' in lines[i] and 'k' in mat_props:
                    k = mat_props_values[prop_index]
                    lines[i] = f'k : {np.round(k, 4)}\n'
                    prop_index += 1

                elif 'crss_s' in lines[i] and 'crss_s' in mat_props:
                    crss_s = mat_props_values[prop_index]
                    lines[i] = f'crss_s : {np.round(crss_s, 4)}\n'
                    prop_index += 1

                elif 'pw_hd' in lines[i] and 'pw_hd' in mat_props:
                    pw_hd = mat_props_values[prop_index]
                    lines[i] = f'pw_hd : {np.round(pw_hd, 4)}\n'
                    prop_index += 1

                elif 'Adir' in lines[i] and 'Adir' in mat_props:
                    Adir = mat_props_values[prop_index]
                    lines[i] = f'Adir : {np.round(Adir, 4)}\n'
                    prop_index += 1

                elif 'Adyn' in lines[i] and 'Adyn' in mat_props:
                    Adyn = mat_props_values[prop_index]
                    lines[i] = f'Adyn : {np.round(Adyn, 4)}\n'
                    prop_index += 1

        f = open(f'{current_sim_dir}/matdata.inp', 'w+')
        for line in lines:
            f.write(line)
        f.close()
        self.num_props.append(num_curr_props)

class Optimize:

    def __init__(self, flag, ex_data, root, name, mat_params, varbound, algorithm_param, sim_flag, n_jobs):
        self.test_flag = flag
        self.ex_data = ex_data
        self.sim_root = root
        self.job_name = name
        self.mat_params = mat_params
        self.varbound = varbound
        self.algorithm_param = algorithm_param
        self.sim_flag = sim_flag
        self.n_jobs = n_jobs



    def init_optimizer(self):
        sim_object = Simulation(sim_root=self.sim_root, ex_data=self.ex_data, mat_params=self.mat_params,
                                job_name=self.job_name, sim_flag=self.sim_flag, n_jobs=self.n_jobs)
        if len(self.mat_params) >= 1:
            blackbox_func = sim_object.blackbox_multiphase
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

if __name__ == '__main__':
    # Get the arguments list
    ###############################
    # -t: test_flag   (boolean)
    # -r: restart_flag (boolean)
    # -j: job_name     (string)
    # -d: sim_root   (string)
    # -s: sim_flag   (string) cyclic cyclic loading tensile tensile test
    # -n: n_jobs (int) number of simulations running cocurrently
    # --ex_data: 'filename.csv' (string) Must be located in sample_files
    ###############################
    myopts, args = getopt.getopt(sys.argv[2:], "t:r:s:n:j:d:", ['ex_data='])
    ###############################
    # o == option
    # a == argument passed to the o
    ###############################
    if len(myopts) > 0:
        for o, a in myopts:
            if o == '-t':
                test_flag = json.loads(a.lower())
            elif o == '-r':
                restart_flag = json.loads(a.lower())
            elif o == '-j':
                job_name = str(a)
            elif o == '-d':
                sim_root = str(a)
            elif o == '-s':
                sim_flag = str(a)
                # if sim_flag == 'tensile':
                #     print("run tensile test simulation")
            elif o == '-n':
                n_jobs = int(a)
            elif o == '--ex_data':
                ex_data = str(a)
            else:
                print("Usage: %s -t test_flag -r restart_flag -j job_name -d sim_root --ex_data=filename.csv" % sys.argv[0])
                sys.exit(1)
    else:
        sys.exit()

    #process mat params input file
    mat_params = MatParams(root=sim_root)
    varbound = mat_params.get_varbounds()

    os.system('echo python script initialized with follwing input:')
    os.system(f'echo restart: {restart_flag}')
    os.system(f'echo job_name: {job_name}')
    os.system(f'echo sim_root: {sim_root}')
    os.system(f'echo sim_type: {sim_flag}')

    algorithm_param = {'max_num_iteration': 100, \
                       'population_size': 50, \
                       'mutation_probability': 0.1, \
                       'elit_ratio': 0.1, \
                       'parents_portion': 0.3, \
                       'max_iteration_without_improv': None}

    opt = Optimize(flag=test_flag, ex_data=ex_data, root=sim_root, name=job_name, mat_params=mat_params,
                   varbound=varbound, algorithm_param=algorithm_param, sim_flag=sim_flag, n_jobs = n_jobs)

    model, func = opt.init_optimizer()
    os.system('echo optimizer initialized starting simulations now')
    if restart_flag == False:
        model.run(no_plot=True,
                  progress_bar_stream = None,
                  save_last_generation_as = f'{sim_root}/logs/lastgeneration.npz',
                  set_function=ga.set_function_multiprocess(func, n_jobs=n_jobs))
    else:
        model.run(no_plot=True,
                progress_bar_stream = None,
                start_generation=f'{sim_root}/logs/lastgeneration.npz',
                set_function=ga.set_function_multiprocess(func, n_jobs=n_jobs))
    f = open(sim_root + '/logs/results.txt', 'w+')
    f.write(model.output_dict)
    f.close()
