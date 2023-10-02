import pandas as pd
import matplotlib.pyplot as plt
import os
def plot_data2(fig_name:str, x_label:str, y_label:str, sim_data:pd.DataFrame, ex_data:pd.DataFrame,
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
    plt.savefig(f'{fig_name}.png')
    plt.close()

def create_batch_job_script(sample_files, sim_root, job_name, job_index):
        
    with open(f'{sample_files}/simulation_job.sh','r') as f:
        lines = f.readlines()

    new_lines = [line.replace('%ROOT%', f'{sim_root}') for line in lines]
    new_lines = [line.replace('%JOBNAME%', f'{job_name}_{job_index}') for line in new_lines]
    new_lines = [line.replace('simulation', f'simulation_{job_index}') for line in new_lines]
    current_simulation_dir = f'test'
    with open(f'{current_simulation_dir}/simulation_job_{job_index}.sh','w+') as f:
        f.writelines(new_lines)

if __name__ == "__main__":

    # data = pd.read_csv('/home/rwth1393/single/sample_files/ex_data_single.csv')
    # ex_data = data.copy()

    # ex_data['stress_t'] += 100
    # ex_data['stress_alpha'] += 100
    # ex_data['stress_beta'] += 100

    # sim_y_cols = ['stress_t', 'stress_alpha', 'stress_beta']
    # sim_x_cols = ['strain_t'] * len(sim_y_cols)
    # sim_labels = ['sim_1', 'sim_2', 'sim_3']

    # ex_x_cols = sim_x_cols
    # ex_y_cols = ['stress_t', 'stress_alpha', 'stress_beta']
    # ex_labels = ['ex_1', 'ex_2','ex_3']

    # plot_data2(fig_name='test', x_label='Strain', y_label='Stress(MPa)',
    #            sim_data=data, ex_data=ex_data,
    #            sim_x_cols=sim_x_cols, sim_y_cols=sim_y_cols, sim_labels=sim_labels,
    #            ex_x_cols=ex_x_cols, ex_y_cols=ex_y_cols, ex_labels=ex_labels)

    job_name = "CP_Tensile_5489"
    with open(f'simulation_{job_name[-4:]}/{job_name}.sta','r') as f:
        content = f.read()
        if "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" in content:
            os.system("echo yes")