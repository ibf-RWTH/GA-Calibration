import pandas as pd
import numpy as np

if __name__ == "__main__":
    sample_files = '/home/bq407628/CP-Calibration-Genentic/sample_files'
    experimental_df = pd.read_csv(f'{sample_files}/trial_ex_data.csv', sep=',')
    simulation_df = pd.read_csv('sim_results.csv')

    exp_total_stress_interp = np.interp(simulation_df['Strain'], experimental_df['strain_t'], experimental_df['stress_t'])
    exp_alpha_stress_interp = np.interp(simulation_df['Strain'], experimental_df['strain_t'], experimental_df['stress_alpha'])
    exp_beta_stress_interp = np.interp(simulation_df['Strain'], experimental_df['strain_t'], experimental_df['stress_beta'])
    
    mad_stress_total = np.mean((exp_total_stress_interp - simulation_df['Stress']) / exp_total_stress_interp) * 100
    mad_stress_alpha = np.mean((exp_alpha_stress_interp - simulation_df['Stress']) / exp_alpha_stress_interp) * 100
    mad_stress_beta = np.mean((exp_beta_stress_interp - simulation_df['Stress']) / exp_beta_stress_interp) * 100
    