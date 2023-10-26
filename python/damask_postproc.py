"""
Skript for Damask3 Postprocessing, including macro - micro strain and defgrad comparison for RVE with phase field Damage
For 1 - 3 Phases

Niklas Fehlemann, IMS
"""

import damask
import os
import glob
import sys
import numpy as np
import pandas as pd
# import seaborn as sns
# import pyvista as pv
import yaml
import matplotlib.pyplot as plt

if __name__ == "__main__":

    result_file = '/home/p0021070/damask/test/grid_load.hdf5'
    images = '/home/p0021070/damask/test/evaluation_images'
    phases = ['Alpha', 'Beta']

    result = damask.Result(result_file)
    # result.add_strain(F='F')
    # result.add_stress_Cauchy(F='F')
    # result.add_equivalent_Mises('sigma')
    # result.add_equivalent_Mises('epsilon_V^0.0(F)')


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

    #plot results
    plt.plot(strain_vector, np.asarray(stress_vector)/1e6, label='True Stress') #MPa
    plt.xlabel("True Strain")
    plt.ylabel("True Stress (MPa)")

    for phase in phases:
        plt.plot(phase2strain[phase], np.asarray(phase2stress[phase])/1e6, label=f'Sim_{phase}')#MPa

    plt.legend()
    plt.savefig(f"{images}/test_result.png")




