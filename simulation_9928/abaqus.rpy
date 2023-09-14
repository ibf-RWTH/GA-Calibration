# -*- coding: mbcs -*-
#
# Abaqus/CAE Release 2022.HF6 replay file
# Internal Version: 2022_11_16-21.42.28 RELr424 176889
# Run by bq407628 on Wed Sep 13 16:08:18 2023
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
from abaqus import *
from abaqusConstants import *
session.Viewport(name='Viewport: 1', origin=(1.36719, 1.36719), width=201.25, 
    height=135.625)
session.viewports['Viewport: 1'].makeCurrent()
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
execfile('/home/rwth1393/GA-Calibration/python/readOdb.py', __main__.__dict__)
#: /home/rwth1393/GA-Calibration/simulation_9928
#: CP_Calibration_9928
#: Model: /home/rwth1393/GA-Calibration/simulation_9928/CP_Calibration_9928.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       155
#: Number of Node Sets:          738
#: Number of Steps:              1
#: data_Stress.csv has been successfully created
#: data_Strain.csv has been successfully created
#: data_el_vol.csv has been successfully created
#: data_Stress.csv has been successfully created
#: data_Strain.csv has been successfully created
#: data_el_vol.csv has been successfully created
#: data_Stress.csv has been successfully created
#: data_Strain.csv has been successfully created
#: data_el_vol.csv has been successfully created
#: data_Stress.csv has been successfully created
#: data_Strain.csv has been successfully created
#: data_el_vol.csv has been successfully created
#: data_Stress.csv has been successfully created
#: data_Strain.csv has been successfully created
#: data_el_vol.csv has been successfully created
#: data_Stress.csv has been successfully created
#: data_Strain.csv has been successfully created
#: data_el_vol.csv has been successfully created
#: data_Stress.csv has been successfully created
#: data_Strain.csv has been successfully created
#: data_el_vol.csv has been successfully created
#: data_Stress.csv has been successfully created
#: data_Strain.csv has been successfully created
#: data_el_vol.csv has been successfully created
#: data_Stress.csv has been successfully created
#: data_Strain.csv has been successfully created
#: data_el_vol.csv has been successfully created
#: data_Stress.csv has been successfully created
#: data_Strain.csv has been successfully created
#: data_el_vol.csv has been successfully created
#: data_Stress.csv has been successfully created
#: data_Strain.csv has been successfully created
#: data_el_vol.csv has been successfully created
#: data_Stress.csv has been successfully created
#: data_Strain.csv has been successfully created
#: data_el_vol.csv has been successfully created
#: data_Stress.csv has been successfully created
#: data_Strain.csv has been successfully created
#: data_el_vol.csv has been successfully created
#: data_Stress.csv has been successfully created
#: data_Strain.csv has been successfully created
#: data_el_vol.csv has been successfully created
#: data_Stress.csv has been successfully created
#: data_Strain.csv has been successfully created
#: data_el_vol.csv has been successfully created
#: data_Stress.csv has been successfully created
#: data_Strain.csv has been successfully created
#: data_el_vol.csv has been successfully created
#: data_Stress.csv has been successfully created
#: data_Strain.csv has been successfully created
#: data_el_vol.csv has been successfully created
#: data_Stress.csv has been successfully created
#: data_Strain.csv has been successfully created
#: data_el_vol.csv has been successfully created
#: data_Stress.csv has been successfully created
#: data_Strain.csv has been successfully created
#: data_el_vol.csv has been successfully created
#: data_Stress.csv has been successfully created
#: data_Strain.csv has been successfully created
#: data_el_vol.csv has been successfully created
#: data_Stress.csv has been successfully created
#: data_Strain.csv has been successfully created
#: data_el_vol.csv has been successfully created
#: data_Stress.csv has been successfully created
#: data_Strain.csv has been successfully created
#: data_el_vol.csv has been successfully created
#: data_Stress.csv has been successfully created
#: data_Strain.csv has been successfully created
#: data_el_vol.csv has been successfully created
#: data_Stress.csv has been successfully created
#: data_Strain.csv has been successfully created
#: data_el_vol.csv has been successfully created
#: data_Stress.csv has been successfully created
#: data_Strain.csv has been successfully created
#: data_el_vol.csv has been successfully created
#: data_Stress.csv has been successfully created
#: data_Strain.csv has been successfully created
#: data_el_vol.csv has been successfully created
#: data_Stress.csv has been successfully created
#: data_Strain.csv has been successfully created
#: data_el_vol.csv has been successfully created
#: data_Stress.csv has been successfully created
#: data_Strain.csv has been successfully created
#: data_el_vol.csv has been successfully created
#: data_Stress.csv has been successfully created
#: data_Strain.csv has been successfully created
#: data_el_vol.csv has been successfully created
#: data_Stress.csv has been successfully created
#: data_Strain.csv has been successfully created
#: data_el_vol.csv has been successfully created
#: data_Stress.csv has been successfully created
#: data_Strain.csv has been successfully created
#: data_el_vol.csv has been successfully created
#: data_Stress.csv has been successfully created
#: data_Strain.csv has been successfully created
#: data_el_vol.csv has been successfully created
#: data_Stress.csv has been successfully created
#: data_Strain.csv has been successfully created
#: data_el_vol.csv has been successfully created
#: data_Stress.csv has been successfully created
#: data_Strain.csv has been successfully created
#: data_el_vol.csv has been successfully created
#: data_Stress.csv has been successfully created
#: data_Strain.csv has been successfully created
#: data_el_vol.csv has been successfully created
#: data_Stress.csv has been successfully created
#: data_Strain.csv has been successfully created
#: data_el_vol.csv has been successfully created
#: data_Stress.csv has been successfully created
#: data_Strain.csv has been successfully created
#: data_el_vol.csv has been successfully created
#: data_Stress.csv has been successfully created
#: data_Strain.csv has been successfully created
#: data_el_vol.csv has been successfully created
#: data_Stress.csv has been successfully created
#: data_Strain.csv has been successfully created
#: data_el_vol.csv has been successfully created
#: data_Stress.csv has been successfully created
#: data_Strain.csv has been successfully created
#: data_el_vol.csv has been successfully created
#: data_Stress.csv has been successfully created
#: data_Strain.csv has been successfully created
#: data_el_vol.csv has been successfully created
#: data_Stress.csv has been successfully created
#: data_Strain.csv has been successfully created
#: data_el_vol.csv has been successfully created
#: data_Stress.csv has been successfully created
#: data_Strain.csv has been successfully created
#: data_el_vol.csv has been successfully created
#: data_Stress.csv has been successfully created
#: data_Strain.csv has been successfully created
#: data_el_vol.csv has been successfully created
#: data_Stress.csv has been successfully created
#: data_Strain.csv has been successfully created
#: data_el_vol.csv has been successfully created
#: data_Stress.csv has been successfully created
#: data_Strain.csv has been successfully created
#: data_el_vol.csv has been successfully created
#: data_Stress.csv has been successfully created
#: data_Strain.csv has been successfully created
#: data_el_vol.csv has been successfully created
#: data_Stress.csv has been successfully created
#: data_Strain.csv has been successfully created
#: data_el_vol.csv has been successfully created
#: data_Stress.csv has been successfully created
#: data_Strain.csv has been successfully created
#: data_el_vol.csv has been successfully created
#: data_Stress.csv has been successfully created
#: data_Strain.csv has been successfully created
#: data_el_vol.csv has been successfully created
#: data_Stress.csv has been successfully created
#: data_Strain.csv has been successfully created
#: data_el_vol.csv has been successfully created
#: data_Stress.csv has been successfully created
#: data_Strain.csv has been successfully created
#: data_el_vol.csv has been successfully created
#: /home/rwth1393/GA-Calibration/simulation_9928/00_Data/data_grain_ID_{}.csv
#: data_grain_ID.csv has been successfully created
#: /home/rwth1393/GA-Calibration/simulation_9928/00_Data/data_mat_ID_{}.csv
#: data_grain_ID.csv has been successfully created
print 'RT script done'
#: RT script done
