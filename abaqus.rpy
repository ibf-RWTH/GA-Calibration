# -*- coding: mbcs -*-
#
# Abaqus/CAE Release 2022.HF6 replay file
# Internal Version: 2022_11_16-21.42.28 RELr424 176889
# Run by bq407628 on Tue Sep 26 12:52:34 2023
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
execfile('/home/rwth1393/GA_Calibration_Test/python/readOdb.py', 
    __main__.__dict__)
#: /rwthfs/rz/cluster/home/rwth1393/single
#: CP_Calibration_7743
#* OdbError: Cannot open file 
#* /rwthfs/rz/cluster/home/rwth1393/single/CP_Calibration_7743.odb. *** ERROR: 
#* No such file: 
#* /rwthfs/rz/cluster/home/rwth1393/single/CP_Calibration_7743.odb.
#* File "/home/rwth1393/GA_Calibration_Test/python/readOdb.py", line 174, in 
#* <module>
#*     readOdb(odbpath, odbname, key_Stress, key_Strain, key_EVOL, key_Part, 
#* storename)
#* File "/home/rwth1393/GA_Calibration_Test/python/readOdb.py", line 15, in 
#* readOdb
#*     myOdb = odbAccess.openOdb(path=odbfilepath)
