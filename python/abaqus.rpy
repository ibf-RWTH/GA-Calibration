# -*- coding: mbcs -*-
#
# Abaqus/CAE Release 2022.HF6 replay file
# Internal Version: 2022_11_16-21.42.28 RELr424 176889
# Run by bq407628 on Mon Sep 30 10:32:55 2024
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
execfile('readOdb_CPPF.py', __main__.__dict__)
#: Model: /hpcwork/rwth1393/pfcp_simulation/nrb_r3v1/Submodel_r3v1.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       275
#: Number of Node Sets:          3
#: Number of Steps:              1
print 'RT script done'
#: RT script done
