from abaqus import *
from abaqusConstants import *
import os
import numpy as np
import sys


def readOdb(odbpath, odbname):
    odb = session.openOdb(odbpath+'/'+odbname+'.odb')
    
    dataRT2 = [frame.fieldOutputs['RT'].getSubset(region=odb.rootAssembly.instances['WERKZEUG_OBEN'].nodeSets['RP']).values[0].data[1]
               for frame in odb.steps['Bewegung'].frames]
    dataUT2 = [frame.fieldOutputs['UT'].getSubset(region=odb.rootAssembly.instances['WERKZEUG_OBEN'].nodeSets['RP']).values[0].data[1]
              for frame in odb.steps['Bewegung'].frames]
    dataTime = [frame.frameValue for frame in odb.steps['Bewegung'].frames]

    simulated_data = np.array(tuple(zip(dataTime,dataUT2,dataRT2)))
    simulated_data[:,2] = simulated_data[:,2]*20*2/1000
    
    f = open("RF_data.txt","w")
    for value in simulated_data:
        string = str(value[0])+","+str(value[1])+ "," + str(value[2]) + "\n"
        f.write(string)
    f.close() 
    odb.close()    

if __name__ == "__main__":
    odbpath = sys.argv[-2]
    odbname = sys.argv[-1]
    storename = odbname

    print(odbpath)
    print(odbname)
    readOdb(odbpath, odbname)
