from odbAccess import *
import numpy as np
import csv
import sys

if __name__ == "__main__": 
    odbpath = sys.argv[-2]
    odbname = sys.argv[-1]
    odbfilepath = odbpath+'/'+odbname+".odb"
    odb = openOdb(path=odbfilepath)

    step = odb.steps['Step-1']  

    csv_filename = 'results_df.csv'
    with open(odbpath+'/'+csv_filename, mode='w') as csvfile:  
        csvwriter = csv.writer(csvfile)
        
        csvwriter.writerow(['Frame', 'Time', 'PEEQ', 'Mises Stress'])
        
        for i, frame in enumerate(step.frames):

            sdv_field = frame.fieldOutputs['SDV54']  
            element_sdv54 = []
            
            for value in sdv_field.values:
                sdv54 = value.data  
                element_sdv54.append(sdv54)
            
            average_sdv54 = np.mean(element_sdv54)
            
            stress_field = frame.fieldOutputs['S']
            mises_stress_values = []
            
            for stress_value in stress_field.values:
                mises_stress = stress_value.mises  
                mises_stress_values.append(mises_stress)
            
            average_mises = np.mean(mises_stress_values)
            
            time = frame.frameValue
            
            csvwriter.writerow([i+1, time, average_sdv54, average_mises])

    odb.close()

