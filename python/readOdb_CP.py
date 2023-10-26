from odbAccess import openOdb
import odbAccess
import math
import numpy as np
import os
import sys

def readOdb(odbpath, odbname, key_Stress, key_Strain, key_EVOL, key_Part, storename):
    stress_type = 'S22'
    strain_type = 'LE22'
    odbfilepath = odbpath+'/'+odbname+".odb"
    dataDir = odbpath+"/00_Data/"
    if not os.path.isdir(dataDir):
        os.system("mkdir -p "+dataDir)
    myOdb = odbAccess.openOdb(path=odbfilepath)
    myinstance = myOdb.rootAssembly.instances.values()[0]
    mystep = myOdb.steps.values()[-1]
    start_frame = 0
    last_frame = len(mystep.frames)

    time = list()
    for frame in mystep.frames:
        time.append(frame.frameValue)
    with open(dataDir + 'data_time_{}.csv'.format(storename), 'wb+') as f:
        f.write('frame_id,sim_time\n')
        for i, timestep in enumerate(time):
            f.write(str(i)+','+str(timestep).strip())
            f.write('\n')

    for i in range(start_frame, last_frame):
        if i < 10:
            counter = '00'+str(i)
        elif i > 9 and i < 100:
            counter = '0'+str(i)
        else:
            counter = str(i)

        try:
            el_IDs_Stress = list()
            Stress = list()
            for value in mystep.frames[i].fieldOutputs[key_Stress].values:
                el_IDs_Stress.append(value.elementLabel)
                if stress_type == 'maxps':
                    Stress.append(value.maxPrincipal)
                elif stress_type == 'S11':
                    Stress.append(value.data[0])
                elif stress_type == 'S22':
                    Stress.append(value.data[1])
                elif stress_type == 'S33':
                    Stress.append(value.data[2])

            data_Stress = [el_IDs_Stress,Stress]
            data_Stress = list(zip(*data_Stress))

            with open(dataDir +'data_Stress_{}_frame_{}.csv'.format(storename,counter), 'wb+') as f:
                f.write('element_id,Stress\n')
                for value in data_Stress:
                    f.write(str(value)[1:-1])
                    f.write('\n')
            print('data_Stress.csv has been successfully created')
        except Exception as e:
            print('Whoops something went wrong: '+ str(e))

        try:
            el_IDs_Strain = list()
            Strain = list()
            for value in mystep.frames[i].fieldOutputs[key_Strain].values:
                el_IDs_Strain.append(value.elementLabel)
                if strain_type == 'maxps':
                    Strain.append(value.maxPrincipal)
                if strain_type == 'LE11':
                    Strain.append(value.data[0])
                if strain_type == 'LE22':
                    Strain.append(value.data[1])
                if strain_type == 'LE33':
                    Strain.append(value.data[2])

            data_Strain= [el_IDs_Strain,Strain]
            data_Strain= list(zip(*data_Strain))

            with open(dataDir +'data_Strain_{}_frame_{}.csv'.format(storename,counter), 'wb+') as f:
                f.write('element_id,Strain\n')
                for value in data_Strain:
                    f.write(str(value)[1:-1])
                    f.write('\n')
            print('data_Strain.csv has been successfully created')
        except Exception as e:
            print('Whoops something went wrong: '+ str(e))

        try:
            el_IDs_vol = list()
            el_vol = list()
            for value in mystep.frames[i].fieldOutputs[key_EVOL].values:
                el_IDs_vol.append(value.elementLabel)
                el_vol.append(value.data)

            data_el_vol = [el_IDs_vol,el_vol]
            data_el_vol = list(zip(*data_el_vol))

            with open(dataDir + 'data_el_vol_{}_frame_{}.csv'.format(storename,counter), 'wb+') as f:
                f.write('element_id,el_vol\n')
                for line in data_el_vol:
                    f.write(str(line)[1:-1])
                    f.write('\n')
            print('data_el_vol.csv has been successfully created')
        except Exception as e:
            print('Whoops something went wrong: '+ str(e))

    try:
        grainIDs = list()
        elementIDs = list()
        grainID = 0
        elem_set_keys = myinstance.elementSets.keys()
        for set_key in elem_set_keys:
            elem_set = myinstance.elementSets[set_key]
            for element in elem_set.elements:
                elementIDs.append(element.label)
                grainIDs.append(grainID)
            grainID +=1

        data_grainID = [elementIDs, grainIDs]
        data_grainID = list(zip(*data_grainID))
        print(dataDir + 'data_grain_ID_{}.csv')
        with open(dataDir + 'data_grain_ID_{}.csv'.format(storename), 'wb+') as f:
            f.write('element_id,grainID\n')
            for line in data_grainID:
                f.write(str(line)[1:-1])
                f.write('\n')
        print('data_grain_ID.csv has been successfully created')
    except Exception as e:
        print('Whoops something went wrong: '+ str(e))
    try:
        grainIDs = list()
        matIDs = list()
        grainID = 0
        section_keys = myOdb.sections.keys()
        for section_key in section_keys:
            material = myOdb.sections[section_key].material
            if 'FERRITE' in material:
                matID = 1
            elif 'MARTENSITE' in material:
                matID = 2
            elif 'PEARLITE' in material:
                matID = 3
            elif 'BAINITE' in material:
                matID = 4
            grainIDs.append(grainID)
            matIDs.append(matID)
            grainID += 1
        data_matID = [grainIDs, matIDs]
        data_matID = list(zip(*data_matID))
        print(dataDir + 'data_mat_ID_{}.csv')
        with open(dataDir + 'data_mat_ID_{}.csv'.format(storename), 'wb+') as f:
            f.write('grainID,matID\n')
            for line in data_matID:
                f.write(str(line)[1:-1])
                f.write('\n')
        print('data_grain_ID.csv has been successfully created')
    except Exception as e:
        print('Whoops something went wrong: '+ str(e))



if __name__ == "__main__":
    odbpath = sys.argv[-2]
    odbname = sys.argv[-1]
    storename = odbname
    key_Stress = 'S'
    key_Strain = 'LE'
    key_EVOL = 'EVOL'
    key_Part = 'PART-1-1'
    print(odbpath)
    print(odbname)
    readOdb(odbpath, odbname, key_Stress, key_Strain, key_EVOL, key_Part, storename)
