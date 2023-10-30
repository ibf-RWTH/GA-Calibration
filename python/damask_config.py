import damask
m = damask.ConfigMaterial()
m = m.load('/home/p0021070/damask/GA-Calibration-Damask/sample_files_damask_simulation/material.yaml')
# print(m)
alpha = m.get('phase')['Alpha']
print(alpha['mechanical']['plastic']['a_sl'])

