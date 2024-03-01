import re
import numpy as np
import os
import configparser
import shutil

if __name__ == "__main__":

  sim_root = os.getcwd()
  config = configparser.ConfigParser(allow_no_value=True)
  config.optionxform = lambda option : option # Preserve case of the keys
  config.read(f'{sim_root}/configs/configs.ini')

  software = config.get('MainProcessSettings','software')
  JobSettings = 'AbaqusJobSettings' if software == 'abaqus' else 'DamaskJobSettings'
  name = config.get(JobSettings,'sim_job_base_name')
  logs_path = f'{sim_root}/logs_{name}'

  min_error = 1e9
  min_index = 0

  with open(f'{logs_path}/results.txt','r') as f:
    lines = f.readlines()
    line_index = 0
    for line in lines:
      if 'Error: ' in line:
        error_index = line.index('Error: ')
        error = line[error_index + len('Error: '):]
        if float(error) < min_error:
          min_error = float(error)
          min_index = line_index

      line_index += 1

    time_stamp_index = lines[min_index - 3].index('time_Stamp: ')
 
    time_stamp = lines[min_index - 3][time_stamp_index + len('time_Stamp: '): -1]
    img1 = f'{sim_root}/evaluation_images_{name}/stress_strain_{time_stamp}.png'
    img2 = f'{sim_root}/evaluation_images_{name}/stress_strain_phase1_{time_stamp}.png'
    img3 = f'{sim_root}/evaluation_images_{name}/stress_strain_phase2_{time_stamp}.png'

    shutil.copy(img1, f'{sim_root}/postproc/stress_strain_{time_stamp}.png')

    if os.path.exists(img2):
      shutil.copy(img2, f'{sim_root}/postproc/stress_strain_phase1_{time_stamp}.png')

    if os.path.exists(img3):
      shutil.copy(img3, f'{sim_root}/postproc/stress_strain_phase2_{time_stamp}.png')
    params = []
    global_param = lines[min_index + 1].split(':')[2].strip().replace(',','')
    params.append(float(global_param))

    phase1_param_line = re.split(r'[:,]', lines[min_index + 2])

    pattern = r'^\d+(\.\d+)?$'
    for param in phase1_param_line:
      param = param.strip()
      if bool(re.match(pattern, param)):
        params.append(float(param))

    if not lines[min_index + 2].startswith('#'):
      phase2_param_line = re.split(r'[:,]', lines[min_index + 3])
      for param in phase2_param_line:
        param = param.strip()
        if bool(re.match(pattern, param)):
          params.append(float(param))

    params = np.array(params).reshape(1, -1)
    np.save(f'{logs_path}/lastgeneration.npy', params)








