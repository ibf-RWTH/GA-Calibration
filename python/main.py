import os
import configparser
import sys
import ast 
import numpy as np
from calibration_genetic import Simulation, Optimize

def parse_matparams(config:configparser.ConfigParser):
  global_params = config.options('GlobalParams')
  phases = config.get('JobSettings','PHASES').split(',')

  varbound = []
  matparams = []
  for phase in phases:
    options = config.options(phase)
    # add global params to varbound
    for global_param in global_params:
      global_param_value = ast.literal_eval(config.get('GlobalParams',global_param))
      
      if global_param != "abaqus_id":
        # check params type
        if not isinstance(global_param_value, list):
          os.system("echo config error: global params to be calibrated must be list")
          sys.exit()
        
        matparams.append(global_param)
        # add global params
        varbound.append(global_param_value)

    #add phase params to varbound
    for option in options:
      value = config.get(phase, option)
      if value.startswith('['):
        value = ast.literal_eval(value)
        matparams.append(option)
        varbound.append(value)

    varbound = np.array(varbound)

    return matparams, varbound

def gen_cmd(config):

  batch_options = config.options('BatchSettings')
  cmd = 'sbatch'
  for opt in batch_options:
    cmd += f" --{opt}={config.get('BatchSettings', opt)}"

  conda_options = config.options('CondaSettings')
  with open("./configs/jobConfig.sh",'w+') as jc:
    jc.write("export")
    for opt in conda_options:
      jc.write(f' {opt}={config.get("CondaSettings", opt)}')

    jc.write("\n")

  job_options = config.options('JobSettings')
  with open("./configs/jobConfig.sh",'a+') as jc:
    jc.write("export")
    for opt in job_options:
      jc.write(f' {opt}={config.get("JobSettings", opt)}')

  return cmd

if __name__ == "__main__":
  config = configparser.ConfigParser()
  config.read('./configs/configs.ini')

  cmd = gen_cmd(config)
  print(cmd)

