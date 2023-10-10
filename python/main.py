import os
import configparser
import sys
import ast
import numpy as np
from calibration_genetic import Simulation, Optimize
import subprocess
import json
from geneticalgorithm2 import geneticalgorithm2 as ga

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

  batch_cmd = {}
  batch_options = config.options('BatchSettings')
  cmd = 'sbatch'
  for opt in batch_options:
    cmd += f" --{opt}={config.get('BatchSettings', opt)}"
    batch_cmd[opt] = config.get('BatchSettings', opt)


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
      batch_cmd[opt] = config.get('BatchSettings', opt)

  return batch_cmd

if __name__ == "__main__":
  #get current path
  sim_root = os.getcwd()

  #config simulation
  config = configparser.ConfigParser()
  config.read(f'{sim_root}/configs/configs.ini')
  mat_params, varbound = parse_matparams(config)

  #initialize Optimizer
  algorithm_param = {'max_num_iteration': 100, \
                       'population_size': 50, \
                       'mutation_probability': 0.1, \
                       'elit_ratio': 0.1, \
                       'parents_portion': 0.3, \
                       'max_iteration_without_improv': None}

  opt = Optimize(flag=config.get('JobSettings','test_flag'),
                 ex_data=config.get('JobSettings','ex_data'),
                 root=sim_root,
                 name=config.get('JobSettings','sim_job_base_name'),
                 mat_params=mat_params,
                 varbound=varbound,
                 algorithm_param=algorithm_param,
                 sim_flag=config.get('JobSettings','sim_type'),
                 n_jobs = config.get('BatchSettings','ntasks'))

  model, func = opt.init_optimizer()
  os.system('echo optimizer initialized starting simulations now')
  if not ast.literal_eval(config.get('JobSettings','restart_flag')):
      model.run(no_plot=True,
                progress_bar_stream = None,
                save_last_generation_as = f'{sim_root}/logs/lastgeneration.npz',
                set_function=ga.set_function_multiprocess(func, n_jobs=ast.literal_eval(config.get('BatchSettings','ntasks'))))
  else:
      model.run(no_plot=True,
              progress_bar_stream = None,
              start_generation=f'{sim_root}/logs/lastgeneration.npz',
              set_function=ga.set_function_multiprocess(func, n_jobs=ast.literal_eval(config.get('BatchSettings','ntasks'))))

  f = open(sim_root + '/logs/final_results.txt', 'w')
  json.dump(model.output_dict, f, indent=4)






