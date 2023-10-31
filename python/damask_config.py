import damask
import configparser
import ast
import numpy as np

class Parser:
    def __init__(self, JobSettings) -> None:
        self.JobSettings = JobSettings

    def parse_matparams(self):
        if self.JobSettings == 'AbaqusJobSettings':
          return self.parse_matparams_abaqus()
        elif self.JobSettings == 'DamaskJobSettings':
          return self.parse_matparams_damask()

    def parse_matparams_abaqus(self):
        config = Utils.CONFIG
        varbound = []
        constantParams = {}
        matparams = {}
        params_names = []
        phases = config.get(self.JobSettings,'PHASES').split(',')
        sim_type = config.get(JobSettings,'sim_type')
        if 'CP' in sim_type:
            global_params = config.options('GlobalParams')

            # add global params to varbound
            constantParams[0] = {}
            for global_param in global_params:
                global_param_value = config.get('GlobalParams',global_param)
                if global_param != "phase_id":
                    # add global params
                    if global_param_value.startswith('['):
                        global_param_value = ast.literal_eval(global_param_value)
                        varbound.append(global_param_value)
                        params_names.append(global_param)
                    else:
                        constantParams[0][global_param] = global_param_value
            matparams[0] = params_names
            params_names = []

            for phase in phases:
                phase_id = config.get(phase, 'phase_id')
                options = config.options(phase)
                constantParams[phase_id] = {}
                #add phase params to varbound
                for option in options:
                    value = config.get(phase, option)
                    if value.startswith('['):
                        value = ast.literal_eval(value)
                        params_names.append(option)
                        varbound.append(value)
                    else:
                        if 'd' in value:
                            value = value.replace('d','e') # subroutine input file sytax

                        if option != "lattice" and option != "plastic":
                            value = ast.literal_eval(value)
                        constantParams[phase_id][option] = value

                matparams[phase_id] = params_names
                params_names = []

            varbound = np.array(varbound)

        elif 'Chaboche' in self.sim_type:
            phase_id = config.get(phases[0], 'phase_id')
            constantParams[phase_id] = {}
            options = config.options(phases[0])
            #add phase params to varbound
            for option in options:
                value = config.get(phases[0], option)
                if value.startswith('['):
                    value = ast.literal_eval(value)
                    params_names.append(option)
                    varbound.append(value)
                else:
                    if 'd' in value:
                        value = value.replace('d','e') # subroutine input file sytax
                    value = ast.literal_eval(value)
                    constantParams[phase_id][option] = value


            matparams[phase_id] = params_names
            varbound = np.array(varbound)

        return constantParams, matparams, varbound

    def parse_matparams_damask(self):
      config = Utils.CONFIG
      varbound = []
      constantParams = {}
      matparams = {}
      params_names = []
      phases = config.get(self.JobSettings,'PHASES').split(',')

      global_params = config.options('GlobalParams')
      # add global params to varbound
      constantParams[0] = {}
      for global_param in global_params:
          global_param_value = config.get('GlobalParams',global_param)
          if global_param != "phase_id":
              # add global params
              if global_param_value.startswith('['):
                  global_param_value = ast.literal_eval(global_param_value)
                  varbound.append(global_param_value)
                  params_names.append(global_param)
              else:
                  constantParams[0][global_param] = global_param_value
      matparams[0] = params_names
      params_names = []

      # add phase settings
      for phase in phases:
        phase_id = config.get(phase, 'phase_id')
        params_names.extend(config.get(self.JobSettings, f'{phase}.variables').split(','))
        matparams[phase_id] = params_names
        constantParams[phase_id] = {}
        options = config.options(phase)
        for option in options:
            value = config.get(phase, option)
            if option in params_names:
                value = ast.literal_eval(value)
                varbound.append(value)
            else:
                if option != "lattice" and option != "type":
                    value = ast.literal_eval(value)
                constantParams[phase_id][option] = value

        params_names = []

      varbound = np.array(varbound)

      return constantParams, matparams, varbound

class Utils:
    CONFIG = None
    CONSTANTPARAMS = None
    EVALUATINGPARAMS = None
if __name__ == "__main__":
  m = damask.ConfigMaterial()
  m = m.load('/home/p0021070/damask/GA-Calibration-Damask/damask/sample_files_simulation/material.yaml')
  phase = m.get('phase')
  alpha = phase['Alpha']
  alpha['mechanical']['plastic']['a_sl'] = float(np.float64(2.0112197953175355))
  m.save('test.yaml')
  # config = configparser.ConfigParser(allow_no_value=True)
  # config.optionxform = lambda option : option # Preserve case of the keys
  # sim_root = '/home/p0021070/damask/GA-Calibration-Damask'
  # config.read(f'{sim_root}/configs/configs.ini')
  # Utils.CONFIG = config

  # software = config.get('MainProcessSettings','software')
  # JobSettings = 'AbaqusJobSettings' if software == 'abaqus' else 'DamaskJobSettings'
  # parser = Parser(JobSettings)
  # constantParams, matparams, varbound = parser.parse_matparams()

  # print(constantParams)
  # print(matparams)
  # print(varbound)

  # phases = config.get('DamaskJobSettings','PHASES').split(',')
  # for phase in phases:
  #   phase_data = m.get('phase')[phase]
  #   phase_id = config.get(phase, 'phase_id')
  #   # config global constant params
  #   for key, item in constantParams[0].items():
  #     if key in phase_data['mechanical']['plastic'].keys():
  #       phase_data['mechanical']['plastic'][key] = item
  #     elif key in phase_data['mechanical']['elastic'].keys():
  #       phase_data['mechanical']['elastic'][key] = item
  #   #config phase constant params
  #   for key, item in constantParams[phase_id].items():
  #     if key in phase_data['mechanical']['plastic'].keys():
  #       phase_data['mechanical']['plastic'][key] = item
  #     elif key in phase_data['mechanical']['elastic'].keys():
  #       phase_data['mechanical']['elastic'][key] = item
  #     elif key == 'lattice':
  #       phase_data[key] = item



