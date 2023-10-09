import os
import configparser

if __name__ == "__main__":
  config = configparser.ConfigParser()
  config.read('./configs/configs.ini')

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

  cmd += " runBatch.sh"
  print(cmd)
  os.system(cmd)

