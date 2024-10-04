import os
import configparser
import glob
import sys
import shutil

config = configparser.ConfigParser(allow_no_value=True)
config.optionxform = lambda option : option # Preserve case of the keys
config.read('./configs/configs.ini')
#call batch settings
batchSettings = config['MainProcessSettings']
jobName = batchSettings['job-name']
time = batchSettings['time']
memPerCpu = batchSettings['mem-per-cpu']
nodes = batchSettings['nodes']
ntasks = batchSettings['ntasks']
software = batchSettings['software']
partition = batchSettings['partition']

#call conda settings
condaSettings = config['CondaSettings']
condaRoot = condaSettings['CONDA_ROOT']
condaEnv = condaSettings['CONDA_ENV']
if software == "abaqus":
        pythoncode = config.get('AbaqusJobSettings','PYTHONCODE')
        sim_job_base_name = config.get('AbaqusJobSettings','sim_job_base_name')
elif software == "damask":
        pythoncode = config.get('DamaskJobSettings','PYTHONCODE')
        sim_job_base_name = config.get('DamaskJobSettings','sim_job_base_name')
else:
        os.system("echo error: wrong software name!")
        sys.exit()



f = open('./configs/JobConfig.sh', 'w+')
f.write(f'CONDA_ROOT={condaRoot}\n')
f.write(f'CONDA_ENV={condaEnv}\n')
f.write(f'PYTHONCODE={pythoncode}\n')
f.write(f'ROOT={os.getcwd()}')
f.close()

image_path = f'{os.getcwd()}/evaluation_images_{sim_job_base_name}'
logs_path = f'{os.getcwd()}/logs_{sim_job_base_name}'
output = logs_path +'/MainProcess-log.%J'

# prepare data_tree for upcoming simulations
if not os.path.isdir(image_path):
        os.mkdir(f'{image_path}')
elif os.path.isdir(image_path) and len(os.listdir(image_path)) > 0:
        os.system(f'rm -r {image_path}/*.png')

if not os.path.isdir(logs_path):
        os.mkdir(f'{logs_path}')
elif os.path.isdir(logs_path) and len(os.listdir(logs_path)) > 0:
        os.system(f'rm -r {logs_path}/*log*')
        files = glob.glob(f"{logs_path}/*.txt")
        if files:
                os.system(f'rm -r {logs_path}/*txt')

# remove previous simulation files
# Directory to search in, assuming current directory for this example
search_dir = "."

# Find and delete folders starting with 'simulation_' followed by a date
for root, dirs, files in os.walk(search_dir):
    for dir_name in dirs:
        if dir_name.startswith("simulation_") and len(dir_name) == 15:  # e.g., simulation_0203
            shutil.rmtree(os.path.join(root, dir_name))

# run batch file
os.system(f"sbatch --account=rwth1694 --job-name={jobName} --partition={partition} --output={output} --time={time} --mem-per-cpu={memPerCpu} --nodes={nodes} --ntasks={ntasks} runBatch.sh")
