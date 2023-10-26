import os
import configparser
import glob

config = configparser.ConfigParser()
config.read('./configs/configs.ini')
#call batch settings
batchSettings = config['MainProcessSettings']
jobName = batchSettings['job-name']
output = batchSettings['output']+'/MainProcess-log.%J'
time = batchSettings['time']
memPerCpu = batchSettings['mem-per-cpu']
nodes = batchSettings['nodes']
ntasks = batchSettings['ntasks']

#call conda settings
condaSettings = config['CondaSettings']
condaRoot = condaSettings['CONDA_ROOT']
condaEnv = condaSettings['CONDA_ENV']
pythoncode = config.get('JobSettings','PYTHONCODE')

sim_job_base_name = config.get('JobSettings','SIM_JOB_BASE_NAME')

f = open('./configs/JobConfig.sh', 'w+')
f.write(f'CONDA_ROOT={condaRoot}\n')
f.write(f'CONDA_ENV={condaEnv}\n')
f.write(f'PYTHONCODE={pythoncode}\n')
f.write(f'ROOT={os.getcwd()}')
f.close()

image_path = f'{os.getcwd()}/evaluation_images_{sim_job_base_name}'
logs_path = f'{os.getcwd()}/logs_{sim_job_base_name}'

if os.path.isdir(image_path) and len(os.listdir(image_path)) > 0:
        os.system(f'rm -r {image_path}/*.png')

if not os.path.isdir(logs_path):
        os.mkdir(f'{logs_path}')
elif os.path.isdir(logs_path) and len(os.listdir(logs_path)) > 0:
        os.system(f'rm -r {logs_path}/*log*')
        files = glob.glob(f"{logs_path}/*.txt")
        if files:
                os.system(f'rm -r {logs_path}/*txt')
# run batch file
os.system(f"sbatch --job-name={jobName} --output={output} --time={time} --mem-per-cpu={memPerCpu} --nodes={nodes} --ntasks={ntasks} runBatch.sh")
