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

f = open('./configs/SimJobConfig.sh', 'w+')
f.write(f'CONDA_ROOT={condaRoot}\n')
f.write(f'CONDA_ENV={condaEnv}\n')
f.write(f'PYTHONCODE={pythoncode}\n')
f.write(f'ROOT={os.getcwd()}')
f.close()

image_path = f'{os.getcwd()}/evaluation_images'
logs_path = f'{os.getcwd()}/logs'

if os.path.isdir(image_path) and len(os.listdir(image_path)) > 0:
        os.system(f'rm -r {os.getcwd()}/evaluation_images/*.png')

if not os.path.isdir(logs_path):
        os.mkdir(f'{os.getcwd()}/logs')
elif os.path.isdir(logs_path) and len(os.listdir(logs_path)) > 0:
        os.system(f'rm -r {os.getcwd()}/logs/*log*')
        files = glob.glob(f"{os.getcwd()}/logs/*.txt")
        if files:
                os.system(f'rm -r {os.getcwd()}/logs/*txt')
# run batch file
os.system(f"sbatch --job-name={jobName} --output={output} --time={time} --mem-per-cpu={memPerCpu} --nodes={nodes} --ntasks={ntasks} runBatch.sh")
