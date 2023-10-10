import os
import configparser

config = configparser.ConfigParser()
config.read('./configs/configs.ini')
#call batch settings
batchSettings = config['BatchSettings']
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

f = open('./configs/jobConfig.sh', 'w+')
f.write(f'CONDA_ROOT={condaRoot}\n')
f.write(f'CONDA_ENV={condaEnv}\n')
f.close()

# run batch file
os.system(f"sbatch --job-name={jobName} --output={output} --time={time} --mem-per-cpu={memPerCpu} --nodes={nodes} --ntasks={ntasks} runBatch.sh")
