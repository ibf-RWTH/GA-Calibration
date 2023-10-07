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

#call job settings
jobSettings = config['JobSettings']
simJobBaseName = jobSettings['SIM_JOB_BASE_NAME']
pythonCode = jobSettings['PYTHONCODE']
testFlag = jobSettings['TEST_FLAG']
restartFlag = jobSettings['RESTART_FLAG']
simType = jobSettings['SIM_TYPE']
exDataPath = jobSettings['EX_DATA']
phases = jobSettings['PHASES'] 

f = open('./configs/jobConfig.sh', 'w+')
f.write(f'CONDA_ROOT={condaRoot}\n')
f.write(f'CONDA_ENV={condaEnv}\n')
f.write(f'SIM_JOB_BASE_NAME={simJobBaseName}\n')
f.write(f'PYTHONCODE={pythonCode}\n')
f.write(f'TEST_FLAG={testFlag}\n')
f.write(f'RESTART_FLAG={restartFlag}\n')
f.write(f'SIM_TYPE={simType}\n')
f.write(f'EX_DATA={exDataPath}\n')
f.write(f'PHASES={phases}')
f.close()

# run batch file
os.system(f"sbatch --job-name={jobName} --output={output} --time={time} --mem-per-cpu={memPerCpu} --nodes={nodes} --ntasks={ntasks} runBatch.sh")
