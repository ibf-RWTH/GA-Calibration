#!/usr/local_rwth/bin/zsh

### Job name
#SBATCH --job-name=CP_Cal_Single
 
### File/Path where STDOUT will be written to, %J is the job id
#SBATCH --output /home/rwth1393/single/logs/MainProcess-log.%J
 
### Request the time you need for execution. The full format is D-HH:MM:SS
### You must at least specify minutes or days and hours and may add or
### leave out any other parameters
#SBATCH --time=7-00:00:00
 
### Request the memory you need for your job. You can specify this
### in either MB (1024M) or GB (4G). BEWARE: This is a per-cpu limit,
### and will be multiplied with cpus-per-task for the total requested memory
#SBATCH --mem-per-cpu=1024M
### Request one host
#SBATCH --nodes=1
 
### Request number of CPUs/MPI Ranks
#SBATCH --ntasks=8

# Insert this AFTER the #SLURM argument section of your job script
## chose anaconda or miniconda as CONDA_ROOT depending on what is installed for you
## export CONDA_ROOT=$HOME/anaconda3
export CONDA_ROOT=/home/rwth1393/anaconda3
. $CONDA_ROOT/etc/profile.d/conda.sh
export PATH="$CONDA_ROOT/bin:$PATH"

# Now you can activate your configured conda environments
conda activate calibration
export SIM_JOB_BASE_NAME=CP_Single
### Execute your application
### Please remember, to adjust the memory, it must be less than requested above
export PYTHONCODE=$PWD/python/calibration_genetic.py
export TEST_FLAG=False
export RESTART_FLAG=False
export SIM_TYP=tensile
export EX_DATA=ex_data_single.csv

python $PYTHONCODE -t $TEST_FLAG -r $RESTART_FLAG -s $SIM_TYP -n $SLURM_NTASKS -j $SIM_JOB_BASE_NAME -d $PWD --ex_data=$EX_DATA
