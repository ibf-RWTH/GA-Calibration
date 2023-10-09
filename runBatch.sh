#!/usr/local_rwth/bin/zsh

# Read user variables
source configs/jobConfig.sh
. $conda_root/etc/profile.d/conda.sh
export PATH="$conda_root/bin:$PATH"
# Now you can activate your configured conda environments
conda activate $conda_env
### Execute your application
### Please remember, to adjust the memory, it must be less than requested above

### For debugging everything above you can simply echo the command instead of running the python script
#echo $PYTHONCODE -t $TEST_FLAG -r $RESTART_FLAG -s $SIM_TYP -n $SLURM_NTASKS -j $SIM_JOB_BASE_NAME -d $PWD --ex_data=$EX_DATA
python $pythoncode -t $test_flag -r $restart_flag -s $sim_type -n $slurm_ntasks -j $sim_job_base_name -d $PWD --ex_data=$ex_data --phases=$phases
