#!/usr/local_rwth/bin/zsh

# Read user variables
source configs/SimJobConfig.sh
module load $ABAQUS
### Create ABAQUS environment file for current job, you can set/add your own options (Python syntax)
env_file=abaqus_v6.env
cat << EOF > ${env_file}
#verbose = 3
#ask_delete = OFF
#usub_lib_dir=os.getcwd()
mp_file_system = (SHARED, LOCAL)
mp_host_list = $R_WLM_ABAQUSHOSTLIST
mp_mode=MPI
EOF
unset SLURM_GTIDS
unsetopt -o NOMATCH
rm -f *$JOBNAME*.* 2>/dev/null
rm -f abaqus.r* 2>/dev/null
sleep 5
setopt -o NOMATCH
### Execute your application
### Please remember, to adjust the memory, it must be less than requested above
abaqus interactive job=$JOBNAME input=$INPUTFILE cpus=$SLURM_NTASKS threads_per_mpi_process=$THREADS_PER_MPI memory="$ABAQUS_MEM_ARG" user=$SUBROUTINE_PATH
abaqus cae noGUI=$PYTHON_PATH -- $PWD $JOBNAME