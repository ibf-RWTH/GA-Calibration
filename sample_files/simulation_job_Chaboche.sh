#!/usr/local_rwth/bin/zsh
# Read user variables
source SimJobConfig.sh
module load $ABAQUS
### Create ABAQUS environment file
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
#
abaqus interactive job=$JOBNAME input=$INPUTFILE cpus=$SLURM_NTASKS double=both
abaqus cae noGUI=$ROOT/$PYTHON_PATH -- $SIM_DIR $JOBNAME