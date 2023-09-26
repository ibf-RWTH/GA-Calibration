#!/usr/local_rwth/bin/zsh
### Job name
#SBATCH --job-name=%JOBNAME%
### Account
#SBATCH --account=rwth1393
### File/Path where STDOUT will be written to, %J is the job id
#SBATCH --output /home/rwth1393/single/logs/%JOBNAME%-log.%J
### Request the time you need for execution. The full format is D-HH:MM:SS
### You must at least specify minutes or days and hours and may add or
### leave out any other parameters
#SBATCH --time=5:00:00
### Request the memory you need for your job. You can specify this
### in either MB (1024M) or GB (4G). BEWARE: This is a per-cpu limit,
### and will be multiplied with cpus-per-task for the total requested memory
#SBATCH --mem-per-cpu=1G
### Request one host
#SBATCH --nodes=1
### Request number of CPUs/MPI Ranks
#SBATCH --ntasks=16
### Initialization of the software
module load ABAQUS/2022
### Set the amount of memory to be passed to Abaqus as a command line argument
### Beware: This HAS to be lower than the value you requested via --mem-per-cpu*--ntasks
export ABAQUS_MEM_ARG="16 Gb"
### number of threads on each MPI RANK needed since ABAQUS 2019 doesnt work with ealier versions
export THREADS_PER_MPI=1
# directory variables
export ROOT=%ROOT%
export PYTHON_PATH="$ROOT/python/readOdb.py"
export SUBROUTINE_PATH="$ROOT/subroutine/Umat_CP.for"
### name your job HERE, name it DIFFERENT from your input file!
export JOBNAME=%JOBNAME%
export INPUTFILE=DRAGen_RVE.inp
### Change (!) to your desired work directory
cd $ROOT/simulation
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