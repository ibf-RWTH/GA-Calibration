#!/usr/local_rwth/bin/zsh
#SBATCH --partition=c23ml
#SBATCH --account=rwth1393
### Job name
#SBATCH --job-name=Submodel_r3

### File/Path where STDOUT will be written to, %J is the job id
#SBATCH --output submodel_r3-log.%J

### Request the time you need for execution. The full format is D-HH:MM:SS
### You must at least specify minutes or days and hours and may add or
### leave out any other parameters
#SBATCH --time=48:00:00

### Request the memory you need for your job. You can specify this
### in either MB (1024M) or GB (4G). BEWARE: This is a per-cpu limit,
### and will be multiplied with cpus-per-task for the total requested memory
#SBATCH --mem-per-cpu=24G

### Request number of hosts
#SBATCH --nodes=1

### Request number of CPUs/MPI Ranks
#SBATCH --ntasks=8
# Read user variables
module load ABAQUS/2022
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
export THREADS_PER_MPI=1
unset SLURM_GTIDS
unsetopt -o NOMATCH
rm -f abaqus.r* 2>/dev/null
rm -f Submodel*
sleep 5
setopt -o NOMATCH
### Execute your application
### Please remember, to adjust the memory, it must be less than requested above
abaqus interactive job=Submodel_r3v1 input=DRAGen_RVE.inp cpus=8 threads_per_mpi_process=$THREADS_PER_MPI global=nrb_r3 memory=16Gb user=/hpcwork/rwth1393/pfcp_cleavage/Umat_PFCP.for
