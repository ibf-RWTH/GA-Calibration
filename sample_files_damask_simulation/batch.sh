#!/bin/bash
#SBATCH --job-name=DAMASK_Test
#SBATCH --output=DAMASK_Test.txt
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --account=rwth1393
#SBATCH --mem-per-cpu=2G
#SBATCH --cpus-per-task=33

cd /rwthfs/rz/cluster/home/p0021070/damask/test

chmod +x damask.sh

apptainer exec /rwthfs/rz/SW/UTIL.common/singularity/damask-grid-alpha7 ./damask.sh

