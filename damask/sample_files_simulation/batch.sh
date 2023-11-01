#!/bin/bash
cd /rwthfs/rz/cluster/home/p0021070/damask/example

chmod +x damask.sh

apptainer exec /rwthfs/rz/SW/UTIL.common/singularity/damask-grid-alpha7 ./damask.sh

