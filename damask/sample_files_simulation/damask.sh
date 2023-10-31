#!/bin/bash

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

DAMASK_grid --load load.yaml --geom grid.vti