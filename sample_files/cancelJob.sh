#!/usr/local_rwth/bin/zsh 
echo mad:1.1332381952115689 
cd /rwthfs/rz/cluster/home/rwth0925/DFG/Calibration/simulation 
module load ABAQUS/2022 
abaqus terminate job=DFG_Calibration 
cd /rwthfs/rz/cluster/home/rwth0925/DFG/Calibration