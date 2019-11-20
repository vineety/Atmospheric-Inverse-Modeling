#!/bin/sh
#PBS -o /home/yshiga/Documents/ASCENDS/Inversions/Output/July/3hrly/wgt2micron/Total_ff_par_test/
#PBS -j oe 
#PBS -l nodes=8:ppn=1,walltime=10:00:00,mem=350gb
#PBS -M yshiga@stanford.edu 
#PBS -m a 
#PBS -N Total_JULY_ff_par_test 
#PBS -r n 
#PBS -V 
cd $PBS_O_WORKDIR
cd /home/yshiga/Documents/Inversion/Parallelv2/
source /opt/torque/etc/openmpi-setup.sh
mpiexec --hostfile $PBS_NODEFILE -loadbalance inverse_parallel.exe config_July_Total_xFF.txt
