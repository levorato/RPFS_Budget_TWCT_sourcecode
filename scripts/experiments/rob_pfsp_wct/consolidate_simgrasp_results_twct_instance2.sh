#!/bin/bash

#SBATCH --job-name=TWCT.SG.I2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=cpuonly
#SBATCH --mem=32G
#SBATCH --time=20000

echo "Starting job on node $SLURMD_NODENAME"
cd /users/ubuntu/RPFS_Budget_TWCT/src/julia
sh ./install_gurobi_license.sh $SLURMD_NODENAME

echo "RobPFSP Julia TWCT consolidate simgrasp results experiment script start."

cd /users/ubuntu/grasp/RPFS_Budget_TWCT/scripts/experiments/rob_pfsp_wct
/users/ubuntu/julia-1.6.0/bin/julia consolidate_simgrasp_results_twct.jl --machine=$SLURMD_NODENAME --instances=instance_2
