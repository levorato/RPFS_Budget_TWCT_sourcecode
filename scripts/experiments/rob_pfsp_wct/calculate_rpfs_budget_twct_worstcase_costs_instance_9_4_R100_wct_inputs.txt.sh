#!/bin/bash

#SBATCH --job-name=WC.I2.R100
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=cpuonly
#SBATCH --mem=32G
#SBATCH --time=20000

echo "Starting job on node $SLURMD_NODENAME"
cd /users/ubuntu/RPFS_Budget_TWCT/src/julia
sh ./install_gurobi_license.sh $SLURMD_NODENAME

echo "RobPFSP Julia TWCT worstcase experiment script start."

cd /users/ubuntu/grasp/RPFS_Budget_TWCT/scripts/experiments/rob_pfsp_wct
/users/ubuntu/julia-1.6.0/bin/julia calculate_rpfs_budget_twct_worstcase_costs.jl --machine=$SLURMD_NODENAME --instances=instance_9_4_R100_wct_inputs.txt
