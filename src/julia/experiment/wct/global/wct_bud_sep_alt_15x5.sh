#!/bin/bash

#SBATCH --job-name=CGW15x5Manne
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=cpuonly
#SBATCH --mem=32G
#SBATCH --time=20000

echo "Starting job on node $SLURMD_NODENAME"
cd /users/ubuntu/RPFS_Budget_TWCT/src/julia
sh ./install_gurobi_license.sh $SLURMD_NODENAME

echo "RobPFSP Julia CCG WCT experiment script start."

cd /users/ubuntu/RPFS_Budget_TWCT/src/julia/experiment/wct/global
/users/ubuntu/julia-1.6.0/bin/julia exp_robust_pfsp_wct_budget_separation_main.jl --machine=$SLURMD_NODENAME --time-limit=14400 --instances=15x5 --model=manne --single-cut=true --fractional-solution=true --solver=Gurobi --budget-type=global --sp-bigm-type=7 --warm-start=true --max-cores=8

