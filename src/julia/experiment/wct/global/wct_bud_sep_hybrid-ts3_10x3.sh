#!/bin/bash

#SBATCH --job-name=CH10x3TS3
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=cpuonly
#SBATCH --mem=16G
#SBATCH --time=20000

echo "Starting job on node $SLURMD_NODENAME"
cd /users/ubuntu/RPFS_Budget_TWCT/src/julia
sh ./install_gurobi_license.sh $SLURMD_NODENAME

echo "RobPFSP Julia CCG WCT experiment script start."

cd /users/ubuntu/RPFS_Budget_TWCT/src/julia/experiment/wct/global
/users/ubuntu/julia-1.6.0/bin/julia exp_robust_pfsp_wct_budget_separation_main.jl --machine=$SLURMD_NODENAME --gap=1e-7 --time-limit=7200 --instances=10x3 --model=hybrid-ts3 --single-cut=false --fractional-solution=false --solver=Gurobi --budget-type=global --sp-bigm-type=7 --max-cores=16 --hybrid=true
