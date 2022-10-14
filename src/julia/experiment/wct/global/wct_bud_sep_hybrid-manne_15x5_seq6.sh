#!/bin/bash

#SBATCH --job-name=CW15-HM6
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=cpuonly
#SBATCH --mem=32G
#SBATCH --time=50000

echo "Starting job on node $SLURMD_NODENAME"
cd /users/ubuntu/RPFS_Budget_TWCT/src/julia
sh ./install_gurobi_license.sh $SLURMD_NODENAME

echo "RobPFSP Julia CCG WCT experiment script start."

cd /users/ubuntu/grasp/RPFS_Budget_TWCT/src/julia/experiment/wct/global
/users/ubuntu/julia-1.6.0/bin/julia exp_robust_pfsp_wct_budget_separation_main.jl --max-cores=16 --machine=$SLURMD_NODENAME --time-limit=14400 --instances=15x5 --model=hybrid-manne --single-cut=true --fractional-solution=false --solver=Gurobi --budget-type=global --sp-bigm-type=7 --hybrid=true --gap=1e-7 --seq=06

