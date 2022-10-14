# config.jl for Flowshop Scheduling Problem
# Prepares environment to run:

#using JuMP, MathProgBase, Gurobi
solver = "CPLEX"
if solver == "Gurobi"
  using JuMP, Gurobi
elseif solver == "GLPKLP" || solver == "GLPKMIP"
  using JuMP, GLPKMathProgInterface
elseif solver == "CPLEX"
  using  JuMP, CPLEX  #,JuMPeR

else
  println("NO SOLVER DEFINED")
end

MACHINE = "laptop"
USER_NAME = "ubuntu"
if MACHINE == "laptop"
  if Sys.isapple()
    home_prefix = joinpath("/Users", USER_NAME)
    git_prefix = "RPFS_Budget"
  else
    home_prefix = joinpath("/home", USER_NAME)
    git_prefix = "RPFS_Budget"
  end
else
  home_prefix = joinpath("/users", USER_NAME)
  git_prefix = "RPFS_Budget"
end
robust_ying_instances_folder = joinpath(home_prefix, git_prefix, "instances/robust/ying")
robust_tail_instances_folder = joinpath(home_prefix, git_prefix, "instances/robust/taillard")
instances_folder = joinpath(home_prefix, git_prefix, "instances")
project_folder = joinpath(home_prefix, git_prefix)
output_folder = joinpath(home_prefix, "pfsp_experiments")
println("*** Project folder is $(project_folder)")
println("*** Instances folder is $(instances_folder)")
println("*** Output folder is $(output_folder)")

searchdir(path,key) = filter(x->occursin(key,x), readdir(path))
