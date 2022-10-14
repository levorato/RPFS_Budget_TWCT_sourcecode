using ArgParse

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--machine"
            help = "machine name in the cluster"
		"--instances"
            help = "name of the instance group to solve"
		"--model"
            help = "robust model used to solve the CCG procedure"
		"--solver"
			help = "MILP solver to use in Master Problem : 'Gurobi' or 'CPLEX' (default)."
			default = "CPLEX"
		"--sp-solver"
			help = "MILP solver to use in Sub Problem : 'Gurobi' or 'CPLEX'. Default = solver."
			default = nothing
		"--budget-type"
			help = "Budgeted uncertainty type to use: 'global' (default) or 'machine'"
			default = "global"
		"--fractional-solution"
			help="Allow fractional budgeted uncertainty values in worst-case calculation."
			arg_type = Bool
			default = false
		"--single-cut"
			help = "If true, generates a single cut per CCG iteration."
			arg_type = Bool
			default = false
		"--time-limit"
			help = "Procedure time limit. Default = 7200 (seconds)."
			arg_type = Float64
			default = 7200.0
		"--max-cores"
			help = "Number of processor cores to use when solving the MILP problems. Default=16."
			arg_type = Int
			default = 16
		"--mp-bigm-type"
			help = "Type of bigM calculation to be used inside master MP MILPs: 1 or 2. Default=2."
			arg_type = Int
			default = 2
		"--sp-bigm-type"
			help = "Type of bigM calculation to be used inside sub SP MILPs: 1..7. Default=7."
			arg_type = Int
			default = 7
		"--log-cuts"
			help = "If true, write all generated cuts to disk, to allow recovery. Default=false."
			arg_type = Bool
			default = false
		"--sp-type"
			help = "CCG-Cmax SubProblem type: 'dp' for dynamic programming, 'mip' for worst-case MILP."
			default = "dp"
		"--alpha"
			help = "Instance alpha values to process. Default=-1 (process all files from instance group)."
			default = "-1"
		"--warm-start"
			help = "If true, the first CCG iteration will use GRASP solution. Default=false."
			arg_type = Bool
			default = false
		"--gamma"
			help = "Gamma budget value used to solve the problem. Default=-1 (solve for all gamma values)."
			arg_type = Int
			default = -1
		"--num-runs"
			help = "Number of replications (used in heuristics experiments, like GRASP). Default=20."
			arg_type = Int
			default = 20
		"--suffix"
			help = "Result filename suffix (used in parallel experiment executions)."
			default = ""
		"--hybrid"
			help = "Use Hybrid Combinatorial + MILP Branch-and-Bound algorithm for the PFSP-TWCT problem."
			arg_type = Bool
			default = false
		"--gap"
			help = "Relative gap tolerance to use inside C&CG procedure. Default = 1e-8."
			arg_type = Float64
			default = 1e-8
		"--seq"
			help = "Instance seq values to process. Default=-1 (process all files from instance group)."
			default = "-1"
    end
	parsed_args = parse_args(s)
	println("Parsed command-line args:")
	for (arg,val) in parsed_args
		println("  $arg  =>  $val")
	end
    return parsed_args
end

function get_machine_name()
    parsed_args = parse_commandline()
    println("Parsed args:")
    for (arg,val) in parsed_args
        println("  $arg  =>  $val")
    end
    return parsed_args["machine"]
end

function setup_gurobi_license()
	println("Setup Gurobi License...")
	if MACHINE != "laptop"
		println("Cluster environment detected, setting up Gurobi license...")
	    machine_name = get_machine_name()
		ENV["GRB_LICENSE_FILE"] = joinpath(home_prefix, "$(machine_name)", "gurobi.lic")
		println("GRB_LICENSE_FILE = $(ENV["GRB_LICENSE_FILE"])")
	else
		println("Local environment detected, ignoring Gurobi license...")
	end
	flush(stdout)
end
