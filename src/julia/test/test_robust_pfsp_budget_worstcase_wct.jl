# TEST robust_pfsp_budget_worstcase_wct.jl
# Worst-case calculation for the robust budget problem
# Weighed Completion Time Objective Function

include("../config.jl")
include("../wct/robust_pfsp_budget_worstcase_wct_global.jl")
include("../wct/robust_pfsp_alt_worstcase_wct_mip.jl")
include("../wct/robust_pfsp_wilson_worstcase_wct_mip.jl")
include("../wct/robust_pfsp_budget_worstcase_wct_brute.jl")


using Gurobi
using CPLEX

EXPERIMENT_NAME = "robust_pfsp_wct_worstcase"
# If possible, do not modify the lines below
EXPERIMENT_OUTPUT_FOLDER = joinpath(normpath(output_folder), EXPERIMENT_NAME)
mkpath(EXPERIMENT_OUTPUT_FOLDER)

using ArgParse

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
		"--model"
			help = "robust model used to solve the worst-case WCT procedure: wilson or alt."
			default = "alt"
		"--machine"
			help = "machine name in the cluster"
		"--instances"
            help = "name of the instance group to solve"
			default = "10x5"
		"--solver"
			help = "MILP solver to use : 'Gurobi' (default) or 'CPLEX'"
			default = "Gurobi"
		"--fractional-solution"
			help="Allow fractional budgeted uncertainty values in worst-case calculation (default = true)."
			arg_type = Bool
			default = true
		"--time-limit"
			help = "Procedure time limit. Default = 1800 (seconds)."
			arg_type = Float64
			default = 1800.0
		"--max-cores"
			help = "Number of processor cores to use when solving the MILP problems. Default=8."
			arg_type = Int
			default = 2
		"--bigm-type"
			help = "Type of bigM calculation to be used inside sub SP MILPs: 1,...,7. Default=7."
			arg_type = Int
			default = 7
    end
	parsed_args = parse_args(s)
	println("Parsed command-line args:")
	for (arg,val) in parsed_args
		println("  $arg  =>  $val")
	end
    return parsed_args
end

function setup_gurobi_license(machine_name)
	println("Setup Gurobi License...")
	if MACHINE != "laptop"
		println("Cluster environment detected, setting up Gurobi license...")
	    ENV["GRB_LICENSE_FILE"] = joinpath(home_prefix, "$(machine_name)", "gurobi.lic")
		println("GRB_LICENSE_FILE = $(ENV["GRB_LICENSE_FILE"])")
	else
		if Sys.isapple()
			println("Local environment detected, ignoring Gurobi license...")
		else
			ENV["GRB_LICENSE_FILE"] = "/home/ubuntu/gurobi.lic"
		end
	end
	flush(stdout)
end

function test_m_machine_budget_worst_case_dp_3machines()
    searchdir(path,key) = filter(x->occursin(key,x), readdir(path))
    basedir = "/home/ubuntu/RPFS_Budget_TWCT/instances/robust/ying/data"
    for filename in searchdir(basedir, ".txt")
        if filename == "RB0104010_20perc_3m.txt"
            inputfile = "$(basedir)/$(filename)"
            m, n, P_bar, P_hat = read_robust_ying_input_file(inputfile)
            # φ_star is the optimal Johnson solution obtained over the nominal processing time array (P_bar)
            #φ_star, pi, x = solve_deterministic_pfsp(m, n, P_bar)
            #println("Based on the nominal processing times, the optimal makespan φ_star is $(φ_star).")
            #println("x: $(x)")
            Jobs = 1:n
            N = Jobs
            Machines = 1:m
            M = Machines
            # pi = [ 1, 5, 10, 2, 3, 7, 4, 6, 9, 8 ]  # Tested up to [3,0,0]
            pi = [ 10, 1, 7, 5, 3, 2, 4, 6, 9, 8]
            w = [1.0 for x in 1:n]  # For now, assume all job weights are equal to 1
            println("Calculating worst-case makespan for sequence pi = $(pi)...")
            error_count = 0
            for Γ1 in 0:n
                for Γ2 in 0:n
                    for Γ3 in 0:n
                        println("*** p_Γ1 = $(Γ1); p_Γ2 = $(Γ2); p_Γ3 = $(Γ3)")
                        #Γ1 = (p_Γ1 * n) / 100.0
                        #Γ2 = (p_Γ2 * n) / 100.0
                        wct_brute, worst_Γ_scenario_brute, worst_j = worst_case_wct_brute_force_3_machines(pi, m, n, [Γ1, Γ2, Γ3], P_bar, P_hat, w, false)
                        wct_dp, worst_Γ_scenario_dp = worst_case_wct_dp_m_machines_2(pi, m, n, [Γ1, Γ2, Γ3], P_bar, P_hat, w, false)
                        #println("worst_Γ_scenario_brute = $(worst_Γ_scenario_brute)")
                        if abs(wct_dp - wct_brute) > EPS
                            println("ERROR : wct_dp incorrect : wct_brute = $(wct_brute) x wct_dp = $(wct_dp)")
                            error_count += 1
                            break
                        end
                    end
                    if error_count > 0
                        break
                    end
                end
                if error_count > 0
                    break
                end
            end
        end
    end
end

function test_time_mip_worstcase()
	println("Testing RPFS-TWCT MIP worstcase parameters...")
	parsed_args = parse_commandline()
	instance_group = parsed_args["instances"]
	machine_name = parsed_args["machine"]
	solver = parsed_args["solver"]
	time_limit = parsed_args["time-limit"]
	max_cores = parsed_args["max-cores"]
	fractional_solution = parsed_args["fractional-solution"]
	bigM_type = parsed_args["bigm-type"]
	model = parsed_args["model"]
	if (!haskey(parsed_args, "instances")) | isnothing(instance_group)
		println("ERROR: instances argument is mandatory!")
		return
	end
	if (!haskey(parsed_args, "model")) | isnothing(model)
		println("ERROR: model argument is mandatory!")
		return
	end
	setup_gurobi_license(machine_name)
	# Open output CSV file for results
    result_file_path = joinpath(EXPERIMENT_OUTPUT_FOLDER, "test_rpfs_wct_worstcase_time_spent_$(model)_$(instance_group)_bigm-$(bigM_type)_frac-$(fractional_solution).csv")
    result_file = open(result_file_path, "a+")
    instance_dict = Dict()
    file_list = ["RB0101004_10_5_10_wct_inputs.txt", "RB0101007_10_5_10_wct_inputs.txt", "RB0102004_10_5_20_wct_inputs.txt"]
    basedir = joinpath(robust_ying_instances_folder, "rob-pfsp-wct", instance_group)
    for filename in searchdir(basedir, ".txt")
        if !(filename in file_list) | true
            inputfile = joinpath("$(basedir)", "$(filename)")
            instance_dict[filename] = inputfile
        end
    end
    println(result_file, "instance, gamma_perc, gamma, fractional_solution, bigM_type, solver, time_spent, wct_mip, wct_validation")
	total_time = 0.0
	warm_start = true
    for (instance, filepath) in instance_dict
        m, n, P_bar, P_hat, w = read_robust_input_file(filepath)
        Jobs = 1:n
        N = Jobs
        Machines = 1:m
        M = Machines
        pi = [ x for x in 1:n ]
        error = false
        for p_Γ in reverse([10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
            Γ = (p_Γ * m * n) / 100.0
            println("Calculating worst-case makespan for Γ = $Γ and sequence pi = $(pi)...")
            flush(stdout)
            #time_spent = @elapsed x = worst_case_wct_brute_force_global_budget(pi, m, n, Γ, P_bar, P_hat, w, true)
            time_spent_dict = Dict()
			if model == "alt"
				initial_scenario = nothing
				if warm_start
					println("Using DP procedure to warm start WCL MILP model.")
					max_wct, worst_Γ_scenario = worst_case_wct_dp_m_machines(pi, m, n, Γ, P_bar, P_hat, w)
					initial_scenario = convert_scenario_job_deviation_list_to_matrix(m, n, worst_Γ_scenario)
				end
				wct_mip, worst_Γ_scenario_mip, time_spent = solve_robust_model_worst_case_alt(pi, m, n, Γ, P_bar, P_hat, w, 1;
                        solver=solver, budget_type = "global", solver_time_limit = time_limit, single_cut_per_iteration = true,
                        fractional_solution = fractional_solution, max_cores=max_cores, bigM_type=bigM_type, initial_scenario=initial_scenario)
			else  # wilson
				wct_mip, worst_Γ_scenario_mip, time_spent = solve_robust_model_worst_case_wilson(pi, m, n, Γ, P_bar, P_hat, w, 1;
                        solver=solver, budget_type = "global", solver_time_limit = time_limit, single_cut_per_iteration = true,
                        fractional_solution = fractional_solution, max_cores=max_cores, bigM_type=bigM_type)
			end
            wct_validation = 0.0
            if size(worst_Γ_scenario_mip, 1) > 0
                P = calculate_P_given_deviation_matrix(m, n, P_bar, P_hat, worst_Γ_scenario_mip[1])
                wct_validation = calculate_wct(m, n, P, w, pi)
            end
            time_spent_dict[(fractional_solution, bigM_type, solver)] = (time_spent, wct_mip, wct_validation)
            println("[($fractional_solution, $bigM_type, $solver)] = ($time_spent, $wct_mip, $wct_validation)")
            println("wct_mip = $(wct_mip) ; wct_validation = $(wct_validation)")
			total_time += time_spent
            #if abs(wct_mip - wct_validation) > 1e-5
            #    error = true
            #    break
            #end
            #if size(worst_Γ_scenario_mip, 1) == 0
            #    error = true
            #    break
            #end
            #if error  break  end
            for (key, value) in time_spent_dict
                key_str = "$(key)"
                key_str = replace(key_str, "(" => "")
                key_str = replace(key_str, ")" => "")
                value_str = "$(value)"
                value_str = replace(value_str, "(" => "")
                value_str = replace(value_str, ")" => "")
                println(result_file, "$(instance), $(p_Γ), $(Γ), $(key_str), $(value_str)")
                flush(result_file)
            end
            println(time_spent_dict)
        end
        if error
            break
        end
    end
    close(result_file)
	println("Total time spent: $(total_time) s")
    println("DONE.")
    flush(stdout)
end

test_time_mip_worstcase()

# test_m_machine_budget_worst_case_dp_3machines()
#test_time_mip_worstcase([false, true], [7], ["Gurobi"], ["10x5"])  # "20x5"
#test_time_mip_worstcase([false, true], [5, 6], ["Gurobi"], ["15x5"])  # "20x5"
#test_time_mip_worstcase([false, true], [1, 2, 3], ["Gurobi"], ["15x5"])  # "20x5"
# test_time_mip_worstcase([true], [2], ["Gurobi"], ["20x5"])
