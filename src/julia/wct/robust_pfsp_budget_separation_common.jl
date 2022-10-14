# ====================================================================
# robust_pfsp_budget_separation_common.jl
# Common functions for Rob PFSP BUdget Separation procedures.
# ====================================================================

include("../pfsp_util.jl")
include("robust_pfsp_budget_worstcase_wct_brute.jl")
include("robust_pfsp_budget_worstcase_wct_global.jl")

using CSV
using DataFrames
using Dates

# Minimal time limit for the execution of the SubProblem MILP model
MIN_TIME_SP_MIP = 1800
# Folder containing warm start solutions to use in C&CG procedures
WARM_START_SOLUTION_FOLDER = abspath(create_full_dir(normpath(output_folder),
	["run_bb_rpfs_wct_global"]))


function solve_separation_problem(m, n, Γ, P_bar, P_hat, w, permutation, wct_star, verbose = true; max_cores = 16,
            budget_type = "machine", time_limit = 7200, single_cut_per_iteration = false, fractional_solution = false,
            use_brute_force = true, solver = "CPLEX", bigM_type = 7, sp_warm_start=false)
    wct = 0.0
    list_worst_Γ_scenario_z_sep = []
    time_spent = 0.0
    if budget_type == "global"
		if worst_case_brute_force_is_worth(m, n, Int64(ceil(Γ))) && use_brute_force
			max_solutions = 5
			if single_cut_per_iteration
				max_solutions = 1
			end
	        time_spent = @elapsed wct, list_worst_Γ_scenario = worst_case_wct_brute_force_global_budget(permutation, m, n, Γ, P_bar, P_hat, w, true;
																		max_solutions=max_solutions)
	        println("Separation time: $(time_spent) s")
	        flush(stdout)
	        list_worst_Γ_scenario_z_sep = [convert_scenario_job_deviation_list_to_matrix(m, n, s) for s in list_worst_Γ_scenario]
			is_opt = size(list_worst_Γ_scenario_z_sep, 1) > 0
		    return wct, is_opt, list_worst_Γ_scenario_z_sep, time_spent
		end
	end
	initial_scenario = nothing
	if sp_warm_start
		println("Using DP procedure to warm start WCL MILP model.")
		max_wct, worst_Γ_scenario = worst_case_wct_dp_m_machines(permutation, m, n, Γ, P_bar, P_hat, w)
		initial_scenario = convert_scenario_job_deviation_list_to_matrix(m, n, worst_Γ_scenario)
	end
    wct, list_worst_Γ_scenario_z_sep, time_spent = solve_robust_model_worst_case_alt(permutation, m, n, Γ, P_bar, P_hat, w;
        budget_type=budget_type, solver_time_limit=time_limit, single_cut_per_iteration=single_cut_per_iteration,
        fractional_solution=fractional_solution, solver=solver, max_cores=max_cores, bigM_type=bigM_type,
		initial_scenario=initial_scenario)
    is_opt = size(list_worst_Γ_scenario_z_sep, 1) > 0
    return wct, is_opt, list_worst_Γ_scenario_z_sep, time_spent
end

function get_current_model_solution_assignment(model, n)
    Jobs = 1:n
    Z_star = zeros(Float64, (n, n))
    if JuMP.termination_status(model) == MOI.OPTIMAL
        #println("\n===========  R O B U S T    S O L U T I O N  ===========\n")
        #println("Optimal Objective Function value: ", objective_value(model))
        #println("Solve time : ", solve_time(model))
        permutation = zeros(Int64, n)
        for k in Jobs  # Position k
            for j in Jobs  # Job j
                Z_jk = variable_by_name(model, "Z[$(j),$(k)]")
                if value(Z_jk) >= 0.9  # Z[j, k]
                    Z_star[j, k] = 1
                    if permutation[k] == 0
                        permutation[k] = j
                    else
                        println("ERROR assembling permutation from model solution!")
                    end
                end
            end
        end
        return objective_value(model), permutation, Z_star, true
    else
        println("ERROR: [get_current_model_solution] Optimal solution not found! Best bound: ", objective_bound(model))
        println("Optimal solution not found! Best bound: ", objective_bound(model))
        println("Model solve time: $(solve_time(model)); time_limit = $(time_limit)")
        flush(stdout)
        return objective_bound(model), zeros(Int64, n), Z_star, false
    end
end

function get_current_model_solution_dichotomous(model, n)
    Jobs = 1:n
    D_star = zeros(Float64, (n, n))
    if JuMP.termination_status(model) == MOI.OPTIMAL
        #println("\n===========  R O B U S T    S O L U T I O N  ===========\n")
        #println("Optimal Objective Function value: ", objective_value(model))
        #println("Solve time : ", solve_time(model))
        for i in Jobs  # Job i
            for k in Jobs  # Job k
                D_ik = variable_by_name(model, "D[$(i),$(k)]")
                if value(D_ik) >= 0.9  # D[i, k]
                    D_star[i, k] = 1
                end
            end
        end
        permutation = topological_sort_by_dfs(D_star, n)
        return objective_value(model), permutation, D_star, true
    else
        println("ERROR: [get_current_model_solution] Optimal solution not found! Best bound: ", objective_bound(model))
        println("Optimal solution not found! Best bound: ", objective_bound(model))
        println("Model solve time: $(solve_time(model)); time_limit = $(time_limit)")
        flush(stdout)
        return objective_bound(model), zeros(Int64, n), D_star, false
    end
end

function extract_bb_solution_info(sol_str, hybrid)
	result_dict = Dict()
	try
		r = findlast("Objective value:", sol_str)
		s = sol_str[r[2]:length(sol_str)]
		result_dict["solution_value"] = parse(Float64, s[findfirst(':', s) + 2:findfirst('\n', s) - 1])

		r = findlast("Permutation:", sol_str)
		s = sol_str[r[2]:length(sol_str)]
		result_dict["permutation"] = s[findfirst('[', s) + 2:findfirst(']', s) - 1]

		r = findlast("Time Spent (s):", sol_str)
		s = sol_str[r[2]:length(sol_str)]
		result_dict["time_spent"] = parse(Float64, s[findfirst(':', s) + 2:findfirst('\n', s) - 1])

		if !hybrid
			r = findlast("Number of processed nodes (NN):", sol_str)
			s = sol_str[r[2]:length(sol_str)]
			result_dict["iterations"] = parse(Int64, s[findfirst(':', s) + 2:findfirst('\n', s) - 1])

			r = findlast("Total number of phatomed nodes:", sol_str)
			s = sol_str[r[2]:length(sol_str)]
			result_dict["num_phatomed"] = parse(Int64, s[findfirst(':', s) + 2:findfirst('\n', s) - 1])
		else
			r = findlast("Iterations:", sol_str)
			s = sol_str[r[2]:length(sol_str)]
			result_dict["iterations"] = parse(Int64, s[findfirst(':', s) + 2:findfirst('\n', s) - 1])
			result_dict["num_phatomed"] = 0
		end

		r = findlast("Number of improvements:", sol_str)
		s = sol_str[r[2]:length(sol_str)]
		result_dict["num_improvements"] = parse(Int64, s[findfirst(':', s) + 2:findfirst('\n', s) - 1])

		r = findlast("execution_id:", sol_str)
		s = sol_str[r[2]:length(sol_str)]
		result_dict["execution_id"] = s[findfirst(':', s) + 2:findfirst('\n', s) - 1]

		r = findlast("seed:", sol_str)
		s = sol_str[r[2]:length(sol_str)]
		result_dict["seed"] = s[findfirst(':', s) + 2:findfirst('\n', s) - 1]

		r = findlast("gap:", sol_str)
		s = sol_str[r[2]:length(sol_str)]
		result_dict["gap"] = parse(Float64, strip(s[findfirst(':', s) + 2:findfirst('\n', s) - 1]))

		r = findlast("is_optimal:", sol_str)
		s = sol_str[r[2]:length(sol_str)]
		result_dict["is_optimal"] = parse(Bool, strip(s[findfirst(':', s) + 2:findfirst('\n', s) - 1]))

		r = findlast("validated:", sol_str)
		s = sol_str[r[2]:length(sol_str)]
		result_dict["validated"] = parse(Bool, strip(s[findfirst(':', s) + 2:findfirst('\n', s) - 1]))
		return result_dict
	catch e
		bt = catch_backtrace()
		msg = sprint(showerror, e, bt)
		println("[extract_bb_solution_info] ERROR extracting B & B C++ solution info.")
		println(msg)
		# solver error
		result_dict["solution_value"] = 0
		result_dict["permutation"] = "[]"
		result_dict["time_spent"] = 0
		result_dict["time_to_best_sol"] = 0
		result_dict["iterations"] = 0
		result_dict["num_phatomed"] = 0
		result_dict["num_improvements"] = 0
		if length(sol_str) > 0 && occursin("execution_id:", sol_str)
			r = findlast("execution_id:", sol_str)
			s = sol_str[r[2]:length(sol_str)]
			result_dict["execution_id"] = s[findfirst(':', s) + 2:findfirst('\n', s) - 1]
		else
			result_dict["execution_id"] = "unknown"
		end
		result_dict["seed"] = 0
		result_dict["gap"] = 0
		result_dict["is_optimal"] = false
		result_dict["validated"] = false
		return result_dict
	end
end

function extract_grasp_solution_info(sol_str, success)
	result_dict = Dict()
	try
		r = findlast("Objective function value:", sol_str)
		s = sol_str[r[2]:length(sol_str)]
		result_dict["solution_value"] = parse(Float64, s[findfirst(':', s) + 2:findfirst('\n', s) - 1])

		r = findlast("Permutation:", sol_str)
		s = sol_str[r[2]:length(sol_str)]
		result_dict["permutation"] = s[findfirst(':', s) + 2:findfirst('\n', s) - 1]

		r = findlast("Time Spent (s):", sol_str)
		s = sol_str[r[2]:length(sol_str)]
		result_dict["time_spent"] = parse(Float64, s[findfirst(':', s) + 2:findfirst('\n', s) - 1])

		r = findlast("ime to best solution (s):", sol_str)
		s = sol_str[r[2]:length(sol_str)]
		result_dict["time_to_best_sol"] = parse(Float64, s[findfirst(':', s) + 2:findfirst('\n', s) - 1])

		r = findlast("Iterations:", sol_str)
		s = sol_str[r[2]:length(sol_str)]
		result_dict["iterations"] = parse(Int64, s[findfirst(':', s) + 2:findfirst('\n', s) - 1])

		r = findlast("Number of visited solutions:", sol_str)
		s = sol_str[r[2]:length(sol_str)]
		result_dict["num_visited_solutions"] = parse(Int64, s[findfirst(':', s) + 2:findfirst('\n', s) - 1])

		r = findlast("Number of improvements:", sol_str)
		s = sol_str[r[2]:length(sol_str)]
		result_dict["num_improvements"] = parse(Int64, s[findfirst(':', s) + 2:findfirst('\n', s) - 1])

		r = findlast("execution_id:", sol_str)
		s = sol_str[r[2]:length(sol_str)]
		result_dict["execution_id"] = s[findfirst(':', s) + 2:findfirst('\n', s) - 1]

		r = findlast("seed:", sol_str)
		s = sol_str[r[2]:length(sol_str)]
		result_dict["seed"] = s[findfirst(':', s) + 2:findfirst('\n', s) - 1]
		return result_dict
	catch e
		bt = catch_backtrace()
		msg = sprint(showerror, e, bt)
		println("[extract_grasp_solution_info] ERROR extracting GRASP C++ solution info.")
		println(msg)
		# solver error
		result_dict["solution_value"] = 0
		result_dict["permutation"] = "[]"
		result_dict["time_spent"] = 0
		result_dict["time_to_best_sol"] = 0
		result_dict["iterations"] = 0
		result_dict["num_visited_solutions"] = 0
		result_dict["num_improvements"] = 0
		if length(sol_str) > 0 && occursin("execution_id:", sol_str)
			r = findlast("execution_id:", sol_str)
			s = sol_str[r[2]:length(sol_str)]
			result_dict["execution_id"] = s[findfirst(':', s) + 2:findfirst('\n', s) - 1]
		else
			result_dict["execution_id"] = "unknown"
		end
		result_dict["seed"] = 0
		return result_dict
	end
end

function string_to_permutation(permutation)
	s = "$(permutation)"
	s = replace(s, "Any[" => "")
	s = replace(s, "[" => "")
	s = replace(s, "]" => "")
	s = replace(s, "," => "")
	s = strip(s)
	return [parse(Int, ss) for ss in split(s)]
end

"""
Run `cmd`, and collect the output such that `stdout` and `stderr` are captured
independently, but with the time of each line recorded such that they can be
stored/analyzed independently, but replayed synchronously.
https://discourse.julialang.org/t/collecting-all-output-from-shell-commands/15592
https://github.com/mbaz/Gaston.jl/blob/a272e736f93920861249e5bb8364ff8c63e3be29/src/gaston_llplot.jl#L7
"""
function communicate(cmd::Base.AbstractCmd, stdout_file_path, stderr_file_path)
    #out_pipe = Pipe()
    #err_pipe = Pipe()
	process = try
        run(pipeline(cmd, stdin=devnull, stdout=stdout_file_path, stderr=stderr_file_path); wait=false)
    catch
        @warn("Could not spawn $(cmd)")
        rethrow()
    end
    #close(out_pipe.in)
    #close(err_pipe.in)
    #stdout = @async String(read(out_pipe))
    #stderr = @async String(read(err_pipe))
    wait(process)
    return (
        #stdout = fetch(stdout),
        #stderr = fetch(stderr),
        code = process.exitcode
    )
end

function generate_hybrid_master_output_path(exp_output_folder, optimization_type, instance_filename)
	this_output_folder = abspath(create_full_dir(exp_output_folder, ["pfsp_wct", optimization_type, "$(instance_filename)"]))
	output_folder_stdout = abspath(create_full_dir(this_output_folder, ["stdout"]))
	return this_output_folder, output_folder_stdout
end

# GRASP obj-update-freq parameter values:
# 0 == always recalculate, i.e., use worst-case MIP (exact solution);
# 1 == no recalculation (use DP v1.0 all the time). Only use worst-case MIP (exact solution) to give final result;
# 3 == obj recalculated at the end of each MH iteration (use DP v1.0 otherwise);
# 5 == obj recalculated at the end of each MH iteration and also at each VND improvement (use DP v1.0 otherwise);
function solve_pfsp_wct_cpp(filepath, Γ, budget_type; optimization_type = "grasp", time_limit = 200,
		cut_list = [], max_cores = 16)
	filename = basename(filepath)
	output_folder2, output_folder_stdout = generate_hybrid_master_output_path(EXPERIMENT_OUTPUT_FOLDER, optimization_type, filename)
	filename = basename(filepath)
	println("=========================================================")
	println(" Solving instance $(filename)...")
	println(" (optimization_type = $(optimization_type))")
	println("=========================================================")
	datetimenow = Dates.now()
	println("Date: $(datetimenow)")
	flush(stdout)
	m, n, P_bar, P_hat, w = read_robust_input_file(filepath)
	# Create an output text file to save the stdout and sterr process output
	iteration = size(cut_list, 1)
	pid = getpid()
	stdout_file_path = joinpath(output_folder_stdout, "pfsp_cpp_$(pid)_$(filename)_gamma-$(Γ)_$(budget_type)_$(iteration)_stdout.txt")
	stderr_file_path = joinpath(output_folder_stdout, "pfsp_cpp_$(pid)_$(filename)_gamma-$(Γ)_$(budget_type)_$(iteration)_stderr.txt")
	# Hybrid BB - Create a temp file with the list of cuts
	cuts_file_path = joinpath(output_folder_stdout, "pfsp_cpp_$(pid)_$(filename)_gamma-$(Γ)_$(budget_type)_cuts-list.txt")
	cuts_file = open(cuts_file_path, "w+")
	# write cuts to cuts_filepath
	for cut in cut_list
		_v = collect(Iterators.flatten(cut'))
		println(cuts_file, "$(_v)")
	end
	close(cuts_file)
	println("[solve_pfsp_wct_cpp] Cuts file written to: $(cuts_file_path)")
	println("*** Gamma = $(Γ)")
	Γ_perc = Int64(floor((Γ / (m * n)) * 100.0))
	println("*** Gamma % = $(Γ_perc)")
	scenario_processing_start_time = time_ns()
	flush(stdout)
	result_dict = Dict()
	first_improvement = true
	vnd_size = 4
	random_vnd = true
	adaptive = false
	beta1 = 0.4
	beta2 = 1.0
	if optimization_type == "grasp"
		model_version = 604
		time_factor = 3000
		obj_update_freq = 4
	else   # WCT Hybrid combinatorial branch-and-bound / C&CG Master Problem MP
		model_version = 607
		time_factor = 300
		obj_update_freq = 3
		Γ_perc = 0  # The MP does not use a budget parameter, only a list of cuts
	end
	cmd = ``
	try
		execpath = abspath(joinpath(project_folder, "bin", "pfsp"))
		cmd = `$(execpath) --input-file="$(filepath)"
			--output-folder $(output_folder2)
			--model-version=$(model_version)
			--budget-gamma="$(Γ_perc)"
			--rob-ub-type=grasp
			--first-improvement=$(first_improvement)
			--vnd-size=$(vnd_size)
			--random-vnd=$(random_vnd)
			--adaptive-construction=$(adaptive)
			--grasp-tf=$(time_factor)
			--const-beta1=$(beta1)
			--const-beta2=$(beta2)
			--budget-type=$(budget_type)
			--obj-update-freq=$(obj_update_freq)
			--time-limit=$(time_limit)
			--cuts-file="$(cuts_file_path)"
			--sel-strategy=cplex
			--mip-use-c-bounds=true
			--max-cores=$(max_cores)
			--mp-model-name=$(optimization_type)
			--dominance=false
			--seed=42`
		println("[solve_pfsp_wct_cpp] cmd: $(cmd)")
		flush(stdout)
		# Invoke process via command-line
		@show exitcode = communicate(cmd, stdout_file_path, stderr_file_path)
		#oc = OutputCollector(cmd; verbose = true)
		#stdout_str = collect_stdout(oc)
		#stderr_str = collect_stderr(oc)
		#println(stdout_file, "*** stdout: ", stdout_str, "\n")
		#println(stdout_file, "*** stderr: ", stderr_str, "\n")
		#println(stdout_file, "*** exit status: ", proc.exitcode, "\n")
		println("[solve_pfsp_wct_cpp] C++ Invocation finished.")
		datetimenow = Dates.now()
		println("Date: $(datetimenow)")
		flush(stdout)
		#flush(stdout_file)
		#exit_code = proc.exitcode
		# Read stdout file to solution_string
		solution_string = read(stdout_file_path, String)
		if optimization_type == "grasp"
			result_dict = extract_grasp_solution_info(solution_string)
			result_dict["gap"] = -100.0
			result_dict["is_optimal"] = false
		else
			result_dict = extract_bb_solution_info(solution_string, true)
		end
	catch e
		elapsed_time = (time_ns() - scenario_processing_start_time) * 1e-9
		bt = catch_backtrace()
		msg = sprint(showerror, e, bt)
		println(msg)
		stdout_file = open(stdout_file_path, "a+")
		println(stdout_file, "[solve_pfsp_wct_cpp] cmd: $(cmd)")
		println(stdout_file, "*** [solve_pfsp_wct_cpp] EXCEPTION OCCURED: $(msg)")
		println(stdout_file, "*** [solve_pfsp_wct_cpp] EXCEPTION OCCURED: $(msg)")
		flush(stdout)
		flush(stdout_file)
		close(stdout_file)
		exit(1)
	end
	if result_dict["permutation"] == "[]"
		stdout_file = open(stdout_file_path, "a+")
		println("stdout_file, [solve_pfsp_wct_cpp] cmd: $(cmd)")
		println(stdout_file, "*** [solve_pfsp_wct_cpp] PFSP C++ RETURNED ERROR.")
		println("ERROR: *** [solve_pfsp_wct_cpp] PFSP C++ RETURNED ERROR.")
		flush(stdout)
		flush(stdout_file)
		close(stdout_file)
		exit(1)
	end
	elapsed_time = (time_ns() - scenario_processing_start_time) * 1e-9
    return result_dict["solution_value"], string_to_permutation(result_dict["permutation"]),
		elapsed_time, result_dict["gap"], result_dict["is_optimal"]
end

function create_jump_model(solver)
	if solver == "Gurobi"
		println("Using Gurobi solver.")
		model = Model(Gurobi.Optimizer)
	elseif solver == "CPLEX"
		println("Using CPLEX solver.")
		model = Model(CPLEX.Optimizer)
	else
	  println("No solver defined")
	  model = nothing
	end
	return model
end

function setup_solver_params(model, solver, solver_time_limit, max_cores, ub = typemax(Float64))
	if solver == "Gurobi"
		MOI.set(model, MOI.RawOptimizerAttribute("TimeLimit"), solver_time_limit)
		MOI.set(model, MOI.RawOptimizerAttribute("Threads"), max_cores)
		# https://www.gurobi.com/documentation/9.1/refman/mipgapabs.html
		MOI.set(model, MOI.RawOptimizerAttribute("MIPGapAbs"), 1e-5)   # default is 1e-10
		#MOI.set(model, MOI.RawOptimizerAttribute("Cutoff"), ub)
		# https://www.gurobi.com/documentation/9.1/refman/integralityfocus.html#parameter:IntegralityFocus
		#MOI.set(model, MOI.RawOptimizerAttribute("IntegralityFocus"), 1)
		# https://www.gurobi.com/documentation/9.1/refman/presolve.html#parameter:Presolve
		#MOI.set(model, MOI.RawOptimizerAttribute("Presolve"), 2)  # agressive = 2, default = -1
		# https://www.gurobi.com/documentation/9.1/refman/mipfocus.html#parameter:MIPFocus
		#MOI.set(model, MOI.RawOptimizerAttribute("MIPFocus"), 2)  # proving optimality = 2, default = 0
		# https://www.gurobi.com/documentation/9.1/refman/heuristics.html#parameter:Heuristics
		#MOI.set(model, MOI.RawOptimizerAttribute("Heuristics"), 0.10)  # 10% of time spent on heuristics, default = 5%
	elseif solver == "CPLEX"
		MOI.set(model, MOI.RawOptimizerAttribute("CPX_PARAM_TILIM"), solver_time_limit)
		MOI.set(model, MOI.RawOptimizerAttribute("CPX_PARAM_THREADS"), max_cores)
		MOI.set(model, MOI.RawOptimizerAttribute("CPX_PARAM_EPAGAP"), 1e-5)   # default is 1e-6
		# Feasibility Pump Heuristic: 1 or 2
		MOI.set(model, MOI.RawOptimizerAttribute("CPXPARAM_MIP_Strategy_FPHeur"), 2)
		#MOI.set(model, MOI.RawOptimizerAttribute("CPXPARAM_MIP_Tolerances_UpperCutoff"), ub)
	else
	  println("No solver defined")
	end
end

function read_initial_results_from_csv(result_filename)
	println("Reading result file: $(result_filename)")
	if isfile(result_filename)
	    println("Result file exists.")
	    df = DataFrame(CSV.File(result_filename, delim=','))
		#println("First 6 results: ")
		#println(first(df, 6))
		flush(stdout)
	    return df
	else
		println("No existing result file found!")
		flush(stdout)
		return DataFrame()
	end
end

function read_warm_start_solutions(instance_group, instance_name)
	result_file = joinpath(WARM_START_SOLUTION_FOLDER, "bb_rpfs_wct_$(instance_group)_fi_true_vs_4_rv_true.csv")
	df = read_initial_results_from_csv(result_file)
	nrows, ncols = size(df)
	solution_list = []
	elapsed_time = 0.0
	if nrows > 0
		df[!, :instance_name] = strip.(df[!, :instance_name])
		if (typeof(df[!,:time_spent][1]) != Float64) & (typeof(df[!,:time_spent][1]) != Int64)
			df[!,:time_spent] = strip.(df[!,:time_spent])
			df[!,:time_spent] = map(x->(v = tryparse(Float64,x); v == nothing ? -1.0 : v), df[!,:time_spent])
		end
		df = df[df[!, :instance_name] .== instance_name, :]
		nrows, ncols = size(df)
		for row in 1:nrows
			permutation_str = strip(df[row, :permutation])
			permutation = string_to_permutation(permutation_str)
			if length(permutation) > 1
				push!(solution_list, permutation)
				elapsed_time += df[row, :time_spent]
				# Return only 1 solution if instance_group == 15x5
				if occursin("15", instance_group)
					println("[Warm-start] Returning only one warm-start solution, since instance_group == 15x5.")
					break
				end
			end
		end
	end
	if length(solution_list) == 0
		println("WARN: no warm-start solutions found!")
	end
	return solution_list, elapsed_time
end

function generate_initial_solutions(filepath, m, n, Γ, budget_type)
	println("filepath = $(filepath)")
	instance_type = "tail"
	if occursin("ying", filepath)
		instance_type = "ying"
	end
	instance_group = instance_type * "_$(n)x$(m)"
	instance_name = basename(filepath)
	initial_solutions, elapsed_time = read_warm_start_solutions(instance_group, instance_name)
	if length(initial_solutions) == 0
		println("WARN: generating initial solution with GRASP.")
		obj, permutation_current, elapsed_time = solve_pfsp_wct_cpp(filepath, Γ, budget_type; optimization_type="grasp")
		initial_solutions = [permutation_current]
		println("\n===========  R O B U S T    S O L U T I O N  -  G R A S P ===========\n")
		println("[GRASP] Optimal Objective Function value: ", obj)
		println("[GRASP] Solve time : ", elapsed_time)
	else
		println("\n===========  R O B U S T    S O L U T I O N  -  W A R M   S T A R T ===========\n")
		println("[WARM] Solve time : ", elapsed_time)
	end
	println("Initial solutions generated: $(initial_solutions)")
	flush(stdout)
	return initial_solutions, elapsed_time
end

function generate_permutation_from_Z_matrix(Z, n)
	permutation = zeros(Int64, n)
	Z_star = zeros(Float64, (n, n))
	for k in 1:n
		for j in 1:n
			if value(Z[j, k]) >= 0.9
				Z_star[j, k] = 1
				permutation[k] = j
			end
		end
	end
	return permutation, Z_star
end

function set_warm_start_assignment_model(Z, n, permutation)
	for k in 1:n
		for j in 1:n
			set_start_value(Z[k, j], 0)
		end
	end
	for k in 1:n
		set_start_value(Z[permutation[k], k], 1)
	end
end

# D_ij = 1, if job i is scheduled any time before job j; 0, otherwise.
function set_warm_start_dichotomous_model(D, n, permutation)
	for k in 1:n
		for j in 1:n
			set_start_value(D[k, j], 0)
		end
	end
	for k1 in 1:(n-1)
		i = permutation[k1]
		for k2 in (k1+1):n
			j = permutation[k2]
			# Job i is scheduled any time before job j
			set_start_value(D[i, j], 1)
			set_start_value(D[j, i], 0)
		end
	end
end

function read_twct_cuts_from_file(filename)
	cut_list = []
	open(filename) do io
		while !eof(io)
			line = readline(io)
			# remove open and closing brackets, if exists
			line = replace(line, "[" => "")
			line = replace(line, "]" => "")
			line = replace(line, "," => "")
			cut = readdlm(IOBuffer(line))
			push!(cut_list, cut)
		end
	end
	return cut_list
end

function load_cut_list_and_partial_solution_twct(partial_results_full_filename, log_cuts_dir, instance_name, Γ)
    existing_cut_list = []
	max_iter = 1
	UB = typemax(Float64)
	best_pi = []
	gamma_str = "_gamma-$(Γ)_"
    println("Reading existing cuts from dir: $(log_cuts_dir)")
	for filename in readdir(log_cuts_dir)
		# Filename pattern: pfsp_cpp_4400_instance_15_5_R100_wct_inputs.txt_gamma-11.25_global_cuts-list.txt
		if (startswith(filename, "pfsp_cpp_") && endswith(filename, "_global_cuts-list.txt")
				&& occursin(gamma_str, filename) && occursin(instance_name, filename))
			println("Processing cut file: $(filename)")
			filepath = joinpath(log_cuts_dir, filename)
			this_cut_list = read_twct_cuts_from_file(filepath)
			existing_cut_list = [existing_cut_list ; this_cut_list]
			println("Added $(size(this_cut_list, 1)) cuts.")
			flush(stdout)
		end
    end
	# Remove duplicate cuts from list
	existing_cut_list = collect(Set(existing_cut_list))
	count = size(existing_cut_list, 1)
	max_iter = count + 1
	println("Total number of cuts/iterations is $(size(existing_cut_list, 1)) cuts.")
	flush(stdout)
	# Read partial_results_full_filename as dataframe, to extract UB and best_pi
	if count > 0
		println("Reading result file: $(partial_results_full_filename)")
		if isfile(partial_results_full_filename)
		    println("Result file exists.")
		    df = DataFrame(CSV.File(partial_results_full_filename, delim=',', openquotechar='[', closequotechar=']'))
			nrows, ncols = size(df)
			if nrows > 0  # Pre-process df
				df = df[completecases(df), :]
				nrows, ncols = size(df)
				df[!, :instance_name] = strip.(df[!, :instance_name])
				# Filter results only for this instance_name and gamma
				df[ ( df.instance_name .== instance_name ), :]  # .& ( df.budget_T .== Γ ), :]
				println("Last 6 results: ")
				println(last(df, 6))
				for row in reverse(1:nrows)  # Process results, from newest to oldest
					gamma = df[row, :budget_T]
					if (typeof(df[row,:budget_T]) != Float64) & (typeof(df[row,:budget_T]) != Int64)
						df[row,:budget_T] = strip.(df[row,:budget_T])
						gamma = parse(Float64, df[row,:budget_T])
					end
					if gamma != Γ
						continue
					end
					println("Processing row: $(df[row, :])")
					permutation_str = strip(df[row, :permutation])
					best_pi = string_to_permutation(permutation_str)
					UB = df[row, :wct_validation]
					if (typeof(df[row,:wct_validation]) != Float64) & (typeof(df[row,:wct_validation]) != Int64)
						df[row,:wct_validation] = strip.(df[row,:wct_validation])
						UB = parse(Float64, df[row,:wct_validation])
					end
					break
				end
				println("Sucessfully read $(count) cuts with max_iter = $(max_iter) and UB = $(UB).")
			else
				println("No existing cuts found, returning empty cut list.")
			end
		else
			println("No existing cuts found, returning empty cut list.")
		end
	else
		println("No existing cuts found, returning empty cut list.")
	end
	flush(stdout)
    return existing_cut_list, max_iter, UB, best_pi
end
