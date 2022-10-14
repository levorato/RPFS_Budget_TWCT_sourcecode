# ===================================================================================================
# EXPERIMENT exp_robust_pfsp_wct_budget_separation_main.jl
# Solve the robust problem via separation approach, with a specified robust counterpart.
# Robust Permutation Flow Shop Problem to Minimize the Weighted Sum of Job Completion Times (WCT)
# ===================================================================================================

include("../../../config.jl")
include("../../../wct/robust_pfsp_alt_worstcase_wct_mip.jl")
include("../../../wct/robust_pfsp_budget_wct_separation_wilson.jl")
include("../../../wct/robust_pfsp_budget_wct_separation_tba.jl")
include("../../../wct/robust_pfsp_budget_wct_separation_ts2.jl")
include("../../../wct/robust_pfsp_budget_wct_separation_ts3.jl")
include("../../../wct/robust_pfsp_budget_wct_separation_wagner_wst2.jl")
include("../../../wct/robust_pfsp_budget_wct_separation_manne.jl")
include("../../../wct/robust_pfsp_budget_wct_separation_liao_you.jl")
include("../../../wct/robust_pfsp_budget_wct_separation_hybrid.jl")
include("../../../cli_arguments.jl")

using CSV
using DataFrames
using Dates

# To call this script, use, for example:
# julia exp_robust_pfsp_wct_budget_separation_main.jl --instances 10jobs --model=ts3

EXPERIMENT_NAME = "run_ccg_rpfs_wct_global"
# If possible, do not modify the lines below
EXPERIMENT_OUTPUT_FOLDER = abspath(create_full_dir(normpath(output_folder),
	[EXPERIMENT_NAME]))

function get_instance_list(instance_group, alpha_value, seq_value)
	basedir = ""
	instance_list = []
	if alpha_value != "-1"
		println("Filter instances based on alpha value = $(alpha_value)...")
	end
	if instance_group == "tail"  # Taillard instances 20x5
		println("Processing Taillard instances...")
		t_instance_list = ["tail00$(x)_" for x in 1:9]
		push!(t_instance_list, "tail010_")
		println("Instance list is $(t_instance_list)")
	    basedir = joinpath(robust_tail_instances_folder, "rob-pfsp-wct")
		for filename in searchdir(basedir, ".txt")
			prefix = filename[1:findnext("_", filename, 1)]
			alpha_idx = findlast("_wct_inputs", filename)[1]
			alpha_str = filename[begin:alpha_idx-1]
			alpha = alpha_str[findlast("_", alpha_str)[1]+1:end]
			seq = prefix[5:end]
			if prefix in t_instance_list
			    inputfile = joinpath("$(basedir)", "$(filename)")
				if (alpha_value == "-1") && (seq_value == "-1")
	                push!(instance_list, (filename, inputfile))
		        elseif (seq_value == "-1") && (alpha == alpha_value)
			        push!(instance_list, (filename, inputfile))
				elseif (alpha_value == "-1") && (seq == seq_value)
					push!(instance_list, (filename, inputfile))
				elseif (alpha == alpha_value) && (seq == seq_value)
					push!(instance_list, (filename, inputfile))
				end
			end
		end
	elseif occursin("petro", instance_group)
		println("Processing Petro instances...")
		if instance_group == "petro"
			basedir = joinpath(instances_folder, "$(instance_group)", "v3")
		else
			basedir = joinpath(instances_folder, "petro", "v3", "random_alpha")
		end
		println("Instance folder = $(basedir)")
		for filename in searchdir(basedir, ".txt")
			println("$(filename)")
			if occursin("wct_inputs", filename)
				inputfile = joinpath("$(basedir)", "$(filename)")
				alpha_idx = findlast("_wct_inputs", filename)[1]
				alpha_str = filename[begin:alpha_idx-1]
				alpha = alpha_str[findlast("_", alpha_str)[1]+1:end]
				seq_idx = findfirst("_", filename)[1] + 1
				seq_str = filename[seq_idx:end]
				seq = seq_str[begin:findfirst("_", seq_str)[1] - 1]
				if (alpha_value == "-1") && (seq_value == "-1")
					push!(instance_list, (filename, inputfile))
				elseif (seq_value == "-1") && (alpha == alpha_value)
					push!(instance_list, (filename, inputfile))
				elseif (alpha_value == "-1") && (seq == seq_value)
					push!(instance_list, (filename, inputfile))
				elseif (alpha == alpha_value) && (seq == seq_value)
					push!(instance_list, (filename, inputfile))
				end
			end
		end
	else  # Ying-based instances
		println("Processing Ying-based instances...")
		basedir = joinpath(robust_ying_instances_folder, "rob-pfsp-wct", "$(instance_group)")
		for filename in searchdir(basedir, ".txt")
			inputfile = joinpath("$(basedir)", "$(filename)")
			alpha_idx = findlast("_wct_inputs", filename)[1]
			alpha_str = filename[begin:alpha_idx-1]
			alpha = alpha_str[findlast("_", alpha_str)[1]+1:end]
			seq_idx = findfirst("_", filename)[1] - 1
			seq = filename[seq_idx-1:seq_idx]
			if (alpha_value == "-1") && (seq_value == "-1")
				push!(instance_list, (filename, inputfile))
			elseif (seq_value == "-1") && (alpha == alpha_value)
				push!(instance_list, (filename, inputfile))
			elseif (alpha_value == "-1") && (seq == seq_value)
				push!(instance_list, (filename, inputfile))
			elseif (alpha == alpha_value) && (seq == seq_value)
				push!(instance_list, (filename, inputfile))
			end
		end
	end
	num_instances = size(instance_list, 1)
	println("# instances to process: $(num_instances)")
	# println("- Instances to process: $(instance_list)")
	return instance_list
	flush(stdout)
end

function read_rpfs_results_from_csv(result_filename)
	println("Reading result file: $(result_filename)")
	if isfile(result_filename)
	    println("Result file exists.")
	    df = DataFrame(CSV.File(result_filename, delim=','))
		println("First 6 results: ")
		println(first(df, 6))
		flush(stdout)
	    return df
	else
		println("No existing result file found!")
		flush(stdout)
		return DataFrame()
	end
end

function get_list_of_missing_results(instance_list, df, gamma_values, rel_gap_tol, is_gamma_perc)
	run_list = []
	nrows, ncols = size(df)
	#df = df[completecases(df), :]
	# Pre-process df
	if nrows > 0
		df = df[completecases(df), :]
		nrows, ncols = size(df)
		df[!, :instance_name] = strip.(df[!, :instance_name])
	end
	processed_instance_list = []
	for (instance, inputfile) in instance_list
		m, n, P_bar, P_hat, w = read_robust_input_file(inputfile)
		existing_list = []
		if nrows > 0
			for row in reverse(1:nrows)  # Process results, from newest to oldest
				instance_name = strip(df[row, :instance_name])
				is_opt = df[row, :is_optimal]
				validated = df[row, :validated]
				gap = df[row, :gap]
				mp_total_time = df[row, :mp_total_time]
				sp_total_time = df[row, :sp_total_time]
				gamma = df[row,:budget_T]
				if !(occursin(".txt", instance_name))  # Skip invalid row
					continue
				end
				if (typeof(df[row,:budget_T]) != Float64) & (typeof(df[row,:budget_T]) != Int64)
					df[row,:budget_T] = strip.(df[row,:budget_T])
					gamma = parse(Float64, df[row,:budget_T])
				end
				if (typeof(df[row,:gap]) != Float64) & (typeof(df[row,:gap]) != Int64)
					df[row,:gap] = strip.(df[row,:gap])
					gap = parse(Float64, df[row,:gap])
				end
				if (typeof(df[row,:is_optimal]) != Bool)
					df[row,:is_optimal] = strip.(df[row,:is_optimal])
					is_opt = parse(Bool, df[row,:is_optimal])
				end
				if (typeof(df[row,:validated]) != Bool)
					df[row,:validated] = strip.(df[row,:validated])
					validated = parse(Bool, df[row,:validated])
				end
				if (typeof(df[row,:mp_total_time]) != Float64)
					df[row,:mp_total_time] = strip.(df[row,:mp_total_time])
					mp_total_time = parse(Float64, df[row,:mp_total_time])
				end
				if (typeof(df[row,:sp_total_time]) != Float64)
					df[row,:sp_total_time] = strip.(df[row,:sp_total_time])
					sp_total_time = parse(Float64, df[row,:sp_total_time])
				end
				if (instance == instance_name)
					gamma = Float64(gamma)
					if (instance, gamma) in processed_instance_list  # if we've already processed results for this instance and gamma value
						continue
					end
					println("*** $validated, $gap, $is_opt, $mp_total_time, $sp_total_time")
					if ((validated === Missing) || (gap === Missing) || (is_opt === Missing) || ismissing(validated) || ismissing(gap) || ismissing(is_opt)) || ismissing(mp_total_time) || ismissing(sp_total_time)
						println("Skipping invalid line.")
						continue
					end
					if ((!validated) || ((abs(gap) > rel_gap_tol) && is_opt))
						println("Reprocessing instance $(instance) and budget $(gamma), since (is_optimal and rel_gap > $(rel_gap_tol)) or (validated == false).")
					elseif (!is_opt) && (mp_total_time + sp_total_time < 7000.0)
						println("Reprocessing instance $(instance) and budget $(gamma), since (is_optimal == false) and (mp_total_time + sp_total_time < 7000.0).")
					else
						push!(existing_list, gamma)
						println("Reusing existing result for instance $(instance) and budget $(gamma).")
					end
					push!(processed_instance_list, (instance, gamma))
				end
		    end
		end
		println("existing_list = $(existing_list)")
		if is_gamma_perc
			for p_Γ in gamma_values
				Γ = (p_Γ * m * n) / 100.0
				println("*** p_Γ = $(p_Γ) => Γ = $(Γ)")
				if !(Float64(Γ) in existing_list)
					println("Adding missing result for instance $(instance) and budget $(p_Γ)% = $(Γ).")
					push!(run_list, (instance, inputfile, p_Γ, Γ))
				else
					println("Reusing result.")
				end
			end
		else  # process all integer gamma values
			gamma_values = [x for x in (1:(m*n))]
			for Γ in gamma_values
				p_Γ = 100.0 * Γ / (m * n)
				println("*** p_Γ = $(p_Γ) => Γ = $(Γ)")
				if !(Float64(Γ) in existing_list)
					println("Adding missing result for instance $(instance) and budget Γ = $(Γ).")
					push!(run_list, (instance, inputfile, p_Γ, Γ))
				else
					println("Reusing result.")
				end
			end
		end
	end
	flush(stdout)
	return run_list
end

function get_instance_label(instance_group)
	if instance_group == "10jobs"
		return "10x2"
	elseif instance_group == "tail"
		return "tail_20x5"
	else
		return instance_group
	end
end

function solve_robust_model_flowshop_wct_budget_separation(m, n, P_bar, P_hat, w, Γ, time_limit;
			mp_bigM_type = 2, sp_bigM_type = 3, max_cores = 16, rel_gap_tol=1e-8,
            solver = "CPLEX", budget_type = "machine", single_cut_per_iteration = false, fractional_solution = false,
			use_brute_force = false, model = "wilson", warm_start = false, inputfile = nothing, sp_solver = nothing)
	if model == "wilson"
		return solve_robust_model_flowshop_wct_budget_separation_wilson(m, n, P_bar, P_hat, w, Γ, time_limit; max_cores=max_cores,
		            solver=solver, budget_type=budget_type, single_cut_per_iteration=single_cut_per_iteration, fractional_solution=fractional_solution,
					use_brute_force=use_brute_force, mp_bigM_type=mp_bigM_type, sp_bigM_type=sp_bigM_type, warm_start=warm_start,
					instance_filename=inputfile, sp_solver=sp_solver)
	elseif model == "wagner-wst2"
		return solve_robust_model_flowshop_wct_budget_separation_wagner_wst2(m, n, P_bar, P_hat, w, Γ, time_limit; max_cores=max_cores,
		            solver=solver, budget_type=budget_type, single_cut_per_iteration=single_cut_per_iteration, fractional_solution=fractional_solution,
					use_brute_force=use_brute_force, mp_bigM_type=mp_bigM_type, sp_bigM_type=sp_bigM_type, warm_start=warm_start,
					instance_filename=inputfile, sp_solver=sp_solver)
	elseif model == "ts2"
		return solve_robust_model_flowshop_wct_budget_separation_ts2(m, n, P_bar, P_hat, w, Γ, time_limit; max_cores=max_cores,
		            solver=solver, budget_type=budget_type, single_cut_per_iteration=single_cut_per_iteration, fractional_solution=fractional_solution,
					use_brute_force=use_brute_force, mp_bigM_type=mp_bigM_type, sp_bigM_type=sp_bigM_type, warm_start=warm_start,
					instance_filename=inputfile, sp_solver=sp_solver)
	elseif model == "ts3"
		return solve_robust_model_flowshop_wct_budget_separation_ts3(m, n, P_bar, P_hat, w, Γ, time_limit; max_cores=max_cores,
		            solver=solver, budget_type=budget_type, single_cut_per_iteration=single_cut_per_iteration, fractional_solution=fractional_solution,
					use_brute_force=use_brute_force, mp_bigM_type=mp_bigM_type, sp_bigM_type=sp_bigM_type, warm_start=warm_start,
					instance_filename=inputfile, sp_solver=sp_solver)
	elseif model == "tba"
		return solve_robust_model_flowshop_wct_budget_separation_tba(m, n, P_bar, P_hat, w, Γ, time_limit; max_cores=max_cores,
		            solver=solver, budget_type=budget_type, single_cut_per_iteration=single_cut_per_iteration, fractional_solution=fractional_solution,
					use_brute_force=use_brute_force, mp_bigM_type=mp_bigM_type, sp_bigM_type=sp_bigM_type, warm_start=warm_start,
					instance_filename=inputfile, sp_solver=sp_solver)
	elseif model == "manne"
		return solve_robust_model_flowshop_wct_budget_separation_manne(m, n, P_bar, P_hat, w, Γ, time_limit; max_cores=max_cores,
		            solver=solver, budget_type=budget_type, single_cut_per_iteration=single_cut_per_iteration, fractional_solution=fractional_solution,
					use_brute_force=use_brute_force, mp_bigM_type=mp_bigM_type, sp_bigM_type=sp_bigM_type, warm_start=warm_start,
					instance_filename=inputfile, sp_solver=sp_solver)
	elseif model == "liao-you"
		return solve_robust_model_flowshop_wct_budget_separation_liao_you(m, n, P_bar, P_hat, w, Γ, time_limit; max_cores=max_cores,
		            solver=solver, budget_type=budget_type, single_cut_per_iteration=single_cut_per_iteration, fractional_solution=fractional_solution,
					use_brute_force=use_brute_force, mp_bigM_type=mp_bigM_type, sp_bigM_type=sp_bigM_type, warm_start=warm_start,
					instance_filename=inputfile, sp_solver=sp_solver)
	elseif occursin("hybrid", model)
		return solve_robust_model_flowshop_wct_budget_separation_hybrid(m, n, P_bar, P_hat, w, Γ, time_limit; max_cores=max_cores,
		            solver=solver, budget_type=budget_type, single_cut_per_iteration=single_cut_per_iteration, fractional_solution=fractional_solution,
					use_brute_force=use_brute_force, mp_bigM_type=mp_bigM_type, sp_bigM_type=sp_bigM_type, warm_start=warm_start,
					instance_filename=inputfile, sp_solver=sp_solver, rel_gap_tol=rel_gap_tol, mp_model_name=model)
	else
		println("ERROR: CCG model not found: $(model)")
		flush(stdout)
		return Nothing
	end
end

function run_experiment()
	setup_gurobi_license()
	parsed_args = parse_commandline()
	if !haskey(parsed_args, "instances") || isnothing(parsed_args["instances"])
		println("ERROR: instances argument is mandatory!")
		return
	elseif !haskey(parsed_args, "model") || isnothing(parsed_args["model"])
		println("ERROR: model argument is mandatory!")
		return
	end
	instance_group = parsed_args["instances"]
	model = parsed_args["model"]
	solver = parsed_args["solver"]
	sp_solver = parsed_args["sp-solver"]
	budget_type = parsed_args["budget-type"]
	fractional_solution = parsed_args["fractional-solution"]
	single_cut_per_iteration = parsed_args["single-cut"]
	time_limit = parsed_args["time-limit"]
	max_cores = parsed_args["max-cores"]
	mp_bigm_type = parsed_args["mp-bigm-type"]
	sp_bigm_type = parsed_args["sp-bigm-type"]
	warm_start = parsed_args["warm-start"]
	alpha_value = parsed_args["alpha"]
	seq_value = parsed_args["seq"]
	gamma_value = parsed_args["gamma"]
	rel_gap_tol = parsed_args["gap"]
	println("Running Rob-PFSP-WCT C&CG Experiment for instance_group = $(instance_group), model = $(model) and rel_gap_tol = $(rel_gap_tol)...")
	instance_list = get_instance_list(instance_group, alpha_value, seq_value)
	if gamma_value < 0
		gamma_values = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
		is_gamma_perc = true
	elseif gamma_value == 200
		gamma_values = []
		is_gamma_perc = false
	elseif gamma_value > 100
		gamma_values = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 90, 100]
		is_gamma_perc = true
	else
		gamma_values = [gamma_value]
		is_gamma_perc = true
	end
	if isnothing(sp_solver)
		sp_solver = solver
	end
	println("Main MIP solver is $(solver), MIP SP solver is $(sp_solver).")
	println("Processing experiment for gamma_values = $(gamma_values), is_gamma_perc = $(is_gamma_perc).")

    # Open output CSV file for results
	alpha_filename_suffix = ""
	seq_filename_suffix = ""
	gamma_filename_suffix = ""
	if alpha_value != "-1"
		alpha_filename_suffix = "_alpha-$(alpha_value)"
	end
	if seq_value != "-1"
		seq_filename_suffix = "_seq-$(seq_value)"
	end
	if (gamma_value >= 0) && (gamma_value <= 100)
		gamma_filename_suffix = "_gamma-$(gamma_value)"
	end
	full_filename = "separation_wct_$(model)_randomweights_$(get_instance_label(instance_group))_warmstart-$(warm_start)_singlecut-$(single_cut_per_iteration)$(alpha_filename_suffix)$(seq_filename_suffix)$(gamma_filename_suffix).csv"
    result_file_path = joinpath(EXPERIMENT_OUTPUT_FOLDER, full_filename)
	result_file_exists = isfile(result_file_path)
	result_df = read_rpfs_results_from_csv(result_file_path)
	run_list = get_list_of_missing_results(instance_list, result_df, gamma_values, rel_gap_tol, is_gamma_perc)
    result_file = open(result_file_path, "a+")
	if !result_file_exists
    	println(result_file, "executionId,ub_name,instance_name,alpha,n,m,budget_T,wct,permutation,time_spent,time_to_best_sol,mp_total_time,sp_total_time,iterations,num_visited_solutions,num_improvements,is_optimal,validated,gap,lb,cost,wct_validation,cut_count")
	end
	prev_inputfile = Nothing
	m, n, P_bar, P_hat, w = 0, 0, zeros(Float64, (2, 2)), zeros(Float64, (2, 2)), zeros(Float64, (2))
    for (filename, inputfile, p_Γ, Γ) in run_list
		if inputfile != prev_inputfile
			m, n, P_bar, P_hat, w = read_robust_input_file(inputfile)
		end
		alpha_tuple = findlast("_wct_inputs", filename)
		if isnothing(alpha_tuple)
			alpha = "NA"
		else
			alpha_idx = alpha_tuple[1]
			alpha = filename[alpha_idx-2:alpha_idx-1]
		end
		datetimenow = Dates.now()
		println("Date: $(datetimenow)")
		println("Calculating robust PFSP wct budget solution (instance = $(filename), Γ = $(Γ))...")
		flush(stdout)
		scenario_processing_start_time = time_ns()
		cost, pi, x, run_time, is_opt, iter, gap, lb, cut_count, mp_total_time, sp_total_time = solve_robust_model_flowshop_wct_budget_separation(
			m, n, P_bar, P_hat, w, Γ, time_limit; solver=solver, budget_type=budget_type, fractional_solution=fractional_solution,
			model=model, max_cores=max_cores, mp_bigM_type=mp_bigm_type, sp_bigM_type=sp_bigm_type, use_brute_force = true,
			single_cut_per_iteration=single_cut_per_iteration, warm_start=warm_start, inputfile=inputfile, sp_solver=sp_solver, rel_gap_tol=rel_gap_tol)
		pi_str = join(["$(x)" for x in pi], " ")
		if (!is_opt) && (size(pi, 1) == 0)
			elapsed_time = (time_ns() - scenario_processing_start_time) * 1e-9  # ********************** stop time
			println(result_file, "none, mip_separation, $(filename), $(alpha), $(n), $(m), $(Γ), $(cost), $(pi_str), $(elapsed_time), $(elapsed_time), $(mp_total_time), $(sp_total_time), $(iter), $(iter), $(iter-1), $(is_opt), false, $(gap), $(lb), $(cost), 0, $(cut_count)")
			flush(result_file)
			continue
		end
		wct, Γ_scenario, e_time = solve_robust_model_worst_case_alt(pi, m, n, Γ, P_bar, P_hat, w, 1;
			solver=sp_solver, budget_type=budget_type, fractional_solution=fractional_solution, max_cores=max_cores,
			single_cut_per_iteration=single_cut_per_iteration, bigM_type=sp_bigm_type, solver_time_limit=time_limit)
		validated = true
		if abs(wct - cost) > 0.0001
			println("ERROR : wct incorrect : cost = $(cost) x wct = $(wct)")
			validated = false
		end
		elapsed_time = (time_ns() - scenario_processing_start_time) * 1e-9  # ********************** stop time
		println(result_file, "none, mip_separation, $(filename), $(alpha), $(n), $(m), $(Γ), $(wct), $(pi_str), $(elapsed_time), $(elapsed_time), $(mp_total_time), $(sp_total_time), $(iter), $(iter), $(iter-1), $(is_opt), $(validated), $(gap), $(lb), $(cost), $(wct), $(cut_count)")
		flush(result_file)
		prev_inputfile = inputfile
    end
    close(result_file)
    println("DONE.")
	flush(stdout)
end

run_experiment()
