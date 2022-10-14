# =============================================================================
# calculate_rpfs_budget_twct_worstcase_costs.jl
# =============================================================================
# Julia script to calculate worstcase costs for the RPFS Budgted Uncertainty
# (TWCT objective), based on an existing solution (permutation \sigma).
# =============================================================================

include("../../../src/julia/cli_arguments.jl")
include("../../../src/julia/config.jl")
include("../../../src/julia/pfsp_file_reader.jl")
include("../../../src/julia/wct/robust_pfsp_alt_worstcase_wct_mip.jl")

using CSV
using DataFrames
using Dates

EXPERIMENT_NAME = "calculate_rpfs_twct_worstcase"
EXPERIMENT_INSTANCE_FOLDER = abspath(joinpath(instances_folder))
CONSOLIDATED_RESULT_DIR = abspath(joinpath(home_prefix, "pfsp_experiments", "run_ccg_rpfs_wct_global"))
SCRIPT_OUTPUT_FOLDER = abspath(create_full_dir(joinpath(home_prefix, "pfsp_experiments"),
	[EXPERIMENT_NAME]))
SOLVER_TIME_LIMIT = 14400

function get_instance_list(csv_filepath, instances = "")
	basedir = ""
	instance_list = []
	println("Processing robust results for Petro instances...")
	# Read deterministic results (det0 and det100) into dataframe
	df = DataFrame(CSV.File(csv_filepath; delim=","))
	for row in eachrow(df)
		instance_name = strip(row["instance_name"])
		executionId = strip(row["executionId"])
		model = strip(row["ub_name"])
		n = row["n"]
		m = row["m"]
		alpha = row["alpha"]
		permutation = strip(row["permutation"])
		perm_pi = parse_permutation_from_string(permutation)
		Gamma = -1
		cost = 0
		Gamma = row["budget_T"]
		cost = row["wct"]
		if (length(instances) > 0) && (!occursin(instances, instance_name))
			println("Skipping robust instance result, since instances = $(instances) NOT IN $(instance_name).")
			continue
		end
		push!(instance_list, (instance_name, executionId, model, n, m, alpha, permutation, perm_pi, Gamma, cost))
	end
	num_instances = size(instance_list, 1)
	println("# instances to process: $(num_instances)")
	flush(stdout)
	return instance_list
end

function read_experiment_results_from_csv(result_filename)
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

function get_list_of_missing_results(instance_list, experiment_result_file_path, gamma_values)
	run_list = []
	df = read_experiment_results_from_csv(experiment_result_file_path)
	nrows, ncols = size(df)
	#df = df[completecases(df), :]
	# Pre-process df
	if nrows > 0
		df = df[completecases(df), :]
		nrows, ncols = size(df)
		df[!, :instance_name] = strip.(df[!, :instance_name])
	end
	processed_instance_list = []
	for (instance, executionId, model, n, m, alpha, permutation, perm_pi, GammaRob, cost) in instance_list
		existing_list = []
		if nrows > 0
			for row in reverse(1:nrows)  # Process results, from newest to oldest
				executionId = strip(df[row, :executionId])
				model = strip(df[row, :model])
				instance_name = strip(df[row, :instance_name])
				time_spent = df[row, :time_spent]
				Gamma = df[row,:Gamma]
				GammaPercRobustCost = df[row,:GammaPercRobustCost]
				if !(occursin(".txt", instance_name))  # Skip invalid row
					continue
				end
				if (typeof(df[row,:Gamma]) != Float64) & (typeof(df[row,:Gamma]) != Int64)
					df[row,:Gamma] = strip.(df[row,:Gamma])
					GammaPerc = parse(Float64, df[row,:Gamma])
				end
				if (typeof(df[row,:GammaPercRobustCost]) != Float64) & (typeof(df[row,:GammaPercRobustCost]) != Int64)
					df[row,:GammaPercRobustCost] = strip.(df[row,:GammaPercRobustCost])
					GammaPercRobustCost = parse(Float64, df[row,:GammaPercRobustCost])
				end
				if (typeof(df[row,:time_spent]) != Float64)
					df[row,:time_spent] = strip.(df[row,:time_spent])
					time_spent = parse(Float64, df[row,:time_spent])
				end
				if (instance == instance_name) && (GammaRob == Gamma)
					GammaPercRobustCost = Float64(GammaPercRobustCost)
					#if (instance, GammaRob) in processed_instance_list  # if we've already processed results for this instance and gamma value
					#	continue
					#end
					push!(existing_list, GammaPercRobustCost)
					println("Reusing existing result for instance $(instance) and budget $(GammaRob).")
					push!(processed_instance_list, (instance, GammaRob))
				end
		    end
		end
		println("existing_list = $(existing_list)")
		for p_Γ in gamma_values
			Γ = (p_Γ * m * n) / 100.0
			println("*** p_Γ = $(p_Γ) => Γ = $(Γ)")
			if !(p_Γ in existing_list)
				println("Adding missing result for instance ($(instance), $(GammaRob)) and GammaPercRobustCost $(p_Γ)% = $(Γ).")
				push!(run_list, (instance, executionId, model, n, m, alpha, permutation, perm_pi, GammaRob, cost, p_Γ, Γ))
			else
				println("Reusing result.")
			end
		end
	end
	num_experiments = size(run_list, 1)
	println("# experiments to run: $(num_experiments)")
	flush(stdout)
	return run_list
end

function permutation_to_string(permutation)
	s = "$(permutation)"
	s = replace(s, "Any[" => "")
	s = replace(s, "]" => "")
	s = replace(s, "," => "")
	return s
end

function parse_permutation_from_string(s)
	return [parse(Int64, ss) for ss in split(s)]
end

function calculate_worstcase_cost(m, n, P_bar, P_hat, w, permutation, GammaPercRobustCost)
	p_Γ = GammaPercRobustCost
	Γ = (p_Γ * m * n) / 100.0
	worstcase_cost, list_worst_Γ_scenario_z_sep, time_spent = solve_robust_model_worst_case_alt(permutation, m, n, Γ, P_bar, P_hat, w;
		budget_type="global", solver_time_limit=SOLVER_TIME_LIMIT, single_cut_per_iteration=true,
		fractional_solution=true, solver="Gurobi", max_cores=16, bigM_type=7)
	return worstcase_cost, time_spent
end

# Result Type can be: "deterministic" or "robust".
function process_consolidated_csv_result(csv_filename, exact = true)
	setup_gurobi_license()
	parsed_args = parse_commandline()
	gamma_value = parsed_args["gamma"]
	instances = parsed_args["instances"]
	# Open output CSV file for results
    result_file_path = joinpath(SCRIPT_OUTPUT_FOLDER, EXPERIMENT_NAME)
	if !haskey(parsed_args, "instances") || isnothing(parsed_args["instances"])
		result_file_path = result_file_path * ".csv"
		instances = ""
	else
		result_file_path = result_file_path * "-" * instances * ".csv"
	end
	result_file_exists = isfile(result_file_path)
    result_file = open(result_file_path, "a+")
	println("Saving worstcase costs result file to " * result_file_path)
	if !result_file_exists
		print(result_file, "executionId,model,instance_name,alpha,n,m,")
		println(result_file, "Gamma,wct_star,permutation,GammaPercRobustCost,RobustCost,time_spent")
		flush(result_file)
	end

	instance_list = get_instance_list(joinpath(normpath(CONSOLIDATED_RESULT_DIR), csv_filename), instances)
	if gamma_value < 0
		gamma_values = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 90, 100]
	else
		gamma_values = [gamma_value]
	end
	println("Processing experiment for gamma_values = $(gamma_values).")
	run_list = get_list_of_missing_results(instance_list, result_file_path, gamma_values)

	instance_dict = Dict()
	for (instance_name, executionId, model, n, m, alpha, permutation, perm_pi, Gamma, cost, GammaPercRobustCost, Γ) in run_list
		println("================================================================================================================")
		println("    Processing instance ($(instance_name), Gamma=$(Gamma)) for GammaPercRobustCost = $(GammaPercRobustCost) %...")
		println("================================================================================================================")
		fullfilepath = joinpath(instances_folder, "petro", "v3", instance_name)
		if !haskey(instance_dict, fullfilepath)
			m, n, P_bar, P_hat, w = read_robust_input_file(fullfilepath)
			instance_dict[fullfilepath] = m, n, P_bar, P_hat, w
		else
			m, n, P_bar, P_hat, w = instance_dict[fullfilepath]
		end
		println("Calculating PFSP worst-case budget scenarios...")
		print(result_file, "$(executionId),$(model),$(instance_name),$(alpha),$(n),$(m),")
		print(result_file, "$(Gamma),$(cost),")
		print(result_file, "$(permutation),")
		worstcase_cost, time_spent = calculate_worstcase_cost(m, n, P_bar, P_hat, w, perm_pi, GammaPercRobustCost)
		println(result_file, "$(GammaPercRobustCost),$(worstcase_cost),$(time_spent)")
		flush(result_file)
		flush(stdout)
	end
    close(result_file)
	println("DONE.")
end

#process_consolidated_csv_result("separation_wct_hybrid-manne_randomweights_petro_warmstart-false_singlecut-true.csv")
#process_consolidated_csv_result("separation_wct_hybrid-wilson_randomweights_petro-random_warmstart-false_singlecut-true_alpha-R100_seq-15.csv")
process_consolidated_csv_result("separation_wct_hybrid-wilson_randomweights_petro-random_warmstart-false_singlecut-true_alpha-R100_seq-9.csv")
