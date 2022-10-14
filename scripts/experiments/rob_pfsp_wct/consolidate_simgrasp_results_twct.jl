# =============================================================================
# consolidate_simgrasp_results_twct.jl
# =============================================================================
# Julia script to consolidate SimGRASP output files for the Stochastic PFSP
# (TWCT objective), based on the GRASP simheuristic.
# =============================================================================

include("../../../src/julia/cli_arguments.jl")
include("../../../src/julia/config.jl")
include("../../../src/julia/pfsp_file_reader.jl")
include("../../../src/julia/wct/robust_pfsp_alt_worstcase_wct_mip.jl")

using CSV
using DataFrames
using Dates

EXPERIMENT_NAME = "run_simgrasp_twct"
SIMGRASP_OUTPUTDIR = abspath(joinpath(home_prefix, "pfsp_experiments", "run_simgrasp_twct"))
SCRIPT_OUTPUT_FOLDER = abspath(create_full_dir(joinpath(home_prefix, "pfsp_experiments"),
	[EXPERIMENT_NAME]))
SOLVER_TIME_LIMIT = 14400

function get_instance_list(simgrasp_resultdir, instances = "")
	println("Processing SimGRASP results...")
	outputfile_list = ["$(file)" for file in searchdir_recursive(simgrasp_resultdir, "outputs.txt")]
	instance_list = []
	for result_filepath in outputfile_list
		result_filename = basename(result_filepath)
		println("Result filename: $(result_filename)")
		partial_instance_name, n, m, nehtdist, beta1, beta2, seed = interpret_filename(result_filename)
		println("Instance: $(partial_instance_name), n=$(n), m=$(m), nehtdist=$(nehtdist), beta1=$(beta1), beta2=$(beta2), seed=$(seed)")
		alpha = Int64(ceil(beta2 * 100))
		# Process only simgrasp results with alpha == 100%
		if alpha != 100
			println("Skipping SimGRASP result, since alpha = $(alpha) != 100.")
			continue
		end
		seq = result_filename[findfirst("_", result_filename)[1]+1:end]
		seq = seq[begin:findfirst("_", seq)[1]-1]
		instance_filename = result_filename[begin:findfirst("_t_", result_filename)[1]-1] * "_wct_inputs.txt"
		if (length(instances) > 0) && (!occursin(instances, instance_filename))
			println("Skipping SimGRASP result, since instances = $(instances) NOT IN $(instance_filename).")
			continue
		end
		if instance_filename == "instance_2_9_4_wct_inputs.txt"
			instance_name_list = ["instance_2_9_4_wct_inputs.txt", "instance_9_4_R30_wct_inputs.txt", "instance_9_4_R50_wct_inputs.txt", "instance_9_4_R100_wct_inputs.txt"]
		else  # "instance_1_15_5_wct_inputs.txt"
			instance_name_list = ["instance_1_15_5_wct_inputs.txt", "instance_15_5_R30_wct_inputs.txt", "instance_15_5_R50_wct_inputs.txt", "instance_15_5_R100_wct_inputs.txt"]
		end
		for instance_name in instance_name_list
			push!(instance_list, (result_filepath, result_filename, instance_name, n, m, nehtdist, beta1, beta2, seed, alpha, seq))
		end
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
	# Pre-process df
	if nrows > 0
		df = df[completecases(df), :]
		nrows, ncols = size(df)
		df[!, :instance_filename] = strip.(df[!, :instance_filename])
	end
	processed_instance_list = []
	for (result_filepath, result_filename, instance_name, n, m, nehtdist, beta1, beta2, seed, alpha, seq) in instance_list
		existing_list = []
		if nrows > 0
			for row in reverse(1:nrows)  # Process results, from newest to oldest
				this_result_filename = strip(df[row, :result_filename])
				this_instance_name = strip(df[row, :instance_filename])
				GammaPercRobustCost = df[row,:GammaPercRobustCost]
				if (typeof(df[row,:GammaPercRobustCost]) != Float64) & (typeof(df[row,:GammaPercRobustCost]) != Int64)
					df[row,:GammaPercRobustCost] = strip.(df[row,:GammaPercRobustCost])
					GammaPercRobustCost = parse(Float64, df[row,:GammaPercRobustCost])
				end
				if (this_result_filename == result_filename) && (this_instance_name == instance_name)
					GammaPercRobustCost = Float64(GammaPercRobustCost)
					push!(existing_list, GammaPercRobustCost)
					println("Reusing existing result for result_filename $(result_filename) and instance_name $(instance_name).")
					push!(processed_instance_list, instance_name)
				end
		    end
		end
		println("existing_list = $(existing_list)")
		for p_Γ in gamma_values
			Γ = (p_Γ * m * n) / 100.0
			println("*** p_Γ = $(p_Γ) => Γ = $(Γ)")
			if !(p_Γ in existing_list)
				println("Adding missing result for result_filename $(result_filename), instance_name $(instance_name) and GammaPercRobustCost $(p_Γ)% = $(Γ).")
				push!(run_list, (result_filepath, result_filename, instance_name, n, m, nehtdist, beta1, beta2, seed, alpha, seq, p_Γ, Γ))
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

# interpret the filename of simgrasp result
# e.g. instance_1_15_5_t_1.0_1.0_441454_outputs.txt
function interpret_filename(filename)
	items = split(filename, "_")
	instance_name = items[1] * "_" * items[2]
	n = parse(Int64, items[3])
	m = parse(Int64, items[4])
	nehtdist = items[5]
	beta1 = parse(Float64, items[6])
	beta2 = parse(Float64, items[7])
	seed = parse(Int64, items[8])
	return instance_name, n, m, nehtdist, beta1, beta2, seed
end

# Order NEH solution, Out best solution, Our best stoch-sol
# Sol ID : 644
# Sol costs: 276
# Sol expCosts: 276.1002685393885
# Sol time: 0h 0m 0s (0.001984315 sec.)
# List of jobs:
# 1
# (...)
function read_simgrasp_solution_file(filepath)
	neh_sol = Dict()  # id, cost, expcost, time, permutation
	best_sol = Dict()  #
	stoch_sol = Dict()  #
	content = open(filepath) do file  # read whole file to str
	    read(file, String)
	end
	i1 = findfirst("NEH", content)[1]
	i2 = findfirst("Our best solution", content)[1]
	i3 = findfirst("Our best stoch-sol", content)[1]
    neh_sol = read_individual_solution(content[i1:i2-1])
    best_sol = read_individual_solution(content[i2:i3-1])
	stoch_sol = read_individual_solution(content[i3:length(content)])
	println(stoch_sol)
	return neh_sol, best_sol, stoch_sol
end

function read_individual_solution(lines)
	solution = Dict()
	lines = split(lines, "\n")
	count = 1
	while count <= size(lines, 1)
		ln = strip(lines[count])
        if length(ln) > 0 && (ln[1] != '*') && (ln[1] != '-')
			if occursin("Sol ID", ln)
				solution["sol_id"] = parse(Int64, strip(ln[findfirst(':', ln)+1:length(ln)]))
			elseif occursin("Sol costs", ln)
				solution["cost"] = parse(Float64, strip(ln[findfirst(':', ln)+1:length(ln)]))
			elseif occursin("Sol expCosts", ln)
				solution["exp_cost"] = parse(Float64, strip(ln[findfirst(':', ln)+1:length(ln)]))
			elseif occursin("Sol time", ln)
				time = ln[findfirst(':', ln)+1:length(ln)]
				time = time[findfirst('(', time)+1:findfirst("sec", time)[1]-1]
				solution["time"] = parse(Float64, strip(time))
			elseif occursin("List of jobs", ln)
				permutation = []
				count += 1
				while count <= size(lines, 1)
					ln = strip(lines[count])
					if length(ln) > 0 && ln[1] != '-'
						push!(permutation, parse(Int64, strip(ln)))
					end
					count += 1
				end
				solution["permutation"] = permutation
			end
		end
		count += 1
	end
	return solution
end

function permutation_to_string(permutation)
	s = "$(permutation)"
	s = replace(s, "Any[" => "")
	s = replace(s, "]" => "")
	s = replace(s, "," => "")
	return s
end

function calculate_worstcase_cost(m, n, P_bar, P_hat, w, permutation, GammaPercRobustCost)
	p_Γ = GammaPercRobustCost
	Γ = (p_Γ * m * n) / 100.0
	worstcase_cost, list_worst_Γ_scenario_z_sep, time_spent = solve_robust_model_worst_case_alt(permutation, m, n, Γ, P_bar, P_hat, w;
		budget_type="global", solver_time_limit=SOLVER_TIME_LIMIT, single_cut_per_iteration=true,
		fractional_solution=true, solver="Gurobi", max_cores=16, bigM_type=7)
	return worstcase_cost, time_spent
end

function searchdir_recursive(directory, filename)
	all_files_list = []
	for (dir_name, subdir_list, file_list) in walkdir(directory)
		file_list = [joinpath(dir_name, x) for x in file_list if occursin(filename, x)]
		all_files_list = [all_files_list ; file_list]
	end
	return all_files_list
end

function process_simgrasp_outputdir(base_instance_filepath)
	setup_gurobi_license()
	parsed_args = parse_commandline()
	gamma_value = parsed_args["gamma"]
	instances = parsed_args["instances"]
    # Open output CSV file for results
	result_file_path_stochgrasp = joinpath(SCRIPT_OUTPUT_FOLDER, EXPERIMENT_NAME)
	if !haskey(parsed_args, "instances") || isnothing(parsed_args["instances"])
		result_file_path_stochgrasp = result_file_path_stochgrasp * "_stochgrasp_results.csv"
		instances = ""
	else
		result_file_path_stochgrasp = result_file_path_stochgrasp * "-" * instances * "_stochgrasp_results.csv"
	end
	result_file_exists = isfile(result_file_path_stochgrasp)
    result_file_stochgrasp = open(result_file_path_stochgrasp, "a+")
	println("Saving consolidated result file to " * result_file_path_stochgrasp)
	if !result_file_exists
		println(result_file_stochgrasp, "result_filename,instance_filename,n,m,nehtdist,beta1,beta2,seed,alpha,seq,"
			* "stochsol_sol_id,stochsol_cost,stochsol_exp_cost,stochsol_time,stochsol_permutation,GammaPercRobustCost,RobustCost,time_spent")
		flush(result_file_stochgrasp)
	end

	instance_list = get_instance_list(SIMGRASP_OUTPUTDIR, instances)
	if gamma_value < 0
		gamma_values = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 90, 100]
	else
		gamma_values = [gamma_value]
	end
	println("Processing experiment for gamma_values = $(gamma_values).")
	run_list = get_list_of_missing_results(instance_list, result_file_path_stochgrasp, gamma_values)

	for (result_filepath, result_filename, instance_filename, n, m, nehtdist, beta1, beta2, seed, alpha, seq, GammaPercRobustCost, Γ) in run_list
		println("=======================================================================================================================================")
		println(" Processing Instance: $(instance_filename), seq=$(seq), n=$(n), m=$(m), nehtdist=$(nehtdist), beta1=$(beta1), beta2=$(beta2), seed=$(seed)")
		println("    For GammaPercRobustCost = $(GammaPercRobustCost) %...")
		println("=======================================================================================================================================")
		instance_filepath = abspath(joinpath(base_instance_filepath, instance_filename))
		m, n, P_bar, P_hat, w = read_robust_input_file(instance_filepath)
		# println("$m, $n, $P_bar, $P_hat, $w")
		neh_sol, best_sol, stoch_sol = read_simgrasp_solution_file(result_filepath)

		print(result_file_stochgrasp, "$(result_filename),$(instance_filename),$(n),$(m),$(nehtdist),$(beta1),$(beta2),$(seed),"
			* "$(alpha),$(seq)")
		print(result_file_stochgrasp, ",$(stoch_sol["sol_id"])")
		print(result_file_stochgrasp, ",$(stoch_sol["cost"])")
		print(result_file_stochgrasp, ",$(stoch_sol["exp_cost"])")
		print(result_file_stochgrasp, ",$(stoch_sol["time"])")
		print(result_file_stochgrasp, ",$(permutation_to_string(stoch_sol["permutation"]))")
		worstcase_cost, time_spent = calculate_worstcase_cost(m, n, P_bar, P_hat, w, stoch_sol["permutation"], GammaPercRobustCost)
		println(result_file_stochgrasp, ",$(GammaPercRobustCost),$(worstcase_cost),$(time_spent)")

		flush(result_file_stochgrasp)
		flush(stdout)
	end
    close(result_file_stochgrasp)
	println("DONE.")
end

base_instance_filepath = joinpath(instances_folder, "petro", "v3")
process_simgrasp_outputdir(base_instance_filepath)
