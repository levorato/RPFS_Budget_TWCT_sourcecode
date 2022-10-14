# =============================================================================
# run_bb_rpfs_wct.jl
# =============================================================================
# Julia script to run the branch-and-bound procedure for the Robust PFS Problem
# (TWCT objective), based on the budget uncertainty set.
# =============================================================================
# To call this script, use, for example:
# julia run_bb_rpfs_wct.jl --instances 10x15 --model=ts3

using Random
using UUIDs
using CSV
using DataFrames
using Dates

include("../../src/julia/cli_arguments.jl")
include("../../src/julia/config.jl")
include("../../src/julia/pfsp_file_reader.jl")

EXPERIMENT_NAME = "run_bb_rpfs_wct_global"

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
			inputfile = joinpath("$(basedir)", "$(filename)")
			alpha = "RX"
			seq_idx = findfirst("_", filename)[1] + 1
			seq = filename[seq_idx:seq_idx+1]
			if seq_value == "-1"
				push!(instance_list, (filename, inputfile))
			elseif seq == seq_value
				push!(instance_list, (filename, inputfile))
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

function get_instance_label(instance_group)
	if instance_group == "10jobs"
		return "10x2"
	elseif instance_group == "tail"
		return "tail_20x5"
	elseif occursin("petro", instance_group)
		return instance_group
	else
		return "ying_" * instance_group
	end
end

function get_list_of_missing_results(instance_list, df, gamma_values)
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
				exit_code = df[row, :exit_code]
				gamma = df[row, :budget_gamma]
				if !(occursin(".txt", instance_name))  # Skip invalid row
					continue
				end
				if (typeof(df[row,:budget_gamma]) != Float64) & (typeof(df[row,:budget_gamma]) != Int64)
					df[row,:budget_gamma] = strip.(df[row,:budget_gamma])
					gamma = parse(Float64, df[row,:budget_gamma])
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
				if (instance == instance_name)
					println("*** $validated, $gap, $is_opt, $exit_code")
					if ((validated === Missing) || (gap === Missing) || (is_opt === Missing) || ismissing(validated) || ismissing(gap) || ismissing(is_opt)) || ismissing(exit_code)
						println("Skipping invalid line.")
						continue
					end
					if (exit_code != 0)
						println("Reprocessing instance $(instance) and budget $(gamma), since exit_code != 0.")
					elseif ((is_opt & (!validated)) == false) && (((abs(gap) > 1e-8) & is_opt) == false)
						println("Reusing existing result for instance $(instance) and budget $(gamma).")
						push!(existing_list, gamma)
					else
						println("Reprocessing instance $(instance) and budget $(gamma), since is_optimal == true and (gap > 1e-8 or validated == false).")
					end
					push!(processed_instance_list, (instance, gamma))
				end
		    end
		end
		println("existing_list = $(existing_list)")
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
	end
	flush(stdout)
	return run_list
end

function extract_solution_info(sol_str, success, hybrid)
	result_dict = Dict()
	if !success  # solver error
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
		result_dict["is_optimal"] = "false"
		result_dict["validated"] = "false"
	else
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
		result_dict["gap"] = parse(Float64, s[findfirst(':', s) + 2:findfirst('\n', s) - 1])

		r = findlast("is_optimal:", sol_str)
		s = sol_str[r[2]:length(sol_str)]
		result_dict["is_optimal"] = s[findfirst(':', s) + 2:findfirst('\n', s) - 1]

		r = findlast("validated:", sol_str)
		s = sol_str[r[2]:length(sol_str)]
		result_dict["validated"] = s[findfirst(':', s) + 2:findfirst('\n', s) - 1]
	end
	return result_dict
end

# first_improvement, vnd_size, vnd_permutation, random_vnd
function run_bb_pfsp_budget()
	parsed_args = parse_commandline()
	if !haskey(parsed_args, "instances") || isnothing(parsed_args["instances"])
		println("ERROR: instances argument is mandatory!")
		return
	end
	instance_group = parsed_args["instances"]
	budget_type = parsed_args["budget-type"]
	time_limit = parsed_args["time-limit"]
	alpha_value = parsed_args["alpha"]
	seq_value = parsed_args["seq"]
	gamma_value = parsed_args["gamma"]
	hybrid = parsed_args["hybrid"]
	# GRASP - best parameter set for timefactor = 300
	first_improvement = true
	vnd_size = 4
	random_vnd = true
	adaptive = false
	beta1 = 0.4
	beta2 = 1.0
	time_factor = 300
	println("Running Rob-PFSP-TWCT B&B Experiment for instance_group = $(instance_group)...")
	instance_list = get_instance_list(instance_group, alpha_value, seq_value)
	if gamma_value < 0
		gamma_values = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
	else
		gamma_values = [gamma_value]
	end
	println("Processing experiment for gamma_values = $(gamma_values).")

    # Open output CSV file for results
	prefix = ""
	if hybrid
		prefix = "hybrid_"
	end
	alpha_filename_suffix = ""
	seq_filename_suffix = ""
	if alpha_value != "-1"
		alpha_filename_suffix = "_alpha-$(alpha_value)"
	end
	if seq_value != "-1"
		seq_filename_suffix = "_seq-$(seq_value)"
	end
	full_filename = "bb_rpfs_wct_$(prefix)$(get_instance_label(instance_group))_fi_$(first_improvement)_vs_$(vnd_size)_rv_$(random_vnd)$(alpha_filename_suffix)$(seq_filename_suffix).csv"
	result_file_path = joinpath(EXPERIMENT_OUTPUT_FOLDER, full_filename)
	result_df = read_rpfs_results_from_csv(result_file_path)
	run_list = get_list_of_missing_results(instance_list, result_df, gamma_values)
	file_exists = isfile(result_file_path)
    result_file = open(result_file_path, "a+")
	if !file_exists
    	println(result_file, "execution_id,seed,ub_name,instance_name,alpha,n,m,budget_gamma,time_spent,exit_code"
			* ",solution_value,permutation,time_spent,iterations,num_phatomed,num_improvements"
			* ",first_improvement,vnd_size,vnd_permutation,random_vnd,adaptive,const_beta1,const_beta2,time_factor"
			* ",gap,is_optimal,validated")
	end
    output_folder2 = abspath(create_full_dir(EXPERIMENT_OUTPUT_FOLDER, ["bb-output"]))
	output_folder_stdout = abspath(create_full_dir(output_folder, ["stdout"]))
	prev_inputfile = Nothing
	m, n, P_bar, P_hat, w = 0, 0, zeros(Float64, (2, 2)), zeros(Float64, (2, 2)), zeros(Float64, (2))
	seed = rand(Int, 2)[1]
    rng = MersenneTwister(abs(seed))
	uid = uuid4(rng)
	batchId = Dates.format(Dates.now(), "yyyy_mm_dd-HH_MM_SS") * "-$(uid)"
	println("batchId = $(batchId)")
	for (filename, inputfile, p_Γ, Γ) in run_list
        println("=========================================================")
        println(" Solving instance $(filename) for Γ=$(Γ) (p_Γ=$(p_Γ)%)...")
        println("=========================================================")
		flush(stdout)
		if inputfile != prev_inputfile
        	m, n, P_bar, P_hat, w = read_robust_input_file(inputfile)
		end
        # Detect the value of alpha from the file name
		alpha_tuple = findlast("_wct_inputs", filename)
		if isnothing(alpha_tuple)
			alpha = "NA"
		else
			alpha_idx = alpha_tuple[1]
			alpha = filename[alpha_idx-2:alpha_idx-1]
		end
        # Create an output text file to save the stdout and sterr process output
		stdout_file_path = joinpath(output_folder_stdout, "a$(alpha)_$(filename)_gamma$(Γ)_stdout.txt")
	    stdout_file = open(stdout_file_path, "a+")
		# The experiment will be conducted such that \Gamma has the same values for all machines
        println(stdout_file, "*** Gamma = $(Γ)")
        Γ = "$(Γ)"
        scenario_processing_start_time = time_ns()
		exit_code = -500
		result_dict = Dict()
		if hybrid
			model_version = 607
		else
			model_version = 606
		end
		try
			execpath = abspath(joinpath(project_folder, "bin", "pfsp"))
            cmd = `$(execpath) --input-file="$(inputfile)"
                --output-folder $(output_folder2)
                --model-version=$(model_version)
				--time-limit=$(time_limit)
                --budget-gamma="$(Γ)"
                --sel-strategy=dfs
                --rob-ub-type=grasp
				--first-improvement=$(first_improvement)
				--vnd-size=$(vnd_size)
				--random-vnd=$(random_vnd)
				--adaptive-construction=$(adaptive)
				--grasp-tf=$(time_factor)
				--const-beta1=$(beta1)
				--const-beta2=$(beta2)
				--budget-type=global
				--obj-update-freq=3
				--seed=42`
			println("cmd: $(cmd)")
			flush(stdout)
			err = Pipe()
		    out = Pipe()
		    proc = run(pipeline(ignorestatus(cmd),stdout=out,stderr=err))
		    close(err.in)
		    close(out.in)
			stdout_str = String(read(out))
		    println(stdout_file, "*** stdout: ",stdout_str, "\n")
		    println(stdout_file, "*** stderr: ",String(read(err)), "\n")
		    println(stdout_file, "*** exit status: ", proc.exitcode, "\n")
			flush(stdout_file)
			exit_code = proc.exitcode
			if exit_code == 0
				result_dict = extract_solution_info(stdout_str, true, hybrid)
			else
				result_dict = extract_solution_info(stdout_str, false, hybrid)
			end
        catch e
			elapsed_time = (time_ns() - scenario_processing_start_time) * 1e-9
			bt = catch_backtrace()
			msg = sprint(showerror, e, bt)
			println(msg)
			println(stdout_file, "*** EXCEPTION OCCURED: $(msg)")
			flush(stdout_file)
			result_dict = extract_solution_info("", false, hybrid)
        end
		elapsed_time = (time_ns() - scenario_processing_start_time) * 1e-9
		print(result_file, "$(result_dict["execution_id"]),$(result_dict["seed"]),bb-grasp,$(filename),$(alpha),$(n),$(m),$(Γ),$(elapsed_time),$(exit_code)")
		print(result_file, ",$(result_dict["solution_value"]),$(result_dict["permutation"]),$(result_dict["time_spent"])")
		print(result_file, ",$(result_dict["iterations"]),$(result_dict["num_phatomed"]),$(result_dict["num_improvements"])")
		vnd_perm = [x for x in 1:vnd_size]
		vnd_perm_str = "$(vnd_perm)"
		vnd_perm_str = replace(vnd_perm_str, "," => "")
		print(result_file, ",$(first_improvement),$(vnd_size),$(vnd_perm_str),$(random_vnd),$(adaptive),$(beta1),$(beta2),$(time_factor)")
		println(result_file, ",$(result_dict["gap"]),$(result_dict["is_optimal"]),$(result_dict["validated"])")
        flush(result_file)
		close(stdout_file)
		prev_inputfile = inputfile
    end
    close(result_file)
    println("DONE.")
	flush(stdout)
end

run_bb_pfsp_budget()
