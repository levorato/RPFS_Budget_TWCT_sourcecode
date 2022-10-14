# =============================================================================
# run_grasp_rpfs_wct.jl
# =============================================================================
# Julia script to run the GRASP metaheuristic for the Robust PFS Problem
# (WCT objective), based on the budget uncertainty set.
# =============================================================================

using Random
using UUIDs
using Dates
include("../../../src/julia/config.jl")
include("../../../src/julia/pfsp_file_reader.jl")

EXPERIMENT_NAME = "run_grasp_rpfs_wct"

# If possible, do not modify the lines below
EXPERIMENT_OUTPUT_FOLDER = abspath(create_full_dir(normpath(output_folder),
	[EXPERIMENT_NAME]))

function extract_solution_info(sol_str, success)
	result_dict = Dict()
	if !success  # solver error
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
	else
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
	end
	return result_dict
end

# first_improvement, vnd_size, vnd_permutation, random_vnd
function run_grasp_pfsp_wct_budget(file_list, num_runs, first_improvement, vnd_size, random_vnd, adaptive, beta1, beta2, time_factor, obj_update_freq)
    # Open output CSV file for results
    result_file_path = joinpath(EXPERIMENT_OUTPUT_FOLDER, "ying10_grasp_wct_fi_$(first_improvement)_vs_$(vnd_size)_rv_$(random_vnd).csv")
    result_file = open(result_file_path, "a+")
    println(result_file, "batch_id, run_id, execution_id, seed, ub_name, instance_name, alpha, n, m, budget_gamma, time_spent, exit_code"
		* ", solution_value, permutation, time_spent, time_to_best_sol, iterations, num_visited_solutions, num_improvements"
		* ", first_improvement, vnd_size, vnd_permutation, random_vnd, adaptive, const_beta1, const_beta2, time_factor")
    rng = MersenneTwister(1234)
    for run_id in 1:num_runs
        uid = uuid4(rng)
        println("$(uid)")
        batchId = Dates.format(Dates.now(), "yyyy_mm_dd-HH_MM_SS") * "-$(uid)"
        output_folder2 = abspath(create_full_dir(EXPERIMENT_OUTPUT_FOLDER, [batchId]))
		output_folder_stdout = abspath(create_full_dir(output_folder2, ["stdout"]))
        for filepath in file_list
			filename = filepath[findlast('/', filepath)+1:length(filepath)]
            println("============================================================================================")
            println("  [GRASP-RobPFSP-WCT] Solving instance $(filename)...")
            println("============================================================================================")
			fullfilepath = joinpath(instances_folder, filepath)
            m, n, P_bar, P_hat, w = read_robust_input_file(fullfilepath)
            # Detect the value of alpha from the parent folder name
            alpha_str = filepath[1:findlast("_wct_inputs", filepath)[1] - 1]
            alpha_str = alpha_str[findlast('_', alpha_str) + 1:length(alpha_str)]
            alpha_str = strip(alpha_str)
            alpha = parse(Int64, alpha_str)
			# Create an output text file to save the stdout and sterr process output
			stdout_file_path = joinpath(output_folder_stdout, "r$(run_id)_a$(alpha)_$(filename)_stdout.txt")
		    stdout_file = open(stdout_file_path, "a+")
			# The experiment will be conducted such that \Gamma has the same values for all machines
            for Γ1 in [20, 40, 60, 80, 100]
                println(stdout_file, "*** Gamma = $(Γ1)")
                Γ = ""
                for x in 1:m
					Γ *= "$(Γ1)"
                    if x < m
                        Γ *= " "
                    end
                end
                scenario_processing_start_time = time_ns()
				exit_code = -500
				result_dict = Dict()
				try
					execpath = abspath(joinpath(project_folder, "bin", "pfsp"))
                    cmd = `$(execpath) --input-file="$(fullfilepath)"
                        --output-folder $(output_folder2)
                        --model-version=604
                        --budget-T="$(Γ)"
                        --rob-ub-type=grasp
						--first-improvement=$(first_improvement)
						--vnd-size=$(vnd_size)
						--random-vnd=$(random_vnd)
						--adaptive-construction=$(adaptive)
						--grasp-tf=$(time_factor)
						--const-beta1=$(beta1)
						--const-beta2=$(beta2)
						--obj-update-freq=$(obj_update_freq)`
					println("cmd: $(cmd)")
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
						result_dict = extract_solution_info(stdout_str, true)
					else
						result_dict = extract_solution_info(stdout_str, false)
					end
                catch e
					elapsed_time = (time_ns() - scenario_processing_start_time) * 1e-9
					bt = catch_backtrace()
					msg = sprint(showerror, e, bt)
					println(msg)
					println(stdout_file, "*** EXCEPTION OCCURED: $(msg)")
					flush(stdout_file)
					result_dict = extract_solution_info("", false)
                end
				elapsed_time = (time_ns() - scenario_processing_start_time) * 1e-9
				print(result_file, "$(batchId), $(run_id), $(result_dict["execution_id"]), $(result_dict["seed"]), grasp, $(filepath), $(alpha), $(n), $(m), $(Γ), $(elapsed_time), $(exit_code)")
				print(result_file, ", $(result_dict["solution_value"]), $(result_dict["permutation"]), $(result_dict["time_spent"]), $(result_dict["time_to_best_sol"])")
				print(result_file, ", $(result_dict["iterations"]), $(result_dict["num_visited_solutions"]), $(result_dict["num_improvements"])")
				vnd_perm = [x for x in 1:vnd_size]
				vnd_perm_str = "$(vnd_perm)"
				vnd_perm_str = replace(vnd_perm_str, "," => "")
				println(result_file, ", $(first_improvement), $(vnd_size), $(vnd_perm_str), $(random_vnd), $(adaptive), $(beta1), $(beta2), $(time_factor)")
                flush(result_file)
            end
			close(stdout_file)
        end
    end
    close(result_file)
    println("DONE.")
end

function generate_instance_list(instance_set, size_list)
    instance_folder_list = []
	if instance_set == "ying"
	    for folder in ["robust/ying/rob-pfsp-wct/$(num_jobs)jobs" for num_jobs in size_list]
	        instance_folder_list = [ instance_folder_list ; [folder] ]
	    end
	end
	if instance_set == "tail"
    	instance_folder_list = [instance_folder_list ; ["robust/taillard/rob-pfsp-wct"]]
	end
    println("Instances folders to process: $(instance_folder_list)")
    # Randomly extract 20% of instances from each instance folder
    instances_list = []
    instance_count = 0
    for folder in instance_folder_list
        instances = ["$(joinpath(folder, file))" for file in searchdir(joinpath(instances_folder, folder), ".txt")]
        instance_count += length(instances)
        instances_list = [ instances_list; instances ]
    end
    println("All instances count: $(instance_count)")
    return instances_list
end

# Ying
file_list = generate_instance_list("ying", ["10"])
#file_list = generate_instance_list("ying", ["100"])
#file_list = generate_instance_list("ying", ["150"])
#file_list = generate_instance_list("ying", ["200"])
# Taillard
#file_list = generate_instance_list("tail", [])

# best parameter set for timefactor = 30
num_runs = 100
first_improvement = true
vnd_size = 4
random_vnd = true
adaptive = false
beta1 = 1.0
beta2 = 1.0
time_factor = 30
###run_grasp_pfsp_wct_budget(file_list, num_runs, first_improvement, vnd_size, random_vnd, adaptive, beta1, beta2, time_factor, 4)

# best parameter set for timefactor = 300
num_runs = 1
first_improvement = true
vnd_size = 4
random_vnd = true
adaptive = false
beta1 = 1.0
beta2 = 1.0
time_factor = 3000
obj_update_freq = 4
run_grasp_pfsp_wct_budget(file_list, num_runs, first_improvement, vnd_size, random_vnd, adaptive, beta1, beta2, time_factor, obj_update_freq)
