# ================================================================================
#  test_deterministic_pfsp_wct.jl
# ================================================================================
# Tests the deterministic models for the permutation flow shop problem
# with weighted sum of completion times (wct) objective.
# ================================================================================

include("../../wct/deterministic_pfsp_wagner_wct_mip.jl")
include("../../wct/deterministic_pfsp_wilson_wct_mip.jl")
include("../../wct/deterministic_pfsp_alt_wct_mip.jl")
include("../../wct/deterministic_pfsp_manne_wct_mip.jl")

include("../../config.jl")
include("../../cli_arguments.jl")

EXPERIMENT_NAME = "run_rpfs_wct_milp"
# If possible, do not modify the lines below
EXPERIMENT_OUTPUT_FOLDER = abspath(create_full_dir(normpath(output_folder),
	[EXPERIMENT_NAME]))

# The criterion under consideration is sum_i (w_i * C_i), which is the
# the sum of the weighted job completion times.
function my_calculate_wct(m, n, p, w, pi)
    C = zeros(Float64, (m, n))
    k = pi[1]
    C[1, 1] = p[1, k]
    for i in 2:m
        C[i, 1] = C[i - 1, 1] + p[i, k]
    end
    for j in 2:n
        k = pi[j]
        C[1, j] = C[1, j - 1] + p[1, k]
        for i in 2:m
            C[i, j] = max(C[i, j - 1], C[i - 1, j]) + p[i, k]
        end
    end
    wct = 0
    for j in 1:n  # C[m, j] is the final completion time of job in position j
        k = pi[j]
        wct += C[m, j] * w[k]
    end
    return wct
end

function test_deterministic_model(model = "wagner"; m = 2, n = 10, instance_name = "", time_limit = 1800, solver="CPLEX")
    basedir = joinpath(robust_ying_instances_folder, "rob-pfsp-wct", "$(n)x$(m)")
    count = 0
    # Open output CSV file for results
    my_output_folder = create_full_dir(EXPERIMENT_OUTPUT_FOLDER, [model])
    result_file_path = joinpath(my_output_folder, model * "_$(n)x$(m).csv")
    result_file = open(result_file_path, "a+")
    println(result_file, "executionId,model_name,instance_name,n,m,wct,permutation,time_spent,is_optimal,validated,best_bound,wct_validation")
    if length(instance_name) > 0
        println("Filtering by instance_name = '$(instance_name)'.")
    end
    file_list = []
    for filename in searchdir(basedir, ".txt")
        if length(instance_name) > 0
            if filename != instance_name
                continue
            end
        else  # Let's solve only the instances with alpha = 10%, since we only use T = P_bar matrix
            _alpha = filename[begin:findlast("_wct_inputs", filename)[1]-1]
            _alpha = _alpha[findlast("_", _alpha)[1]+1:end]
            if _alpha != "10"
                continue
            end
        end
        inputfile = "$(basedir)/$(filename)"
        push!(file_list, (inputfile, filename))
    end
    println("List of instances to process: $(file_list)")
    for (inputfile, instance_name) in file_list
        println("Processing file $(inputfile)...")
        m, n, P_bar, P_hat, w = read_robust_input_file(inputfile)
        # wct_star is the optimal solution value obtained over the nominal processing time array (P_bar)
        if model == "wagner"
            wct_star, pi, time, is_opt, Z, best_bound = solve_deterministic_wagner_model_wct(m, n, P_bar, w, time_limit; solver=solver)
        elseif model == "wilson"
            wct_star, pi, time, is_opt, Z, best_bound = solve_deterministic_wilson_model_wct(m, n, P_bar, w, time_limit; solver=solver)
		elseif model == "manne"
			wct_star, pi, time, is_opt, Z, best_bound = solve_deterministic_manne_wct_model(m, n, P_bar, w, time_limit; solver=solver)
        else  # "alt"
            wct_star, pi, time, is_opt, Z, best_bound = solve_deterministic_alt_model_wct(m, n, P_bar, w, time_limit)
        end
        println("Based on the nominal processing times, the optimal WCT wct_star is $(wct_star).")
        println("Time spent: $(time)")
        wct_validation = my_calculate_wct(m, n, P_bar, w, pi)
        println("WCT Validation = $(wct_validation)")
        pi_str = join(["$(x)" for x in pi], " ")
        validated = true
        if abs(wct_validation - wct_star) > 0.0001
            println("ERROR : wct incorrect : wct_star = $(wct_star) x wct_validation = $(wct_validation)")
            validated = false
        end
        println(result_file, "none, $(model), $(instance_name), $(n), $(m), $(wct_star), $(pi_str), $(time), $(is_opt), $(validated), $(best_bound), $(wct_validation)")
        flush(result_file)
    end
    close(result_file)
    println("DONE.")
end

function run_experiment()
	parsed_args = parse_commandline()
	if !haskey(parsed_args, "instances") || (parsed_args["instances"] == nothing)
		println("ERROR: instances argument is mandatory!")
		return
	elseif !haskey(parsed_args, "model") || (parsed_args["model"] == nothing)
		println("ERROR: model argument is mandatory!")
		return
	end
	instance_group = parsed_args["instances"]
	model = parsed_args["model"]
	solver = parsed_args["solver"]
	time_limit = parsed_args["time-limit"]
	max_cores = parsed_args["max-cores"]
	n = tryparse(Int64, instance_group[begin:findlast("x", instance_group)[1]-1])
	m = tryparse(Int64, instance_group[findlast("x", instance_group)[1]+1:end])
	println("m = $(m), n = $(n)")
	test_deterministic_model(model; m=m, n=n, solver=solver, time_limit=time_limit)
end

run_experiment()
