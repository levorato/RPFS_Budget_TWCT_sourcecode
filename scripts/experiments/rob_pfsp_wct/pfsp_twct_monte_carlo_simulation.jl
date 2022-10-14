# =============================================================================
# pfsp_twct_monte_carlo_simulation.jl
# =============================================================================
# Julia script to execute a Monte Carlo simulation for the Permutation Flowshop
# (TWCT objective), based on a given sequence of jobs (permutation).
# =============================================================================

# Apply simulation to all PFS Solutions (det0, det100, robpfs and stochpfs) using the same random seed.
using Random
using UUIDs
using Dates
using CSV
using DataFrames
using Distributions
using GZip
include("../../../src/julia/config.jl")
include("../../../src/julia/pfsp_file_reader.jl")
include("../../../src/julia/deterministic_pfsp.jl")

# If possible, do not modify the lines below
EXPERIMENT_NAME = "montecarlo_sim_rpfs_twct"
EXPERIMENT_OUTPUT_FOLDER = abspath(create_full_dir(normpath(output_folder),
	[EXPERIMENT_NAME]))
RANDOM_SEED = 809132

# https://stats.stackexchange.com/questions/12232/calculating-the-parameters-of-a-beta-distribution-using-the-mean-and-variance
# Code to estimate the parameters of the Beta distribution from a given mean, mu, and variance, var:
function estBetaParams(mu, var)
	alpha = ((1 - mu) / var - 1 / mu) * mu ^ 2
	beta = alpha * (1 / mu - 1)
	return alpha, beta
end

# Returns a list of nIterations+1 MCS in which the first entry is the average of the
# nIterations following entries.
function get_MCS_as_list(m, n, pi, w, beta2, p_time, distribution, nIterations = 10000)
    Random.seed!(RANDOM_SEED)
    exp_p_time = zeros(m, n)
    results = []
	beta = Beta(5, 2)
	triangular = TriangularDist(0.0, 1.0)
	uniform = Uniform(0.0, 1.0)
    for i in 1:nIterations+1
		# job j has the same variation of processing time on both machines (M1 and M2)
		beta_v = rand(beta, m, n)
		triangular_v = rand(triangular, m, n)
		uniform_v = rand(uniform, m, n)
		for j in 1:n
            # generate random processing times p_{ij}
            for k in 1:m
				mean = Float64(p_time[k, pi[j]])  # P_bar is the expected processing time value
				variance = beta2 * mean
				min_value = floor(max(0.0, mean - variance))
				max_value = ceil(mean + variance)
				delta = max_value - min_value
				random_processing_time = 0.0
				# default distribution is uniform
				if distribution == "uniform"
					random_processing_time = min_value + uniform_v[k, j] * delta
				elseif distribution == "lognormal" # Log-normal distribution with log-mean mu and scale sigma
					sigma = sqrt(log(1 + variance / (mean * mean)))
	                mu = log(mean) - sigma * sigma / 2
					d = Truncated(LogNormal(mu, sigma), min_value, max_value)
					random_processing_time = rand(d, 1)[1]  # FIXME invocar tudo de uma vez : rand(d, 10000)
					previous_dev = (random_processing_time - min_value) / delta
				elseif distribution == "normal"
					d = Truncated(Normal(mean, sqrt(variance)), min_value, max_value)
					random_processing_time = rand(d, 1)[1]
				elseif distribution == "triangular"  # symmetric triangular
					random_processing_time = min_value + triangular_v[k, j] * delta
				elseif distribution == "beta"
					random_processing_time = min_value + beta_v[k, j] * delta
				end
                exp_p_time[k, pi[j]] = random_processing_time
            end
        end
       push!(results, calculate_wct(m, n, exp_p_time, w, pi))
   end
   max_twct = maximum(results)
   println("Max TWCT in this simulation: $(max_twct)")
   return results
end

# Export each simulation to a gzipped text file
function save_simulations_to_file(simulation_result_list, output_folder, filename)
	# Open output text file for results
	result_file_path = joinpath(output_folder, filename)
    result_file = GZip.gzopen(result_file_path, "a9")
	for value in simulation_result_list
		println(result_file, "$(value)")
	end
	close(result_file)
    println("DONE.")
end

function parse_permutation_from_string(s)
	return [parse(Int64, ss) for ss in split(s)]
end

function run_monte_carlo_simulation(nIterations, model_list, distribution_list, instance;
		robust_result_filename = "separation_wct_hybrid-manne_randomweights_petro_warmstart-false_singlecut-true.csv")
	rng = MersenneTwister(RANDOM_SEED)
	uid = uuid1(rng)
	println("Simulation execution id: $(uid)")
	batchId = Dates.format(Dates.now(), "yyyy_mm_dd-HH_MM_SS") * "-$(uid)"
	instance_dict = Dict()
	for distribution in distribution_list
		if "robust" in model_list
			previous_filepath = ""
			# Read robust results into dataframe
			df_rob = DataFrame(CSV.File(joinpath(normpath(output_folder), "run_ccg_rpfs_wct_global",
					robust_result_filename); delim=","))
			println(df_rob)
			#df_rob = df_rob[!, (df_rob["model"] .== "Wagner")]
			#df_rob = filter(row -> row.model ∈ ["Wagner"], df_rob)
			for row in eachrow(df_rob)
				instance_name = strip(row["instance_name"])
				if !(instance_name == instance)
					continue
				end
				println("=========================================================")
				println("    (rob)  Processing instance $(instance_name)...")
				println("=========================================================")
				n = row["n"]
				alpha = row["alpha"]
				Γ = row["budget_T"]
				perm = parse_permutation_from_string(row["permutation"])
				fullfilepath = joinpath(instances_folder, "petro", "v3", instance_name)
				if !haskey(instance_dict, fullfilepath)
					m, n, P_bar, P_hat, w = read_robust_input_file(fullfilepath)
					instance_dict[fullfilepath] = m, n, P_bar, P_hat, w
				else
					m, n, P_bar, P_hat, w = instance_dict[fullfilepath]
				end
				p_time_nominal = P_bar
				# The alpha parameter is equivalent to the beta2 parameter of the MCS
				alpha = 100  # row["alpha"]
				beta2 = Float64(alpha) / 100.0
				rob_simulations = get_MCS_as_list(m, n, perm, w, beta2, p_time_nominal, distribution, nIterations)
				filename = "MCS_rob_$(Γ)_$(instance_name)_$(alpha)_$(distribution)_iter$(nIterations).txt.gz"
				sim_output_folder = abspath(create_full_dir(EXPERIMENT_OUTPUT_FOLDER,
					["robust_pfsp", "$(distribution)", "alpha$(alpha)%"]))  # [batchId,
				save_simulations_to_file(rob_simulations, sim_output_folder, filename)
				previous_filepath = fullfilepath
			end
		end
		if "simgrasp" in model_list
			println("Processing sim_grasp for distribution $(distribution)...")
			# Read deterministic results (det0 and det100) into dataframe
			df_sim = DataFrame(CSV.File(joinpath(normpath(output_folder), "run_simgrasp_twct",
					"run_simgrasp_twct_stochgrasp_results.csv"); delim=","))
			#df_ = by(df_sim, [:simgrasp_instance, :n, :m, :nehtdist, :beta1, :beta2, :rob_pfsp_instance, :alpha, :seq],
			#	df_sim -> DataFrame(cost = [df_sim[:stochsol_exp_cost]], pi = [df_sim[:stochsol_permutation]]))
			#print("by: $(df_)")
			for row in eachrow(df_sim)  # Process all 25 SimGRASP solutions
				instance_name = row["instance_filename"]
				if !(instance_name == instance)
					continue
				end
				println("=========================================================")
				println("    (SimGRASP)  Processing instance $(instance_name)...")
				println("=========================================================")
				n = row["n"]
				alpha = row["alpha"]
				sol_index = argmin(row["stochsol_exp_cost"])
				pi = parse_permutation_from_string(row["stochsol_permutation"])
				fullfilepath = joinpath(instances_folder, "petro", "v3", instance_name)
				if !haskey(instance_dict, fullfilepath)
					m, n, P_bar, P_hat, w = read_robust_input_file(fullfilepath)
					instance_dict[fullfilepath] = m, n, P_bar, P_hat, w
				else
					m, n, P_bar, P_hat, w = instance_dict[fullfilepath]
				end
				p_time_nominal = P_bar
				# The alpha parameter is equivalent to the beta2 parameter of the MCS
				alpha = row["alpha"]
				beta2 = Float64(alpha) / 100.0
				simgrasp_simulations = get_MCS_as_list(m, n, pi, w, beta2, p_time_nominal, distribution, 10000)
				filename = "MCS_SimGRASP_$(instance_name)_$(alpha)_$(distribution)_iter$(nIterations).txt.gz"
				sim_output_folder = abspath(create_full_dir(EXPERIMENT_OUTPUT_FOLDER,
					["simgrasp", "$(distribution)", "alpha$(alpha)%"]))  # [batchId,
				save_simulations_to_file(simgrasp_simulations, sim_output_folder, filename)
				previous_filepath = fullfilepath
			end
		end
		flush(stdout)
	end
end

#run_monte_carlo_simulation(10000, ["deterministic"], ["normal"], "RB1505006.txt")# "RB1505007.txt")
#run_monte_carlo_simulation(10000, ["robust", "simgrasp"], ["triangular", "lognormal", "uniform", "normal"], "instance_1_15_5_wct_inputs.txt")
###run_monte_carlo_simulation(10000, ["robust"], ["triangular", "lognormal", "uniform", "normal"], "instance_15_5_R100_wct_inputs.txt";
###	robust_result_filename = "separation_wct_hybrid-wilson_randomweights_petro-random_warmstart-false_singlecut-true_alpha-R100_seq-15.csv")
#run_monte_carlo_simulation(10000, ["robust"], ["triangular", "lognormal", "uniform", "normal"], "instance_15_5_R30_wct_inputs.txt")
#run_monte_carlo_simulation(10000, ["robust"], ["triangular", "lognormal", "uniform", "normal"], "instance_15_5_R50_wct_inputs.txt")
#run_monte_carlo_simulation(10000, ["robust", "simgrasp"], ["triangular", "lognormal", "uniform", "normal"], "instance_2_9_4_wct_inputs.txt")
#run_monte_carlo_simulation(10000, ["robust"], ["triangular", "lognormal", "uniform", "normal"], "instance_9_4_R100_wct_inputs.txt")
#run_monte_carlo_simulation(10000, ["robust"], ["triangular", "lognormal", "uniform", "normal"], "instance_9_4_R30_wct_inputs.txt")
#run_monte_carlo_simulation(10000, ["robust"], ["triangular", "lognormal", "uniform", "normal"], "instance_9_4_R50_wct_inputs.txt")
run_monte_carlo_simulation(10000, ["robust"], ["triangular", "lognormal", "uniform", "normal"], "instance_9_4_R100_wct_inputs.txt";
	robust_result_filename = "separation_wct_hybrid-wilson_randomweights_petro-random_warmstart-false_singlecut-true_alpha-R100_seq-9.csv")
