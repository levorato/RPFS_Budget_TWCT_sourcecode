# =============================================================================
# pfsp_monte_carlo_simulation.jl
# =============================================================================
# Julia script to execute a Monte Carlo simulation for the Permutation Flowshop
# (Cmax objective), based on a given sequence of jobs (permutation).
# =============================================================================

# TODO Apply simulation to all PFS Solutions (det0, det100, robpfs and stochpfs) using the same random seed.
using Random
using UUIDs
using Dates
using CSV
using DataFrames
using Distributions
using GZip
include("../../src/julia/config.jl")
include("../../src/julia/pfsp_file_reader.jl")
include("../../src/julia/deterministic_pfsp.jl")

# If possible, do not modify the lines below
EXPERIMENT_NAME = "montecarlo_sim_2rpfs"
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
function get_MCS_as_list(m, n, pi, beta2, p_time, distribution, nIterations = 10000)
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
       push!(results, calculate_cmax(m, n, exp_p_time, pi))
   end
   max_cmax = maximum(results)
   println("Max makespan in this simulation: $(max_cmax)")
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

function run_monte_carlo_simulation(nIterations, model_list, distribution_list, instance)
	rng = MersenneTwister(RANDOM_SEED)
	uid = uuid1(rng)
	println("Simulation execution id: $(uid)")
	batchId = Dates.format(Dates.now(), "yyyy_mm_dd-HH_MM_SS") * "-$(uid)"
	instance_dict = Dict()
	for distribution in distribution_list
		if "deterministic" in model_list
			# Read deterministic results (det0 and det100) into dataframe
			df_det = DataFrame(CSV.File(joinpath(normpath(output_folder), "consolidated",
					"PFSP_Cmax_deterministic_all_results.csv"); delim=";"))
			#println("$(df_det[1:4, :])")
			#df_det = df_det[(df_det["model"] .== "Wagner") .& (df_det["Gamma1"] .== 20) .& (df_det["Gamma2"] .== 20), :]
			df_det = filter(row -> row.model ∈ ["Wagner"], df_det)
			df_det = filter(row -> row.Gamma1 ∈ [20], df_det)
			df_det = filter(row -> row.Gamma2 ∈ [20], df_det)
			for row in eachrow(df_det)
				instance_name = row["instance_name"]
				if !(instance_name == instance)
					continue
				end
				println("=========================================================")
				println("    (det)  Processing instance $(instance_name)...")
				println("=========================================================")
				n = row["n"]
				alpha = row["alpha"]
				pi = parse_permutation_from_string(row["permutation"])
				perc_deviation_p_bar = row["perc_deviation_p_bar"]
				fullfilepath = joinpath(robust_ying_instances_folder, "data", "$(n)jobs", "$(alpha)%", instance_name)
				if !haskey(instance_dict, fullfilepath)
					m, n, P_bar, P_hat = read_robust_ying_input_file(fullfilepath)
					instance_dict[fullfilepath] = m, n, P_bar
				else
					m, n, P_bar = instance_dict[fullfilepath]
				end
				p_time_nominal = P_bar
				# The alpha parameter is equivalent to the beta2 parameter of the MCS
				beta2 = Float64(row["alpha"]) / 100.0
				det_simulations = get_MCS_as_list(m, n, pi, beta2, p_time_nominal, distribution, nIterations)
				filename = "MCS_det$(perc_deviation_p_bar)_$(instance_name)_$(alpha)_$(distribution)_iter$(nIterations).txt.gz"
				sim_output_folder = abspath(create_full_dir(EXPERIMENT_OUTPUT_FOLDER,
					[batchId, "deterministic_pfsp", "$(distribution)", "alpha$(alpha)%"]))
				save_simulations_to_file(det_simulations, sim_output_folder, filename)
				previous_filepath = fullfilepath
			end
		end
		if "robust" in model_list
			previous_filepath = ""
			# Read robust results into dataframe
			df_rob = DataFrame(CSV.File(joinpath(normpath(output_folder), "run_ccg_rpfs_cmax",
					"2RPFS_Cmax_all_results.csv"); delim=";"))
			#df_rob = df_rob[!, (df_rob["model"] .== "Wagner")]
			df_rob = filter(row -> row.model ∈ ["Wagner"], df_rob)
			for row in eachrow(df_rob)
				instance_name = row["instance_name"]
				if !(instance_name == instance)
					continue
				end
				println("=========================================================")
				println("    (rob)  Processing instance $(instance_name)...")
				println("=========================================================")
				n = row["n"]
				alpha = row["alpha"]
				Γ1 = row["Gamma1"]
				Γ2 = row["Gamma2"]
				pi = parse_permutation_from_string(row["permutation"])
				fullfilepath = joinpath(robust_ying_instances_folder, "data", "$(n)jobs", "$(alpha)%", instance_name)
				if !haskey(instance_dict, fullfilepath)
					m, n, P_bar, P_hat = read_robust_ying_input_file(fullfilepath)
					instance_dict[fullfilepath] = m, n, P_bar
				else
					m, n, P_bar = instance_dict[fullfilepath]
				end
				p_time_nominal = P_bar
				# The alpha parameter is equivalent to the beta2 parameter of the MCS
				beta2 = Float64(row["alpha"]) / 100.0
				rob_simulations = get_MCS_as_list(m, n, pi, beta2, p_time_nominal, distribution, nIterations)
				filename = "MCS_rob_$(Γ1)_$(Γ2)_$(instance_name)_$(alpha)_$(distribution)_iter$(nIterations).txt.gz"
				sim_output_folder = abspath(create_full_dir(EXPERIMENT_OUTPUT_FOLDER,
					[batchId, "robust_pfsp", "$(distribution)", "alpha$(alpha)%"]))
				save_simulations_to_file(rob_simulations, sim_output_folder, filename)
				previous_filepath = fullfilepath
			end
		end
		if "simgrasp" in model_list
			println("Processing sim_grasp for distribution $(distribution)...")
			# Read deterministic results (det0 and det100) into dataframe
			df_sim = DataFrame(CSV.File(joinpath(normpath(output_folder), "consolidated",
					"simgrasp_cmax_ying_stochgrasp_results.csv"); delim=","))
			#df_ = by(df_sim, [:simgrasp_instance, :n, :m, :nehtdist, :beta1, :beta2, :rob_pfsp_instance, :alpha, :seq],
			#	df_sim -> DataFrame(cost = [df_sim[:stochsol_exp_cost]], pi = [df_sim[:stochsol_permutation]]))
			#print("by: $(df_)")
			for row in eachrow(df_sim)  # Process all 25 SimGRASP solutions
				instance_name = row["rob_pfsp_instance"]
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
				fullfilepath = joinpath(robust_ying_instances_folder, "data", "$(n)jobs", "$(alpha)%", instance_name)
				if !haskey(instance_dict, fullfilepath)
					m, n, P_bar, P_hat = read_robust_ying_input_file(fullfilepath)
					instance_dict[fullfilepath] = m, n, P_bar
				else
					m, n, P_bar = instance_dict[fullfilepath]
				end
				p_time_nominal = P_bar
				# The alpha parameter is equivalent to the beta2 parameter of the MCS
				beta2 = Float64(row["alpha"]) / 100.0
				simgrasp_simulations = get_MCS_as_list(m, n, pi, beta2, p_time_nominal, distribution, 10000)
				filename = "MCS_SimGRASP_$(instance_name)_$(alpha)_$(distribution)_iter$(nIterations).txt.gz"
				sim_output_folder = abspath(create_full_dir(EXPERIMENT_OUTPUT_FOLDER,
					[batchId, "simgrasp", "$(distribution)", "alpha$(alpha)%"]))
				save_simulations_to_file(simgrasp_simulations, sim_output_folder, filename)
				previous_filepath = fullfilepath
			end
		end
	end
end

#run_monte_carlo_simulation(10000, ["deterministic"], ["normal"], "RB1505006.txt")# "RB1505007.txt")
#run_monte_carlo_simulation(10000, ["deterministic"], ["normal"], "RB1505007.txt")# "RB1505007.txt")
run_monte_carlo_simulation(10000, ["deterministic", "robust"], ["triangular", "lognormal", "uniform", "normal"], "RB1501008.txt")# "RB1505007.txt")
run_monte_carlo_simulation(10000, ["deterministic", "robust"], ["triangular", "lognormal", "uniform", "normal"], "RB1502008.txt")# "RB1505007.txt")
run_monte_carlo_simulation(10000, ["deterministic", "robust"], ["triangular", "lognormal", "uniform", "normal"], "RB1503008.txt")# "RB1505007.txt")
run_monte_carlo_simulation(10000, ["deterministic", "robust"], ["triangular", "lognormal", "uniform", "normal"], "RB1504008.txt")# "RB1505007.txt")
run_monte_carlo_simulation(10000, ["deterministic", "robust"], ["triangular", "lognormal", "uniform", "normal"], "RB1505008.txt")# "RB1505007.txt")

# run_monte_carlo_simulation(10000, ["robust"], ["triangular", "lognormal", "uniform", "normal"], "RB2005010.txt")# "RB1505007.txt")
#run_monte_carlo_simulation(10000, ["deterministic", "robust", "simgrasp"], ["triangular", "lognormal", "uniform"])
#run_monte_carlo_simulation(10000, ["deterministic", "robust", "simgrasp"], ["triangular", "lognormal", "uniform", "normal"], "RB2005010.txt")# "RB1505007.txt")
