# =============================================================================
# calculate_rpfs_budget_worstcase_costs.jl
# =============================================================================
# Julia script to calculate worstcase costs for the RPFS Budgted Uncertainty
# (Cmax objective), based on an existing solution (permutation \sigma).
# =============================================================================

include("../../src/julia/config.jl")
include("../../src/julia/pfsp_file_reader.jl")
include("../../src/julia/robust_pfsp_budget_worstcase_m_mach.jl")

using CSV
using DataFrames


EXPERIMENT_NAME = "calculate_rpfs_cmax_worstcase"
EXPERIMENT_INSTANCE_FOLDER = abspath(joinpath(instances_folder, "robust"))
CONSOLIDATED_RESULT_DIR = abspath(joinpath(home_prefix, "pfsp_experiments", "consolidated"))
SCRIPT_OUTPUT_FOLDER = abspath(create_full_dir(joinpath(home_prefix, "pfsp_experiments"),
	["results", EXPERIMENT_NAME]))


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

function calculate_worstcase_costs(m, n, P_bar, P_hat, permutation)
	s = ""
	for p_Γ1 in [20, 40, 60, 80, 100]
		for p_Γ2 in [20, 40, 60, 80, 100]
			Γ1 = (p_Γ1 * n) / 100.0
			Γ2 = (p_Γ2 * n) / 100.0
			worstcase_cost, worst_Γ_scenario = worst_case_cmax_dp_m_machines(permutation, m, n, [Γ1, Γ2], P_bar, P_hat)
			s = s * ";$(worstcase_cost)"
		end
	end
	return s
end

# Result Type can be: "deterministic" or "robust".
function process_consolidated_csv_result(csv_filename, result_type, exact = true)
	# Open output CSV file for results
    result_file_path = joinpath(SCRIPT_OUTPUT_FOLDER, EXPERIMENT_NAME * "_" * csv_filename)
    result_file = open(result_file_path, "w+")
	println("Saving worstcase costs result file to " * result_file_path)
	print(result_file, "executionId;model;instance_name;alpha;n;m;perc_deviation_p_bar;")
	print(result_file, "Gamma1;Gamma2;cmax_star;permutation")
	s = ""
	for p_Γ1 in [20, 40, 60, 80, 100]
		for p_Γ2 in [20, 40, 60, 80, 100]
			s = s * ";$(p_Γ1) $(p_Γ2)"
		end
	end
	println(result_file, "$(s)")
	# Read deterministic results (det0 and det100) into dataframe
	df = DataFrame(CSV.File(joinpath(normpath(CONSOLIDATED_RESULT_DIR),
			csv_filename); delim=";"))
	#println("$(df[1:4, :])")
	if exact
		df = df[(df["model"] .== "Wagner"), :]
	end
	instance_dict = Dict()
	for row in eachrow(df)
		instance_name = row["instance_name"]
		println("=========================================================")
		println("    Processing instance $(instance_name)...")
		println("=========================================================")
		executionId = row["executionId"]
		model = row["model"]
		n = row["n"]
		m = row["m"]
		alpha = row["alpha"]
		permutation = row["permutation"]
		pi = parse_permutation_from_string(permutation)
		perc_deviation_p_bar = -1
		Gamma1 = -1
		Gamma2 = -1
		cost = 0
		if result_type == "determinisic"
			perc_deviation_p_bar = row["perc_deviation_p_bar"]
			cost = row["cmax_star"]
		elseif result_type == "robust"
			Gamma1 = row["Gamma1"]
			Gamma2 = row["Gamma2"]
			cost = row["cmax"]
		end
		fullfilepath = joinpath(robust_ying_instances_folder, "$(n)jobs", "$(alpha)%", instance_name)
		if !haskey(instance_dict, fullfilepath)
			m, n, P_bar, P_hat = read_robust_ying_input_file(fullfilepath)
			instance_dict[fullfilepath] = m, n, P_bar, P_hat
		else
			m, n, P_bar, P_hat = instance_dict[fullfilepath]
		end
		println("Calculating PFSP worst-case budget scenarios...")
		print(result_file, "$(executionId);$(model);$(instance_name);$(alpha);$(n);$(m);")
		print(result_file, "$(perc_deviation_p_bar);$(Gamma1);$(Gamma2);$(cost);")
		print(result_file, "$(permutation)")
		println(result_file, "$(calculate_worstcase_costs(m, n, P_bar, P_hat, pi))")
	end
    close(result_file)
	println("DONE.")
end

process_consolidated_csv_result("2RPFS_Cmax_all_results.csv", "robust")
# process_consolidated_csv_result("PFSP_Cmax_deterministic_all_results.csv", "deterministic")
