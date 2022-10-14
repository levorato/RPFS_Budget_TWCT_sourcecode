# =============================================================================
# consolidate_simgrasp_results.jl
# =============================================================================
# Julia script to consolidate SimGRASP output files for the Stochastic PFSP
# (Cmax objective), based on the GRASP simheuristic.
# =============================================================================

include("../../src/julia/config.jl")
include("../../src/julia/pfsp_file_reader.jl")
include("../../src/julia/robust_pfsp_budget_worstcase_m_mach.jl")


EXPERIMENT_NAME = "simgrasp_cmax_ying"
EXPERIMENT_INSTANCE_FOLDER = abspath(joinpath(instances_folder, "robust", "ying", "data"))
SIMGRASP_OUTPUTDIR = abspath(joinpath(home_prefix, "pfsp_experiments", "simgrasp_outputs"))
SCRIPT_OUTPUT_FOLDER = abspath(create_full_dir(joinpath(home_prefix, "pfsp_experiments"),
	["results", EXPERIMENT_NAME]))


# interpret the filename of simgrasp result
# e.g. RB0105001_10_2_t_1.0_0.1_124341_outputs.txt
function interpret_filename(filename)
	items = split(filename, "_")
	instance_name = items[1]
	n = parse(Int64, items[2])
	m = parse(Int64, items[3])
	nehtdist = items[4]
	beta1 = parse(Float64, items[5])
	beta2 = parse(Float64, items[6])
	seed = parse(Int64, items[7])
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
function read_simgrasp_solution_file(path, filename)
	neh_sol = Dict()  # id, cost, expcost, time, permutation
	best_sol = Dict()  #
	stoch_sol = Dict()  #
	content = open(joinpath(path, filename)) do file  # read whole file to str
	    read(file, String)
	end
	i1 = findfirst("NEH", content)[1]
	i2 = findfirst("Our best solution", content)[1]
	i3 = findfirst("Our best stoch-sol", content)[1]
    neh_sol = read_individual_solution(content[i1:i2-1])
    best_sol = read_individual_solution(content[i2:i3-1])
	stoch_sol = read_individual_solution(content[i3:length(content)])
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

function calculate_worstcase_costs(m, n, P_bar, P_hat, permutation)
	s = ""
	for p_Γ1 in [20, 40, 60, 80, 100]
		for p_Γ2 in [20, 40, 60, 80, 100]
			Γ1 = (p_Γ1 * n) / 100.0
			Γ2 = (p_Γ2 * n) / 100.0
			worstcase_cost, worst_Γ_scenario = worst_case_cmax_dp_m_machines(permutation, m, n, [Γ1, Γ2], P_bar, P_hat)
			s = s * ",$(worstcase_cost)"
		end
	end
	return s
end

function process_simgrasp_outputdir()
    outputfile_list = ["$(file)" for file in searchdir(SIMGRASP_OUTPUTDIR, "outputs.txt")]
    output_count = length(outputfile_list)
    println("Outputs count: $(output_count)")
	Γ_header = ""
	for p_Γ1 in [20, 40, 60, 80, 100]
		for p_Γ2 in [20, 40, 60, 80, 100]
			Γ_header = Γ_header * ",$(p_Γ1) $(p_Γ2)"
		end
	end
	# Open output CSV file for results
    result_file_path_neh = joinpath(SCRIPT_OUTPUT_FOLDER, EXPERIMENT_NAME * "_neh_results.csv")
    result_file_neh = open(result_file_path_neh, "w+")
	println("Saving consolidated result file to " * result_file_path_neh)
	println(result_file_neh, "filename,simgrasp_instance,n,m,nehtdist,beta1,beta2,seed,rob_pfsp_instance,alpha,seq,"
		* "neh_sol_id,neh_cost,neh_exp_cost,neh_time,neh_permutation$(Γ_header)")
	result_file_path_detgrasp = joinpath(SCRIPT_OUTPUT_FOLDER, EXPERIMENT_NAME * "_detgrasp_results.csv")
    result_file_detgrasp = open(result_file_path_detgrasp, "w+")
	println("Saving consolidated result file to " * result_file_path_detgrasp)
	println(result_file_detgrasp, "filename,simgrasp_instance,n,m,nehtdist,beta1,beta2,seed,rob_pfsp_instance,alpha,seq,"
		* "bestsol_sol_id,bestsol_cost,bestsol_exp_cost,bestsol_time,bestsol_permutation$(Γ_header)")
	result_file_path_stochgrasp = joinpath(SCRIPT_OUTPUT_FOLDER, EXPERIMENT_NAME * "_stochgrasp_results.csv")
    result_file_stochgrasp = open(result_file_path_stochgrasp, "w+")
	println("Saving consolidated result file to " * result_file_path_stochgrasp)
	println(result_file_stochgrasp, "filename,simgrasp_instance,n,m,nehtdist,beta1,beta2,seed,rob_pfsp_instance,alpha,seq,"
		* "stochsol_sol_id,stochsol_cost,stochsol_exp_cost,stochsol_time,stochsol_permutation$(Γ_header)")

	for filename in outputfile_list
		instance_name, n, m, nehtdist, beta1, beta2, seed = interpret_filename(filename)
		println("Instance: $(instance_name), n=$(n), m=$(m), nehtdist=$(nehtdist), beta1=$(beta1), beta2=$(beta2), seed=$(seed)")
		alpha = Int64(ceil(beta2 * 100))
		seq = instance_name[8:9]
		ying_instance_filename = "RB$(lpad(n,3,"0"))$(alpha)$(lpad(seq,2,"0")).txt"
		println("$(ying_instance_filename),$(alpha),$(seq)")
		instance_filepath = abspath(joinpath(EXPERIMENT_INSTANCE_FOLDER, "$(n)jobs", "$(alpha)%", ying_instance_filename))
		m, n, P_bar, P_hat = read_robust_ying_input_file(instance_filepath)
		neh_sol, best_sol, stoch_sol = read_simgrasp_solution_file(SIMGRASP_OUTPUTDIR, filename)

		print(result_file_neh, "$(filename),$(instance_name),$(n),$(m),$(nehtdist),$(beta1),$(beta2),$(seed),"
			* "$(ying_instance_filename),$(alpha),$(seq)")
		print(result_file_detgrasp, "$(filename),$(instance_name),$(n),$(m),$(nehtdist),$(beta1),$(beta2),$(seed),"
			* "$(ying_instance_filename),$(alpha),$(seq)")
		print(result_file_stochgrasp, "$(filename),$(instance_name),$(n),$(m),$(nehtdist),$(beta1),$(beta2),$(seed),"
			* "$(ying_instance_filename),$(alpha),$(seq)")
		print(result_file_neh, ",$(neh_sol["sol_id"])")
		print(result_file_neh, ",$(neh_sol["cost"])")
		print(result_file_neh, ",$(neh_sol["exp_cost"])")
		print(result_file_neh, ",$(neh_sol["time"])")
		print(result_file_neh, ",$(permutation_to_string(neh_sol["permutation"]))")
		println(result_file_neh, "$(calculate_worstcase_costs(m, n, P_bar, P_hat, neh_sol["permutation"]))")

		print(result_file_detgrasp, ",$(best_sol["sol_id"])")
		print(result_file_detgrasp, ",$(best_sol["cost"])")
		print(result_file_detgrasp, ",$(best_sol["exp_cost"])")
		print(result_file_detgrasp, ",$(best_sol["time"])")
		print(result_file_detgrasp, ",$(permutation_to_string(best_sol["permutation"]))")
		println(result_file_detgrasp, "$(calculate_worstcase_costs(m, n, P_bar, P_hat, best_sol["permutation"]))")

		print(result_file_stochgrasp, ",$(stoch_sol["sol_id"])")
		print(result_file_stochgrasp, ",$(stoch_sol["cost"])")
		print(result_file_stochgrasp, ",$(stoch_sol["exp_cost"])")
		print(result_file_stochgrasp, ",$(stoch_sol["time"])")
		print(result_file_stochgrasp, ",$(permutation_to_string(stoch_sol["permutation"]))")
		println(result_file_stochgrasp, "$(calculate_worstcase_costs(m, n, P_bar, P_hat, stoch_sol["permutation"]))")

		flush(result_file_neh)
		flush(result_file_detgrasp)
		flush(result_file_stochgrasp)
	end
    close(result_file_neh)
	close(result_file_detgrasp)
	close(result_file_stochgrasp)
	println("DONE.")
end

process_simgrasp_outputdir()
