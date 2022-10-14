# ====================================================================
# robust_pfsp_cmax_budget_separation_common.jl
# Common functions for Rob PFSP BUdget Separation procedures (Cmax).
# ====================================================================

include("robust_pfsp_budget_worstcase_m_mach.jl")
include("robust_pfsp_budget_worstcase.jl")
include("robust_pfsp_budget_worstcase_global_multi_scenario.jl")
include("robust_pfsp_budget_worstcase_global_single_scenario.jl")
include("robust_pfsp_alt_worstcase_cmax_mip.jl")
include("robust_pfsp_budget_worstcase_cmax_brute.jl")


function solve_separation_problem_cmax(m, n, Γ, P_bar, P_hat, permutation, verbose = true;
			single_solution = true, budget_type = "machine", sp_type = "dp",
			solver = "CPLEX", max_cores = 16, solver_time_limit = 7200,
			fractional_solution = false, sp_bigM_type=7)
    #cmax_dp, (m1, m2) = worst_case_cmax_brute_force(permutation, n, Γ1, Γ2, P_bar, P_hat)
	worst_Γ_scenario_dp_list = []
	cmax_dp = 0.0
	z_sep_list = []
	if sp_type == "dp"
		if budget_type == "machine"
			if single_solution
				cmax_dp, worst_Γ_scenario_dp = worst_case_cmax_dp_m_machines_v1(permutation, m, n, Γ, P_bar, P_hat)
				push!(worst_Γ_scenario_dp_list, worst_Γ_scenario_dp)
			else
		    	cmax_dp, worst_Γ_scenario_dp_list = worst_case_cmax_dp_m_machines(permutation, m, n, Γ, P_bar, P_hat)
			end
			for worst_Γ_scenario_dp in worst_Γ_scenario_dp_list
				z_sep = zeros(Int64, (m, n))
				for q in 1:m
					for v in worst_Γ_scenario_dp[q]
						z_sep[q, v] = 1
					end
				end
				push!(z_sep_list, z_sep)
			end
		else  # "global"
			if single_solution
				cmax_dp, worst_Γ_scenario_dp_list = worst_case_cmax_dp_m_machines_global_single(permutation, m, n, Γ, P_bar, P_hat)
			else
				cmax_dp, worst_Γ_scenario_dp_list = worst_case_cmax_dp_m_machines_global_multi(permutation, m, n, Γ, P_bar, P_hat)
			end
			for worst_Γ_scenario_dp in worst_Γ_scenario_dp_list
				z_sep = zeros(Int64, (m, n))
				for (r, i) in worst_Γ_scenario_dp
					z_sep[r, i] = 1
				end
				push!(z_sep_list, z_sep)
			end
		end
		if size(z_sep_list, 1) > 0
			return cmax_dp, true, z_sep_list
		else
			println("*** ERROR! No DP scenario found! Will fallback to MILP worst-case solution.\n")
		end
	end
	# FALLBACK or "mip"
	if budget_type == "global"
		if worst_case_brute_force_is_worth(m, n, Int64(ceil(Γ)))
			max_solutions = 5
			if single_solution
				max_solutions = 1
			end
			time_spent = @elapsed cmax, list_worst_Γ_scenario = worst_case_cmax_brute_force_global_budget(permutation, m, n, Γ, P_bar, P_hat, true;
																	max_solutions=max_solutions)
			println("Separation time: $(time_spent) s")
			flush(stdout)
			list_worst_Γ_scenario_z_sep = [convert_scenario_job_deviation_list_to_matrix(m, n, s) for s in list_worst_Γ_scenario]
			is_opt = size(list_worst_Γ_scenario_z_sep, 1) > 0
			return cmax, is_opt, list_worst_Γ_scenario_z_sep
		end
	end
	cmax_mip, z_sep_list, time_spent = solve_robust_model_worst_case_alt(permutation, m, n, Γ, P_bar, P_hat;
			max_cores=max_cores, single_cut_per_iteration=single_solution, budget_type=budget_type,
			fractional_solution=fractional_solution, solver_time_limit=solver_time_limit, solver=solver, bigM_type=sp_bigM_type)
	return cmax_mip, true, z_sep_list
end


function get_current_model_solution(model, n)
    Jobs = 1:n
    Z_star = zeros(Float64, (n, n))
    if JuMP.termination_status(model) == MOI.OPTIMAL
        #println("\n===========  R O B U S T    S O L U T I O N  ===========\n")
        #println("Optimal Objective Function value: ", objective_value(model))
        #println("Solve time : ", solve_time(model))
        permutation = zeros(Int64, n)
        for k in Jobs  # Position k
            for j in Jobs  # Job j
                x_jk = variable_by_name(model, "Z[$(j),$(k)]")
                if value(x_jk) >= 0.9  # Z[j, k]
                    Z_star[j, k] = 1
                    if permutation[k] == 0
                        permutation[k] = j
                    else
                        println("ERROR assembling permutation from model solution!")
						flush(stdout)
                        return 0, [], [], false
                    end
                end
            end
        end
		flush(stdout)
        return objective_value(model), permutation, Z_star, true
    else
        println("ERROR: [get_current_model_solution] Optimal solution not found! Best bound: ", objective_bound(model))
        println("Optimal solution not found! Best bound: ", objective_bound(model))
        println("Model solve time: $(solve_time(model)); time_limit = $(time_limit)")
		flush(stdout)
        return objective_bound(model), zeros(Int64, n), Z_star, false  #, solve_time(model), false
    end
end


function load_cut_list_from_csv(log_cuts_dir)
    existing_cut_list = []
	max_iter = 1
	UB = typemax(Float64)
	best_pi = []
    println("Reading cuts from dir: $(log_cuts_dir)")
	for filename in readdir(log_cuts_dir)
		if startswith(filename, "zsep_")
			filepath = joinpath(log_cuts_dir, filename)
			z_sep = readdlm(filepath)
			push!(existing_cut_list, z_sep)
			iter_str = filename[findfirst("zsep_", filename).stop+1:findfirst("-", filename)[1]-1]
			#println("iter_str = $(iter_str)")
			iter = parse(Int, iter_str)
			max_iter = max(max_iter, iter)
		end
    end
	count = size(existing_cut_list, 1)
	if count > 0
		UB_array = readdlm(joinpath(log_cuts_dir, "UB.txt"))
		UB = UB_array[1]
		best_pi = readdlm(joinpath(log_cuts_dir, "best_pi.txt"), '\t', Int, '\n')
	end
	println("Sucessfully read $(count) cuts with max_iter = $(max_iter) and UB = $(UB).")
	flush(stdout)
    return existing_cut_list, max_iter, UB, best_pi
end

function create_jump_model(solver)
	if solver == "Gurobi"
		println("Using Gurobi solver.")
		model = Model(Gurobi.Optimizer)
	elseif solver == "CPLEX"
		println("Using CPLEX solver.")
		model = Model(CPLEX.Optimizer)
	else
	  println("No solver defined")
	  model = nothing
	end
	return model
end

function setup_solver_params(model, solver, solver_time_limit, max_cores)
	if solver == "Gurobi"
		MOI.set(model, MOI.RawOptimizerAttribute("TimeLimit"), solver_time_limit)
		MOI.set(model, MOI.RawOptimizerAttribute("Threads"), max_cores)
		# https://www.gurobi.com/documentation/9.1/refman/mipgapabs.html
		MOI.set(model, MOI.RawOptimizerAttribute("MIPGapAbs"), 1e-5)   # default is 1e-10
	elseif solver == "CPLEX"
		MOI.set(model, MOI.RawOptimizerAttribute("CPX_PARAM_TILIM"), solver_time_limit)
		MOI.set(model, MOI.RawOptimizerAttribute("CPX_PARAM_THREADS"), max_cores)
		MOI.set(model, MOI.RawOptimizerAttribute("CPX_PARAM_EPAGAP"), 1e-5)   # default is 1e-6
	else
	  println("No solver defined")
	end
end
