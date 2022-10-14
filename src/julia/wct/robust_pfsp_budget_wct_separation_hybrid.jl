# ===================================================================================================
# robust_pfsp_budget_wct_separation_hybrid.jl
# Robust Permutation Flow Shop Problem to Minimize the Weighted Sum of Job Completion Times (WCT),
# based on the budgeted uncertainty set (Bertsimas et al.).
# Given all possible scenarios (i.e. for each job i, on which machines it will have its processing times
# deviated to their worst-case values), finds a permutation pi that minimizes the restricted worst-case
# WCT objective function.
# *** Solve the robust problem via separation approach (Column and Constraint Generation procedure),
# based on a Hybrid PFSP Combinatorial+MILP Model, where the Master Problem MP is coded in C++.
# ===================================================================================================

# Worst-case calculation for the robust budget problem
include("robust_pfsp_budget_separation_common.jl")
include("robust_pfsp_alt_worstcase_wct_mip.jl")
include("robust_pfsp_budget_worstcase_wct_brute.jl")

using JuMP
using CPLEX
using Gurobi
using Dates

EXPERIMENT_NAME = "run_ccg_rpfs_wct_global"
# If possible, do not modify the lines below
EXPERIMENT_OUTPUT_FOLDER = abspath(create_full_dir(normpath(output_folder),
	[EXPERIMENT_NAME]))


function solve_master_problem(instance_filepath, Γ, budget_type, cut_list; time_limit = 7200, max_cores = 16,
		model_name = "hybrid-wilson")
	obj_star, permutation_current, elapsed_time, mp_gap, is_mp_opt = solve_pfsp_wct_cpp(instance_filepath, Γ, budget_type;
		optimization_type=model_name, time_limit=time_limit, cut_list=cut_list, max_cores=max_cores)
	return obj_star, permutation_current, elapsed_time, mp_gap, is_mp_opt
end

# According to the separation method results,
# return a set of cuts (new constraints and columns) to be added to the existing model.
# Return a tuple: (separation_done, mip_optimal)
# - separation_done indicates that the separation procedure has found the optimal master problem solution
# - mip_optimal indicates that the newly-generated MILP problem solution is optimal or not
#     ( if it is not optimal, the whole C&CG procedure has to stop and abort! )
# Attention: UB is an in-out function argument (its value if modified)
function separate_hybrid(permutation, m, n, Γ, P_bar, P_hat, w, iteration, LB, UB, best_pi; max_cores = 16,
            single_cut_per_iteration = false, budget_type = "machine", time_limit = 7200, fractional_solution = false,
            use_brute_force = false, solver = "CPLEX", sp_bigM_type = 7, rel_gap_tol = 1e-8)
    println("Generating cuts for permutation = $(permutation)...")
    wct_sep, is_sp_opt, z_sep_list, sp_time = solve_separation_problem(m, n, Γ, P_bar, P_hat, w, permutation, LB, true;
            budget_type=budget_type, time_limit=time_limit, single_cut_per_iteration=single_cut_per_iteration, max_cores=max_cores,
            fractional_solution=fractional_solution, use_brute_force=use_brute_force, solver=solver, bigM_type=sp_bigM_type)
    # Important: round separation values to avoid double-type precision problems
    for z_sep in z_sep_list
        z_sep = round.(z_sep, digits=3)
    end
    if !is_sp_opt
        gap = (UB - LB) / LB
        println("[SP Separation not optimal!] UB (wct_sep) = $(UB)  vs.  LB (wct_current) = $(LB) ; gap = $(gap)")
        flush(stdout)
		# separation_done, sp_optimal, UB, LB, best_pi, gap, sp_time, cut_list
        return true, false, UB, LB, best_pi, gap, sp_time, []
    end
    #   ==>  UB = min(UB, wct_sep)
    if wct_sep < UB
        UB = wct_sep
        best_pi = permutation
    end
    gap = (UB - LB) / LB
    println("UB (wct_sep) = $(UB)  vs.  LB (wct_current) = $(LB) ; gap = $(gap)")
    flush(stdout)
    if gap > rel_gap_tol  # If Gap = (UB - LB) / LB <= EPS, an eps-optimal solution is found, terminate.
        status = false
        # Will return a set of cuts for each worst-case solution / scenario in z_sep_list
    else
        println("SP Separation done.")
        flush(stdout)
		status = true
        z_sep_list = []
    end
    return status, true, UB, LB, best_pi, gap, sp_time, z_sep_list
end # End Function SEPARATE

# ===================================================================================
# Primal cut algorithm
# ===================================================================================
# As described by:
#  1. Zeng, B. (2011). Solving Two-stage Robust Optimization Problems
#                   by A Constraint-and-Column Generation Method. Constraints, 2011.
#  2. Zeng, B., & Zhao, L. (2013). Solving two-stage robust optimization problems
#                   using a column-and- constraint generation method.
#                   Operations Research Letters, 41(5), 457–461.
#                   https://doi.org/10.1016/j.orl.2013.05.003
# ===================================================================================
function solve_robust_model_flowshop_wct_budget_separation_hybrid(m, n, P_bar, P_hat, w, Γ, time_limit; mp_bigM_type = 2, sp_bigM_type = 3, max_cores = 16,
            solver = "CPLEX", budget_type = "machine", single_cut_per_iteration = false, fractional_solution = false, use_brute_force = false,
            warm_start = false, instance_filename = nothing, sp_solver = nothing, rel_gap_tol = 1e-8, mp_model_name = "hybrid-wilson",
			log_cuts = true)
    solver_time_limit = time_limit
    elapsed_time = 0.0
    sp_total_time = 0.0
    mp_total_time = 0.0
    if isnothing(sp_solver)
        sp_solver = solver
    end
	instance_name = basename(instance_filename)
    println("Starting Hybrid C&CG procedure with $(mp_model_name)...")
	println("Instance_name = $(instance_name)")
    println("Procedure time limit = $(time_limit) s")
    println("Solver time limit = $(solver_time_limit) s")
	# write partial results to csv file
	partial_results_full_filename = "partialresults_separation_wct_$(mp_model_name)_$(instance_name)_gamma-$(Γ)_warmstart-$(warm_start)_singlecut-$(single_cut_per_iteration).csv"
	partial_results_file_path = joinpath(EXPERIMENT_OUTPUT_FOLDER, partial_results_full_filename)
	result_file_exists = isfile(partial_results_file_path)
	result_file = open(partial_results_file_path, "a+")
	if !result_file_exists
    	println(result_file, "executionId,ub_name,instance_name,alpha,n,m,budget_T,wct,permutation,time_spent,time_to_best_sol,mp_total_time,sp_total_time,iterations,num_visited_solutions,num_improvements,is_optimal,validated,gap,lb,cost,wct_validation,cut_count")
	end
    iter = 1
    UB = typemax(Float64)
    best_pi = []
    LB = 0.0
    total_cuts = 0
    all_cut_list = []
	if log_cuts  # log the cuts generated by the C&CG procedure in a CSV file
		method_output_folder, output_folder_stdout = generate_hybrid_master_output_path(EXPERIMENT_OUTPUT_FOLDER, mp_model_name, instance_name)
		log_cuts_dir = output_folder_stdout
		mkpath(log_cuts_dir)
		println("Cut dir is $(log_cuts_dir).")
		all_cut_list, iter, UB, best_pi = load_cut_list_and_partial_solution_twct(partial_results_file_path, log_cuts_dir, instance_name, Γ)
		println("=== Resuming problem solve from iteration $(iter), with $(length(all_cut_list)) existing cuts, UB = $(UB) and best_pi = $(best_pi).")
    end
    while true
        # Solve the master problem
        println("====================================================================\n")
        println("Solving WCT budget robust Hybrid model C & CG (Iteration $(iter))...")
        println("====================================================================\n")
        datetimenow = Dates.now()
		time_left = trunc(Int, max(100, solver_time_limit - (mp_total_time + sp_total_time)))
		println("Date: $(datetimenow)")
		println("Time left: $(time_left) s")
		flush(stdout)
		obj_star, permutation_current, elapsed_time, mp_gap, is_mp_opt = solve_master_problem(instance_filename,
			Γ, budget_type, all_cut_list; time_limit=time_left, max_cores=max_cores, model_name=mp_model_name)
        if is_mp_opt
            println("\n===========  $(mp_model_name)  -  R O B U S T    S O L U T I O N  ===========\n")
            println("Optimal Objective Function value: ", obj_star)
            println("Solve time : ", elapsed_time)
            println("Permutation: $(permutation_current)")
            flush(stdout)
            mp_total_time += elapsed_time
			# IMPORTANT : Check for time_limit
            if (mp_total_time + sp_total_time >= time_limit)
                println("Time limit exceeded!")
                LB = obj_star
                gap = (UB - LB) / LB
                println("Model solve time: time_limit = $(time_limit); total_time = $(mp_total_time + sp_total_time); gap = $(gap); UB = $(UB); LB = $(LB)")
                flush(stdout)
                return UB, best_pi, [], mp_total_time + sp_total_time, false, iter, gap, LB, total_cuts, mp_total_time, sp_total_time
            end
			LB = obj_star
            separation_done, sp_optimal, UB, LB, best_pi, gap, sp_time, cut_list = separate_hybrid(permutation_current, m, n, Γ, P_bar, P_hat, w, iter + 1, LB, UB, best_pi;
                    budget_type=budget_type, time_limit=trunc(Int, max(MIN_TIME_SP_MIP, time_limit - (mp_total_time + sp_total_time))),
                    single_cut_per_iteration=single_cut_per_iteration, max_cores=max_cores,
                    fractional_solution=fractional_solution, use_brute_force=use_brute_force, solver=sp_solver, sp_bigM_type=sp_bigM_type,
					rel_gap_tol=rel_gap_tol)
			all_cut_list = unique([all_cut_list ; cut_list])
            sp_total_time += sp_time
            total_cuts = length(all_cut_list)
			# write partial results to csv file
			pi_str = join(["$(x)" for x in best_pi], " ")
			println(result_file, "none, mip_separation, $(instance_name), NA, $(n), $(m), $(Γ), $(UB), $(pi_str), $(mp_total_time + sp_total_time), $(mp_total_time + sp_total_time), $(mp_total_time), $(sp_total_time), $(iter), $(iter), $(iter-1), false, false, $(gap), $(LB), $(UB), $(UB), $(total_cuts)")
			flush(result_file)
            # IMPORTANT : Check for time_limit
            if (mp_total_time + sp_total_time >= time_limit)
                println("Time limit exceeded!")
                println("Model solve time: time_limit = $(time_limit); total_time = $(mp_total_time + sp_total_time); gap = $(gap); UB = $(UB); LB = $(LB)")
                flush(stdout)
				close(result_file)
                return UB, best_pi, [], mp_total_time + sp_total_time, false, iter, gap, LB, total_cuts, mp_total_time, sp_total_time
            end
            if (!sp_optimal)
                println("ERROR: [Sub Problem] Optimal solution not found! Best UB: $(UB) ; LB = $(LB)")
                println("Model solve time: time_limit = $(time_limit); elapsed_time = $(elapsed_time); gap = $(gap); UB = $(UB); LB = $(LB)")
                flush(stdout)
				close(result_file)
                return UB, best_pi, [], mp_total_time + sp_total_time, false, iter, gap, LB, total_cuts, mp_total_time, sp_total_time
            end
            if separation_done
                println(" *** Separation done! obj_star = $(obj_star), best_pi = $(best_pi)")
                flush(stdout)
				close(result_file)
                return UB, best_pi, [], mp_total_time + sp_total_time, true, iter, gap, LB, total_cuts, mp_total_time, sp_total_time
            end  # end if separate
        else
			println("\n===========  $(mp_model_name)  -  R O B U S T    S O L U T I O N  ===========\n")
            println("SubOptimal Objective Function value: ", obj_star)
            println("Solve time : ", elapsed_time)
            println("Permutation: $(permutation_current)")
			println("Relative MP gap: $(mp_gap)")
            println("ERROR: [Master Problem] Optimal solution not found! Best bound: ", obj_star)
            mp_total_time += elapsed_time
            println("Model solve time: $(elapsed_time); time_limit = $(time_limit)")
            flush(stdout)
            if length(best_pi) > 0
                println("Time limit exceeded!")
				LB = obj_star
                gap = (UB - LB) / LB
                println("Model solve time: time_limit = $(time_limit); total_time = $(mp_total_time + sp_total_time); gap = $(gap); UB = $(UB); LB = $(LB)")
                flush(stdout)
				close(result_file)
                return UB, best_pi, [], mp_total_time + sp_total_time, false, iter, gap, LB, total_cuts, mp_total_time, sp_total_time
            end
            gap = (UB - obj_star) / obj_star
			close(result_file)
            return UB, best_pi, zeros(Float64, (n, n)), mp_total_time + sp_total_time, false, iter, gap, xlb, total_cuts, mp_total_time, sp_total_time
        end
        iter += 1
    end
end
