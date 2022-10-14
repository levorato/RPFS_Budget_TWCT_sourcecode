# ===================================================================================================
# robust_pfsp_budget_wct_separation_manne.jl
# Robust Permutation Flow Shop Problem to Minimize the Weighted Sum of Job Completion Times (WCT),
# based on the budgeted uncertainty set (Bertsimas et al.).
# Given all possible scenarios (i.e. for each job i, on which machines it will have its processing times
# deviated to their worst-case values), finds a permutation pi that minimizes the restricted worst-case
# WCT objective function.
# *** Solve the robust problem via separation approach (Column and Constraint Generation procedure),
# based on Manne PFSP MILP Model.
# ===================================================================================================

# Worst-case calculation for the robust budget problem
include("robust_pfsp_alt_worstcase_wct_mip.jl")
include("robust_pfsp_budget_worstcase_wct_brute.jl")
include("robust_pfsp_budget_separation_common.jl")

using JuMP
using CPLEX
using Gurobi
using Dates

# According to the separation method results,
# adds new constraints to the existing model.
# Return a tuple: (separation_done, mip_optimal)
# - separation_done indicates that the separation procedure has found the optimal master problem solution
# - mip_optimal indicates that the newly-generated MILP problem solution is optimal or not
#     ( if it is not optimal, the whole C&CG procedure has to stop and abort! )
# Attention: UB is an in-out function argument (its value if modified)
function separate_manne(model, m, n, Γ, P_bar, P_hat, w, iteration, UB, best_pi, bigM; max_cores = 16,
            single_cut_per_iteration = false, budget_type = "machine", time_limit = 7200, fractional_solution = false,
            use_brute_force = false, solver = "CPLEX", sp_bigM_type = 7, warm_start = false, permutation_current = [],
            sp_warm_start=false)
    EPS = 1e-8
    Jobs = 1:n
    N = Jobs
    Machines = 1:m
    M = Machines
    # Get existing model variables
    D = Matrix{VariableRef}(undef, n, n)
    for i in N
        for k in N
            D[i, k] = variable_by_name(model, "D[$(i),$(k)]")
        end
    end
    obj = variable_by_name(model, "obj")
    Q = bigM
    if warm_start
        LB = 0.0
        wct_current = 0.0
        is_opt = true
    else
        wct_current, permutation, D_current, is_opt = get_current_model_solution_dichotomous(model, n)
        permutation_current = [permutation]
        LB = wct_current
    end
    println("Separation for permutation_current = $(permutation_current)")
    if !is_opt  # Model solution failed (not optimal) : return to main procedure
        println("Model solution failed (not optimal) : return to main procedure")
        xgap = (UB - LB) / LB
        return true, is_opt, UB, LB, best_pi, xgap
    end
    #Solve Separation problem using dynamic programming
    status = true
    gap = (UB - LB) / LB
    sp_total_time = 0.0
    total_count = 0
    for permutation in Set(permutation_current)  # Process all initial solutions (or current CCG solution)
        #cmax_dp, best_Γ = worst_case_cmax_dp_m_machines(pi, m, n, Γ, P_bar, P_hat, w)
        println("Generating cuts for permutation = $(permutation)...")
        #cmax_dp, best_Γ = worst_case_cmax_dp_m_machines(pi, m, n, Γ, P_bar, P_hat, w)
        wct_sep, is_opt, z_sep_list, sp_time = solve_separation_problem(m, n, Γ, P_bar, P_hat, w, permutation, wct_current, true;
                budget_type=budget_type, time_limit=time_limit, single_cut_per_iteration=single_cut_per_iteration, max_cores=max_cores,
                fractional_solution=fractional_solution, use_brute_force=use_brute_force, solver=solver, bigM_type=sp_bigM_type,
                sp_warm_start=sp_warm_start)
        sp_total_time += sp_time
        if !is_opt
            gap = (UB - LB) / LB
            println("[Separation not optimal!] UB (wct_sep) = $(UB)  vs.  LB (wct_current) = $(LB) ; gap = $(gap)")
            flush(stdout)
            return status, false, UB, LB, best_pi, (UB - LB) / LB, sp_total_time, 0
        end
        #wc_validation, xbest_Γ, xworst_dev = worst_case_wct_brute_force_2_machines(pi, n, Γ[1], Γ[2], P_bar, P_hat, w, false)
        #@assert(abs(wc_validation - wct_sep) <= 1e-3)
        #   ==>  UB = min(UB, wct_sep)
        if wct_sep < UB
            UB = wct_sep
            best_pi = permutation
        end
        gap = (UB - LB) / LB
        println("UB (wct_sep) = $(UB)  vs.  LB (wct_current) = $(LB) ; gap = $(gap)")
        flush(stdout)
        count = 0
        if gap > EPS  # If Gap = (UB - LB) / LB <= EPS, an eps-optimal solution is found, terminate.
            status = false
            count = 0
            #### The following cut CANNOT be applied to the Manne variant of the PFSP model
            ### @constraint(model, sum(D[i, k] for i in 1:n for k in 1:n if D_current[i, k] == 1) <= n - 2)
            # Generate a set of cuts for each worst-case solution / scenario in z_sep_list
            for z_sep in z_sep_list
                count += 1
                # C_ik completion time of job k on machine i
                newC = @variable(model, [M, N], base_name = "newC_$(iteration)-$(count)", lower_bound = 0)

                # New objective function constraint
                @constraint(model, sum(w[i] * newC[m, i] for i in Jobs) <= obj)
                # New problem constraints for this specific scenario z_sep
                for i in 1:n
                    @constraint(model, newC[1, i] >= P_bar[1, i] + (P_hat[1, i] * z_sep[1, i]))
                end
                for r in 2:m
                    for i in 1:n
                        @constraint(model, newC[r, i] - newC[r-1, i] >= P_bar[r, i] + (P_hat[r, i] * z_sep[r, i]))
                    end
                end
                for r in 1:m
                    for i in 1:(n-1)
                        for k in (i+1):n
                            @constraint(model, newC[r, i] - newC[r, k] + Q * D[i, k] >= P_bar[r, i] + (P_hat[r, i] * z_sep[r, i]))
                            @constraint(model, newC[r, i] - newC[r, k] + Q * D[i, k] <= Q - (P_bar[r, k] + (P_hat[r, k] * z_sep[r, k])))
                        end
                    end
                end
            end
            total_count += count
        else
            println("Separation done.")
            flush(stdout)
            count = 0
            break
        end
    end
    return status, true, UB, LB, best_pi, gap, sp_total_time, total_count
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
# TODO Adapt to the Manne PFSP formulation
# Assignment / Completion-time formulation based was on:
#    Bektaş, T., Hamzadayı, A., & Ruiz, R. (2020). Benders decomposition for the mixed
#    no-idle permutation flowshop scheduling problem. Journal of Scheduling, 513–523.
#    https://doi.org/10.1007/s10951-020-00637-8
function solve_robust_model_flowshop_wct_budget_separation_manne(m, n, P_bar, P_hat, w, Γ, time_limit; mp_bigM_type = 1, sp_bigM_type = 3, max_cores = 16,
            solver = "CPLEX", budget_type = "machine", single_cut_per_iteration = false, fractional_solution = false, use_brute_force = false,
            warm_start = false, instance_filename = nothing, sp_solver = nothing, sp_warm_start = false)
    sp_warm_start = warm_start
    Jobs = 1:n
    N = Jobs
    Machines = 1:m
    M = Machines
    solver_time_limit = time_limit
    elapsed_time = 0.0
    sp_total_time = 0.0
    mp_total_time = 0.0
    # Setup robust model
    model = create_jump_model(solver)
    setup_solver_params(model, solver, solver_time_limit, max_cores)
    if isnothing(sp_solver)
        sp_solver = solver
    end
    println("Starting Manne C&CG procedure...")
    println("Procedure time limit = $(time_limit) s")
    println("Solver time limit = $(solver_time_limit) s")
    max_cmax = calculate_bigM_naive(m, n, P_bar, P_hat)
    worst_Γ_scenario = [ [x for x in 1:n] for r in 1:m ]
    # IMPORTANT: FIX BIGM TYPE TO 1, SINCE TYPE 2 IS NOT GIVING CORRECT RESULTS FOR SOME INSTANCES
    mp_bigM_type = 1
    println("MP BigM type = $(mp_bigM_type)")
    if mp_bigM_type == 2
        # BigM (Q) : Option 1. Calculate the Makespan (Cmax) for the worst-case scenario
        #                      (on all machines r, all jobs j have their processing times
        #                       deviated to their maximum value : P_bar(r,j) + P_hat(r,j)).
        random_pi = [x for x in 1:n]
        max_cmax = calculate_cmax_given_scenario(random_pi, m, n, P_bar, P_hat, worst_Γ_scenario) + 1
    end
    Q = max_cmax  # Big-M
    bigM = max_cmax
    ## ++++++++++ Model variables:
    # D_ik = 1, if job i is scheduled any time before job k, 0 otherwise
    @variable(model, D[N, N], Bin)
    # C_ri completion time of job i on machine r
    C = @variable(model, [M, N], base_name = "newC_1", lower_bound = 0)

    ## ++++++++++ Model constraints
    for i in 1:n
        @constraint(model, C[1, i] >= P_bar[1, i])
    end
    for r in 2:m
        for i in 1:n
            @constraint(model, C[r, i] - C[r-1, i] >= P_bar[r, i])
        end
    end
    for r in 1:m
        for i in 1:(n-1)
            for k in (i+1):n
                @constraint(model, C[r, i] - C[r, k] + Q * D[i, k] >= P_bar[r, i])
                @constraint(model, C[r, i] - C[r, k] + Q * D[i, k] <= Q - P_bar[r, k])
            end
        end
    end
    # ++++++++++ Objective function  (Constraint (4))
    @variable(model, obj >= 0)
    @objective(model, Min, obj)
    @constraint(model, sum(w[i] * C[m, i] for i in Jobs) <= obj)
    iter = 1
    UB = typemax(Float64)
    LB = 0.0
    best_pi = []
    total_cuts = 0
    while true
        # Solve the problem
        println("====================================================================\n")
        println("Solving WCT budget robust Manne model C & CG (Iteration $(iter))...")
        println("====================================================================\n")
        datetimenow = Dates.now()
		println("Date: $(datetimenow)")
		flush(stdout)
        if warm_start
            status = MOI.OPTIMAL
        else
            time_left = trunc(Int, max(100, solver_time_limit - (mp_total_time + sp_total_time)))
            println("time_left = $(time_left)")
            println("UB = $(UB)")
            flush(stdout)
            setup_solver_params(model, solver, time_left, max_cores, UB)
            if length(best_pi) > 0
                set_warm_start_dichotomous_model(D, n, best_pi)
            end
            status = optimize!(model)
            status = JuMP.termination_status(model)
        end
        D_star = zeros(Float64, (n, n))
        if status == MOI.OPTIMAL
            permutation_current = []
            if warm_start
                permutation_current, elapsed_time = generate_initial_solutions(instance_filename, m, n, Γ, budget_type)
                best_pi = permutation_current[1]
                y_star = 0.0
            else
                println("\n===========  M A N N E    R O B U S T    S O L U T I O N  ===========\n")
                println("Optimal Objective Function value: ", objective_value(model))
                println("Solve time : ", solve_time(model))
                flush(stdout)
                elapsed_time = solve_time(model)
                for i in Jobs  # Job i
                    for k in Jobs  # Job k
                        if value(D[i, k]) >= 0.9
                            D_star[i, k] = 1
                        end
                    end
                end
                permutation_current = [topological_sort_by_dfs(D_star, n)]
                y_star = objective_value(model)
            end
            println("Permutation: $(permutation_current)")
            flush(stdout)
            mp_total_time += elapsed_time
            # IMPORTANT : Check for time_limit
            if (mp_total_time + sp_total_time >= time_limit)
                println("Time limit exceeded!")
                LB = y_star
                gap = (UB - LB) / LB
                println("Model solve time: time_limit = $(time_limit); total_time = $(mp_total_time + sp_total_time); gap = $(gap); UB = $(UB); LB = $(LB)")
                flush(stdout)
                return UB, best_pi, D_star, mp_total_time + sp_total_time, false, iter, gap, LB, total_cuts, mp_total_time, sp_total_time
            end

            separation_done, mip_optimal, UB, LB, best_pi, gap, sp_time, cut_count = separate_manne(model, m, n, Γ, P_bar, P_hat, w, iter + 1, UB, best_pi, bigM;
                    budget_type=budget_type, time_limit=trunc(Int, max(MIN_TIME_SP_MIP, time_limit - (mp_total_time + sp_total_time))),
                    single_cut_per_iteration=single_cut_per_iteration, max_cores=max_cores,
                    fractional_solution=fractional_solution, use_brute_force=use_brute_force, solver=sp_solver, sp_bigM_type=sp_bigM_type,
                    warm_start=warm_start, permutation_current=permutation_current, sp_warm_start=sp_warm_start)
            sp_total_time += sp_time
            total_cuts += cut_count
            if warm_start  # disable warm start after first CCG iteration
                warm_start = false
            end

            # IMPORTANT : Check for time_limit
            if (mp_total_time + sp_total_time >= time_limit)
                println("Time limit exceeded!")
                println("Model solve time: time_limit = $(time_limit); total_time = $(mp_total_time + sp_total_time); gap = $(gap); UB = $(UB); LB = $(LB)")
                flush(stdout)
                return UB, best_pi, D_star, mp_total_time + sp_total_time, false, iter, gap, LB, total_cuts, mp_total_time, sp_total_time
            end
            if (!mip_optimal)
                println("ERROR: [Sub Problem] Optimal solution not found! Best UB: $(UB) ; LB = $(LB)")
                println("Model solve time: time_limit = $(time_limit); elapsed_time = $(elapsed_time); gap = $(gap); UB = $(UB); LB = $(LB)")
                flush(stdout)
                return UB, best_pi, D_star, mp_total_time + sp_total_time, false, iter, gap, LB, total_cuts, mp_total_time, sp_total_time
            end
            if separation_done
                println(" *** Separation done! y_star = $(y_star), best_pi = $(best_pi)")
                flush(stdout)
                return UB, best_pi, D_star, solve_time(model), true, iter, gap, LB, total_cuts, mp_total_time, sp_total_time
            end  # end if separate
        else
            println("ERROR: [Master Problem] Optimal solution not found! Best bound: ", objective_bound(model))
            solvetime = solve_time(model)
            mp_total_time += solvetime
            println("Model solve time: $(solvetime); time_limit = $(time_limit)")
            flush(stdout)
            if length(best_pi) > 0
                println("Time limit exceeded!")
                gap = (UB - LB) / LB
                println("Model solve time: time_limit = $(time_limit); total_time = $(mp_total_time + sp_total_time); gap = $(gap); UB = $(UB); LB = $(LB)")
                flush(stdout)
                return UB, best_pi, D_star, mp_total_time + sp_total_time, false, iter, gap, LB, total_cuts, mp_total_time, sp_total_time
            end
            xlb = objective_bound(model)
            gap = (UB - xlb) / xlb
            return UB, best_pi, zeros(Float64, (n, n)), mp_total_time + sp_total_time, false, iter, gap, xlb, total_cuts, mp_total_time, sp_total_time
        end
        iter += 1
    end
end
