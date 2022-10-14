# ===================================================================================================
# robust_pfsp_budget_wct_separation_wagner_wst2.jl
# Robust Permutation Flow Shop Problem to Minimize the Weighted Sum of Job Completion Times (WCT),
# based on the budgeted uncertainty set (Bertsimas et al.).
# Given all possible scenarios (i.e. for each job i, on which machines it will have its processing times
# deviated to their worst-case values), finds a permutation pi that minimizes the restricted worst-case
# WCT objective function.
# Reference for Wagner WST2 PFSP Model:
#   Tseng, F. T., & Stafford, E. F. (2008). New MILP models for the permutation flowshop problem.
#   Journal of the Operational Research Society, 59(10), 1373–1386.
#   https://doi.org/10.1057/palgrave.jors.2602455
# *** Solve the robust problem via separation approach (Column and Constraint Generation procedure),
# based on Wagner WST2 PFSP MILP Model.
# ===================================================================================================

# Worst-case calculation for the robust budget problem
include("robust_pfsp_budget_separation_common.jl")
include("robust_pfsp_alt_worstcase_wct_mip.jl")
include("robust_pfsp_budget_worstcase_wct_brute.jl")

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
function separate_wagner_wst2(model, m, n, Γ, P_bar, P_hat, w, iteration, UB, best_pi, bigM; max_cores = 16,
            single_cut_per_iteration = false, budget_type = "machine", time_limit = 7200, fractional_solution = false,
            use_brute_force = false, solver = "CPLEX", sp_bigM_type = 7, warm_start = false, permutation_current = [])
    EPS = 1e-8
    Jobs = 1:n
    N = Jobs
    Machines = 1:m
    M = Machines
    # Get existing model variables
    Z = Matrix{VariableRef}(undef, n, n)
    for i in N
        for j in N
            Z[i, j] = variable_by_name(model, "Z[$(i),$(j)]")
        end
    end
    obj = variable_by_name(model, "obj")
    if warm_start
        LB = 0.0
        wct_current = 0.0
        is_opt = true
    else
        wct_current, permutation, Z_current, is_opt = get_current_model_solution_assignment(model, n)
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
                fractional_solution=fractional_solution, use_brute_force=use_brute_force, solver=solver, bigM_type=sp_bigM_type)
        sp_total_time += sp_time
        if !is_opt
            gap = (UB - LB) / LB
            println("[Separation not optimal!] UB (wct_sep) = $(UB)  vs.  LB (wct_current) = $(LB) ; gap = $(gap)")
            flush(stdout)
            return status, false, UB, LB, best_pi, (UB - LB) / LB, sp_total_time, 0
        end
        #global solution = getobjectivevalue(model)
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
            # Additional cut
            # @constraint(model, sum(Z[i, k] for i in 1:n for k in 1:n if Z_current[i, k] == 1) <= n - 2)
            # Generate a set of cuts for each worst-case solution / scenario in z_sep_list
            for z_sep in z_sep_list
                count += 1
                # X_rj : idle time on machine r before the start of job in sequence position j
                newX = @variable(model, [M, N], base_name = "newX_$(iteration)-$(count)", lower_bound = 0)
                # Y_rj idle time of job in sequence position j after it finishes processing on machine r
                newY = @variable(model, [M, N], base_name = "newY_$(iteration)-$(count)", lower_bound = 0)
                # F[i] : completion time of job i on the last machine (M_m)
                newF = @variable(model, [N], base_name = "newF_$(iteration)-$(count)", lower_bound = 0)

                # New objective function constraint
                @constraint(model, sum(w[i] * newF[i] for i in Jobs) <= obj)

                # New problem constraints for this specific scenario z_sep
                for r in 2:m-1  # (A.1)
                    for j in 2:n-1
                        @constraint(model, sum((P_bar[r, i] + P_hat[r, i] * z_sep[r, i]) * Z[i, j+1] for i in Jobs)
                            - sum((P_bar[r+1, i] + P_hat[r+1, i] * z_sep[r+1, i]) * Z[i, j] for i in Jobs)
                            + newX[r, j+1] - newX[r+1, j+1] + newY[r, j+1] - newY[r, j] == 0)
                    end
                end
                for j in 2:n-1  # (A.2)
                    @constraint(model, sum((P_bar[1, i] + P_hat[1, i] * z_sep[1, i]) * Z[i, j+1] for i in Jobs)
                        - sum((P_bar[2, i] + P_hat[2, i] * z_sep[2, i]) * Z[i, j] for i in Jobs) - newX[2, j+1] + newY[1, j+1] - newY[1, j] == 0)
                end
                for r in 2:m-1  # (A.3)
                    @constraint(model, sum((P_bar[r, i] + P_hat[r, i] * z_sep[r, i]) * Z[i, 2] for i in Jobs)
                        - sum((P_bar[r+1, i] + P_hat[r+1, i] * z_sep[r+1, i]) * Z[i, 1] for i in Jobs) + newX[r, 2] - newX[r+1, 2] + newY[r, 2] == 0)
                end
                # (A.4)
                @constraint(model, sum((P_bar[1, i] + P_hat[1, i] * z_sep[1, i]) * Z[i, 2] for i in Jobs)
                    - sum((P_bar[2, i] + P_hat[2, i] * z_sep[2, i]) * Z[i, 1] for i in Jobs) + newY[1, 2] - newX[2, 2] == 0)
                for r in 1:m-1  # (A.5)
                    @constraint(model, sum((P_bar[r, i] + P_hat[r, i] * z_sep[r, i]) * Z[i, 1] for i in Jobs) + newX[r, 1] - newX[r+1, 1] == 0)
                end
                @constraint(model, newX[1, 1] == 0)
                for i in 1:n  # Disjuctive constraints used to determine the completion time of job i on the last machine m
                    for j in 1:n
                        @constraint(model, newF[i] >= sum((P_bar[m, x] + P_hat[m, x] * z_sep[m, x]) * Z[x, p] for p in 1:j, x in 1:n)
                            + sum(newX[m, p] for p in 1:j) - bigM * (1 - Z[i, j]) )
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
function solve_robust_model_flowshop_wct_budget_separation_wagner_wst2(m, n, P_bar, P_hat, w, Γ, time_limit; mp_bigM_type = 2, sp_bigM_type = 3, max_cores = 16,
            solver = "CPLEX", budget_type = "machine", single_cut_per_iteration = false, fractional_solution = false, use_brute_force = false,
            warm_start = false, instance_filename = nothing, sp_solver = nothing)
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
    if isnothing(sp_solver)
        sp_solver = solver
    end
    setup_solver_params(model, solver, solver_time_limit, max_cores)
    max_cmax = calculate_bigM_naive(m, n, P_bar, P_hat)
    worst_Γ_scenario = [ [x for x in 1:n] for r in 1:m ]
    if mp_bigM_type == 2
        # BigM (Q) : Option 1. Calculate the Makespan (Cmax) for the worst-case scenario
        #                      (on all machines r, all jobs j have their processing times
        #                       deviated to their maximum value : P_bar(r,j) + P_hat(r,j)).
        random_pi = [x for x in 1:n]
        max_cmax = calculate_cmax_given_scenario(random_pi, m, n, P_bar, P_hat, worst_Γ_scenario) + 1
    end
    bigM = max_cmax  # Big-M
    println("Starting Wagner WST2 C&CG procedure...")
    println("Procedure time limit = $(time_limit) s")
    println("Solver time limit = $(solver_time_limit) s")
    ## ++++++++++ Model variables:
    # Z_ij = 1, if job i is assigned to sequence position j, 0 otherwise
    @variable(model, Z[N, N], Bin)

    # X_rj : idle time on machine r before the start of job in sequence position j
    X = @variable(model, [M, N], base_name = "newX_1", lower_bound = 0)
    # Y_rj idle time of job in sequence position j after it finishes processing on machine r
    Y = @variable(model, [M, N], base_name = "newY_1", lower_bound = 0)
    # F[i] : end time of job i on the last machine (M_m)
    F = @variable(model, [N], base_name = "newF_1", lower_bound = 0)
    # Delta: auxiliary variable used in the BigM-based disjuction (see last constraint)
    ###Delta = @variable(model, [N, N], base_name = "newDelta_1", Bin)

    ## ++++++++++ Model constraints
    T = P_bar
    # Assignement constraints
    for i in Jobs
        @constraint(model, sum(Z[i, j] for j in Jobs) == 1)
    end
    for j in Jobs
        @constraint(model, sum(Z[i, j] for i in Jobs) == 1)
    end
    for r in 2:m-1  # (A.1)
        for j in 2:n-1
            @constraint(model, sum(T[r, i] * Z[i, j+1] for i in Jobs)
                - sum(T[r+1, i] * Z[i, j] for i in Jobs) + X[r, j+1] - X[r+1, j+1] + Y[r, j+1] - Y[r, j] == 0)
        end
    end
    for j in 2:n-1  # (A.2)
        @constraint(model, sum(T[1, i] * Z[i, j+1] for i in Jobs) - sum(T[2, i] * Z[i, j] for i in Jobs) - X[2, j+1] + Y[1, j+1] - Y[1, j] == 0)
    end
    for r in 2:m-1  # (A.3)
        @constraint(model, sum(T[r, i] * Z[i, 2] for i in Jobs) - sum(T[r+1, i] * Z[i, 1] for i in Jobs) + X[r, 2] - X[r+1, 2] + Y[r, 2] == 0)
    end
    # (A.4)
    @constraint(model, sum(T[1, i] * Z[i, 2] for i in Jobs) - sum(T[2, i] * Z[i, 1] for i in Jobs) + Y[1, 2] - X[2, 2] == 0)
    for r in 1:m-1  # (A.5)
        @constraint(model, sum(T[r, i] * Z[i, 1] for i in Jobs) + X[r, 1] - X[r+1, 1] == 0)
    end
    @constraint(model, X[1, 1] == 0)
    # (A.6)
    ### @constraint(model, Cmax == (sum(T[m, i] for i in Jobs) + sum(X[m, p] for p in Jobs)))
    # The Delta[i, j] is needed because there will be several scheduling scenarios, each one associated with
    # different permutations and distinct start times B. So, we cannot simply use the first-stage variable Z
    # in this case.
    for i in 1:n  # Disjuctive constraints used to determine the completion time of job i on the last machine m
        for j in 1:n
            @constraint(model, F[i] >= sum(T[m, x] * Z[x, p] for p in 1:j, x in 1:n) + sum(X[m, p] for p in 1:j) - bigM * (1 - Z[i, j]) )
        end
    end

    # ++++++++++ Objective function  (Constraint (4))
    @variable(model, obj >= 0)  # Makespan
    @objective(model, Min, obj)
    @constraint(model, sum(w[i] * F[i] for i in Jobs) <= obj)
    iter = 1
    UB = typemax(Float64)
    best_pi = []
    LB = 0.0
    total_cuts = 0
    while true
        # Solve the problem
        println("=========================================================================\n")
        println("Solving WCT budget robust Wagner WST2 model C & CG (Iteration $(iter))...")
        println("=========================================================================\n")
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
                set_warm_start_assignment_model(Z, n, best_pi)
            end
            status = optimize!(model)
            status = JuMP.termination_status(model)
        end
        Z_star = zeros(Float64, (n, n))
        if status == MOI.OPTIMAL
            permutation_current = []
            if warm_start
                permutation_current, elapsed_time = generate_initial_solutions(instance_filename, m, n, Γ, budget_type)
                best_pi = permutation_current[1]
                y_star = 0.0
            else
                println("\n===========  W A G N E R    W S T 2    R O B U S T    S O L U T I O N  ===========\n")
                println("Optimal Objective Function value: ", objective_value(model))
                println("Solve time : ", solve_time(model))
                flush(stdout)
                elapsed_time = solve_time(model)
                permutation, Z_star = generate_permutation_from_Z_matrix(Z, n)
                permutation_current = [permutation]
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
                return UB, best_pi, Z_star, mp_total_time + sp_total_time, false, iter, gap, LB, total_cuts, mp_total_time, sp_total_time
            end
            separation_done, mip_optimal, UB, LB, best_pi, gap, sp_time, cut_count = separate_wagner_wst2(model, m, n, Γ, P_bar, P_hat, w, iter + 1, UB, best_pi, bigM;
                    budget_type=budget_type, time_limit=trunc(Int, max(MIN_TIME_SP_MIP, time_limit - (mp_total_time + sp_total_time))),
                    single_cut_per_iteration=single_cut_per_iteration, max_cores=max_cores,
                    fractional_solution=fractional_solution, use_brute_force=use_brute_force, solver=sp_solver, sp_bigM_type=sp_bigM_type,
                    warm_start=warm_start, permutation_current=permutation_current)
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
                return UB, best_pi, Z_star, mp_total_time + sp_total_time, false, iter, gap, LB, total_cuts, mp_total_time, sp_total_time
            end
            if (!mip_optimal)
                println("ERROR: [Sub Problem] Optimal solution not found! Best UB: $(UB) ; LB = $(LB)")
                println("Model solve time: time_limit = $(time_limit); elapsed_time = $(elapsed_time); gap = $(gap); UB = $(UB); LB = $(LB)")
                flush(stdout)
                return UB, best_pi, Z_star, mp_total_time + sp_total_time, false, iter, gap, LB, total_cuts, mp_total_time, sp_total_time
            end
            if separation_done
                println(" *** Separation done! y_star = $(y_star), best_pi = $(best_pi)")
                flush(stdout)
                return UB, best_pi, Z_star, solve_time(model), true, iter, gap, LB, total_cuts, mp_total_time, sp_total_time
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
                return UB, best_pi, Z_star, mp_total_time + sp_total_time, false, iter, gap, LB, total_cuts, mp_total_time, sp_total_time
            end
            xlb = objective_bound(model)
            gap = (UB - xlb) / xlb
            return UB, best_pi, zeros(Float64, (n, n)), mp_total_time + sp_total_time, false, iter, gap, xlb, total_cuts, mp_total_time, sp_total_time
        end
        iter += 1
    end
end
