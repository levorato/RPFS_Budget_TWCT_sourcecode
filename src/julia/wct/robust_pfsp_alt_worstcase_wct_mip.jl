# ===================================================================================================
# robust_pfsp_alt_worstcase_wct_mip.jl
# Worst-case determination procedure for the Robust Permutation Flow Shop Problem to Minimize the
# Weighted Sum of Job Completion Times (WCT), based on the budgeted uncertainty set (Bertsimas et al.).
# Given a permutation pi, finds the restricted worst-case scenario (i.e. for each job i,
# on which machines it will have its processing times deviated to their worst-case values), so that the
# objective function (weighted completion time - wct) reaches the worst-case (maximum) value.
# *** MIP approach for calculating the worst-case scenario.
# ===================================================================================================

using JuMP
using CPLEX
using Gurobi

#include("../config.jl")
include("../deterministic_pfsp.jl")
include("../pfsp_file_reader.jl")
include("robust_pfsp_budget_worstcase_wct_brute.jl")
include("deterministic_pfsp_alt_wct_mip.jl")

function test_model()
    # basedir = joinpath(robust_ying_instances_folder, "rob-pfsp-wct/20jobs")  filename == "RB0205010_20_2_50_wct_inputs.txt"
    basedir = joinpath(robust_ying_instances_folder, "rob-pfsp-wct/10x3")
    time_limit = 1800
    time_spent = [0.0, 0.0, 0.0]
    count = 0
    for bigM_type in [3, 2, 1]
        println("BigM type: $(bigM_type)\n====================================")
        for filename in searchdir(basedir, ".txt")
            if filename == "RB0101001_010_003_10_wct_inputs.txt"
                println("Instance: $(filename)")
                inputfile = "$(basedir)/$(filename)"
                m, n, P_bar, P_hat, w = read_robust_input_file(inputfile)
                # φ_star is the optimal Johnson solution obtained over the nominal processing time array (P_bar)
                # w = [x for x in 1:n]
                if n < 10
                    wct_star, permutation, x = solve_deterministic_pfsp_wct(m, n, P_bar, w)
                    println("Based on the nominal processing times, the optimal WCT wct_star is $(wct_star).")
                    println("x: $(x)")
                    wct_validation = calculate_wct(m, n, P_bar, w, permutation, false)
                    println("WCT Validation = $(wct_validation)")
                    @assert abs(wct_star - wct_validation) < 1e-5
                else
                    permutation = [x for x in 1:n]  # [9, 8, 7, 6, 10, 5, 4, 3, 2, 1] #
                end
                println("pi: $(permutation)")
                println("w: $(w)")
                wct_validation = calculate_wct(m, n, P_bar, w, permutation, false)
                println("WCT Validation = $(wct_validation)")
                Jobs = 1:n
                N = Jobs
                Machines = 1:m
                M = Machines
                println("Calculating robust PFSP budget solution (Alternative model)...")
                for Γ1 in 1:n
                    for Γ2 in 1:n
                        for Γ3 in 1:n
                        #Γ3 = 5
                            cost, scenarios, elapsed_time = solve_robust_model_worst_case_alt(permutation, m, n, [Γ1, Γ2, Γ3], P_bar, P_hat, w, 1;
                                solver="CPLEX", bigM_type=bigM_type, max_cores=2)
                            time_spent[bigM_type] += elapsed_time
                            count += 1
                            scenario = scenarios[1]
                            deviated_jobs = [[] for x in 1:m]
                            for r in 1:m
                                for k in 1:n
                                    if scenario[r, k] >= 0.9
                                        push!(deviated_jobs[r], k)
                                    end
                                end
                            end
                            println("Worst-case scenario: $(deviated_jobs)")
                            println("The worst-case cost from MILP is $(cost).")
                        end
                    end
                end
                #wct = calculate_wct_given_scenario(pi, m, n, P_bar, P_hat, w, scenarios[1], true )     # true)
                #println("[Validation MILP] Given scenario = $(deviated_jobs), the WCT is $(wct).")
                #if n <= 10
                #    wct_brute, worst_Γ_scenario_brute = worst_case_wct_brute_force_2_machines(pi, n, Γ1, Γ2, P_bar, P_hat, w, true)
                #    @assert abs(wct_brute - cost) < 1e-5
                #else
                #    println("Skipping worst-case validation, since n > 10.")
                #end
            end
        end
    end
    println("Total count is $(count)")
    for bigM_type in [3, 2, 1]
        println("bigM_type = $(bigM_type) => Total time spent: $(time_spent[bigM_type])")
    end
end

function test_model_2()
    # basedir = joinpath(robust_ying_instances_folder, "rob-pfsp-wct/20jobs")  filename == "RB0205010_20_2_50_wct_inputs.txt"
    basedir = joinpath(robust_ying_instances_folder, "rob-pfsp-wct", "10x5", "50%")
    time_limit = 4000
    time_spent = [0.0 for x in 1:10]
    count = 0
    for bigM_type in [7]
        println("BigM type: $(bigM_type)\n====================================")
        for filename in searchdir(basedir, ".txt")
            if filename == "RB0105007_10_5_50_wct_inputs.txt"
                println("Instance: $(filename)")
                inputfile = joinpath("$(basedir)", "$(filename)")
                m, n, P_bar, P_hat, w = read_robust_input_file(inputfile)
                # φ_star is the optimal Johnson solution obtained over the nominal processing time array (P_bar)
                # w = [x for x in 1:n]
                if n < 10
                    wct_star, permutation, x = solve_deterministic_pfsp_wct(m, n, P_bar, w)
                    println("Based on the nominal processing times, the optimal WCT wct_star is $(wct_star).")
                    println("x: $(x)")
                    wct_validation = calculate_wct(m, n, P_bar, w, permutation, false)
                    println("WCT Validation = $(wct_validation)")
                    @assert abs(wct_star - wct_validation) < 1e-5
                else
                    permutation = [ 9, 10, 5, 6, 2, 4, 1, 3, 8, 7 ]
                end
                println("pi: $(permutation)")
                println("w: $(w)")
                wct_validation = calculate_wct(m, n, P_bar, w, permutation, false)
                println("WCT Validation = $(wct_validation)")
                Jobs = 1:n
                N = Jobs
                Machines = 1:m
                M = Machines
                println("Calculating robust PFSP budget solution (Alternative model)...")
                Γ = 2
                cost, scenarios, elapsed_time = solve_robust_model_worst_case_alt(permutation, m, n, Γ, P_bar, P_hat, w, 1;
                    solver="CPLEX", bigM_type=bigM_type, max_cores=2, budget_type="global")
                time_spent[bigM_type] += elapsed_time
                count += 1
                scenario = scenarios[1]
                deviated_jobs = [[] for x in 1:m]
                for r in 1:m
                    for k in 1:n
                        if scenario[r, k] >= 0.9
                            push!(deviated_jobs[r], k)
                        end
                    end
                end
                println("Worst-case scenario: $(deviated_jobs)")
                println("The worst-case cost from MILP is $(cost).")

                wct = calculate_wct_given_scenario(permutation, m, n, P_bar, P_hat, w, scenarios[1], true )     # true)
                println("[Validation MILP] Given scenario = $(deviated_jobs), the WCT is $(wct).")
                #wct = calculate_wct_given_scenario(pi, m, n, P_bar, P_hat, w, scenarios[1], true )     # true)
                #println("[Validation MILP] Given scenario = $(deviated_jobs), the WCT is $(wct).")
                #if n <= 10
                #    wct_brute, worst_Γ_scenario_brute = worst_case_wct_brute_force_2_machines(pi, n, Γ1, Γ2, P_bar, P_hat, w, true)
                #    @assert abs(wct_brute - cost) < 1e-5
                #else
                #    println("Skipping worst-case validation, since n > 10.")
                #end
            end
        end
    end
    println("Total count is $(count)")
    for bigM_type in [3, 2, 1]
        println("bigM_type = $(bigM_type) => Total time spent: $(time_spent[bigM_type])")
    end
end

function calculate_bigM_naive(m, n, P_bar, P_hat)
    #bigM = 1
    op_list = []
    for r in 1:m
        for i in 1:n
            # bigM += P_bar[r, i] + P_hat[r, i]
            push!(op_list, P_bar[r, i] + P_hat[r, i])
        end
    end
    # Get the (m + n - 1) largest elements from op_list
    largest_indices = partialsortperm(op_list, 1:(m + n - 1), rev=true)
    # Sum the (m + n - 1) largest elements from op_list
    bigM = 0
    for idx in largest_indices
        bigM += op_list[idx]
    end
    return bigM + 1
end

# Assignment / Completion-time formulation based was on:
# 1. TS2 MILP Model by Tseng and Stafford (2001).
# 2. Bektaş, T., Hamzadayı, A., & Ruiz, R. (2020). Benders decomposition for the mixed
#    no-idle permutation flowshop scheduling problem. Journal of Scheduling, 513–523.
#    https://doi.org/10.1007/s10951-020-00637-8
function solve_robust_model_worst_case_alt(permutation, m, n, Γ, P_bar, P_hat, w, max_solutions = 5; bigM_type = 3, max_cores = 8,
            solver = "CPLEX", budget_type = "machine", solver_time_limit = 7200, single_cut_per_iteration = false,
            fractional_solution = false, initial_scenario = nothing)
    Jobs = 1:n
    N = Jobs
    Machines = 1:m
    M = Machines
    elapsed_time = 0.0
    # Setup robust model
    if solver == "Gurobi"
        model = Model(Gurobi.Optimizer)  # MIPGap=0.005
        MOI.set(model, MOI.RawOptimizerAttribute("TimeLimit"), solver_time_limit)
        MOI.set(model, MOI.RawOptimizerAttribute("Threads"), max_cores)
        # https://support.gurobi.com/hc/en-us/articles/360013195772-What-can-be-done-to-avoid-an-out-of-memory-condition-
        MOI.set(model, MOI.RawOptimizerAttribute("NodefileStart"), 0.5)
        #MOI.set(model, MOI.RawOptimizerAttribute("Method"), 1)  # use dual simplex to save memory
        #MOI.set(model, MOI.RawOptimizerAttribute("IntFeasTol"), 1e-3)
        # https://www.gurobi.com/documentation/9.1/refman/mipgapabs.html
        MOI.set(model, MOI.RawOptimizerAttribute("MIPGapAbs"), 1e-5)
        # https://www.gurobi.com/documentation/9.1/refman/integralityfocus.html
        #MOI.set(model, MOI.RawOptimizerAttribute("IntegralityFocus"), 1)
    elseif solver == "CPLEX"
        model = Model(CPLEX.Optimizer)  # MIPGap=0.005
        MOI.set(model, MOI.RawOptimizerAttribute("CPX_PARAM_TILIM"), solver_time_limit)
        MOI.set(model, MOI.RawOptimizerAttribute("CPX_PARAM_THREADS"), max_cores)
    else
      println("No solver defined")
    end
    num_constraints = 0
    Z = zeros(Float64, (n, n))
    j = 1
    for i in 1:n
        Z[permutation[i], j] = 1.0
        j += 1
    end
    println("z = $(Z)")
    bigM = zeros(Float64, (m, n))
    fill!(bigM, calculate_bigM_naive(m, n, P_bar, P_hat))
    worst_Γ_scenario = [ [x for x in 1:n] for r in 1:m ]
    best_Γ_scenario = [ [] for r in 1:m ]
    max_C = calculate_cmax_given_scenario(permutation, m, n, P_bar, P_hat, worst_Γ_scenario, true)
    min_C = calculate_cmax_given_scenario(permutation, m, n, P_bar, P_hat, best_Γ_scenario, true)
    max_cmax = calculate_cmax_given_scenario(permutation, m, n, P_bar, P_hat, worst_Γ_scenario) + 1
    # https://spartanideas.msu.edu/2018/09/17/choosing-big-m-values/
    if bigM_type == 2
        # BigM (Q) : Option 1. Calculate the Makespan (Cmax) for the worst-case scenario
        #                      (on all machines r, all jobs j have their processing times
        #                       deviated to their maximum value : P_bar(r,j) + P_hat(r,j)).
        fill!(bigM, max_cmax)
    elseif bigM_type == 3
        # BigM (Q) : Option 2. Given all machines r: calculate the maximum difference
        #            between job completion times (C_{r, j}) of a job j and any other
        #            job k, on the same machine r.
        max_diff = 0
        for r in 1:m
            for j in 1:n
                for k in 1:n
                    if j != k
                        max_diff = max(max_diff, abs(max_C[r, j] - max_C[r, k]))
                    end
                end
            end
        end
        fill!(bigM, 2*max_diff + 1)
    elseif bigM_type == 4
        fill!(bigM, max_cmax)
        for i in 2:m  # [i, k-1] x [i-1, k]
            for k in 2:n
                bigM[i, k] = max_C[i, k]  # max(abs(max_C[i-1, k] - min_C[i, k-1]), abs(max_C[i, k-1] - min_C[i-1, k])) + 1
            end
        end
    elseif bigM_type == 5
        fill!(bigM, max_cmax)
        for i in 2:m  # [i, k-1] x [i-1, k]
            for k in 2:n
                if (i == 2) || (k == 2)
                    bigM[i, k] = max_C[i, k]
                else
                    bigM[i, k] = 2 * (max_C[i, k] - min(min_C[i-1, k-2], min_C[i-2, k-1])) + 1  # max(abs(max_C[i-1, k] - min_C[i, k-1]), abs(max_C[i, k-1] - min_C[i-1, k])) + 1
                end
            end
        end
    elseif bigM_type == 6
        fill!(bigM, max_cmax)
        for i in 2:m  # [i, k-1] x [i-1, k]
            for k in 2:n
                #bigM[i, k] = 2 * (max_C[i, k] - min(min_C[i, k-1] - (P_bar[i, k-1] + P_hat[i, k-1]), min_C[i-1, k] - (P_bar[i-1, k] + P_hat[i-1, k]))) + 1
                #bigM[i, k] = max_C[i, k] - min(min_C[i, k-1] - (P_bar[i, k-1]), min_C[i-1, k] - (P_bar[i-1, k])) + 1
                bigM[i, k] = 2 * (max_C[i, k] - min(min_C[i, k-1], min_C[i-1, k])) + 1
            end
        end
    elseif bigM_type == 7
        fill!(bigM, max_cmax)
        for i in 2:m  # [i, k-1] x [i-1, k]
            for k in 2:n
                bigM[i, k] = 2 * (max(max_C[i, k-1] - min_C[i-1, k], max_C[i-1, k] - min_C[i, k-1])) + 1
            end
        end
    elseif bigM_type != 1
        println("ERROR: bigM type not found!")
        return
    end
    ## ++++++++++ Model variables: C (Z is fixed parameter)
    # C_ik completion time of job k on machine i
    @variable(model, abs_dif_d_e[M, N] >= 0)
    @variable(model, C[M, N] >= 0)
    @variable(model, min_d_e[M, N] >= 0)
    @variable(model, disj[M, N], Bin)
    #@variable(model, avg_d_e[M, N] >= 0)
    if fractional_solution  # Allow fractional solutions for dev variable ?
        @variable(model, dev[M, N] >= 0)
        for r in 1:m
            for i in 1:n
                @constraint(model, dev[r, i] <= 1)
            end
        end
    else
        @variable(model, dev[M, N], Bin)
    end
    # +++++++++++ Warm start parameter
    if !isnothing(initial_scenario)
        println("[INFO] Warmstart: using initial_scenario provided.")
        for r in 1:m
            for i in 1:n
                set_start_value(dev[r, i], initial_scenario[r, i])
            end
        end
    end
    ## ++++++++++ Model constraints
    T = P_bar
    # Caso especial (1): quando i = 1 (maquina 1)
    @constraint(model, C[1, 1] == sum(T[1, j] * Z[j, 1] for j in Jobs) + sum(P_hat[1, j] * dev[1, j] * Z[j, 1] for j in Jobs))
    for k in 2:n
        @constraint(model, C[1, k] == C[1, k - 1] + sum(T[1, j] * Z[j, k] for j in Jobs) + sum(P_hat[1, j] * dev[1, j] * Z[j, k] for j in Jobs))
    end
    # Caso especial (2): quando k == 1 (primeira posicao de cada maquina)
    for i in 2:m
        @constraint(model, C[i, 1] == C[i - 1, 1] + sum(T[i, j] * Z[j, 1] for j in Jobs) + sum(P_hat[i, j] * dev[i, j] * Z[j, 1] for j in Jobs))
    end
    for k in 2:n
        # Caso geral (k > 2, i > 2)
        for i in 2:m
            @constraint(model, min_d_e[i, k] <= C[i, k - 1])
            @constraint(model, min_d_e[i, k] <= C[i - 1, k])

            # https://math.stackexchange.com/questions/1774275/mixed-integer-linear-programming-absolute-value-of-a-variable-not-involved-n-t
            @constraint(model, abs_dif_d_e[i, k] >= (C[i, k - 1] - C[i - 1, k]))
            @constraint(model, abs_dif_d_e[i, k] >= -(C[i, k - 1] - C[i - 1, k]))

            @constraint(model, abs_dif_d_e[i, k] <= (C[i, k - 1] - C[i - 1, k]) + bigM[i, k] * disj[i, k])
            @constraint(model, abs_dif_d_e[i, k] <= -(C[i, k - 1] - C[i - 1, k]) + bigM[i, k] * (1 - disj[i, k]))

            # Strengthen option 1:
            #@constraint(model, avg_d_e[i, k] == (C[i, k - 1] + C[i - 1, k]) / 2)
            #@constraint(model, abs_dif_d_e[i, k] == 2 * (avg_d_e[i, k] - min_d_e[i, k]))
            #@constraint(model, C[i, k] <= min_d_e[i, k] + abs_dif_d_e[i, k] + sum(T[i, j] * Z[j, k] for j in Jobs) + sum(P_hat[i, j] * dev[i, j] * Z[j, k] for j in Jobs))
            # Strengthen option 2:
            @constraint(model, abs_dif_d_e[i, k] == C[i, k - 1] + C[i - 1, k] - 2 * min_d_e[i, k])
            @constraint(model, C[i, k] <= min_d_e[i, k] + abs_dif_d_e[i, k] + sum(T[i, j] * Z[j, k] for j in Jobs) + sum(P_hat[i, j] * dev[i, j] * Z[j, k] for j in Jobs))
            # Strengthen option 3:
            #@constraint(model, C[i, k] <= C[i, k - 1] + C[i - 1, k] - min_d_e[i, k] + sum(T[i, j] * Z[j, k] for j in Jobs) + sum(P_hat[i, j] * dev[i, j] * Z[j, k] for j in Jobs))
        end
    end
    if budget_type == "machine"
        # for each machine r, which jobs i will have their processing times deviated?
        for r in M
            @constraint(model, sum(dev[r, i] for i in Jobs) <= Γ[r])  # c10
        end
    elseif budget_type == "global"  # a single budget, given all operations on all machines
        @constraint(model, sum(dev[r, i] for r in M, i in Jobs) <= Γ)  # c10
    end
    # ++++++++++ Objective function  (Constraint (4))
    @variable(model, obj)
    @objective(model, Max, obj)
    c11 = @constraint(model, obj == sum(w[j] * sum(C[m, k] * Z[j, k] for k in Jobs) for j in Jobs))
    # Solve the problem
    println("[Alt MIP] Calculating worst-case WCT given a permutation pi = $(permutation), m = $(m), n = $(n), w = $(w) and budget parameter Γ = $(Γ)...")
    println("====================================================================\n")
    println("[Alt] Solving robust flow shop budget model worst-case WCT...")
    println("====================================================================\n")
    println("single_cut_per_iteration = $(single_cut_per_iteration)")
    println("fractional_solution = $(fractional_solution)")
    println("bigM_type = $(bigM_type)")
    flush(stdout)
    status = JuMP.optimize!(model)
    solutions = []
    if JuMP.termination_status(model) == MOI.OPTIMAL
        println("\n===========  W O R S T - C A S E    A L T    W C T   S O L U T I O N  ===========\n")
        println("Optimal Objective Function value: ", objective_value(model))
        flush(stdout)
        obj_star = objective_value(model)
        # Generate other solutions with the same obj function value
        sol_count = 0
        previous_obj = 0.0
        dev_m = zeros(Float64, m, n)
        while JuMP.termination_status(model) == MOI.OPTIMAL
            sol_count += 1
            println("Solution #$(sol_count)")
            println("Solve time : ", solve_time(model))
            flush(stdout)
            elapsed_time += solve_time(model)
            if sol_count > 1
                previous_dev_m = deepcopy(dev_m)
                previous_obj = obj_star
            end
            obj_star = objective_value(model)
            dev_m = zeros(Float64, m, n)
            C_star = zeros(Float64, m, n)
            # IMPORTANT: return the whole deviation matrix, since the solution containing
            # deviation values may be fractional
            for r in 1:m
                 for k in 1:n
                     dev_m[r, k] = value(dev[r, k])
                     dev_m[r, k] = min(dev_m[r, k], 1.0)
                     dev_m[r, k] = max(dev_m[r, k], 0.0)
                     C_star[r, k] = value(C[r, k])
                 end
            end
            println("C_star = $(C_star)")
            println("WCT = $(sum(w[permutation[j]] * C_star[m, j] for j in 1:n))")
            println("Dev: $(dev_m); obj_star = $(obj_star)")
            flush(stdout)
            if single_cut_per_iteration || (max_solutions == 1)
                # Add new solution to the list
                push!(solutions, dev_m)
                println("Finished (single_cut_per_iteration == true).")
                flush(stdout)
                break
            end
            if fractional_solution
                # Add new solution to the list
                push!(solutions, dev_m)
                println("Finished (fractional_solution == true).")
                flush(stdout)
                break
            end
            if sol_count > 1
                if abs(obj_star - previous_obj) >= 1.e-4
                    sol_count -= 1
                    obj_star = previous_obj
                    println("Finished (different obj): $(obj_star)")
                    flush(stdout)
                    break
                end
                if dev_m in solutions
                    sol_count -= 1
                    obj_star = previous_obj
                    println("Finished (repeated solution).")
                    flush(stdout)
                    break
                end
                if sol_count >= max_solutions
                    # Add new solution to the list
                    push!(solutions, dev_m)
                    obj_star = previous_obj
                    println("Finished (max number of solutions == $(max_solutions)).")
                    flush(stdout)
                    break
                end
            end
            # Add new solution to the list
            push!(solutions, dev_m)
            if elapsed_time >= solver_time_limit
                println("Finished (time limit exceeded).")
                flush(stdout)
                break
            end
            # Add cuts to forbid the current solution
            # https://stackoverflow.com/questions/38317040/jump-how-to-get-multiple-solutions-from-getvaluex
            # http://yetanothermathprogrammingconsultant.blogspot.com/2011/10/integer-cuts.html
            # NOTE: Only valid for non-fractional (integer) solutions for dev variable !!!
            eps = 1.e-9
            @constraint(model, sum(dev[r, k] for r in 1:m for k in 1:n if dev_m[r, k] >= 0.9) - sum(dev[r, k] for r in 1:m for k in 1:n if dev_m[r, k] <= 0.1) <= sum(dev_m[r, k] for r in 1:m for k in 1:n) - 1)
            status = optimize!(model)
        end
        println("Found $(sol_count) solution(s).")
        flush(stdout)
        return obj_star, solutions, elapsed_time
    else
        println("ERROR: [Worst-case WCT MIP Problem] Optimal solution not found! Best bound: ", objective_bound(model))
        elapsed_time += solve_time(model)
        println("Model solve time: $(elapsed_time); solver_time_limit = $(solver_time_limit)")
        flush(stdout)
        return objective_bound(model), [], elapsed_time
    end
end

#test_model()
#test_model_2()
