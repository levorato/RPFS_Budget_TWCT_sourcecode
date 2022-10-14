# ===================================================================================================
# robust_pfsp_alt_worstcase_wct_non_linear.jl
# Worst-case determination procedure for the Robust Permutation Flow Shop Problem to Minimize the
# Weighted Sum of Job Completion Times (WCT), based on the budgeted uncertainty set (Bertsimas et al.).
# Given a permutation pi, finds the restricted worst-case scenario (i.e. for each job i,
# on which machines it will have its processing times deviated to their worst-case values), so that the
# objective function (weighted completion time - wct) reaches the worst-case (maximum) value.
# *** Non-linear approach for calculating the worst-case scenario.
# ===================================================================================================

using JuMP #, MathOptFormat
using CPLEX
using Ipopt
using MosekTools, Mosek
using Gurobi

include("../config.jl")
include("../deterministic_pfsp.jl")
include("../pfsp_file_reader.jl")
include("robust_pfsp_budget_worstcase_wct_brute.jl")
include("robust_pfsp_alt_worstcase_wct_mip.jl")
include("robust_pfsp_budget_worstcase_wct_global.jl")
#include("deterministic_pfsp_alt_wct_mip.jl")
include("../robust_pfsp_budget_worstcase_cmax_brute.jl")

function test_model_2()
    # basedir = joinpath(robust_ying_instances_folder, "rob-pfsp-wct/20jobs")  filename == "RB0205010_20_2_50_wct_inputs.txt"
    #basedir = joinpath(robust_ying_instances_folder, "rob-pfsp-wct", "10x5", "50%")
    basedir = joinpath(robust_ying_instances_folder, "rob-pfsp-wct", "15x5")
    time_limit = 4000
    time_spent = [0.0 for x in 1:8]
    count = 0
    for bigM_type in [7]
        println("BigM type: $(bigM_type)\n====================================")
        for filename in searchdir(basedir, ".txt")
            if filename == "RB0151001_15_5_R400_wct_inputs.txt"
                Γ = 4  # 75 - 4  # 30
            #if filename == "RB0155007_15_5_50_wct_inputs.txt"
            #if filename == "RB0105007_10_5_50_wct_inputs.txt"
                #Γ = 2
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
                    #permutation = [ 9, 10, 5, 6, 2, 4, 1, 3, 8, 7 ]
                    #permutation = [ x for x in 1:n ]
                    permutation = [1, 11, 4, 2, 12, 8, 15, 13, 10, 3, 9, 7, 6, 14, 5]
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
                cost_nlp, scenarios_nlp, elapsed_time_nlp = solve_robust_model_worst_case_alt_non_linear(permutation, m, n, Γ, P_bar, P_hat, w, 1;
                    solver="Ipopt", bigM_type=bigM_type, max_cores = 2, budget_type = "global")
                scenario = scenarios_nlp[1]
                deviated_jobs = [[] for x in 1:m]
                for r in 1:m
                    for k in 1:n
                        if scenario[r, k] >= 0.9
                            push!(deviated_jobs[r], k)
                        end
                    end
                end
                println("[Non-linear] Worst-case scenario: $(deviated_jobs) ; $(scenario)")
                println("[Non-linear] The worst-case cost from model is $(cost_nlp).")
                wct = calculate_wct_given_scenario(permutation, m, n, P_bar, P_hat, w, scenarios_nlp[1], true )     # true)
                println("[Non-linear] [Validation] Given scenario = $(deviated_jobs), the WCT is $(wct).")


                #max_wct, scenario_deviated_jobs = worst_case_wct_dp_m_machines(permutation, m, n, Γ, P_bar, P_hat, w)
                #dev_dp = zeros(Float64, m, n)
                #for (r, i) in scenario_deviated_jobs
                #    dev_dp[r, i] = 1.0
                #end
                if worst_case_brute_force_is_worth(m, n, Int64(ceil(Γ)))
                    @time max_fo, worst_Γ_scenario = worst_case_cmax_brute_force_global_budget(permutation, m, n, Γ, P_bar, P_hat, true)
                end

                if worst_case_brute_force_is_worth(m, n, Int64(ceil(Γ)))
                    @time max_fo, worst_Γ_scenario = worst_case_wct_brute_force_global_budget(permutation, m, n, Γ, P_bar, P_hat, w, true)
                end


                cost, scenarios, elapsed_time = solve_robust_model_worst_case_alt(permutation, m, n, Γ, P_bar, P_hat, w, 1;
                    solver="Gurobi", bigM_type=bigM_type, max_cores = 2, budget_type="global")  #, initial_scenario=[dev_dp, max_wct])  # =scenarios_nlp[1])
                # TODO Fazer algoritmo de enumeracao (brute force), pensar no complemento
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
            end
        end
    end
    println("Total count is $(count)")
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
function solve_robust_model_worst_case_alt_non_linear(permutation, m, n, Γ, P_bar, P_hat, w, max_solutions = 5; bigM_type = 3, max_cores = 16,
            solver = "CPLEX", budget_type = "machine", solver_time_limit = 7200.0, single_cut_per_iteration = false, fractional_solution = false)
    Jobs = 1:n
    N = Jobs
    Machines = 1:m
    M = Machines
    elapsed_time = 0.0
    println("*** Solver = $(solver)")
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
    elseif solver == "Ipopt"
        model = Model(Ipopt.Optimizer)
        MOI.set(model, MOI.RawOptimizerAttribute("max_cpu_time"), Float64(solver_time_limit))
        #MOI.set(model, MOI.RawOptimizerAttribute("print_level"), 1)
        #set_optimizer_attribute(model, "max_cpu_time", solver_time_limit)
        #set_optimizer_attribute(model, "print_level", 0)
    elseif solver == "NLopt"
        model = Model(NLopt.Optimizer)
        set_optimizer_attribute(model, "verbose", 1)
        set_optimizer_attribute(model, "algorithm", :LD_SLSQP)  # :LN_COBYLA), LD_TNEWTON_PRECOND_RESTART
    elseif solver == "Mosek"
        model = Model(optimizer_with_attributes(Mosek.Optimizer, "QUIET" => false, "INTPNT_CO_TOL_DFEAS" => 1e-7))
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
    ## ++++++++++ Model variables: C (Z is fixed parameter)
    # C_ik completion time of job k on machine i
    @variable(model, dif_d_e[M, N])
    @variable(model, abs_dif_d_e[M, N] >= 0)
    @variable(model, C[M, N] >= 0)
    @variable(model, min_d_e[M, N] >= 0)
    #@variable(model, disj[M, N], Bin)
    #@variable(model, avg_d_e[M, N] >= 0)
    @variable(model, dev[M, N] >= 0)
    for r in 1:m
        for i in 1:n
            @constraint(model, dev[r, i] <= 1)
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
            #@constraint(model, abs_dif_d_e[i, k] >= (C[i, k - 1] - C[i - 1, k]))
            #@constraint(model, abs_dif_d_e[i, k] >= -(C[i, k - 1] - C[i - 1, k]))
            #@constraint(model, abs_dif_d_e[i, k] <= (C[i, k - 1] - C[i - 1, k]) + bigM[i, k] * disj[i, k])
            #@constraint(model, abs_dif_d_e[i, k] <= -(C[i, k - 1] - C[i - 1, k]) + bigM[i, k] * (1 - disj[i, k]))
            @constraint(model, dif_d_e[i, k] == (C[i, k - 1] - C[i - 1, k]))
            @NLconstraint(model, abs_dif_d_e[i, k] == abs(dif_d_e[i, k]))
            ###@constraint(model, avg_d_e[i, k] == (C[i, k - 1] + C[i - 1, k]) / 2)
            ###@constraint(model, abs_dif_d_e[i, k] == 2 * (avg_d_e[i, k] - min_d_e[i, k]))

            @constraint(model, C[i, k] <= min_d_e[i, k] + abs_dif_d_e[i, k] + sum(T[i, j] * Z[j, k] for j in Jobs) + sum(P_hat[i, j] * dev[i, j] * Z[j, k] for j in Jobs))
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
    JuMP.optimize!(model)
    solutions = []
    status = MOI.get(model, MOI.TerminationStatus())
    println("status = $(status)")
    if status == MOI.LOCALLY_SOLVED
        println("\n===========  W O R S T - C A S E    A L T    W C T   S O L U T I O N  ===========\n")
        println("Optimal Objective Function value: ", objective_value(model))
        flush(stdout)
        obj_star = objective_value(model)
        # Generate other solutions with the same obj function value
        sol_count = 0
        previous_obj = 0.0
        dev_m = zeros(Float64, m, n)
        sol_count += 1
        println("Solution #$(sol_count)")
        println("Solve time : ", solve_time(model))  # MOI.get(model, MOI.SolveTime())
        flush(stdout)
        elapsed_time += solve_time(model)
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
        if elapsed_time >= solver_time_limit
            println("Finished (time limit exceeded).")
            flush(stdout)
            return obj_star, solutions, elapsed_time
        end
        # Add new solution to the list
        push!(solutions, dev_m)
        println("Finished (single_cut_per_iteration == true).")
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

test_model_2()
