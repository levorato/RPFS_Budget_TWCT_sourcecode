# ===================================================================================================
# robust_pfsp_wilson_worstcase_wct_mip.jl
# Worst-case determination procedure for the Robust Permutation Flow Shop Problem to Minimize the
# Weighted Sum of Job Completion Times (WCT), based on the budgeted uncertainty set (Bertsimas et al.).
# Given a permutation pi, finds the restricted worst-case scenario (i.e. for each job i,
# on which machines it will have its processing times deviated to their worst-case values), so that the
# objective function (weighted completion time - wct) reaches the worst-case (maximum) value.
# *** MIP approach for calculating the worst-case scenario.
# ===================================================================================================

using JuMP
using CPLEX

include("../deterministic_pfsp.jl")
include("../pfsp_file_reader.jl")
#include("robust_pfsp_budget_worstcase_wct_brute.jl")
#include("deterministic_pfsp_wilson_wct_mip.jl")

function test_wilson()
    basedir = joinpath(robust_ying_instances_folder, "data/10jobs/10%")
    time_limit = 1800
    for filename in searchdir(basedir, ".txt")
        if filename == "RB0101010.txt"
            inputfile = "$(basedir)/$(filename)"
            m, n, P_bar, P_hat = read_robust_ying_input_file(inputfile)
            # φ_star is the optimal Johnson solution obtained over the nominal processing time array (P_bar)
            w = [x for x in 1:n]
            if n < 10
                wct_star, pi, x = solve_deterministic_pfsp_wct(m, n, P_bar, w)
                println("Based on the nominal processing times, the optimal WCT wct_star is $(wct_star).")
                println("x: $(x)")
                wct_validation = calculate_wct(m, n, P_bar, w, pi, false)
                println("WCT Validation = $(wct_validation)")
                @assert abs(wct_star - wct_validation) < 1e-5
            else
                pi = [ 2, 3, 1, 4, 5, 6, 7, 8, 9, 10 ] #[x for x in 1:n]  # [9, 8, 7, 6, 10, 5, 4, 3, 2, 1] #
            end
            println("pi: $(pi)")
            println("w: $(w)")
            wct_validation = calculate_wct(m, n, P_bar, w, pi, false)
            println("WCT Validation = $(wct_validation)")
            Jobs = 1:n
            N = Jobs
            Machines = 1:m
            M = Machines
            println("Calculating robust PFSP budget solution (Wilson)...")
            Γ1 = (60 * n) / 100.0
            Γ2 = (60 * n) / 100.0
            Γ3 = (60 * n) / 100.0
            cost, scenarios, elapsed_time = solve_robust_model_worst_case_wilson(pi, m, n, [Γ1, Γ2, Γ3], P_bar, P_hat, w)
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
            if n <= 10
                wct_brute, worst_Γ_scenario_brute = worst_case_wct_brute_force_2_machines(pi, n, Γ1, Γ2, P_bar, P_hat, w, true)
                @assert abs(wct_brute - cost) < 1e-5
            else
                println("Skipping worst-case validation, since n > 10.")
            end
        end
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

# TODO implementar o retorno de multiplas solucoes para o pior caso, vide:
# http://yetanothermathprogrammingconsultant.blogspot.com/2011/10/integer-cuts.html
# https://stackoverflow.com/questions/38317040/jump-how-to-get-multiple-solutions-from-getvaluex
function solve_robust_model_worst_case_wilson(pi, m, n, Γ, P_bar, P_hat, w, max_solutions = 5; bigM_type = 3, max_cores = 16,
            solver = "CPLEX", budget_type = "machine", solver_time_limit = 7200, single_cut_per_iteration = false, fractional_solution = false)
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
        Z[pi[i], j] = 1.0
        j += 1
    end
    bigM = zeros(Float64, (m, n))
    fill!(bigM, calculate_bigM_naive(m, n, P_bar, P_hat))
    worst_Γ_scenario = [ [x for x in 1:n] for r in 1:m ]
    best_Γ_scenario = [ [] for r in 1:m ]
    max_C = calculate_cmax_given_scenario(pi, m, n, P_bar, P_hat, worst_Γ_scenario, true)
    min_C = calculate_cmax_given_scenario(pi, m, n, P_bar, P_hat, best_Γ_scenario, true)
    # https://spartanideas.msu.edu/2018/09/17/choosing-big-m-values/
    if bigM_type == 2
        # BigM (Q) : Option 1. Calculate the Makespan (Cmax) for the worst-case scenario
        #                      (on all machines r, all jobs j have their processing times
        #                       deviated to their maximum value : P_bar(r,j) + P_hat(r,j)).
        max_cmax = calculate_cmax_given_scenario(pi, m, n, P_bar, P_hat, worst_Γ_scenario) + 1
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
        for i in 2:m  # [i, k-1] x [i-1, k]
            for k in 2:n
                bigM[i, k] = max_C[i, k]  # max(abs(max_C[i-1, k] - min_C[i, k-1]), abs(max_C[i, k-1] - min_C[i-1, k])) + 1
            end
        end
    elseif bigM_type == 5
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
        for i in 2:m  # [i, k-1] x [i-1, k]
            for k in 2:n
                #bigM[i, k] = 2 * (max_C[i, k] - min(min_C[i, k-1] - (P_bar[i, k-1] + P_hat[i, k-1]), min_C[i-1, k] - (P_bar[i-1, k] + P_hat[i-1, k]))) + 1
                #bigM[i, k] = max_C[i, k] - min(min_C[i, k-1] - (P_bar[i, k-1]), min_C[i-1, k] - (P_bar[i-1, k])) + 1
                bigM[i, k] = 2 * (max_C[i, k] - min(min_C[i, k-1], min_C[i-1, k])) + 1
            end
        end
    elseif bigM_type == 7
        for i in 2:m  # [i, k-1] x [i-1, k]
            for k in 2:n
                bigM[i, k] = 2 * (max(max_C[i, k-1] - min_C[i-1, k], max_C[i-1, k] - min_C[i, k-1])) + 1
            end
        end
    else
        println("ERROR: bigM type not found!")
        return
    end
    ## ++++++++++ Model variables: B and Z
    # B_rj : start (begin) time of job in sequence position j on machine r
    @variable(model, B[M, N] >= 0)
    @variable(model, min_d_e[M, N] >= 0)
    @variable(model, disj[M, N], Bin)
    @variable(model, abs_dif_d_e[M, N] >= 0)
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
    ## ++++++++++ Model constraints
    T = P_bar
    c1 = Array{Any}(undef, n-1)
    for j in 1:n-1
        c1[j] = @constraint(model, B[1, j] + sum((P_hat[1, i] * dev[1, i]) * Z[i, j] for i in Jobs) + sum(T[1, i] * Z[i, j] for i in Jobs) == B[1, j+1])
    end
    @constraint(model, B[1, 1] == 0)
    c2 = Array{Any}(undef, m-1)
    for r in 1:m-1
        c2[r] = @constraint(model, B[r, 1] + sum((T[r, i] + (P_hat[r, i] * dev[r, i])) * Z[i, 1] for i in Jobs) == B[r+1, 1])
    end
    for r in 2:m
        for j in 2:n
            @constraint(model, min_d_e[r, j] <= B[r-1, j] + sum((T[r-1, i] + (P_hat[r-1, i] * dev[r-1, i])) * Z[i, j] for i in Jobs) )  # c5
            @constraint(model, min_d_e[r, j] <= B[r, j-1] + sum((T[r, i] + (P_hat[r, i] * dev[r, i])) * Z[i, j-1] for i in Jobs) )  # c6

            @constraint(model, abs_dif_d_e[r, j] <= (B[r-1, j] + sum((T[r-1, i] + (P_hat[r-1, i] * dev[r-1, i])) * Z[i, j] for i in Jobs))
                            - (B[r, j-1] + sum((T[r, i] + (P_hat[r, i] * dev[r, i])) * Z[i, j-1] for i in Jobs)) + bigM[r, j] * disj[r, j])  # c7
            @constraint(model, abs_dif_d_e[r, j] <= - ( (B[r-1, j] + sum((T[r-1, i] + (P_hat[r-1, i] * dev[r-1, i])) * Z[i, j] for i in Jobs))
                            - (B[r, j-1] + sum((T[r, i] + (P_hat[r, i] * dev[r, i])) * Z[i, j-1] for i in Jobs)) ) + bigM[r, j] * (1 - disj[r, j]))  # c8

            @constraint(model, B[r, j] <= min_d_e[r, j] + abs_dif_d_e[r, j])  #c10
            # max[r, j] == min[r, j] + abs(d[r,j] - e[r,j])
            #@constraint(model, dif_a_b[r, j] == d[r, j] - e[r, j])
            # https://math.stackexchange.com/questions/1954992/linear-programming-minimizing-absolute-values-and-formulate-in-lp
            # Given that the we are maximizing the objective function : max max(d-e, e-d)
            ###@constraint(model, B[r, j] <= min_d_e[r, j] + (d[r, j] - e[r, j]) + bigM * disj1[r, j])  # c7
            ###@constraint(model, B[r, j] <= min_d_e[r, j] - (d[r, j] - e[r, j]) + bigM * disj2[r, j])  # c8
            ###@constraint(model, disj1[r, j] + disj2[r, j] == 1)  # c9
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
    c12 = @constraint(model, obj == sum(w[i] * sum((B[m, j] + T[m, i] + (P_hat[m, i] * dev[m, i])) * Z[i, j] for j in Jobs) for i in Jobs) )
    # Solve the problem
    println("[Wilson] Calculating worst-case WCT given a permutation pi = $(pi), m = $(m), n = $(n), w = $(w) and budget parameter Γ = $(Γ)...")
    println("====================================================================\n")
    println("[Wilson] Solving robust flow shop budget model worst-case WCT...")
    println("====================================================================\n")
    println("single_cut_per_iteration = $(single_cut_per_iteration)")
    println("fractional_solution = $(fractional_solution)")
    println("bigM_type = $(bigM_type)")
    flush(stdout)
    status = optimize!(model)
    solutions = []
    if JuMP.termination_status(model) == MOI.OPTIMAL
        println("\n===========  W O R S T - C A S E   W I L S O N   W C T   S O L U T I O N  ===========\n")
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
            # IMPORTANT: return the whole deviation matrix, since the solution containing
            # deviation values may be fractional
            for r in 1:m
                 for k in 1:n
                     dev_m[r, k] = value(dev[r, k])
                     dev_m[r, k] = min(dev_m[r, k], 1.0)
                     dev_m[r, k] = max(dev_m[r, k], 0.0)
                 end
            end
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

# TODO O MODELO AINDA ESTA COM ERROS: O RESULTADO NAO BATE COM O CALCULO DO PIOR CASO (FORCA BRUTA)
# test_wilson()
