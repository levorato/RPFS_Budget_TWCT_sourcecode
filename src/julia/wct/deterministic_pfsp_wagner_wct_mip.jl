# ===================================================================================================
# deterministic_pfsp_wagner_wct_mip.jl
# Deterministic Permutation Flow Shop Problem to Minimize the Weighted Sum of Job Completion Times (WCT).
# *** MILP approach based on Wagner PFSP Cmax model.
# ===================================================================================================

using JuMP
using CPLEX

include("../pfsp_file_reader.jl")
include("robust_pfsp_budget_worstcase_wct_brute.jl")

function solve_deterministic_pfsp_wct(m, n, T, w, time_limit = 1800)
    cmax, permutation, time, is_opt, Z = solve_deterministic_wagner_model_wct(m, n, T, w, time_limit)
    return cmax, permutation, Z
end

# Solves the Wagner PFSP MILP Model to minimize the weighted sum of completion times (Deterministic)
function solve_deterministic_wagner_model_wct(m, n, T, w, time_limit; solver="CPLEX")
    if solver == "Gurobi"
        model = Model(Gurobi.Optimizer)  # MIPGap=0.005
        MOI.set(model, MOI.RawOptimizerAttribute("TimeLimit"), time_limit)
    elseif solver == "CPLEX"
        model = Model(CPLEX.Optimizer)  # MIPGap=0.005
        MOI.set(model, MOI.RawOptimizerAttribute("CPX_PARAM_TILIM"), time_limit)
    else
      println("No solver defined")
    end
    set_optimizer_attribute(model, "CPXPARAM_TimeLimit", time_limit)
    Jobs = 1:n
    N = Jobs
    Machines = 1:m
    M = Machines
    println("Solving Wagner PFSP MILP Model - WCT (Deterministic)")
    bigM = 0
    for r in M
        for i in Jobs
            bigM += T[r, i]
        end
    end
    ## ++++++++++ Model variables: X, Y and Z
    # X_rj : idle time on machine r before the start of job in sequence position j
    @variable(model, X[M, N] >= 0)
    # Y_rj idle time of job in sequence position j after it finishes processing on machine r
    @variable(model, Y[M, N] >= 0)
    # Z_ij = 1, if job i is assigned to sequence position j, 0 otherwise
    #   x x x Originally: @variable(model, Z[N, N], Bin)
    @variable(model, Z[N, N], Bin)
    # Xt[i] : idle time, before the start of job in sequence position j, on the last machine (M_m)
    @variable(model, Xt[N] >= 0)
    @variable(model, St[N] >= 0)
    # b1: auxiliary variable used in the BigM-based disjuction (see last constraint)
    @variable(model, b1[N, N], Bin)
    ## ++++++++++ Model constraints
    # wct : weighted sum of job completion times
    @variable(model, wct >= 0)
    ## ++++++++++ Model constraints
    # Assignement constraints
    for i in Jobs
        @constraint(model, sum(Z[i, j] for j in Jobs) == 1)
    end
    for j in Jobs
        @constraint(model, sum(Z[i, j] for i in Jobs) == 1)
    end
    for r in 1:m-1
        for j in 1:n-1
            @constraint(model, sum(T[r, i] * Z[i, j+1] for i in Jobs) - sum(T[r+1, i] * Z[i, j] for i in Jobs) + X[r, j+1] - X[r+1, j+1] + Y[r, j+1] - Y[r, j] == 0)
        end
    end
    for r in 1:m-1
        @constraint(model, sum(T[r, i] * Z[i, 1] for i in Jobs) + X[r, 1] - X[r+1, 1] + Y[r, 1] == 0)
    end
    for i in 1:n  # Disjuctive constraints used to determine the start time of job i on the last machine m
        for j in 1:n
            @constraint(model, Xt[i] >= sum(X[m, p] for p in 1:j) - bigM * (1 - b1[i, j]) )
            @constraint(model, St[i] >= sum(sum(T[m, p] * Z[p, k] for p in Jobs) for k in 1:j) - bigM * (1 - b1[i, j]) )
            @constraint(model, b1[i, j] <= Z[i, j])
        end
        @constraint(model, sum(b1[i, j] for j in 1:n) == 1)
    end
    # sum(w[i] * sum(T[m, i] * Z[i, j] for j in Jobs) for i in Jobs))
    # sum(w[i] * T[m, i] for i in Jobs)
    # @constraint(model, Cmax == (sum(T[m, i] for i in Jobs) + sum(X[m, p] for p in Jobs)))
    @constraint(model, wct == sum(w[i] * St[i] for i in Jobs) + sum(w[i] * Xt[i] for i in Jobs) )
    # ++++++++++ Objective function
    @objective(model, Min, wct )
    println("\n============== SOLVING PFSP WAGNER MODEL - WCT ==============\n")
    status = optimize!(model)
    Z_values = zeros(Float64, (n, n))
    optimal = false
    if termination_status(model) == MOI.OPTIMAL
        println("Optimal Objective Function value: ", objective_value(model))
        permutation = zeros(Int64, n)
        for j in Jobs
            for i in Jobs
                Z_values[i, j] = value(Z[i, j])
                if Z_values[i, j] >= 0.99  # Fix rounding error problem
                    if j <= n
                        permutation[j] = i
                    end
                end
            end
        end
        println("Permutation: $(permutation)")
        return objective_value(model), permutation, solve_time(model), true, Z_values, objective_bound(model)
    elseif termination_status(model) == MOI.TIME_LIMIT && has_values(model)
        println("Optimal solution not found!")
        println("Suboptimal Objective Function value: ", objective_value(model))
        permutation = zeros(Int64, n)
        for j in Jobs
            for i in Jobs
                Z_values[i, j] = value(Z[i, j])
                if Z_values[i, j] >= 0.99  # Fix rounding error problem
                    if j <= n
                        permutation[j] = i
                    end
                end
            end
        end
        println("Permutation: $(permutation)")
        return objective_value(model), permutation, solve_time(model), false, Z_values, objective_bound(model)
    else
        println("Optimal solution not found! Best bound: ", objective_bound(model))
        return objective_bound(model), zeros(Int64, n), solve_time(model), false, Z_values, objective_bound(model)
    end
end
