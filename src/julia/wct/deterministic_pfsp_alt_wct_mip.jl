# ===================================================================================================
# deterministic_pfsp_alt_wct_mip.jl
# Deterministic Permutation Flow Shop Problem to Minimize the Weighted Sum of Job Completion Times (WCT).
# *** MILP approach based on the Alternative PFSP Cmax model.
# ===================================================================================================

using JuMP
using CPLEX

function solve_deterministic_pfsp_wct(m, n, T, w, time_limit = 1800)
    cmax, permutation, time, is_opt, Z = solve_deterministic_alt_model_wct(m, n, T, w, time_limit)
    return cmax, permutation, Z
end

# Solves the Alternative PFSP MILP Model to minimize the weighted sum of completion times (Deterministic)
# Based on Manne PFSP Model
function solve_deterministic_alt_model_wct(m, n, T, w, time_limit)
    model = Model(CPLEX.Optimizer)  # MIPGap=0.005
    set_optimizer_attribute(model, "CPXPARAM_TimeLimit", time_limit)
    Jobs = 1:n
    N = Jobs
    Machines = 1:m
    M = Machines
    println("Solving Alternative (Manne) PFSP MILP Model - WCT (Deterministic)")
    P = 0
    for r in M
        for i in Jobs
            P += T[r, i]
        end
    end
    ## ++++++++++ Model variables: C and Z
    # # C_ri completion time of job i on machine r
    @variable(model, C[M, N] >= 0)
    # D_ik = 1, if job i is scheduled any time before job k; 0, otherwise;
    @variable(model, D[N, N], Bin)
    ## ++++++++++ Model constraints
    for i in 1:n
        @constraint(model, C[1, i] >= T[1, i])
    end
    for r in 2:m
        for i in 1:n
            @constraint(model, C[r, i] - C[r-1, i] >= T[r, i])
        end
    end
    for r in 1:m
        for i in 1:(n-1)
            for k in (i+1):n
                @constraint(model, C[r, i] - C[r, k] + P * D[i, k] >= T[r, i])
                @constraint(model, C[r, i] - C[r, k] + P * D[i, k] <= P - T[r, k])
            end
        end
    end
    # ++++++++++ Objective function
    @objective(model, Min, sum(w[j] * C[m, j] for j in Jobs))
    println("\n============== SOLVING PFSP ALTERNATIVE MODEL - WCT ==============\n")
    status = optimize!(model)
    Z_values = zeros(Float64, (n, n))
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
