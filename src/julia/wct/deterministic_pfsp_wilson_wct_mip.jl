# ===================================================================================================
# deterministic_pfsp_wilson_wct_mip.jl
# Deterministic Permutation Flow Shop Problem to Minimize the Weighted Sum of Job Completion Times (WCT).
# *** MILP approach based on Wilson PFSP Cmax model.
# ===================================================================================================

using JuMP
using CPLEX

function solve_deterministic_pfsp_wct(m, n, T, w, time_limit = 1800)
    cmax, permutation, time, is_opt, Z = solve_deterministic_wilson_model_wct(m, n, T, w, time_limit)
    return cmax, permutation, Z
end

# Solves the Wilson PFSP MILP Model to minimize the weighted sum of completion times (Deterministic)
function solve_deterministic_wilson_model_wct(m, n, T, w, time_limit; solver="CPLEX")
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
    println("Solving Wilson PFSP MILP Model - WCT (Deterministic)")
    bigM = 0
    for r in M
        for i in Jobs
            bigM += T[r, i]
        end
    end
    ## ++++++++++ Model variables: B and Z
    # B_rj : start (begin) time of job in sequence position j on machine r
    @variable(model, B[M, N])
    for i in M
        for j in N
            @constraint(model, B[i, j] >= 0)
        end
    end
    # Z_ij = 1, if job i is assigned to sequence position j, 0 otherwise
    #   x x x Originally: @variable(model, Z[N, N], Bin)
    @variable(model, Z[N, N], Bin)
    # Bt[i] : start time of job i on the last machine (M_m)
    @variable(model, Bt[N] >= 0)
    # b1: auxiliary variable used in the BigM-based disjuction (see last constraint)
    @variable(model, b1[N, N], Bin)
    ## ++++++++++ Model constraints
    # Assignement constraints
    for i in Jobs
        @constraint(model, sum(Z[i, j] for j in Jobs) == 1)
    end
    for j in Jobs
        @constraint(model, sum(Z[i, j] for i in Jobs) == 1)
    end
    for j in 1:n-1
        @constraint(model, B[1, j] + sum(T[1, i] * Z[i, j] for i in Jobs) == B[1, j+1])
    end
    @constraint(model, B[1, 1] == 0)
    for r in 1:m-1
        @constraint(model, B[r, 1] + sum(T[r, i] * Z[i, 1] for i in Jobs) == B[r+1, 1])
    end
    for r in 1:m-1
        for j in 2:n
            @constraint(model, B[r, j] + sum(T[r, i] * Z[i, j] for i in Jobs) <= B[r+1, j])
        end
    end
    for r in 2:m
        for j in 1:n-1
            @constraint(model, B[r, j] + sum(T[r, i] * Z[i, j] for i in Jobs) <= B[r, j+1])
        end
    end
    for i in 1:n  # Disjuctive constraints used to determine the start time of job i on the last machine m
        for j in 1:n
            @constraint(model, Bt[i] >= B[m, j] - bigM * (1 - b1[i, j]) )
            @constraint(model, b1[i, j] <= Z[i, j])
        end
        @constraint(model, sum(b1[i, j] for j in 1:n) == 1)
    end
    # ++++++++++ Objective function
    #@objective(model, Min, B[m, n] + sum(T[m, i] * Z[i, n] for i in Jobs))
    @objective(model, Min, sum(w[i] * Bt[i] for i in Jobs) + sum(w[i] * sum(T[m, i] * Z[i, j] for j in Jobs) for i in Jobs))
    println("\n============== SOLVING PFSP WILSON MODEL - WCT ==============\n")
    status = optimize!(model)
    Z_values = zeros(Float64, (n, n))
    if termination_status(model) == MOI.OPTIMAL
        println("Optimal Objective Function value: ", objective_value(model))
        permutation = zeros(Int64, n)
        begin_times = zeros(Float64, n)
        for j in Jobs
            begin_times[j] = value(Bt[j])
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
        println("Begin times of each job: $(begin_times)")
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


# Solves the Wilson PFSP MILP Model to minimize the weighted sum of completion times (Deterministic)
function solve_deterministic_wilson_model_wct_v2(m, n, T, w, time_limit)
    model = Model(CPLEX.Optimizer)  # MIPGap=0.005
    set_optimizer_attribute(model, "CPXPARAM_TimeLimit", time_limit)
    Jobs = 1:n
    N = Jobs
    Machines = 1:m
    M = Machines
    println("Solving Wilson PFSP MILP Model - WCT (Deterministic)")
    bigM = 0
    for r in M
        for i in Jobs
            bigM += T[r, i]
        end
    end
    ## ++++++++++ Model variables: B and Z
    # B_rj : start (begin) time of job in sequence position j on machine r
    @variable(model, B[M, N] >= 0)
    # Z_ij = 1, if job i is assigned to sequence position j, 0 otherwise
    #   x x x Originally: @variable(model, Z[N, N], Bin)
    @variable(model, Z[N, N], Bin)
    # Bt[i] : start time of job i on the last machine (M_m)
    @variable(model, Bt[N] >= 0)
    ## ++++++++++ Model constraints
    # Assignement constraints
    for i in Jobs
        @constraint(model, sum(Z[i, j] for j in Jobs) == 1)
    end
    for j in Jobs
        @constraint(model, sum(Z[i, j] for i in Jobs) == 1)
    end
    for j in 1:n-1
        @constraint(model, B[1, j] + sum(T[1, i] * Z[i, j] for i in Jobs) == B[1, j+1])
    end
    @constraint(model, B[1, 1] == 0)
    for r in 1:m-1
        @constraint(model, B[r, 1] + sum(T[r, i] * Z[i, 1] for i in Jobs) == B[r+1, 1])
    end
    for r in 1:m-1
        for j in 2:n
            @constraint(model, B[r, j] + sum(T[r, i] * Z[i, j] for i in Jobs) <= B[r+1, j])
        end
    end
    for r in 2:m
        for j in 1:n-1
            @constraint(model, B[r, j] + sum(T[r, i] * Z[i, j] for i in Jobs) <= B[r, j+1])
        end
    end
    for i in 1:n  # Disjuctive constraints used to determine the start time of job i on the last machine m
        for j in 1:n
            @constraint(model, Bt[i] >= B[m, j] - bigM * (1 - Z[i, j]) )


        end
    end

    for j in 1:n
        for k in 1:n
            if k <= j  # pi_jk
                π[j, j]
            else

            end

            if k >= j  # ro_jk

            else

            end
        end
    end
    # ++++++++++ Objective function
    #@objective(model, Min, B[m, n] + sum(T[m, i] * Z[i, n] for i in Jobs))
    @objective(model, Min, sum(w[i] * Bt[i] for i in Jobs) + sum(w[i] * sum(T[m, i] * Z[i, j] for j in Jobs) for i in Jobs))
    println("\n============== SOLVING PFSP WILSON MODEL - WCT ==============\n")
    status = optimize!(model)
    Z_values = zeros(Float64, (n, n))
    if termination_status(model) == MOI.OPTIMAL
        println("Optimal Objective Function value: ", objective_value(model))
        permutation = zeros(Int64, n)
        begin_times = zeros(Float64, n)
        for j in Jobs
            begin_times[j] = value(Bt[j])
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
        println("Begin times of each job: $(begin_times)")
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
