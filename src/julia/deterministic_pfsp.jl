# ===================================================================================================
# deterministic_pfsp.jl
# Deterministic Permutation Flow Shop Problem to Minimize the Makespan (Cmax).
# *** MILP approach based on Wilson PFSP Cmax model.
# ===================================================================================================

# Attention: run with Julia >= 1.3.0
using JuMP
include("pfsp_util.jl")

function calculate_cmax_given_scenario(pi, m, n, P_bar, P_hat, m_dev, whole_matrix = false; global_budget = false)
    if global_budget
        dev_machine = [[] for x in 1:m]
        for (r, i) in m_dev
            push!(dev_machine[r], i)
        end
        m_dev = dev_machine
    end
    for i in 1:m
        for j in m_dev[i]
            if j > 0
                P_bar[i, j] += P_hat[i, j]
            end
        end
    end
    cmax = calculate_cmax(m, n, P_bar, pi, whole_matrix)
    for i in 1:m
        for j in m_dev[i]
            if j > 0
                P_bar[i, j] -= P_hat[i, j]
            end
        end
    end
    return cmax
end

# Solves the Wilson PFSP MILP Model (Deterministic)
function solve_deterministic_wilson_model(m, n, T, time_limit)
    if solver == "Gurobi"
        model = Model(Gurobi.Optimizer)  # MIPGap=0.005
        MOI.set(model, MOI.RawOptimizerAttribute("TimeLimit"), time_limit)
    elseif solver == "CPLEX"
        model = Model(CPLEX.Optimizer)  # MIPGap=0.005
        MOI.set(model, MOI.RawOptimizerAttribute("CPX_PARAM_TILIM"), time_limit)
        MOI.set(model, MOI.RawOptimizerAttribute("CPX_PARAM_THREADS"), 16)
    else
      println("No solver defined")
    end
    Jobs = 1:n
    N = Jobs
    Machines = 1:m
    M = Machines
    println("Solving Wilson PFSP MILP Model (Deterministic)")
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
    # ++++++++++ Objective function
    @objective(model, Min, B[m, n] + sum(T[m, i] * Z[i, n] for i in Jobs))
    println("\n============== SOLVING PFSP WILSON MODEL ==============\n")
    status = optimize!(model)
    Z_values = zeros(Float64, (n, n))
    if termination_status(model) == MOI.OPTIMAL
        println("Optimal Objective Function value: ", objective_value(model))
        permutation = zeros(Int64, n)
        for j in Jobs
            for i in Jobs
                Z_values[i, j] = value(Z[i, j])
                if Z_values[i, j] >= 0.9  # Fix rounding error problem
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

function solve_deterministic_wagner_model(m, n, T, time_limit)
    if solver == "Gurobi"
        model = Model(Gurobi.Optimizer)  # MIPGap=0.005
        MOI.set(model, MOI.RawOptimizerAttribute("TimeLimit"), time_limit)
    elseif solver == "CPLEX"
        model = Model(CPLEX.Optimizer)  # MIPGap=0.005
        MOI.set(model, MOI.RawOptimizerAttribute("CPX_PARAM_TILIM"), time_limit)
        MOI.set(model, MOI.RawOptimizerAttribute("CPX_PARAM_THREADS"), 16)
    else
      println("No solver defined")
    end
    Jobs = 1:n
    N = Jobs
    Machines = 1:m
    M = Machines
    println("Solving Wagner PFSP MILP Model (Deterministic)")
    ## ++++++++++ Model variables: X, Y and Z
    # X_rj : idle time on machine r before the start of job in sequence position j
    @variable(model, X[M, N] >= 0)
    # Y_rj idle time of job in sequence position j after it finishes processing on machine r
    @variable(model, Y[M, N] >= 0)
    # Z_ij = 1, if job i is assigned to sequence position j, 0 otherwise
    @variable(model, Z[N, N], Bin)
    # Cmax : maximum flowtime (makespan) of the schedule determined by the completion time of
    #        the job in the last sequence position on the last machine.
    @variable(model, Cmax >= 0)
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
    @constraint(model, Cmax == (sum(T[m, i] for i in Jobs) + sum(X[m, p] for p in Jobs)))
    # ++++++++++ Objective function
    @objective(model, Min, Cmax )
    println("\n============== SOLVING WAGNER MODEL ==============\n")
    status = optimize!(model)
    Z_values = zeros(Float64, (n, n))
    if termination_status(model) == MOI.OPTIMAL
        println("Optimal Objective Function value: ", objective_value(model))
        permutation = zeros(Int64, n)
        for j in Jobs
            for i in Jobs
                Z_values[i, j] = value(Z[i, j])
                if Z_values[i, j] >= 0.9  # Fix rounding error problem
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

function solve_deterministic_manne_model(m, n, T, time_limit)
    if solver == "Gurobi"
        model = Model(Gurobi.Optimizer)  # MIPGap=0.005
        MOI.set(model, MOI.RawOptimizerAttribute("TimeLimit"), time_limit)
    elseif solver == "CPLEX"
        model = Model(CPLEX.Optimizer)  # MIPGap=0.005
        MOI.set(model, MOI.RawOptimizerAttribute("CPX_PARAM_TILIM"), time_limit)
        MOI.set(model, MOI.RawOptimizerAttribute("CPX_PARAM_THREADS"), 16)
    else
      println("No solver defined")
    end
    Jobs = 1:n
    N = Jobs
    Machines = 1:m
    M = Machines
    println("Solving Manne PFSP MILP Model (Deterministic)")
    # P is a very large constant (big-M)
    random_pi = [x for x in 1:n]
    P = calculate_cmax(m, n, T, random_pi) + 1
    ## ++++++++++ Model variables: Cmax, C and D
    @variable(model, Cmax >= 0)
    # # C_ri completion time of job i on machine r
    @variable(model, C[M, N])
    for i in M
        for j in N
            @constraint(model, C[i, j] >= 0)
        end
    end
    # D_ik = 1, if job i is scheduled any time before job k; 0, otherwise (i < k)
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
    for i in 1:n
        @constraint(model, Cmax >= C[m, i])
    end
    # ++++++++++ Objective function
    @objective(model, Min, Cmax)
    println("\n============== SOLVING PFSP MANNE MODEL ==============\n")
    status = optimize!(model)
    D_star = zeros(Float64, (n, n))
    if termination_status(model) == MOI.OPTIMAL
        println("Optimal Objective Function value: ", objective_value(model))
        for i in Jobs  # Job i
            for k in Jobs  # Job k
                if value(D[i, k]) >= 0.9
                    D_star[i, k] = 1
                end
            end
        end
        permutation = topological_sort_by_dfs(D_star, n)
        println("Permutation: $(permutation)")
        return objective_value(model), permutation, solve_time(model), true, D_star, objective_bound(model)
    elseif termination_status(model) == MOI.TIME_LIMIT && has_values(model)
        println("Optimal solution not found!")
        println("Suboptimal Objective Function value: ", objective_value(model))
        for i in Jobs  # Job i
            for k in Jobs  # Job k
                if value(D[i, k]) >= 0.9
                    D_star[i, k] = 1
                end
            end
        end
        permutation = topological_sort_by_dfs(D_star, n)
        println("Permutation: $(permutation)")
        return objective_value(model), permutation, solve_time(model), false, D_star, objective_bound(model)
    else
        println("Optimal solution not found! Best bound: ", objective_bound(model))
        return objective_bound(model), zeros(Int64, n), solve_time(model), false, D_star, objective_bound(model)
    end
end

# Assignment / Completion-time formulation based was on:
#    Bektaş, T., Hamzadayı, A., & Ruiz, R. (2020). Benders decomposition for the mixed
#    no-idle permutation flowshop scheduling problem. Journal of Scheduling, 513–523.
#    https://doi.org/10.1007/s10951-020-00637-8
function solve_deterministic_alternative_assignment_model(m, n, T, time_limit)
    if solver == "Gurobi"
        model = Model(Gurobi.Optimizer)  # MIPGap=0.005
        MOI.set(model, MOI.RawOptimizerAttribute("TimeLimit"), time_limit)
    elseif solver == "CPLEX"
        model = Model(CPLEX.Optimizer)  # MIPGap=0.005
        MOI.set(model, MOI.RawOptimizerAttribute("CPX_PARAM_TILIM"), time_limit)
        MOI.set(model, MOI.RawOptimizerAttribute("CPX_PARAM_THREADS"), 16)
    else
      println("No solver defined")
    end
    Jobs = 1:n
    N = Jobs
    Machines = 1:m
    M = Machines
    println("Solving Alternative Assignment PFSP MILP Model (Deterministic)")
    ## ++++++++++ Model variables: C and Z
    # # C_ik completion time of job k on machine i
    @variable(model, C[M, N])
    for i in M
        for k in N
            @constraint(model, C[i, k] >= 0)
        end
    end
    # Z_ij = 1, if job i is assigned to sequence position j, 0 otherwise
    #   x x x Originally: @variable(model, Z[N, N], Bin)
    @variable(model, Z[N, N], Bin)
    ## ++++++++++ Model constraints
    # Assignement constraints
    for i in Jobs
        @constraint(model, sum(Z[i, j] for j in Jobs) == 1)
    end
    for j in Jobs
        @constraint(model, sum(Z[i, j] for i in Jobs) == 1)
    end
    for k in 1:n
        @constraint(model, C[1, k] >= sum(T[1, j] * Z[j, 1] for j in Jobs))
    end
    for k in 1:n  # (A.18)
        for i in 2:m
            @constraint(model, C[i, k] - C[i-1, k] >= sum(T[i, j] * Z[j, k] for j in Jobs))
        end
    end
    for k in 2:n  # (A.17)
        for i in 1:m
            @constraint(model, C[i, k] - C[i, k-1] >= sum(T[i, j] * Z[j, k] for j in Jobs))
        end
    end
    # ++++++++++ Objective function
    @objective(model, Min, C[m, n])
    println("\n============== SOLVING PFSP ALTERNATIVE ASSIGNMENT MODEL ==============\n")
    status = optimize!(model)
    Z_values = zeros(Float64, (n, n))
    if termination_status(model) == MOI.OPTIMAL
        println("Optimal Objective Function value: ", objective_value(model))
        permutation = zeros(Int64, n)
        for j in Jobs
            for i in Jobs
                Z_values[i, j] = value(Z[i, j])
                if Z_values[i, j] >= 0.9  # Fix rounding error problem
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

function solve_deterministic_alternative_assignment_model_2(m, n, T, time_limit)
    if solver == "Gurobi"
        model = Model(Gurobi.Optimizer)  # MIPGap=0.005
        MOI.set(model, MOI.RawOptimizerAttribute("TimeLimit"), time_limit)
    elseif solver == "CPLEX"
        model = Model(CPLEX.Optimizer)  # MIPGap=0.005
        MOI.set(model, MOI.RawOptimizerAttribute("CPX_PARAM_TILIM"), time_limit)
        MOI.set(model, MOI.RawOptimizerAttribute("CPX_PARAM_THREADS"), 16)
    else
      println("No solver defined")
    end
    Jobs = 1:n
    N = Jobs
    Machines = 1:m
    M = Machines
    println("Solving Beta PFSP MILP Model (Deterministic)")
    ## ++++++++++ Model variables: C and Z
    # # C_ik completion time of job k on machine i
    @variable(model, C[M, N])
    for i in M
        for k in N
            @constraint(model, C[i, k] >= 0)
        end
    end
    # Z_ij = 1, if job i is assigned to sequence position j, 0 otherwise
    #   x x x Originally: @variable(model, Z[N, N], Bin)
    @variable(model, Z[N, N], Bin)
    # X[i, k, l] = 1, if job i is part of the critical path in machine interval (k, l) (from machine k to machine l)
    @variable(model, X[N, M, M], Bin)
    ## ++++++++++ Model constraints
    # Assignement constraints
    for i in Jobs
        @constraint(model, sum(Z[i, j] for j in Jobs) == 1)
    end
    for j in Jobs
        @constraint(model, sum(Z[i, j] for i in Jobs) == 1)
    end
    for k in 1:n
        @constraint(model, C[1, k] >= sum(T[1, j] * Z[j, 1] for j in Jobs))
    end
    for k in 1:n
        for i in 2:m
            @constraint(model, C[i, k] - C[i-1, k] >= sum(T[i, j] * Z[j, k] for j in Jobs))
        end
    end
    for k in 2:n
        for i in 1:m
            @constraint(model, C[i, k] - C[i, k-1] >= sum(T[i, j] * Z[j, k] for j in Jobs))
        end
    end
    # ADDITIONAL Problem constraints
    @constraint(model, sum(X[i, 1, k2] for i in Jobs, k2 in 2:m) <= 1)      # caso extremo (k == 1)
    @constraint(model, sum(X[i, k1, m] for i in Jobs, k1 in 1:(m-1)) <= 1)  # caso extremo (k == m)
    for k in 2:(m-1)
        # CASE 1 - OK
        @constraint(model, sum(X[i, k1, k2] for i in Jobs, k1 in 1:(k-1), k2 in (k+1):m) <= 1)
        # CASE 2 (Mario)
        @constraint(model, sum(X[i, k1, k2] for i in Jobs, k1 in 1:(k-1), k2 in (k+1):m)
                # k eh maquina no final do intervalo : (k1, k)
                + sum(X[i, k1, k] for i in Jobs, k1 in 1:(k-1)) <= 1)
        # CASE 3 (Mario)
        @constraint(model, sum(X[i, k1, k2] for i in Jobs, k1 in 1:(k-1), k2 in (k+1):m)
                # k eh maquina no inicio do intervalo: (k, k2)
                + sum(X[i, k, k2] for i in Jobs, k2 in (k+1):m) <= 1)
    end
    # Nova restricao (VI): Caso 1 => ~Caso 4
    # If a machine k is *inside* an interval (k1, k2) of machines in the critical path, then this machine k
    # must have only one job i associated with it and it is *forbidden* to have any other jobs j (including job i),
    # associated in the form X[j, k, k] = 1 => job j in the machine interval (k, k)
    for j in Jobs, k in 2:m-1
        @constraint(model, X[j, k, k] <= 1 - sum(X[i, k1, k2] for i in Jobs, k1 in 1:(k-1), k2 in (k+1):m))
    end
    # Casos particulares
    # Given a job i, avoid intersection between the intervals (k, k) and [ (k1 < k, k) OR (k, k2 > k) ]
    for i in Jobs
        for k in 2:m
            @constraint(model, X[i, k, k] <= 1 - sum(X[i, k1, k] for k1 in 1:(k-1)))
        end
        for k in 1:m-1
            @constraint(model, X[i, k, k] <= 1 - sum(X[i, k, k2] for k2 in (k+1):m))
        end
    end
    # Each job i may belong to at most one machine interval in the form (k1, k2) where k1 <= k2
    for i in Jobs
        @constraint(model, sum(X[i, k1, k2] for k1 in M, k2 in k1:m) <= 1)
    end
    # For each machine interval (k, l) l cannot be less than k
    for i in Jobs
       for k in 2:m
           for l in 1:(k-1)
               @constraint(model, X[i, k, l] == 0)
            end
        end
    end
    # Every machine must be associated to one machine interval that involves changes between machines
    # k == 1
    @constraint(model, sum(X[i, 1, k2] for i in Jobs, k2 in 2:m) >= 1)
    for k in 2:m-1
        @constraint(model, sum(X[i, k1, k2] for i in Jobs, k1 in 1:(k-1), k2 in (k+1):m)
                + sum(X[i, k, k2] for i in Jobs, k2 in (k+1):m) + sum(X[i, k1, k] for i in Jobs, k1 in 1:(k-1)) >= 1)
    end
    # k == m
    @constraint(model, sum(X[i, k1, m] for i in Jobs, k1 in 1:(m-1)) >= 1)
    # the legth of the diagonal of the critical path must be equal to (m + n - 1)
    @constraint(model, sum(X[i, k1, k2] * (k2 - k1 + 1) for i in Jobs, k1 in M, k2 in k1:m) == m + n - 1)
    # Objective Function (from critical path variable X)
    @constraint(model, C[m, n] == sum(X[i, k1, k2] * sum(T[k, i] for k in k1:k2) for i in Jobs, k1 in M, k2 in k1:m))
    # ++++++++++ Objective function
    @objective(model, Min, C[m, n])
    println("\n============== SOLVING PFSP BETA MODEL ==============\n")
    status = optimize!(model)
    Z_values = zeros(Float64, (n, n))
    if termination_status(model) == MOI.OPTIMAL
        println("Optimal Objective Function value: ", objective_value(model))
        permutation = zeros(Int64, n)
        for j in Jobs
            for i in Jobs
                Z_values[i, j] = value(Z[i, j])
                if Z_values[i, j] >= 0.9  # Fix rounding error problem
                    if j <= n
                        permutation[j] = i
                    end
                end
            end
        end
        println("Permutation: $(permutation)")

        permutation = zeros(Int64, n)
        num_diag = 0
        j = 1
        for k1 in M, k2 in M
            for i in Jobs
                if value(X[i, k1, k2]) == 1
                    println("Job $(i) in machines $(k1)-$(k2);")
                    num_diag += (k2 - k1 + 1)
                    if j <= n
                        permutation[j] = i
                    end
                    j += 1
                end
            end
        end
        #num_diag = sum(getvalue(X[i, k, l]) * (l - k) for i in Jobs, k in M, l in k:m)
        println("num_diag: $(num_diag); expected = $(m + n - 1)")
        println("Permutation: $(permutation)")
        cmax = calculate_cmax(m, n, T, permutation)
        println("Calculated Makespan: $(cmax)")

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

function solve_deterministic_pfsp(m, n, T)
    time_limit = 7200
    cmax, permutation, time, is_opt, Z = solve_deterministic_wilson_model(m, n, T, time_limit)
    return cmax, permutation, Z
end

function calculate_cmax(m, n, p, pi, whole_matrix = false)
    C = zeros(Float64, (m, n))
    k = pi[1]
    C[1, 1] = p[1, k]
    for i in 2:m
        C[i, 1] = C[i - 1, 1] + p[i, k]
    end
    for j in 2:n
        k = pi[j]
        C[1, j] = C[1, j - 1] + p[1, k]
        for i in 2:m
            C[i, j] = max(C[i, j - 1], C[i - 1, j]) + p[i, k]
        end
    end
    #println("C = $(C)")
    if whole_matrix
        return C
    else
        return C[m, n]
    end
end

# The criterion under consideration is sum_i (w_i * C_i), which is the
# the sum of the weighted job completion times.
function calculate_wct(m, n, p, w, pi)
    C = zeros(Float64, (m, n))
    k = pi[1]
    C[1, 1] = p[1, k]
    for i in 2:m
        C[i, 1] = C[i - 1, 1] + p[i, k]
    end
    for j in 2:n
        k = pi[j]
        C[1, j] = C[1, j - 1] + p[1, k]
        for i in 2:m
            C[i, j] = max(C[i, j - 1], C[i - 1, j]) + p[i, k]
        end
    end
    wct = 0
    for j in 1:n  # C[m, j] is the final completion time of job j
        k = pi[j]
        wct += C[m, j] * w[k]
    end
    return wct
end

# Solves the TBA PFSP MILP Model (Deterministic)
# Tseng, F. T., & Stafford, E. F. (2008). New MILP models for the permutation flowshop problem.
# Journal of the Operational Research Society, 59(10), 1373–1386. https://doi.org/10.1057/palgrave.jors.2602455
function solve_deterministic_tba_model(m, n, T, time_limit)
    if solver == "Gurobi"
        model = Model(Gurobi.Optimizer)  # MIPGap=0.005
        MOI.set(model, MOI.RawOptimizerAttribute("TimeLimit"), time_limit)
    elseif solver == "CPLEX"
        model = Model(CPLEX.Optimizer)  # MIPGap=0.005
        MOI.set(model, MOI.RawOptimizerAttribute("CPX_PARAM_TILIM"), time_limit)
        MOI.set(model, MOI.RawOptimizerAttribute("CPX_PARAM_THREADS"), 16)
    else
      println("No solver defined")
    end
    Jobs = 1:n
    N = Jobs
    Machines = 1:m
    M = Machines
    println("Solving TBA PFSP MILP Model (Deterministic)")
    ## ++++++++++ Model variables: X and Z
    # X_rj : idle time on machine r before the start of job in sequence position j
    @variable(model, X[M, N] >= 0)
    # Z_ij = 1, if job i is assigned to sequence position j, 0 otherwise
    @variable(model, Z[N, N], Bin)
    ## ++++++++++ Model constraints
    # Assignement constraints - constraints (2)-(3)
    for i in Jobs
        @constraint(model, sum(Z[i, j] for j in Jobs) == 1)
    end
    for j in Jobs
        @constraint(model, sum(Z[i, j] for i in Jobs) == 1)
    end
    for j in 2:n  # Constraints (7)
        @constraint(model, X[1, j] == 0)
    end
    # Constraints (11)
    for r in 2:m
        for j in 2:n
            @constraint(model, sum(T[r - 1, i] * Z[i, 1] for i in 1:n)
              + sum(sum((T[r, i] - T[r-1, i]) * Z[i, q] for i in 1:n) for q in 1:(j-1))
              + sum((X[r, s] - X[r-1, s]) for s in 2:j)
              - sum(T[r-1, i] * Z[i, j] for i in 1:n) >= 0)
        end
    end
    # ++++++++++ Objective function
    # Makespan - Constraint (12)
    @objective(model, Min, sum(sum(T[p, i] * Z[i, 1] for i in 1:n) for p in 1:(m-1)) +
                sum(sum(T[m, i] * Z[i, q] for i in 1:n) for q in 1:n) + sum(X[m, s] for s in 2:n))
    # Mean Flow Time - Constraint (13)
    #@objective(model, Min, sum(sum(T[p, i] * Z[i, 1] for i in 1:n) for p in 1:(m-1))
    #        + ( sum((n-j+1) * sum(T[m, i] * Z[i, j] for i in 1:n) for j in 1:n)
    #            + sum((n-j+1) * X[m, j] for j in 2:n) )/n)
    # Flow time in sequence position j - Constraint (14)
    #for j in 1:n
    #    @constraint(model, F[j] == sum(sum(T[p, i] * Z[i, 1] for i in 1:n) for p in 1:(m-1))
    #         + sum(sum(T[m, i] * Z[i, q] for i in 1:n) for q in 1:j) + sum(X[m, s] for s in 1:j))
    #end
    println("\n============== SOLVING PFSP TBA MODEL ==============\n")
    status = optimize!(model)
    Z_values = zeros(Float64, (n, n))
    if termination_status(model) == MOI.OPTIMAL
        println("Optimal Objective Function value: ", objective_value(model))
        permutation = zeros(Int64, n)
        for j in Jobs
            for i in Jobs
                Z_values[i, j] = value(Z[i, j])
                if Z_values[i, j] >= 0.9  # Fix rounding error problem
                    if j <= n
                        permutation[j] = i
                    end
                end
            end
        end
        println("Permutation: $(permutation)")
        cmax = calculate_cmax(m, n, T, permutation)
        println("Makespan: $(cmax)")
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

# Solves the TS3 PFSP MILP Model (Deterministic)
function solve_deterministic_ts3_model(m, n, T, time_limit)
    if solver == "Gurobi"
        model = Model(Gurobi.Optimizer)  # MIPGap=0.005
        MOI.set(model, MOI.RawOptimizerAttribute("TimeLimit"), time_limit)
    elseif solver == "CPLEX"
        model = Model(CPLEX.Optimizer)  # MIPGap=0.005
        MOI.set(model, MOI.RawOptimizerAttribute("CPX_PARAM_TILIM"), time_limit)
        MOI.set(model, MOI.RawOptimizerAttribute("CPX_PARAM_THREADS"), 16)
    else
      println("No solver defined")
    end
    Jobs = 1:n
    N = Jobs
    Machines = 1:m
    M = Machines
    println("Solving TS3 PFSP MILP Model (Deterministic)")
    ## ++++++++++ Model variables: TT, Y and Z
    # Y_rj idle time of job in sequence position j after it finishes processing on machine r
    @variable(model, Y[M, N] >= 0)
    # Z_ij = 1, if job i is assigned to sequence position j, 0 otherwise
    #   x x x Originally: @variable(model, Z[N, N], Bin)
    @variable(model, Z[N, N], Bin)
    # The term TT rj is used to represent the processing time of the job in sequence position j on machine r. (1)
    #@variable(model, TT[M, N] >= 0)
    ## ++++++++++ Model constraints
    # Assignement constraints - Constraints (2)-(3)
    for i in Jobs
        @constraint(model, sum(Z[i, j] for j in Jobs) == 1)
    end
    for j in Jobs
        @constraint(model, sum(Z[i, j] for i in Jobs) == 1)
    end
    # Constraints (18)
    for r in 2:m
        for j in 2:n
            @constraint(model, sum(T[1, i] * Z[i, j-1] for i in 1:n) - sum(T[r, i] * Z[i, j-1] for i in 1:n)
                + sum(sum(T[q, i] * (Z[i, j] - Z[i, j-1]) for i in 1:n) for q in 1:(r - 1))
                + sum(Y[q, j] - Y[q, j-1] for q in 1:(r - 1)) >= 0)
        end
    end
    # Constraints (20)
    for r in 1:(m - 1)
        @constraint(model, Y[r, 1] == 0)
    end
    # The term TT rj is used to represent the processing time of the job in sequence position j on machine r. (1)
    #for r in 1:m
    #    for j in 1:n
    #        @constraint(model, TT[r, j] == sum(T[r, i] * Z[i, j] for i in 1:n))
    #    end
    #end
    # ++++++++++ Objective function
    # Makespan - Constraint (21)
    @objective(model, Min, sum(T[1, p] for p in 1:n) + sum(sum(T[q, i] * Z[i, n] for i in 1:n) for q in 2:m)
            + sum(Y[q, n] for q in 1:(m - 1)))
    # Mean Flow Time - Constraint (22)
    #@objective(model, Min, ( sum((n-j+1) * sum(T[1, i] * Z[i, j] for i in 1:n) for j in 1:n)
    #        + sum(sum(T[r, i] for r in 2:m) for i in 1:n) + sum(sum(Y[r, j] for r in 1:(m - 1)) for j in 1:n) ) / n)
    # Flow time in sequence position j - Constraint (23)
    #for j in 1:n
    #    @constraint(model, F[j] == sum(TT[1, p] for p in 1:j) + sum(TT[q, j] for q in 2:m) + sum(Y[q, j] for q in 1:(m - 1)))
    #end
    println("\n============== SOLVING PFSP TS3 MODEL ==============\n")
    status = optimize!(model)
    Z_values = zeros(Float64, (n, n))
    if termination_status(model) == MOI.OPTIMAL
        println("Optimal Objective Function value: ", objective_value(model))
        permutation = zeros(Int64, n)
        for j in Jobs
            for i in Jobs
                Z_values[i, j] = value(Z[i, j])
                if Z_values[i, j] >= 0.9  # Fix rounding error problem
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


function solve_deterministic_ts2_model(m, n, T, time_limit)
    if solver == "Gurobi"
        model = Model(Gurobi.Optimizer)  # MIPGap=0.005
        MOI.set(model, MOI.RawOptimizerAttribute("TimeLimit"), time_limit)
    elseif solver == "CPLEX"
        model = Model(CPLEX.Optimizer)  # MIPGap=0.005
        MOI.set(model, MOI.RawOptimizerAttribute("CPX_PARAM_TILIM"), time_limit)
        MOI.set(model, MOI.RawOptimizerAttribute("CPX_PARAM_THREADS"), 16)
    else
      println("No solver defined")
    end
    Jobs = 1:n
    N = Jobs
    Machines = 1:m
    M = Machines
    println("Solving TS2 PFSP MILP Model (Deterministic)")
    ## ++++++++++ Model variables: E and Z
    # # E_ik completion time of job k on machine i
    @variable(model, E[M, N] >= 0)
    # Z_ij = 1, if job i is assigned to sequence position j, 0 otherwise
    #   x x x Originally: @variable(model, Z[N, N], Bin)
    @variable(model, Z[N, N], Bin)
    ## ++++++++++ Model constraints
    # Assignement constraints
    for i in Jobs
        @constraint(model, sum(Z[i, j] for j in Jobs) == 1)
    end
    for j in Jobs
        @constraint(model, sum(Z[i, j] for i in Jobs) == 1)
    end
    for r in 2:m  # (A.17)
        for j in 1:(n-1)
            @constraint(model, E[r, j+1] >= E[r, j] + sum(T[r, i] * Z[i, j+1] for i in Jobs))
        end
    end
    for r in 1:(m-1)  # (A.18)
        for j in 2:n
            @constraint(model, E[r+1, j] >= E[r, j] + sum(T[r+1, i] * Z[i, j] for i in Jobs))
        end
    end
    for j in 1:(n-1)  # (A.19)
        @constraint(model, E[1, j+1] == E[1, j] + sum(T[1, i] * Z[i, j+1] for i in Jobs))
    end
    for r in 1:(m-1)  # (A.20)
        @constraint(model, E[r+1, 1] == E[r, 1] + sum(T[r+1, i] * Z[i, 1] for i in Jobs))
    end
    # (A.21)
    @constraint(model, E[1, 1] == sum(T[1, i] * Z[i, 1] for i in Jobs))
    # ++++++++++ Objective function - (A.22) Makespan
    @objective(model, Min, E[m, n])
    println("\n============== SOLVING PFSP TS2 MODEL ==============\n")
    status = optimize!(model)
    Z_values = zeros(Float64, (n, n))
    if termination_status(model) == MOI.OPTIMAL
        println("Optimal Objective Function value: ", objective_value(model))
        permutation = zeros(Int64, n)
        for j in Jobs
            for i in Jobs
                Z_values[i, j] = value(Z[i, j])
                if Z_values[i, j] >= 0.9  # Fix rounding error problem
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

function u(data, nbj)
    my_u = []
    for j in 1:nbj
        if data[1, j] <= data[2, j]
            push!(my_u, j)
        end
    end
    return sort(my_u, by=x->data[1, x], rev=false)
end

# p_ij should be a two dimension array
function v(data, nbj)
    my_v = []
    for j in 1:nbj
        if data[1, j] > data[2, j]
            push!(my_v, j)
        end
    end
    return sort(my_v, by=x->data[2, x], rev=true)
end

function johnson(p_ij, m, n)
    L, R = u(p_ij, n), v(p_ij, n)
    # println("L = $(L) ; R = $(R)")
    pi = vcat(L, R)
    return pi, calculate_cmax(m, n, p_ij, pi)
end

function calculate_johnson_given_scenario(m, n, P_bar, P_hat, m_dev; global_budget = false)
    if global_budget
        dev_machine = [[] for x in 1:m]
        for (r, i) in m_dev
            push!(dev_machine[r], i)
        end
        m_dev = dev_machine
    end
    for i in 1:m
        for j in m_dev[i]
            if j > 0
                P_bar[i, j] += P_hat[i, j]
            end
        end
    end
    pi, cmax = johnson(P_bar, m, n)
    for i in 1:m
        for j in m_dev[i]
            if j > 0
                P_bar[i, j] -= P_hat[i, j]
            end
        end
    end
    return pi, cmax
end

function calculate_critical_path(m, n, P, C, pi)
	return calculate_critical_path_recursive(m, n, P, C, pi)
end

function calculate_critical_path_recursive(r, i, P, C, pi)
	if i == 1 && r == 1
		return [(1, pi[1])]
	end
	if (r >= 2) && (abs( (C[r, i] - P[r, pi[i]]) - C[r - 1, i]) < 1e-5)
		return [calculate_critical_path_recursive(r - 1, i, P, C, pi) ; [ (r, pi[i]) ] ]
	end
    if (i >= 2) && (abs( (C[r, i] - P[r, pi[i]]) - C[r, i - 1]) < 1e-5)
		return [calculate_critical_path_recursive(r, i - 1, P, C, pi) ; [ (r, pi[i]) ] ]
	end
	println("ERROR !!! r = $r, i = $i.")
end

function calculate_P_given_scenario(m, n, P_bar, P_hat, m_dev; global_budget = false)
    if global_budget
        dev_machine = [[] for x in 1:m]
        for (r, i) in m_dev
            push!(dev_machine[r], i)
        end
        m_dev = dev_machine
    end
	P = deepcopy(P_bar)
	for i in 1:m
        for j in m_dev[i]
            if j > 0
                P[i, j] += P_hat[i, j]
            end
        end
    end
	return P
end

function solve_deterministic_wagner_wst2_model(m, n, T, time_limit)
    if solver == "Gurobi"
        model = Model(Gurobi.Optimizer)  # MIPGap=0.005
        MOI.set(model, MOI.RawOptimizerAttribute("TimeLimit"), time_limit)
    elseif solver == "CPLEX"
        model = Model(CPLEX.Optimizer)  # MIPGap=0.005
        MOI.set(model, MOI.RawOptimizerAttribute("CPX_PARAM_TILIM"), time_limit)
        MOI.set(model, MOI.RawOptimizerAttribute("CPX_PARAM_THREADS"), 16)
    else
      println("No solver defined")
    end
    Jobs = 1:n
    N = Jobs
    Machines = 1:m
    M = Machines
    println("Solving Wagner WST2 PFSP MILP Model (Deterministic)")
    ## ++++++++++ Model variables: X, Y and Z
    # X_rj : idle time on machine r before the start of job in sequence position j
    @variable(model, X[M, N] >= 0)
    # Y_rj idle time of job in sequence position j after it finishes processing on machine r
    @variable(model, Y[M, N] >= 0)
    # Z_ij = 1, if job i is assigned to sequence position j, 0 otherwise
    @variable(model, Z[N, N], Bin)
    # Cmax : maximum flowtime (makespan) of the schedule determined by the completion time of
    #        the job in the last sequence position on the last machine.
    @variable(model, Cmax >= 0)
    ## ++++++++++ Model constraints
    # Assignement constraints
    for i in Jobs
        @constraint(model, sum(Z[i, j] for j in Jobs) == 1)
    end
    for j in Jobs
        @constraint(model, sum(Z[i, j] for i in Jobs) == 1)
    end
    for r in 2:m-1  # (A.1)
        for j in 2:n-1
            @constraint(model, sum(T[r, i] * Z[i, j+1] for i in Jobs) - sum(T[r+1, i] * Z[i, j] for i in Jobs)
                + X[r, j+1] - X[r+1, j+1] + Y[r, j+1] - Y[r, j] == 0)
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
    @constraint(model, Cmax == (sum(T[m, i] for i in Jobs) + sum(X[m, p] for p in Jobs)))
    # ++++++++++ Objective function
    @objective(model, Min, Cmax )
    println("\n============== SOLVING WAGNER WST2 MODEL ==============\n")
    status = optimize!(model)
    Z_values = zeros(Float64, (n, n))
    if termination_status(model) == MOI.OPTIMAL
        println("Optimal Objective Function value: ", objective_value(model))
        permutation = zeros(Int64, n)
        for j in Jobs
            for i in Jobs
                Z_values[i, j] = value(Z[i, j])
                if Z_values[i, j] >= 0.9  # Fix rounding error problem
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

function solve_deterministic_liao_you_model(m, n, T, time_limit)
    if solver == "Gurobi"
        model = Model(Gurobi.Optimizer)  # MIPGap=0.005
        MOI.set(model, MOI.RawOptimizerAttribute("TimeLimit"), time_limit)
    elseif solver == "CPLEX"
        model = Model(CPLEX.Optimizer)  # MIPGap=0.005
        MOI.set(model, MOI.RawOptimizerAttribute("CPX_PARAM_TILIM"), time_limit)
        #MOI.set(model, MOI.RawOptimizerAttribute("CPX_PARAM_THREADS"), 16)
    else
      println("No solver defined")
    end
    Jobs = 1:n
    N = Jobs
    Machines = 1:m
    M = Machines
    println("Solving Liao-You PFSP MILP Model (Deterministic)")
    # P is a very large constant (big-M)
    random_pi = [x for x in 1:n]
    P = calculate_cmax(m, n, T, random_pi) + 1
    ## ++++++++++ Model variables: Cmax, C and D
    @variable(model, Cmax >= 0)
    # S_ri start time of job i on machine r
    @variable(model, S[M, N] >= 0)
    # Surplus variables q[r, i, k]
    @variable(model, q[M, N, N] >= 0)
    # D_ik = 1, if job i is scheduled any time before job k; 0, otherwise (i < k)
    @variable(model, D[N, N], Bin)
    ## ++++++++++ Model constraints
    for r in 1:(m-1)
        for i in 1:n
            @constraint(model, S[r, i] + T[r, i] <= S[r+1, i])
        end
    end
    for r in 1:m
        for i in 1:(n-1)
            for k in (i+1):n
                @constraint(model, S[r, i] - S[r, k] + P * D[i, k] - T[r, k] == q[r, i, k])
                @constraint(model, P - T[r, i] - T[r, k] >= q[r, i, k])
            end
        end
    end
    for i in 1:n
        @constraint(model, Cmax >= S[m, i] + T[m, i])
    end
    # ++++++++++ Objective function
    @objective(model, Min, Cmax)
    println("\n============== SOLVING PFSP LIAO-YOU MODEL ==============\n")
    status = optimize!(model)
    D_star = zeros(Float64, (n, n))
    if termination_status(model) == MOI.OPTIMAL
        println("Optimal Objective Function value: ", objective_value(model))
        for i in Jobs  # Job i
            for k in Jobs  # Job k
                if value(D[i, k]) >= 0.9
                    D_star[i, k] = 1
                end
            end
        end
        permutation = topological_sort_by_dfs(D_star, n)
        println("Permutation: $(permutation)")
        return objective_value(model), permutation, solve_time(model), true, D_star, objective_bound(model)
    elseif termination_status(model) == MOI.TIME_LIMIT && has_values(model)
        println("Optimal solution not found!")
        println("Suboptimal Objective Function value: ", objective_value(model))
        for i in Jobs  # Job i
            for k in Jobs  # Job k
                if value(D[i, k]) >= 0.9
                    D_star[i, k] = 1
                end
            end
        end
        permutation = topological_sort_by_dfs(D_star, n)
        println("Permutation: $(permutation)")
        return objective_value(model), permutation, solve_time(model), false, D_star, objective_bound(model)
    else
        println("Optimal solution not found! Best bound: ", objective_bound(model))
        return objective_bound(model), zeros(Int64, n), solve_time(model), false, D_star, objective_bound(model)
    end
end

function calculate_P_given_deviation_matrix(m, n, P_bar, P_hat, m_dev)
	P = deepcopy(P_bar)
	for i in 1:m
        for j in 1:n
            P[i, j] += P_hat[i, j] * m_dev[i, j]
        end
    end
	return P
end

function solve_deterministic_wagner_2machine_model(m, n, T, time_limit)
    if solver == "Gurobi"
        model = Model(Gurobi.Optimizer)  # MIPGap=0.005
        MOI.set(model, MOI.RawOptimizerAttribute("TimeLimit"), time_limit)
    elseif solver == "CPLEX"
        model = Model(CPLEX.Optimizer)  # MIPGap=0.005
        MOI.set(model, MOI.RawOptimizerAttribute("CPX_PARAM_TILIM"), time_limit)
        MOI.set(model, MOI.RawOptimizerAttribute("CPX_PARAM_THREADS"), 2)
    else
      println("No solver defined")
    end
    Jobs = 1:n
    N = Jobs
    m = 2
    Machines = 1:m
    M = Machines
    println("Solving Wagner 2-machine PFSP MILP Model (Deterministic)")
    ## ++++++++++ Model variables: X, Y and Z
    # X_j : idle time on machine 2 before the start of job in sequence position j
    @variable(model, X[N] >= 0)
    # Y_rj idle time of job in sequence position j after it finishes processing on machine 1
    @variable(model, Y[N] >= 0)
    # Z_ij = 1, if job i is assigned to sequence position j, 0 otherwise
    @variable(model, Z[N, N], Bin)
    # Cmax : maximum flowtime (makespan) of the schedule determined by the completion time of
    #        the job in the last sequence position on the last machine.
    @variable(model, Cmax >= 0)
    ## ++++++++++ Model constraints
    # Assignement constraints
    for i in Jobs
        @constraint(model, sum(Z[i, j] for j in Jobs) == 1)
    end
    for j in Jobs
        @constraint(model, sum(Z[i, j] for i in Jobs) == 1)
    end
    for j in 1:n-1
        @constraint(model, sum(T[1, i] * Z[i, j+1] for i in Jobs) - sum(T[2, i] * Z[i, j] for i in Jobs) - X[j+1] + Y[j+1] - Y[j] == 0)
    end
    @constraint(model, sum(T[1, i] * Z[i, 1] for i in Jobs) - X[1] == 0)
    @constraint(model, Y[1] == 0)
    #for j in 1:n  # Constraints (7) from TBA model
    #    @constraint(model, X[1, j] == 0)
    #end
    # Constraints (20) from TS3 model
    #for r in 1:(m - 1)
    #    @constraint(model, Y[r, 1] == 0)
    #end
    #for j in 1:n  # Constraints derived from X[1, j] == 0
    #    @constraint(model, Y[1, j] == 0)
    #end
    @constraint(model, Cmax == (sum(T[m, i] for i in Jobs) + sum(X[p] for p in Jobs)))
    # ++++++++++ Objective function
    @objective(model, Min, Cmax )
    println("\n============== SOLVING 2-MACHINE WAGNER MODEL ==============\n")
    status = optimize!(model)
    Z_values = zeros(Float64, (n, n))
    if termination_status(model) == MOI.OPTIMAL
        println("Optimal Objective Function value: ", objective_value(model))
        permutation = zeros(Int64, n)
        for j in Jobs
            for i in Jobs
                Z_values[i, j] = value(Z[i, j])
                if Z_values[i, j] >= 0.9  # Fix rounding error problem
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
