# ===================================================================================================
# deterministic_pfsp_alt_wct_mip.jl
# Deterministic Permutation Flow Shop Problem to Minimize the Weighted Sum of Job Completion Times (WCT).
# *** MILP approach based on the Alternative PFSP Cmax model.
# ===================================================================================================

using JuMP
using CPLEX
#using Gurobi
include("../config.jl")
include("../pfsp_file_reader.jl")
include("../deterministic_pfsp.jl")

function test_wct_15_5_instance()
    # basedir = joinpath(robust_ying_instances_folder, "rob-pfsp-wct/20jobs")  filename == "RB0205010_20_2_50_wct_inputs.txt"
    #basedir = joinpath(robust_ying_instances_folder, "rob-pfsp-wct", "10x5", "50%")
    basedir = joinpath(robust_ying_instances_folder, "rob-pfsp-wct", "15x5")
    time_limit = 4000
    for filename in searchdir(basedir, ".txt")
        #if filename == "RB0105007_10_5_50_wct_inputs.txt"
        if filename == "RB0155007_15_5_50_wct_inputs.txt"
            println("Instance: $(filename)")
            inputfile = joinpath("$(basedir)", "$(filename)")
            m, n, P_bar, P_hat, w = read_robust_input_file(inputfile)
            println("w: $(w)")
            solve_deterministic_manne2_wct_model(m, n, P_bar, w, time_limit; solver="Gurobi")
            solve_deterministic_manne_wct_model(m, n, P_bar, w, time_limit; solver="Gurobi")

        end
    end
end

function solve_deterministic_pfsp_wct(m, n, T, w, time_limit = 1800)
    cmax, permutation, time, is_opt, Z = solve_deterministic_alt_model_wct(m, n, T, w, time_limit)
    return cmax, permutation, Z
end

function solve_deterministic_manne_wct_model(m, n, T, w, time_limit; solver="CPLEX")
    if solver == "Gurobi"
        model = Model(Gurobi.Optimizer)  # MIPGap=0.005
        MOI.set(model, MOI.RawOptimizerAttribute("TimeLimit"), time_limit)
    elseif solver == "CPLEX"
        model = Model(CPLEX.Optimizer)  # MIPGap=0.005
        MOI.set(model, MOI.RawOptimizerAttribute("CPX_PARAM_TILIM"), time_limit)
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
    ## ++++++++++ Model variables: C and D
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
        for k in 1:n
            # triangle inequalities, see: Grotschel, Junger and Reinelt (1984, 1985).
            if i < k  # https://www2.isye.gatech.edu/people/faculty/Martin_Savelsbergh/publications/co_nftp.pdf
                @constraint(model, D[i, k] + D[k, i] == 1)
            end
            for j in 1:n
                if (i != k) && (i != j) && (j != k)
                    @constraint(model, D[i, j] + D[j, k] + D[k, i] <= 2)
                end
            end
        end
    end
    # ++++++++++ Objective function
    @variable(model, obj >= 0)
    @objective(model, Min, obj)
    @constraint(model, sum(w[i] * C[m, i] for i in Jobs) <= obj)
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

# Solves the Alternative PFSP MILP Model to minimize the weighted sum of completion times (Deterministic)
# Based on Manne PFSP Model
function solve_deterministic_manne2_wct_model(m, n, T, w, time_limit; solver="CPLEX")
    if solver == "Gurobi"
        model = Model(Gurobi.Optimizer)  # MIPGap=0.005
        MOI.set(model, MOI.RawOptimizerAttribute("TimeLimit"), time_limit)
    elseif solver == "CPLEX"
        model = Model(CPLEX.Optimizer)  # MIPGap=0.005
        MOI.set(model, MOI.RawOptimizerAttribute("CPX_PARAM_TILIM"), time_limit)
    else
      println("No solver defined")
    end
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
    # # C2_ik completion time of job in sequence position k on machine i
    @variable(model, Cb[M, N] >= 0)
    # D_ik = 1, if job i is scheduled any time before job k; 0, otherwise;
    @variable(model, D[N, N], Bin)
    # Z_ij = 1, if job i is assigned to sequence position j, 0 otherwise
    #   x x x Originally: @variable(model, Z[N, N], Bin)
    @variable(model, Z[N, N], Bin)
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
    ## ++++++++++ Model constraints
    # Assignement constraints
    for i in Jobs
        @constraint(model, sum(Z[i, j] for j in Jobs) == 1)
    end
    for j in Jobs
        @constraint(model, sum(Z[i, j] for i in Jobs) == 1)
    end
    for k in 1:n
        @constraint(model, Cb[1, k] >= sum(T[1, j] * Z[j, 1] for j in Jobs))
    end
    for k in 1:n  # (A.18)
        for i in 2:m
            @constraint(model, Cb[i, k] - Cb[i-1, k] >= sum(T[i, j] * Z[j, k] for j in Jobs))
        end
    end
    for k in 2:n  # (A.17)
        for i in 1:m
            @constraint(model, Cb[i, k] - Cb[i, k-1] >= sum(T[i, j] * Z[j, k] for j in Jobs))
        end
    end
    for i in 1:n  # Disjuctive constraints used to determine the start time of job i on the last machine m
        for j in 1:n
            @constraint(model, C[m, i] >= Cb[m, j] - P * (1 - Z[i, j]) )
        end
    end
    # TODO Create a constraint associating decision variables Z and D
    for i in 1:n
        for k in 1:n
            if i != k
                @constraint(model, sum(j2 * Z[k, j2] for j2 in 1:n) - sum(j1 * Z[i, j1] for j1 in 1:n) >= D[i, k] - (n - 1) * (1 - D[i, k]))
                #@constraint(model, sum(j1 * Z[i, j1] for j1 in 1:n) - sum(j2 * Z[k, j2] for j2 in 1:n) >= D[k, i] - (n - 1) * (1 - D[k, i]))
            end
            # triangle inequalities, see: Grotschel, Junger and Reinelt (1984, 1985).
            if i < k  # https://www2.isye.gatech.edu/people/faculty/Martin_Savelsbergh/publications/co_nftp.pdf
                @constraint(model, D[i, k] + D[k, i] == 1)
            end
            for j in 1:n
                if (i != k) && (i != j) && (j != k)
                    @constraint(model, D[i, j] + D[j, k] + D[k, i] <= 2)
                end
            end
        end
    end

    # ++++++++++ Objective function
    @objective(model, Min, sum(w[j] * C[m, j] for j in Jobs))
    println("\n============== SOLVING PFSP MANNE+ALTERNATIVE MODEL - WCT ==============\n")
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

#test_wct_15_5_instance()
