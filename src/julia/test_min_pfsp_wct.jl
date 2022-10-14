using JuMP
#using Cbc
using CSV
using DataFrames
using DataArrays
using Gurobi
using StatsBase

function read_input_file(filepath)
    file = open(filepath)
    text = ""
    count = 1
    m = 0
    n = 0
    for ln in eachline(file)
        if length(ln) > 0 && (! contains(ln, "!"))  # ignore comments and empty lines
            if count == 1
                input = split(ln, r" |\t")
                n = parse(Int, input[1])
                m = parse(Int, input[2])
            else
                ln = lstrip(ln)
                #println("$(count): $(ln)")
                if length(ln) > 0
                    text = text * ln
                    #if ! contains(ln, r"(\w+)\:")
                    text = text * "\n"
                    #end
                end
            end
            #println("$(ln)")
        end
        count += 1
    end
    println("M = $(m), N = $(n)")
    count = 0
    # T = {Trj} is the M × N matrix of job processing times, with T_ri = processing time of job i on machine r.
    T = readdlm(IOBuffer(text))
    #println("Matrix T: $(T)")
    println("Instance file read OK.")
    return m, n, T
end

# Solves a relaxation of Wagner PFSP MILP Model (LB), where, for all jobs j in S, Z_ij = {0, 1}
function solve_wagner_model_fix_S(inputfile, time_limit, pre_solve, S = [])
    m, n, T = read_input_file(inputfile)
    model = Model(solver=GurobiSolver(Presolve=pre_solve,TimeLimit=time_limit))
    Jobs = 1:n
    N = Jobs
    Machines = 1:m
    M = Machines
    if isempty(S)
        S = Jobs
        println("Fixing all Z variables to binary.")
    end
    # TODO Up to now, we are using Taillard instance files, which DO NOT contain job weights
    # For now, job weights will be fixed to 1.0 for every job i [1, ..., 1]
    W = [1.0 for i in Jobs]
    println("Job weights: $(W)")
    ## ++++++++++ Model variables: X, Y and Z
    # X_rj : idle time on machine r before the start of job in sequence position j
    @variable(model, X[M, N] >= 0)
    # Y_rj idle time of job in sequence position j after it finishes processing on machine r
    @variable(model, Y[M, N] >= 0)
    # Z_ij = 1, if job i is assigned to sequence position j, 0 otherwise
    #   x x x Originally: @variable(model, Z[N, N], Bin)
    @variable(model, Z[N, N])
    for i in Jobs
        for j in Jobs
            if j in S
                setcategory(Z[i, j], :Bin)
            else
                setlowerbound(Z[i, j], 0)
                setupperbound(Z[i, j], 1)
                #@variable(model, 0 <= Z[N, j] <= 1)
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
    for r in 1:m-1
        for j in 1:n-1
            @constraint(model, sum(T[r, i] * Z[i, j+1] for i in Jobs) - sum(T[r+1, i] * Z[i, j] for i in Jobs) + X[r, j+1] - X[r+1, j+1] + Y[r, j+1] - Y[r, j] == 0)
        end
    end
    for r in 1:m-1
        @constraint(model, sum(T[r, i] * Z[i, 1] for i in Jobs) + X[r, 1] - X[r+1, 1] + Y[r, 1] == 0)
    end
    # ++++++++++ Objective function
    # C[j] * W[pi[j]] where C[j] is the completion time of job in sequence position j on the last machine (C[m, j])
    @variable(model, C[N] >= 0)
    for j in Jobs  # job i in sequence position j
        # sum(sum(T[m, pi[p]] for p in 1:j) + sum(X[m, p] for p in 1:j)
        @constraint(model, C[j] == sum(sum(T[m, i] * Z[i, p] for i in Jobs) for p in 1:j) + sum(X[m, p] for p in 1:j))
    end
    @variable(model, D[N] >= 0)
    for i in Jobs  # job i in sequence position j
        @constraint(model, D[i] == sum(T[m, i] * Z[i, p] for i in Jobs))
    end
    #@objective(model, Min, sum(C[j] * sum(W[i] * Z[i, j] for i in Jobs) for j in Jobs))
    @objective(model, Min, sum(D[i] * W[i] for j in Jobs))
    # Cmax : maximum flowtime (makespan) of the schedule determined by the completion time of
    #        the job in the last sequence position on the last machine.
    #@variable(model, Cmax >= 0)
    #@constraint(model, Cmax == (sum(T[m, i] for i in Jobs) + sum(X[m, p] for p in Jobs)))
    #@objective(model, Min, Cmax )
    println("\n============== SOLVING WAGNER MODEL ==============\n")
    status = solve(model)
    if status == :Optimal
        println("Optimal Objective Function value: ", getobjectivevalue(model))
        permutation = zeros(Int64, n)
        for j in Jobs
            for i in Jobs
                if getvalue(Z[i, j]) == 1
                    if j <= n
                        permutation[j] = i
                    end
                end
            end
        end
        println("Permutation: $(permutation)")
        return getobjectivevalue(model), permutation, getsolvetime(model), (status == :Optimal ? true : false)
    else
        println("Optimal solution not found! Best bound: ", getobjectivebound(model))
        return getobjectivebound(model), zeros(Int64, n), getsolvetime(model), (status == :Optimal ? true : false)
    end
end

# Solves a relaxation of Wilson PFSP MILP Model (LB), where, for all jobs j in S, Z_ij = {0, 1}
function solve_wilson_model_fix_S(inputfile, time_limit, pre_solve, S = [])
    m, n, T = read_input_file(inputfile)
    model = Model(solver=GurobiSolver(Presolve=pre_solve,TimeLimit=time_limit))
    Jobs = 1:n
    N = Jobs
    Machines = 1:m
    M = Machines
    if isempty(S)
        println("Fixing all Z variables to binary.")
        S = Jobs
    end
    # TODO Up to now, we are using Taillard instance files, which DO NOT contain job weights
    # For now, job weights will be fixed to 1.0 for every job i [1, ..., 1]
    W = [1.0 for i in Jobs]
    println("Job weights: $(W)")
    ## ++++++++++ Model variables: B and Z
    # B_rj : start (begin) time of job in sequence position j on machine r
    @variable(model, B[M, N] >= 0)
    # Z_ij = 1, if job i is assigned to sequence position j, 0 otherwise
    #   x x x Originally: @variable(model, Z[N, N], Bin)
    @variable(model, Z[N, N])
    for i in Jobs
        for j in Jobs
            if j in S
                setcategory(Z[i, j], :Bin)
            else
                setlowerbound(Z[i, j], 0)
                setupperbound(Z[i, j], 1)
                #@variable(model, 0 <= Z[N, j] <= 1)
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
    # C[j] * W[pi[j]] where C[j] is the completion time of job in sequence position j on the last machine (C[m, j])
    @variable(model, C[N] >= 0)
    for j in Jobs  # job in sequence position j
        # C[j] == B[m, j] + sum(T[m, pi[p]] for p in 1:j)
        @constraint(model, C[j] == B[m, j] + sum(sum(T[m, i] * Z[i, p] for i in Jobs) for p in 1:j))
    end
    #@objective(model, Min, sum(C[j] * sum(W[i] * Z[i, j] for i in Jobs) for j in Jobs))
    obj = @objective(model, Min, sum(C[j] * sum(W[i] * Z[i, j] for i in Jobs) for j in n:n))
    println("FO: $(obj)")
    #@objective(model, Min, B[m, n] + sum(T[m, i] * Z[i, n] for i in Jobs))
    println("\n============== SOLVING WILSON MODEL ==============\n")
    status = solve(model)
    if status == :Optimal
        println("Optimal Objective Function value: ", getobjectivevalue(model))
        permutation = zeros(Int64, n)
        for j in Jobs
            for i in Jobs
                if getvalue(Z[i, j]) == 1
                    if j <= n
                        permutation[j] = i
                    end
                end
            end
        end
        println("Permutation: $(permutation)")
        return getobjectivevalue(model), permutation, getsolvetime(model), (status == :Optimal ? true : false)
    else
        println("Optimal solution not found! Best bound: ", getobjectivebound(model))
        return getobjectivebound(model), zeros(Int64, n), getsolvetime(model), (status == :Optimal ? true : false)
    end
end

function process_instance_file(filepath, instance_name, time_limit, test_name, pre_solve)
    df = DataFrame(Instance = String[], TestName = String[], ModelName = String[],
        Value = Float64[], Time = Float64[], IsOptimal = Bool[], Permutation = String[])
    println("Processing instance file $(filepath)...")
    cmax, permutation, time, is_opt = solve_wilson_model_fix_S(filepath, time_limit, pre_solve)
    push!(df, [instance_name, test_name, "Wilson", cmax, time, is_opt, string(permutation)])
    cmax, permutation, time, is_opt = solve_wagner_model_fix_S(filepath, time_limit, pre_solve)
    push!(df, [instance_name, test_name, "Wagner", cmax, time, is_opt, string(permutation)])

    println("\n===================== RESULTS =====================\n")
    showall(df)
    println("")
    return df
end

# ["20x5", "20x10", "20x20", "50x5", "50x10", "50x20", "100x5", "100x10", "100x20", "200x10", "200x20", "500x20"]
time_limit = 3600
for test_set in ["20x5"] #, "20x10", "20x20", "50x5", "50x10", "50x20"]
    folder = "../../instances/taillard/output/" * test_set
    test_name = "Tail_" * test_set
    df = DataFrame(Instance = String[], TestName = String[], ModelName = String[],
        Value = Float64[], Time = Float64[], IsOptimal = Bool[], Permutation = String[])
    pre_solve = 1
    #if test_set in ["200x20", "500x20"]  # For large instances, disable pre_solve to avoid out of memory error
    #        pre_solve = 0
    #end
    for afile in filter!(f -> endswith(f, "in"), readdir(folder))
        df_file = process_instance_file(joinpath(normpath(folder), afile), afile, time_limit, test_name, pre_solve)
        df = vcat(df, df_file)
    end
    # Create an output file to write the dataframe to disk
    output_file = test_name * "_Min_PFSP_WCT_Test.csv"
    println("\nSaving to results file...")
    CSV.write(output_file, df)
end
