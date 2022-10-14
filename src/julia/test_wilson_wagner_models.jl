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

# Solves a relaxation of Wagner PFSP MILP Model (LB), where, for all jobs i in S, Z_ij = {0, 1}
function solve_wagner_model(inputfile, time_limit, relax = false, S = [])
    m, n, T = read_input_file(inputfile)
    model = Model(solver=GurobiSolver(TimeLimit=time_limit))
    Jobs = 1:n
    N = Jobs
    Machines = 1:m
    M = Machines
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
            if relax
                if i in S
                    setcategory(Z[i, j], :Bin)
                else
                    setlowerbound(Z[i, j], 0)
                    setupperbound(Z[i, j], 1)
                    #@variable(model, 0 <= Z[N, j] <= 1)
                end
            else
                setcategory(Z[i, j], :Bin)
            end
        end
    end
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
    status = solve(model)
    if status == :Optimal
        println("Optimal Objective Function value: ", getobjectivevalue(model))
        return getobjectivevalue(model), getsolvetime(model), (status == :Optimal ? true : false)
    else
        println("Optimal solution not found! Best bound: ", getobjectivebound(model))
        return getobjectivebound(model), getsolvetime(model), (status == :Optimal ? true : false)
    end
end

# Solves a relaxation of Wilson PFSP MILP Model (LB), where, for all jobs i in S, Z_ij = {0, 1}
function solve_wilson_model(inputfile, time_limit, relax = false, S = [])
    m, n, T = read_input_file(inputfile)
    model = Model(solver=GurobiSolver(TimeLimit=time_limit))
    Jobs = 1:n
    N = Jobs
    Machines = 1:m
    M = Machines
    ## ++++++++++ Model variables: B and Z
    # B_rj : start (begin) time of job in sequence position j on machine r
    @variable(model, B[M, N] >= 0)
    # Z_ij = 1, if job i is assigned to sequence position j, 0 otherwise
    #   x x x Originally: @variable(model, Z[N, N], Bin)
    @variable(model, Z[N, N])
    for i in Jobs
        for j in Jobs
            if relax
                if i in S
                    setcategory(Z[i, j], :Bin)
                else
                    setlowerbound(Z[i, j], 0)
                    setupperbound(Z[i, j], 1)
                    #@variable(model, 0 <= Z[N, j] <= 1)
                end
            else
                setcategory(Z[i, j], :Bin)
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
    @objective(model, Min, B[m, n] + sum(T[m, i] * Z[i, n] for i in Jobs))
    println("\n============== SOLVING WILSON MODEL ==============\n")
    status = solve(model)
    if status == :Optimal
        println("Optimal Objective Function value: ", getobjectivevalue(model))
        return getobjectivevalue(model), getsolvetime(model), (status == :Optimal ? true : false)
    else
        println("Optimal solution not found! Best bound: ", getobjectivebound(model))
        return getobjectivebound(model), getsolvetime(model), (status == :Optimal ? true : false)
    end
end

function process_instance_folder_LB(folder, time_limit, test_name)
    df_LB = DataFrame(Instance = String[], Measure = String[], k = Int64[], Criteria = String[],
        LB_Type = String[], LB_value = Float64[], Time = Float64[], IsOptimal = Bool[])
    # Create an output file to write the dataframe to disk
    output_file = test_name * "_PFSP_LB_Pure_Wilson_Wagner_Test.csv"
    for filename in readdir(folder)
        # Model parameters (input file)
        filepath = folder * "/" * filename
        println("Processing instance file $(filepath)")

        m, n, T = read_input_file(filepath)
        # Transpoe a matriz T, colocando os jobs nas linhas ao inves das colunas
        T_t = transpose(T)
        # Converte a matriz de tempos de processamento (m x n) para um dataframe
        df_T = convert(DataFrame, T_t)
        names!(df_T, [Symbol("m$i") for i in 1:m])

        # Objetivo: calcular as medidas estatisticas dos tempos de processamento de cada job j,
        #           considerando todas as m maquinas
        #df_T
        #describe(df_T)
        df_stats = DataFrame(j = Int64[], Mean = Float64[], Median = Float64[], Sum = Float64[], Std = Float64[])
        map((row, j) -> push!(df_stats, [
            j,
            mean(convert(Array, row)),
            median(convert(Array, row)),
            sum(convert(Array, row)),
            std(convert(Array, row))]),

            eachrow(df_T), [j for j in 1:n]);

        # Invocando os modelos Wagner e Wilson e armazenando os valores obtidos de Lower Bound
        S = []
        LB_wagner, time_wagner, opt_wagner = solve_wagner_model(filepath, time_limit, true, S)
        LB_wilson, time_wilson, opt_wilson = solve_wilson_model(filepath, time_limit, true, S)
        push!(df_LB, [filename, "None", 0, "None", "Wagner", LB_wagner, time_wagner, opt_wagner])
        push!(df_LB, [filename, "None", 0, "None", "Wilson", LB_wilson, time_wilson, opt_wilson])
        println("Saving to results file...")
        CSV.write(output_file, df_LB)

    end
    println("\n===================== RESULTADOS =====================\n")
    showall(df_LB)
    return df_LB
end

function process_instance_folder_opt(folder, time_limit, test_name)
    df = DataFrame(Instance = String[], Model_Type = String[], Value = Float64[], Time = Float64[], IsOptimal = Bool[])
    # Create an output file to write the dataframe to disk
    output_file = test_name * "_PFSP_Opt_Wilson_Wagner_Test.csv"
    for filename in readdir(folder)
        # Model parameters (input file)
        filepath = folder * "/" * filename
        println("Processing instance file $(filepath)")

        m, n, T = read_input_file(filepath)
        # Transpoe a matriz T, colocando os jobs nas linhas ao inves das colunas
        T_t = transpose(T)
        # Converte a matriz de tempos de processamento (m x n) para um dataframe
        df_T = convert(DataFrame, T_t)
        names!(df_T, [Symbol("m$i") for i in 1:m])

        # Invocando os modelos Wagner e Wilson e armazenando os valores obtidos de solução
        wagner, time_wagner, opt_wagner = solve_wagner_model(filepath, time_limit)
        wilson, time_wilson, opt_wilson = solve_wilson_model(filepath, time_limit)
        push!(df, [filename, "Wagner", wagner, time_wagner, opt_wagner])
        push!(df, [filename, "Wilson", wilson, time_wilson, opt_wilson])
        println("Saving to results file...")
        CSV.write(output_file, df)

    end
    println("\n===================== RESULTADOS =====================\n")
    showall(df)
    return df
end

function process_instance_file(filepath, time_limit, test_name, start_k, final_k, measure_list = [:Median, :Std])
    df_LB = DataFrame(Instance = String[], Measure = String[], k = Int64[], Criteria = String[],
        LB_Type = String[], LB_value = Float64[], Time = Float64[], IsOptimal = Bool[], Job_indices = String[])
    # Create an output file to write the dataframe to disk
    output_file = test_name * "_PFSP_LB_Wilson_Wagner_Test.csv"
    for k in start_k:final_k
        println("Processing instance file $(filepath) for k = $(k)")

        m, n, T = read_input_file(filepath)
        # Transpoe a matriz T, colocando os jobs nas linhas ao inves das colunas
        T_t = transpose(T)
        # Converte a matriz de tempos de processamento (m x n) para um dataframe
        df_T = convert(DataFrame, T_t)
        names!(df_T, [Symbol("m$i") for i in 1:m])

        # Objetivo: calcular as medidas estatisticas dos tempos de processamento de cada job j,
        #           considerando todas as m maquinas
        #df_T
        #describe(df_T)
        df_stats = DataFrame(j = Int64[], Mean = Float64[], Median = Float64[], Sum = Float64[], Std = Float64[])
        map((row, j) -> push!(df_stats, [
            j,
            mean(convert(Array, row)),
            median(convert(Array, row)),
            sum(convert(Array, row)),
            std(convert(Array, row))]),

            eachrow(df_T), [j for j in 1:n]);


        # Selecionando para o subconjunto S os k jobs com maior tempo de processamento, de acordo com os criterios:
        # media, mediana, soma e desvio padrao
        #k = floor(Int, (2 * n) / m)
        for measure in measure_list  # [:Median, :Std] :Mean
            for crit in [true, false]
                # Media - obtem os k jobs com maior media de tempo de processamento
                # Ordenacao decendente -> rev = true
                df = sort(df_stats, [order(measure, rev = crit)])
                println("*** $(k) jobs com $(crit ? "maior" : "menor") tempo de processamento")
                println("Criterio de medida eh $(measure)")
                showall(df[1:k,[:j, :Mean, :Median, :Sum, :Std]])
                S = convert(Array, df[1:k,[:j]])
                println("S = $(S)")

                # Invocando os modelos Wagner e Wilson e armazenando os valores obtidos de Lower Bound
                LB_wilson, time_wilson, opt_wilson = solve_wilson_model(filepath, time_limit, true, S)
                LB_wagner, time_wagner, opt_wagner = solve_wagner_model(filepath, time_limit, true, S)
                push!(df_LB, [test_name, measure, k, crit ? ">" : "<", "Wagner", LB_wagner, time_wagner, opt_wagner, string(S)])
                push!(df_LB, [test_name, measure, k, crit ? ">" : "<", "Wilson", LB_wilson, time_wilson, opt_wilson, string(S)])
                println("Saving to results file...")
                CSV.write(output_file, df_LB)
            end
        end
    end
    println("\n===================== RESULTADOS =====================\n")
    showall(df_LB)
    return df_LB
end

#file = "../../instances/taillard/output/50x20/tail051.in"
#df_LB_50_20 = process_instance_file(file, 1800, "Tail051_", 5, 40)

function test_20x5()
    folder = "../../instances/taillard/output/20x5"
    df_LB_20_5 = process_instance_folder_LB(folder, 1800, "20x5")
    df_20_5_opt = process_instance_folder_opt(folder, 1800, "20x5")
    println("\n===================== RESULTADOS 20x5 =====================\n")
    showall(df_20_5_opt)
end

function test_david()
    filepath = "../../instances/taillard/output/50x20/tail051.in"
    process_instance_file(filepath, 1800, "ta51_dec_sum", 25, 25, [:Sum, :Median, :Std])
    #filepath = "../../instances/taillard/output/20x5/tail01.in"
    #process_instance_file(filepath, 1800, "ta01_dec_sum", 10, 10, [:Sum])
end

# test_20x5()
test_david()
