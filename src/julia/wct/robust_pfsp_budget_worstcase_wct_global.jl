# robust_pfsp_budget_worstcase_wct_global.jl
# Worst-case calculation for the robust budget problem
# Weighed Completion Time Objective Function

# Include file reader and util functions
include("../pfsp_file_reader.jl")
# Deterministic PFSP functions
#include("deterministic_pfsp.jl")
include("../robust_pfsp_budget_worstcase_m_mach.jl")

include("robust_pfsp_budget_worstcase_wct_brute.jl")
include("../config.jl")

using Combinatorics
using OffsetArrays  # Fortran-like arrays with arbitrary, zero or negative starting indices
# using BenchmarkTools

EPS = 1e-5

function IsEqual(A::Float64,B::Float64,epsilon::Float64=1.0e-5)
    return abs( A - B ) <= epsilon
end

function EPSLT(x::Float64,y::Float64,epsilon::Float64=1.0e-5)
    return ((x)-(y) < -(epsilon))
end

function EPSLE(x::Float64,y::Float64,epsilon::Float64=1.0e-5)
    return ((x)-(y) <= (epsilon))
end

function EPSGT(x::Float64,y::Float64,epsilon::Float64=1.0e-5)
    return ((x)-(y) > (epsilon))
end

function EPSGE(x::Float64,y::Float64,epsilon::Float64=1.0e-5)
    return ((x)-(y) >= -(epsilon))
end

function calculate_max_time_deviation(m, n, Γ, P_bar, P_hat)
    max_time_deviation = Int64(0)
    max_gamma = Int64(ceil(Γ))
    #sort process_time_deviation and pick the gamma bigger ones
    for r in 1:m
        process_time_deviation = P_hat[r, :]
        ptd = sort(process_time_deviation, rev=true)
        max_dev_r = 0
        for i in 1:n  # min(max_gamma + 1, m)
            max_dev_r = max_dev_r + Int64(ceil(ptd[i]))
        end
        max_dev_r = Int64(ceil(max_dev_r)) + 1
        max_time_deviation = maximum([max_time_deviation, max_dev_r]) + 1
    end
    return max_time_deviation
end

# f(pi, r, k, C, deviated_jobs, m, n, P_bar, P_hat, w)
# pi                : the sequence / permutation pi
# r                 : the current machine index r in 1:m
# k                 : until which sequence position should the objective function be calculated ?
# C                 : current completion time matrix
# deviated_jobs     : the list of deviated jobs in each machine i in 1:m
# m                 : the number of machines m
# n                 : the number of jobs n
# P_har             : the nominal processing times matrix
# P_hat             : the processing time deviations matrix
# w                 : the weight array for the weighted completion time calculation
function f!(pi, r, k, C, deviated_jobs, m, n, P_bar, P_hat, w, verbose = false)
    #println("Deviated jobs : $(deviated_jobs)")
    # indices r and k (machine and sequence position) are already iterated inside the DP procedure
    # FIXME Fazer novos testes com o recalculo completo de f()
    x = pi[k]
    if r <= 1
        C[1, k] = (k >= 2 ? C[1, k - 1] : 0) + P_bar[1, x]
        if (r, x) in deviated_jobs[1]
            if verbose
                println("[Machine $(r)] Deviating job $(x): $(C[1, k]) + $(P_hat[1, x])")
            end
            C[1, k] += P_hat[1, x]
        end
    else
        for q in (r-1):r
            C[q, k] = max((k >= 2 ? C[q, k - 1] : 0), (q >= 2 ? C[q - 1, k] : 0)) + P_bar[q, x]
            if (q, x) in deviated_jobs[q]
                if verbose
                    println("[Machine $(q)] Deviating job $(x): $(C[q, k]) + $(P_hat[q, x])")
                end
                C[q, k] += P_hat[q, x]
            end
        end
    end
    # New: propagate the makespan calculation (without deviation applied)
    for j in (k+1):n
        x = pi[j]
        C[r, j] = max(C[r, j - 1], (r >= 2 ? C[r - 1, j] : 0)) + P_bar[r, x]
    end
    for j in 1:n
        x = pi[j]
        for i in (r+1):m
            C[i, j] = max((j >= 2 ? C[i, j - 1] : 0), C[i - 1, j]) + P_bar[i, x]
        end
    end
    sum_wct = 0  # Calculate the partial completion time for the current machine, until sequence index k
    for j in 1:n  #sequence_pos  # C[r, j] is the final completion time of job on position j on machine r
        k = pi[j]   # Job k is in position j of the permutation
        sum_wct += C[m, j] * w[k]
    end
    #if verbose
    #    println("[f!] wct = $(sum_wct); C = $(C)")
    #end
    return sum_wct
end     #end function f

# Calculates the worst-case weighted completion time given a sequence / permutation pi,
# global budget parameter Γ, and a weight vector w.
# Uses a dynamic programming method.
function worst_case_wct_dp_m_machines(pi, m, n, Γ, P_bar, P_hat, w)
    # First index: r \in {1..m}  ( for each machine r in 1..m )
    # Second index: k \in {0..n} ( for each position k in the permutation pi )
    # Third index: g1 \in {0..n} ( for each number of processing time deviations in machine [r - 1] ) : g1 \in 0:Γ[r-1]
    # Forth index: g2 \in {0..n} ( for each number of processing time deviations in machine [r] )     : g2 \in 0:Γ[r]
    println("[DP] Calculating worst-case WCT given a permutation pi = $(pi), m = $(m), n = $(n), w = $(w) and budget parameter Γ = $(Γ)...")
    max_time_deviation = calculate_max_time_deviation(m, n, Γ, P_bar, P_hat)
    println("[DP] Max time deviation = $(max_time_deviation)")

    max_gamma = m * n
    #sort process_time_deviation and pick the gamma bigger ones

    # Create dynamic programming matrix to store results
    # First index : which machine r in 1:m
    # Second index : which position k in 1:n
    # Third index : number of deviated processing times g on Machine r in 0:max_gamma
    # Forth index : accumulated processing time deviation (discretized) d in 0:max_time_deviation
    dynmatrix = OffsetArray{Float64}(undef, 1:m, 0:n, 0:max_gamma+1, 0:max_time_deviation+1)
    deviated_jobs = OffsetArray{Array{Any, 1}}(undef, 1:m, 0:n, 0:max_gamma+1, 0:max_time_deviation+1)
    C = OffsetArray{Array{Float64, 2}}(undef, 1:m, 0:n, 0:max_gamma+1, 0:max_time_deviation+1)
    fill!(dynmatrix, -Inf)
    fill!(deviated_jobs, [ [] for r in 1:m ])
    fill!(C, zeros(Float64, m, n))

    ceil_T = Int64(ceil(Γ))
    # NO position will vary only a % for the worst-case processing time in Machine M_1
    r = 1
    for k in 1:n
        jobk = pi[k]
        ptk = Int64(ceil(P_hat[r, jobk]))
        #println("*** jobk = $(jobk)")
        for g in 0:(ceil_T+1)  # was 0:n
            for d in 0:max_time_deviation+1
                # calculate processing time of job in position k
                deviated_jobs[r, k, g, d] = deepcopy(deviated_jobs[r, k - 1, g, d])
                C[r, k, g, d] = deepcopy(C[r, k - 1, g, d])
                dynmatrix[r, k, g, d] = f!(pi, r, k, C[r, k, g, d], deviated_jobs[r, k, g, d], m, n, P_bar, P_hat, w)
                if !(d - ptk < 0 || g == 0)  #  WITH DEVIATION OF OPERATION (1, pi[k])
                    #println("Deviating job $(jobk)")
                    temp_deviated_jobs = deepcopy(deviated_jobs[r, k - 1, g - 1, d - ptk])
                    push!(temp_deviated_jobs[r], (r, jobk))
                    temp_C = deepcopy(C[r, k - 1, g - 1, d - ptk])
                    temp_f = f!(pi, r, k, temp_C, temp_deviated_jobs, m, n, P_bar, P_hat, w)
                    if EPSGT(temp_f, dynmatrix[r, k, g, d])  # >
                        deviated_jobs[r, k, g, d] = temp_deviated_jobs
                        C[r, k, g, d] = temp_C
                        dynmatrix[r, k, g, d] = temp_f
                    end
                end
            end
            #if dynmatrix[k-1,g,d] < dynmatrix[k-1,g-1,d-ptk]
            #  dynmatrix[k,g,d]= ff(k,d-1,njobs,processtime,duedate,jobpos)+ dynmatrix[k-1,g-1,d-ptk]
            #  delayedjobs[k,g,d]=[delayedjobs[k-1,g-1,d-ptk];jobk]
            #else
            #  dynmatrix[k,g,d]= ff(k,d-1,njobs,processtime,duedate,jobpos)+ dynmatrix[k-1,g,d]
            #  delayedjobs[k,g,d]=delayedjobs[k-1,g,d]
            #end
        end
    end
    # Now filling matrix Ca[r] : g1 refers to the budget of Machine [r - 1]
    #   and g2 refers to the budget of uncertainty of Machine [r]
    for r in 2:m
        # NO position will vary only a % for the worst-case processing time in Machine M_q
        for k in 1:n
            # calculate processing time of job in position k
            jobk = pi[k]
            ptk = Int64(ceil(P_hat[r, jobk]))
            #println("*** jobk = $(jobk)")
            for g in 0:(ceil_T+1)  # was 0:n, should be 0:Γ
                for d in 0:max_time_deviation+1
                    #x = max(Ca[q, k, y2, 0], max(Ca[q - 1, k, y2, 0], Ca[q, k - 1, y2, 0]) + P_bar[q, pi[k]])
                    # ---------------- NO DEVIATION ---------------------
                    deviated_jobs[r, k, g, d] = deepcopy(deviated_jobs[r - 1, k, g, d])
                    C[r, k, g, d] = deepcopy(C[r - 1, k, g, d])
                    #println("*** Before: $(C[r, k, g, d])")
                    dynmatrix[r, k, g, d] = f!(pi, r, k, C[r, k, g, d], deviated_jobs[r, k, g, d], m, n, P_bar, P_hat, w)
                    #println("*** After: $(C[r, k, g, d])")
                    if k > 1
                        temp_deviated_jobs = deepcopy(deviated_jobs[r, k - 1, g, d])
                        temp_C = deepcopy(C[r, k - 1, g, d])
                        temp_f = f!(pi, r, k, temp_C, temp_deviated_jobs, m, n, P_bar, P_hat, w)
                        if EPSGT(temp_f, dynmatrix[r, k, g, d])  # >
                            deviated_jobs[r, k, g, d] = temp_deviated_jobs
                            C[r, k, g, d] = temp_C
                            dynmatrix[r, k, g, d] = temp_f
                        end
                    end
                    # ---------------------------------------------------
                    if !(d - ptk < 0 || g == 0)  #  WITH DEVIATION OF OPERATION (r, pi[k])
                        #println("jobk = $(jobk) ; deviated_jobs[r, k - 1, g - 1, d - ptk] = $(deviated_jobs[r, k - 1, g - 1, d - ptk])")
                        #println("Deviating job $(jobk)")
                        temp_deviated_jobs = deepcopy(deviated_jobs[r - 1, k, g - 1, d - ptk])
                        push!(temp_deviated_jobs[r], (r, jobk))
                        temp_C = deepcopy(C[r - 1, k, g - 1, d - ptk])
                        temp_f = f!(pi, r, k, temp_C, temp_deviated_jobs, m, n, P_bar, P_hat, w)
                        if EPSGT(temp_f, dynmatrix[r, k, g, d])  # >
                            deviated_jobs[r, k, g, d] = temp_deviated_jobs
                            C[r, k, g, d] = temp_C
                            dynmatrix[r, k, g, d] = temp_f
                        end
                        #println("*** After: $(C[r, k, g, d])")
                        #println("CHANGE : r = $(r); k = $(k), g = $(g), d = $(d) : $(deviated_jobs[r, k, g, d]) v = $(dynmatrix[r, k, g, d])")
                        if k > 1
                            temp_deviated_jobs = deepcopy(deviated_jobs[r, k - 1, g - 1, d - ptk])
                            push!(temp_deviated_jobs[r], (r, jobk))
                            temp_C = deepcopy(C[r, k - 1, g - 1, d - ptk])
                            temp_f = f!(pi, r, k, temp_C, temp_deviated_jobs, m, n, P_bar, P_hat, w)
                            if EPSGT(temp_f, dynmatrix[r, k, g, d])  # >
                                deviated_jobs[r, k, g, d] = temp_deviated_jobs
                                C[r, k, g, d] = temp_C
                                dynmatrix[r, k, g, d] = temp_f
                            end
                        end
                    end
                end
            end
        end
    end  # for each machine r

    max_wct = dynmatrix[m, n, Int64(ceil(Γ)), max_time_deviation+1]
    println("worst_wct = $(max_wct)")
    max_deviated_jobs = deviated_jobs[m, n, Int64(ceil(Γ)), max_time_deviation+1]
    for d in 0:max_time_deviation
        if EPSLT(max_wct, dynmatrix[m, n, Int64(ceil(Γ)), d])  # <
            max_wct = dynmatrix[m, n, Int64(ceil(Γ)), d]
            max_deviated_jobs = deviated_jobs[m, n, Int64(ceil(Γ)), d]
            println("*** Improved! deviation d = $(d)")
        end
    end
    println("max_wct = $(max_wct); max_deviated_jobs = $(max_deviated_jobs)")
    worst_Γ_scenario = vcat(max_deviated_jobs...)
    println("worst_Γ_scenario = $(worst_Γ_scenario)")
    #worst_Γ_scenario_dev = calculate_deviation_matrix(m, n, pi, Γ, worst_Γ_scenario, [], "machine")
    wct = calculate_wct_given_global_scenario(pi, m, n, P_bar, P_hat, w, worst_Γ_scenario, true)
    println("Budget = $(Γ): Given scenario = $(worst_Γ_scenario), the WCT is $(wct).")
    println("[DP] The m-machine worst-case WCT is $(max_wct); validation = $(wct).\n")
    return max_wct, worst_Γ_scenario  # wct[m, n, Int64(ceil(Γ[m])), 0]
end


function test_reverse_makespan()
    searchdir(path,key) = filter(x->occursin(key,x), readdir(path))
    basedir = joinpath(robust_ying_instances_folder, "data")
    for filename in searchdir(basedir, ".txt")
        if filename == "RB0104010_20perc_3m.txt"
            inputfile = joinpath("$(basedir)", "$(filename)")
            m, n, P_bar, P_hat = read_robust_ying_input_file(inputfile)
            pi = reverse([x for x in 1:n])
            w = [1.0 for x in 1:n]
            deviated_jobs = [ [] for i in 1:m ]
            wct1 = f(pi, 1, 0, 0, deviated_jobs, m, n, P_bar, P_hat, w)
            wct2 = f_reverse(pi, 1, 0, 0, deviated_jobs, m, n, P_bar, P_hat, w)
            println("wct1 = $(wct1) ; wct2 = $(wct2)")
        end
    end
end

function test_m_machine_budget_worst_case_dp_specific()
    searchdir(path,key) = filter(x->occursin(key,x), readdir(path))
    basedir = joinpath(robust_ying_instances_folder, "data", "10jobs", "30%")
    total_error = 0
    error_dict = Dict()
    for filename in searchdir(basedir, ".txt")
        if filename == "RB0103003.txt"
            inputfile = joinpath("$(basedir)", "$(filename)")
            m, n, P_bar, P_hat = read_robust_ying_input_file(inputfile)
            # φ_star is the optimal Johnson solution obtained over the nominal processing time array (P_bar)
            #φ_star, pi, x = solve_deterministic_pfsp(m, n, P_bar)
            #println("Based on the nominal processing times, the optimal makespan φ_star is $(φ_star).")
            pi = [x for x in 1:n]
            println("pi: $(pi)")
            Jobs = 1:n
            N = Jobs
            Machines = 1:m
            M = Machines
            pi = [ 1, 5, 10, 2, 3, 7, 4, 6, 9, 8 ]
            # pi = [ 10, 1, 7, 5, 3, 2, 4, 6, 9, 8]
            w = [1.0 for x in 1:n]  # For now, assume all job weights are equal to 1
            println("Calculating worst-case wct for sequence pi...")
            error_count = 0
            for Γ in 1:(m*n)
                println("*** Γ = $(Γ)")
                wct_brute, best_Γ = worst_case_wct_brute_force_global_budget(pi, m, n, Γ, P_bar, P_hat, w)
                wct_dp, scenario = worst_case_wct_dp_m_machines(pi, m, n, Γ, P_bar, P_hat, w)
                if abs(wct_dp - wct_brute) > EPS
                    println("ERROR : wct_dp incorrect : wct_brute = $(wct_brute) x wct_dp = $(wct_dp)")
                    error_count += 1
                    break
                end
            end
            total_error += error_count
            error_dict[filename] = error_count
            println("The error count for $(filename) is $(error_count).")
        end
    end
    println("The total error count is $(total_error): $(error_dict).")
end

function test_m_machine_budget_worst_case_dp_3machines()
    searchdir(path,key) = filter(x->occursin(key,x), readdir(path))
    basedir = joinpath(robust_ying_instances_folder, "data")
    for filename in searchdir(basedir, ".txt")
        if filename == "RB0104010_20perc_3m.txt"
            inputfile = joinpath("$(basedir)", "$(filename)")
            m, n, P_bar, P_hat = read_robust_ying_input_file(inputfile)
            # φ_star is the optimal Johnson solution obtained over the nominal processing time array (P_bar)
            #φ_star, pi, x = solve_deterministic_pfsp(m, n, P_bar)
            #println("Based on the nominal processing times, the optimal makespan φ_star is $(φ_star).")
            #println("x: $(x)")
            Jobs = 1:n
            N = Jobs
            Machines = 1:m
            M = Machines
            # pi = [ 1, 5, 10, 2, 3, 7, 4, 6, 9, 8 ]  # Tested up to [3,0,0]
            pi = [ 10, 1, 7, 5, 3, 2, 4, 6, 9, 8]
            w = [1.0 for x in 1:n]  # For now, assume all job weights are equal to 1
            println("Calculating worst-case makespan for sequence pi = $(pi)...")
            error_count = 0
            for Γ in 1:n
                println("*** p_Γ = $(Γ)")
                wct_brute, worst_Γ_scenario_brute, worst_j = worst_case_wct_brute_force_global_budget(pi, m, n, Γ, P_bar, P_hat, w, false)
                wct_dp, worst_Γ_scenario_dp = worst_case_wct_dp_m_machines(pi, m, n, Γ, P_bar, P_hat, w)
                #println("worst_Γ_scenario_brute = $(worst_Γ_scenario_brute)")
                if abs(wct_dp - wct_brute) > EPS
                    println("ERROR : wct_dp incorrect : wct_brute = $(wct_brute) x wct_dp = $(wct_dp)")
                    error_count += 1
                    break
                end
            end
        end
    end
end

function test_m_machine_budget_worst_case_dp_toy_instance()
    searchdir(path,key) = filter(x->occursin(key,x), readdir(path))
    basedir = joinpath(robust_ying_instances_folder, "data")
    for filename in searchdir(basedir, ".txt")
        if filename == "test_toy_3m.txt"
            inputfile = joinpath("$(basedir)", "$(filename)")
            m, n, P_bar, P_hat = read_robust_ying_input_file(inputfile)
            # φ_star is the optimal Johnson solution obtained over the nominal processing time array (P_bar)
            #φ_star, pi, x = solve_deterministic_pfsp(m, n, P_bar)
            #println("Based on the nominal processing times, the optimal makespan φ_star is $(φ_star).")
            #println("x: $(x)")
            Jobs = 1:n
            N = Jobs
            Machines = 1:m
            M = Machines
            pi = [ 2, 1, 3, 4 ]
            w = [1.0 for x in 1:n]  # For now, assume all job weights are equal to 1
            println("Calculating worst-case makespan for sequence pi = $(pi)...")
            error_count = 0
            for Γ in 1:n
                println("*** p_Γ = $(Γ)")
                wct_brute, worst_Γ_scenario_brute, worst_j = worst_case_wct_brute_force_global_budget(pi, m, n, Γ, P_bar, P_hat, w, false)
                wct_dp, worst_Γ_scenario_dp = worst_case_wct_dp_m_machines(pi, m, n, Γ, P_bar, P_hat, w)
                #println("worst_Γ_scenario_brute = $(worst_Γ_scenario_brute)")
                if abs(wct_dp - wct_brute) > EPS
                    println("ERROR : wct_dp incorrect : wct_brute = $(wct_brute) x wct_dp = $(wct_dp)")
                    error_count += 1
                    break
                end
            end
        end
    end
end

#@btime test_m_machine_budget_worst_case_dp_10jobs_2machines()
#test_m_machine_budget_worst_case_dp_3machines()
#@btime
#test_m_machine_budget_worst_case_dp_specific()

# test_m_machine_budget_worst_case_dp_toy_instance()

# test_reverse_makespan()
