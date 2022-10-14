    # robust_pfsp_budget_worstcase_wct.jl
# Worst-case calculation for the robust budget problem
# Weighed Completion Time Objective Function

# Include file reader and util functions
include("../pfsp_file_reader.jl")
# Deterministic PFSP functions
#include("deterministic_pfsp.jl")
include("../robust_pfsp_budget_worstcase_m_mach.jl")

include("robust_pfsp_budget_worstcase_wct_brute.jl")

using Combinatorics
using OffsetArrays  # Fortran-like arrays with arbitrary, zero or negative starting indices

EPS = 1e-5

# Calculates the worst-case weighted completion time given a sequence / permutation pi,
# budget parameters Γ[q] (for q in 1:m) and a weight vector w.
# Uses a dynamic programming method.
# Note: The calculation of the worst completion time of each machine M[q]
#     depends only on the previous machine M[q-1] and the machine M[q] itself.
#     So the complexity will be O( m * n^3 ).
function worst_case_wct_dp_m_machines_2(pi, m, n, Γ, P_bar, P_hat, w, verbose = false)
    #Solve Separation problem using dynamic programming
    #println("jobpos= ",jobpos)
    #calculate max_time_deviation for this solution
    max_time_deviation = Int64(0)
    max_gamma = Int64(0)
    #sort process_time_deviation and pick the gamma bigger ones
    for r in 1:m
        gamma = Int64(ceil(Γ[r]))
        max_gamma = maximum([max_gamma, gamma])
    end
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
    #max_time_deviation *= 4
    println("[DP] Max time deviation = $(max_time_deviation)")

    # For each machine r, sort the jobs by non-increasing order of processing times
    proctimes = Array{Array{Tuple{Int64,Float64},1}}(undef, m)
    for r in 1:m
        proctimes_r = []
        for i in 1:n
            push!(proctimes_r, (i, P_bar[r, i] + P_hat[r, i]))
        end
        sort!(proctimes_r, lt = (x, y) -> (x[2] > y[2] || (x[2] == y[2] && x[1] < y[1])) )
        proctimes[r] = proctimes_r
        println("proctimes[$(r)] = $(proctimes_r)")
    end

    # Create dynamic programming matrix to store results
    # First index : which position k in 1:n
    # Second index : which machine r in 1:m
    # Third index : number of deviated processing times g on Machine r in 0:max_gamma
    # Forth index : accumulated processing time deviation (discretized) d in 0:max_time_deviation
    dynmatrix = OffsetArray{Float64}(undef, 1:n, 1:m, 0:max_gamma, 0:max_time_deviation)
    deviated_jobs = OffsetArray{Array{Any, 1}}(undef, 1:n, 1:m, 0:max_gamma, 0:max_time_deviation)
    Ca = OffsetArray{Float64}(undef, 1:n, 1:m, 0:max_gamma, 0:max_time_deviation)
    fill!(dynmatrix, 0)
    fill!(deviated_jobs, [ [] for i in 1:m ])
    fill!(Ca, 0)

    gamma = max_gamma
    max_wct = zeros(m)
    max_delayed_jobs = Array{Array}(undef, n, m)
    max_delayed_jobs_all = [ [] for r in 1:m]
    fill!(max_delayed_jobs, [])
    for k in 1:n
        for r in 1:m
            # calculate processing time of job in position k
            jobk = pi[k]
            #jobk = proctimes[r][k][1]
            ptk = Int64(ceil(P_hat[r, jobk]))
            for g in 1:gamma
                for d in 0:max_time_deviation
                    if d - ptk < 0
                        if k > 1 && r > 1
                            Ca[k, r, g, d] = max(Ca[k, r - 1, g, d], Ca[k - 1, r, g, d]) + P_bar[r, pi[k]]
                            if Ca[k, r - 1, g, d] > Ca[k - 1, r, g, d]
                                dynmatrix[k, r, g, d] = dynmatrix[k, r - 1, g, d]
                            else
                                dynmatrix[k, r, g, d] = dynmatrix[k - 1, r, g, d]
                            end
                        elseif k > 1
                            Ca[k, r, g, d] = Ca[k - 1, r, g, d] + P_bar[r, pi[k]]
                        elseif r > 1
                            Ca[k, r, g, d] = Ca[k, r - 1, g, d] + P_bar[r, pi[k]]
                        else
                            Ca[k, r, g, d] = P_bar[r, pi[k]]
                        end
                        #dynmatrix[k, r, g, d] =
                    else
                        if k > 1 && r > 1
                            Ca[k, r, g, d] = max(Ca[k, r - 1, g - 1, d - ptk], Ca[k - 1, r, g - 1, d - ptk]) + P_bar[r, pi[k]] + P_hat[r, pi[k]]
                        elseif k > 1
                            Ca[k, r, g, d] = Ca[k - 1, r, g - 1, d - ptk] + P_bar[r, pi[k]] + P_hat[r, pi[k]]
                        elseif r > 1
                            Ca[k, r, g, d] = Ca[k, r - 1, g - 1, d - ptk] + P_bar[r, pi[k]] + P_hat[r, pi[k]]
                        else
                            Ca[k, r, g, d] = P_bar[r, pi[k]] + P_hat[r, pi[k]]
                        end
                    end
                end
            end
        end
    end
    r = m
    cmax = 0
    for k in 1:n
        for d in reverse(0:max_time_deviation)
            #if verbose  print("dynmatrix[$(r), $(k), $(Int64(ceil(Γ[r]))), $(d)] = $(dynmatrix[k, r, Int64(ceil(Γ[r])), d]) ; deviated_jobs = $(deviated_jobs[k, r, Int64(ceil(Γ[r])), d])")  end
            if cmax < Ca[k, r, Int64(ceil(Γ[r])), d]
                cmax =  Ca[k, r, Int64(ceil(Γ[r])), d]
            end
            if max_wct[r] < dynmatrix[k, r, Int64(ceil(Γ[r])), d]
                max_wct[r] = dynmatrix[k, r, Int64(ceil(Γ[r])), d]
                max_delayed_jobs[r] = deepcopy(deviated_jobs[k, r, Int64(ceil(Γ[r])), d])
                max_delayed_jobs_all[r] = [ deepcopy(deviated_jobs[k, r, Int64(ceil(Γ[r])), d]) ]
            elseif max_wct[r] <= dynmatrix[k, r, Int64(ceil(Γ[r])), d]  # tie
                if !(deviated_jobs[k, r, Int64(ceil(Γ[r])), d] in max_delayed_jobs_all[r])
                    push!(max_delayed_jobs_all[r], deepcopy(deviated_jobs[k, r, Int64(ceil(Γ[r])), d]))
                end
            end
        end
    end
    println("Cmax = $(cmax)")
    if verbose
        for r in 1:m
            println("max_wct[$(r)] = $(max_wct[r]) ; max_delayed_jobs_all[$(r)] = $(max_delayed_jobs_all[r])")
            println("max_delayed_jobs[$(r)] = $(max_delayed_jobs[r])")
            println("max_delayed_jobs_all[$(r)] = $(max_delayed_jobs_all[r])")
        end
    end
    Γ_scenario_list = []
    max_wct = 0
    which_scenario = []
    println("max_delayed_jobs_all = $(max_delayed_jobs_all)")
    for worst_Γ_scenario in max_delayed_jobs_all[m]
        wct = calculate_wct_given_scenario(pi, m, n, P_bar, P_hat, w, worst_Γ_scenario, false )     # true)
        println("Testing scenario $(worst_Γ_scenario), wct = $(wct).")
        push!(Γ_scenario_list, worst_Γ_scenario)
        max_wct = max(max_wct, wct)
        which_scenario = worst_Γ_scenario
    end
    #for r in 1:m
    #    println("max_wct[$(r)] = $(max_wct[r])")
    #end
    println("[DP]  The m-machine worst-case WCT is $(max_wct) ; which_scenario = $(which_scenario)\n")
    # validating objective function value
    wct_validation = calculate_wct_given_scenario(pi, m, n, P_bar, P_hat, w, which_scenario)
    #@assert abs(wct_validation - max_wct[m]) < 1e-5
    return max_wct, [which_scenario]  ## max_wct[m]
end

function scenario_fits(scenario, Γ, m)
    fits = true
    for i in 1:m
        if size(scenario[i], 1) != Γ[i]
            fits = false
            break
        end
    end
    return fits
end

# f(1, ptk, pi, 1, n, P_bar, P_hat, w)
# pi                : the sequence / permutation pi
# r                 : the machine index r in 1:m
# sequence_pos      : until which sequence position should the objective function be calculated ?
# current_deviation : current processing time deviation for machine r
# deviated_jobs     : the list of deviated jobs in each machine i in 1:m
# n                 : the number of jobs n
# P_har             : the nominal processing times matrix
# P_hat             : the processing time deviations matrix
# w                 : the weight array for the weighted completion time calculation
function f(pi, r, sequence_pos, current_deviation, deviated_jobs, m, n, P_bar, P_hat, w, Γ, proctimes, use_r0 = false)
    r0 = 1
    if use_r0
        r0 = r
    end
    #r = m  DO NOT UNCOMMENT THIS LINE !
    C = zeros(Float64, (m, n))
    k = pi[1]
    C[r0, 1] = P_bar[r0, k]
    if k in deviated_jobs[r0]
        C[r0, 1] += P_hat[r0, k]
    end
    for i in (r0+1):r
        C[i, 1] = C[i - 1, 1] + P_bar[i, k]
        if k in deviated_jobs[i]
            C[i, 1] += P_hat[i, k]
        end
    end
    for j in 2:n
        k = pi[j]
        C[r0, j] = C[r0, j - 1] + P_bar[r0, k]
        if k in deviated_jobs[r0]
            C[r0, j] += P_hat[r0, k]
        end
        for i in (r0+1):r
            C[i, j] = max(C[i, j - 1], C[i - 1, j]) + P_bar[i, k]
            if k in deviated_jobs[i]
                C[i, j] += P_hat[i, k]
            end
        end
    end
    sum_wct = 0
    for j in 1:n  #sequence_pos  # C[r, j] is the final completion time of job on position j on machine r
        k = pi[j]   # Job k is in position j of the permutation
        sum_wct += C[r, j] * w[k]
        #print("C[$(r), $(j)] = $(C[r, j]) ")
    end
    #println()
    return sum_wct
end     #end function f

# f(1, ptk, pi, 1, n, P_bar, P_hat, w)
# pi                : the sequence / permutation pi
# r                 : the machine index r in 1:m
# sequence_pos      : until which sequence position should the objective function be calculated ?
# current_deviation : current processing time deviation for machine r
# deviated_jobs     : the list of deviated jobs in each machine i in 1:m
# n                 : the number of jobs n
# P_har             : the nominal processing times matrix
# P_hat             : the processing time deviations matrix
# w                 : the weight array for the weighted completion time calculation
function f_reverse(pi, r, sequence_pos, current_deviation, deviated_jobs, m, n, P_bar, P_hat, w)
    #sum = current_deviation
    # fixme
    #println("Deviated jobs : $(deviated_jobs)")
    #r = 1
    Cr = zeros(Float64, (m + 1, n + 1))
    for ji in reverse(1:n)
        j = pi[ji]
        for i in reverse(r:m)
            Cr[i, ji] = max(Cr[i+1, ji], Cr[i, ji+1]) + P_bar[i, j]
            if j in deviated_jobs[i]
                Cr[i, ji] += P_hat[i, j]
            end
        end
    end
    sum_wct = 0
    for ji in 1:n  #sequence_pos  # C[r, j] is the final completion time of job on position j on machine r
        j = pi[ji]
        sum_wct += Cr[r, ji] * w[j]
        #print("Cr[$(r), $(ji)] = $(Cr[r, ji]) ")
    end
    #println()
    return sum_wct
end     #end function f

function test_reverse_makespan()
    searchdir(path,key) = filter(x->occursin(key,x), readdir(path))
    basedir = "/home/ubuntu/RPFS_Budget_TWCT/instances/robust/ying/data"
    for filename in searchdir(basedir, ".txt")
        if filename == "RB0104010_20perc_3m.txt"
            inputfile = "$(basedir)/$(filename)"
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

#test_reverse_makespan()
