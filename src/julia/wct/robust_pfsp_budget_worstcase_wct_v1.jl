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
using BenchmarkTools

EPS = 1e-5

function calculate_max_time_deviation(m, n, Γ, P_bar, P_hat)
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
function f(pi, r, k, C, deviated_jobs, m, n, P_bar, P_hat, w)
    #println("Deviated jobs : $(deviated_jobs)")
    # indices r and k (machine and sequence position) are already iterated inside the DP procedure
    x = pi[k]
    if r == 1
        C[1, k] = (k >= 2 ? C[1, k - 1] : 0) + P_bar[1, x]
        if x in deviated_jobs[1]
            C[1, k] += P_hat[1, x]
        end
    else
        C[r, k] = max((k >= 2 ? C[r, k - 1] : 0), C[r - 1, k]) + P_bar[r, x]
        if x in deviated_jobs[r]
            C[r, k] += P_hat[r, x]
        end
    end
    sum_wct = 0  # Calculate the partial completion time for the current machine, until sequence index k
    for j in 1:k  # C[r, j] is the final completion time of job on position j on machine r
        x = pi[j]   # Job x is in position j of the permutation
        sum_wct += C[r, j] * w[x]
    end
    return sum_wct
end     #end function f

# Calculates the worst-case weighted completion time given a sequence / permutation pi,
# budget parameters Γ[r] (for r in 1:m) and a weight vector w.
# Uses a dynamic programming method.
# Note: The calculation of the worst completion time of each machine M[r]
#     depends only on the previous machine M[r-1] and the machine M[r] itself.
#     So the complexity will be O( m * n^3 ).
function worst_case_wct_dp_m_machines(pi, m, n, Γ, P_bar, P_hat, w)
    # First index: r \in {1..m}  ( for each machine r in 1..m )
    # Second index: k \in {0..n} ( for each position k in the permutation pi )
    # Third index: g1 \in {0..n} ( for each number of processing time deviations in machine [r - 1] ) : g1 \in 0:Γ[r-1]
    # Forth index: g2 \in {0..n} ( for each number of processing time deviations in machine [r] )     : g2 \in 0:Γ[r]
    println("[DP] Calculating worst-case WCT given a permutation pi = $(pi), m = $(m), n = $(n), w = $(w) and budget parameters Γ = $(Γ)...")
    max_time_deviation = calculate_max_time_deviation(m, n, Γ, P_bar, P_hat)
    println("[DP] Max time deviation = $(max_time_deviation)")

    # Create dynamic programming matrix to store results
    # First index : which machine r in 1:m
    # Second index : which position k in 1:n
    # Third index : number of deviated processing times g on Machine r in 0:max_gamma
    # Forth index : accumulated processing time deviation (discretized) d in 0:max_time_deviation
    dynmatrix = OffsetArray{Float64}(undef, 1:m, 1:n, 0:max_gamma, 0:max_time_deviation)
    deviated_jobs = OffsetArray{Array{Any, 1}}(undef, 1:m, 1:n, 0:max_gamma, 0:max_time_deviation)
    C = OffsetArray{Array{Float64, 2}}(undef, 1:m, 1:n, 0:max_gamma, 0:max_time_deviation)
    fill!(dynmatrix, -Inf)
    fill!(deviated_jobs, [ [] for i in 1:m ])
    fill!(C, zeroes(m, n))

    ceil_T1 = Int64(ceil(Γ[1]))
    # NO position will vary only a % for the worst-case processing time in Machine M_1
    for k in 1:n
        for g1 in 0:(ceil_T1+1)  # was 0:n
            for d in 0:max_time_deviation+1
                # calculate processing time of job in position k
                jobk = pi[k]
                ptk = Int64(ceil(P_hat[r, jobk]))
                # TODO alterar a logica abaixo para a maquina r == 1
                if d - ptk < 0
                    deviated_jobs[r, k, g, d] = deepcopy(deviated_jobs[r, k - 1, g, d])
                    C[r, k, g, d] = deepcopy(C[r, k - 1, g, d])
                    dynmatrix[r, k, g, d] = f(pi, r, k, C[r, k, g, d], deviated_jobs[r, k, g, d], m, n, P_bar, P_hat, w)
                else
                    #println("jobk = $(jobk) ; deviated_jobs[r, k - 1, g - 1, d - ptk] = $(deviated_jobs[r, k - 1, g - 1, d - ptk])")
                    if (dynmatrix[r, k - 1, g, d] < dynmatrix[r, k - 1, g - 1, d - ptk]) #&& !(jobk in deviated_jobs[r, k - 1, g - 1, d - ptk][r])
                        #println("Deviating job $(jobk)")
                        # deviated_jobs[r, k, g, d] = [ deviated_jobs[r, k - 1, g - 1, d - ptk] ; jobk ]
                        deviated_jobs[r, k, g, d] = deepcopy(deviated_jobs[r, k - 1, g - 1, d - ptk])
                        #println("Before change (k = $(k - 1); g = $(g - 1); d = $(d - ptk)): $(deviated_jobs[r, k - 1, g - 1, d - ptk])")
                        push!(deviated_jobs[r, k, g, d][r], jobk)
                        #println("After change (k = $(k); g = $(g); d = $(d)): $(deviated_jobs[r, k, g, d])")
                        C[r, k, g, d] = deepcopy(C[r, k - 1, g - 1, d - ptk])
                        dynmatrix[r, k, g, d] = f(pi, r, k, C[r, k, g, d], deviated_jobs[r, k, g, d], m, n, P_bar, P_hat, w)
                        #println("CHANGE : r = $(r); k = $(k), g = $(g), d = $(d) : $(deviated_jobs[r, k, g, d]) v = $(dynmatrix[r, k, g, d])")
                    else
                        deviated_jobs[r, k, g, d] = deepcopy(deviated_jobs[r, k - 1, g, d])
                        C[r, k, g, d] = deepcopy(C[r, k - 1, g, d])
                        dynmatrix[r, k, g, d] = f(pi, r, k, C[r, k, g, d], deviated_jobs[r, k, g, d], m, n, P_bar, P_hat, w)
                    end
                end
            end

            if dynmatrix[k-1,g,d] < dynmatrix[k-1,g-1,d-ptk]
              dynmatrix[k,g,d]= ff(k,d-1,njobs,processtime,duedate,jobpos)+ dynmatrix[k-1,g-1,d-ptk]
              delayedjobs[k,g,d]=[delayedjobs[k-1,g-1,d-ptk];jobk]
            else
              dynmatrix[k,g,d]= ff(k,d-1,njobs,processtime,duedate,jobpos)+ dynmatrix[k-1,g,d]
              delayedjobs[k,g,d]=delayedjobs[k-1,g,d]
            end
        end
    end
    # Now filling matrix Ca[r] : g1 refers to the budget of Machine [r - 1]
    #   and g2 refers to the budget of uncertainty of Machine [r]
    for r in 2:m
        for g1 in 0:Int64(ceil(Γ[r - 1]))
            #g1 = Int64(ceil(Γ[r - 1]))
            ceil_Tr = Int64(ceil(Γ[r]))
            # NO position will vary only a % for the worst-case processing time in Machine M_q
            for k in 1:n
                for g2 in 0:(ceil_Tr+1)  # was 0:n, should be 0:Γ[r]
                    for d in 0:max_time_deviation+1
                        # calculate processing time of job in position k
                        jobk = pi[k]
                        ptk = Int64(ceil(P_hat[r, jobk]))
                        if d - ptk < 0
                            deviated_jobs[r, k, g, d] = deepcopy(deviated_jobs[r, k - 1, g, d])
                            C[r, k, g, d] = deepcopy(C[r, k - 1, g, d])
                            dynmatrix[r, k, g, d] = f(pi, r, k, C[r, k, g, d], deviated_jobs[r, k, g, d], m, n, P_bar, P_hat, w)
                        else
                            #println("jobk = $(jobk) ; deviated_jobs[r, k - 1, g - 1, d - ptk] = $(deviated_jobs[r, k - 1, g - 1, d - ptk])")
                            if (dynmatrix[r, k - 1, g, d] < dynmatrix[r, k - 1, g - 1, d - ptk]) #&& !(jobk in deviated_jobs[r, k - 1, g - 1, d - ptk][r])
                                #println("Deviating job $(jobk)")
                                # deviated_jobs[r, k, g, d] = [ deviated_jobs[r, k - 1, g - 1, d - ptk] ; jobk ]
                                deviated_jobs[r, k, g, d] = deepcopy(deviated_jobs[r, k - 1, g - 1, d - ptk])
                                #println("Before change (k = $(k - 1); g = $(g - 1); d = $(d - ptk)): $(deviated_jobs[r, k - 1, g - 1, d - ptk])")
                                push!(deviated_jobs[r, k, g, d][r], jobk)
                                #println("After change (k = $(k); g = $(g); d = $(d)): $(deviated_jobs[r, k, g, d])")
                                C[r, k, g, d] = deepcopy(C[r, k - 1, g - 1, d - ptk])
                                dynmatrix[r, k, g, d] = f(pi, r, k, C[r, k, g, d], deviated_jobs[r, k, g, d], m, n, P_bar, P_hat, w)
                                #println("CHANGE : r = $(r); k = $(k), g = $(g), d = $(d) : $(deviated_jobs[r, k, g, d]) v = $(dynmatrix[r, k, g, d])")
                            else
                                deviated_jobs[r, k, g, d] = deepcopy(deviated_jobs[r, k - 1, g, d])
                                C[r, k, g, d] = deepcopy(C[r, k - 1, g, d])
                                dynmatrix[r, k, g, d] = f(pi, r, k, C[r, k, g, d], deviated_jobs[r, k, g, d], m, n, P_bar, P_hat, w)
                            end
                        end
                    end
                end
            end
        end
    end  # for each machine r

    max_wct = dynmatrix[m, n, Int64(ceil(Γ[r])), max_time_deviation+1]
    max_deviated_jobs = deviated_jobs[m, n, Int64(ceil(Γ[r])), max_time_deviation+1]
    for d in 0:max_time_deviation
        if max_wct < dynmatrix[m, n, Int64(ceil(Γ[r])), d]
            max_wct = dynmatrix[m, n, Int64(ceil(Γ[r])), d]
            max_deviated_jobs = deviated_jobs[m, n, Int64(ceil(Γ[r])), d]
        end
    end
    worst_Γ_scenario = max_deviated_jobs
    wct = calculate_wct_given_scenario(pi, m, n, P_bar, P_hat, w, worst_Γ_scenario, true)
    println("Given scenario = $(worst_Γ_scenario), the WCT is $(wct).")
    println("[DP] The m-machine worst-case WCT is $(max_wct); validation = $(wct).\n")
    return max_wct, max_deviated_jobs  # wct[m, n, Int64(ceil(Γ[m])), 0]
end

function worst_case_wct_dp_m_machines_2(pi, m, n, Γ, P_bar, P_hat, w, expand_search_set = true)
    #Solve Separation problem using dynamic programming
    #println("jobpos= ",jobpos)
    #calculate max_time_deviation for this solution
    max_time_deviation = calculate_max_time_deviation(m, n, Γ, P_bar, P_hat)
    #max_time_deviation *= 2
    println("[DP] Max time deviation = $(max_time_deviation)")
    max_gamma = Int64(0)
    #sort process_time_deviation and pick the gamma bigger ones
    for r in 1:m
        gamma = Int64(ceil(Γ[r]))
        max_gamma = maximum([max_gamma, gamma])
    end

    # Create dynamic programming matrix to store results
    # First index : which machine r in 1:m
    # Second index : which position k in 1:n
    # Third index : number of deviated processing times g in 1:max_gamma+1
    # Forth index : accumulated processing time deviation (discretized) d in 1:max_time_deviation+1
    dynmatrix = zeros(m, n, max_gamma + 1, max_time_deviation + 1)
    deviated_jobs = Array{Array{Any, 1}}(undef, m, n, max_gamma + 1, max_time_deviation + 1)
    fill!(dynmatrix, 0)
    fill!(deviated_jobs, [ [] for i in 1:m ])
    gamma = max_gamma
    max_wct = zeros(m)
    max_delayed_jobs = Array{Array}(undef, m)
    max_delayed_jobs_all = Array{Array}(undef, m)
    delayed_jobs_all = Array{Array}(undef, m)
    fill!(max_delayed_jobs, [])
    fill!(max_delayed_jobs_all, [])
    fill!(delayed_jobs_all, [])
    # Initialize the matrices for the first position (k = 1)
    for r in 1:m    # For each machine r
        for g in 1:(gamma + 1)  # For each gamma value g
            for d in 1:(max_time_deviation + 1)     # For each accumulated proc time deviation
                dynmatrix[r, 1, g, d] = 0
                deviated_jobs[r, 1, g, d] = deepcopy([ [] for i in 1:m ])
            end
        end
    end
    # **************  RULES FOR MACHINE r = 1 *******************************
    # calculate processing time of job in position 1
    jobk = pi[1]
    ptk = P_hat[1, jobk]
    deviated_jobs[1, 1, 1, Int64(ceil(ptk)) + 1] = deepcopy([ [] for i in 1:m ])
    for g in 2:(gamma + 1)
        deviated_jobs[1, 1, g, Int64(ceil(ptk)) + 1] = deepcopy([ [] for i in 1:m ])
        deviated_jobs[1, 1, g, Int64(ceil(ptk)) + 1][1] = [ jobk ]
        dynmatrix[1, 1, g, Int64(ceil(ptk)) + 1] = ff(pi, 1, 1, ptk, deviated_jobs[1, 1, g, Int64(ceil(ptk)) + 1], m, n, P_bar, P_hat, w)
        #println("dynmatrix[1, 1, $(g), $(Int64(ceil(ptk)) + 1)] =  $(dynmatrix[1, 1, g, Int64(ceil(ptk)) + 1]) : deviated = $(deviated_jobs[1, 1, g, Int64(ceil(ptk)) + 1])")
    end

    deviated_jobs[1, 1, 1, 1] = deepcopy([ [] for i in 1:m ])
    dynmatrix[1, 1, 1, 1] = ff(pi, 1, 1, 0, deviated_jobs, m, n, P_bar, P_hat, w)
    for k in 2:n
        deviated_jobs[1, k, 1, 1] = deepcopy(deviated_jobs[1, k - 1, 1, 1])
        dynmatrix[1, k, 1, 1] = ff(pi, 1, k, 0, deviated_jobs[1, k, 1, 1], m, n, P_bar, P_hat, w) # + dynmatrix[1, k - 1, 1, 1]
        #println("dynmatrix[1, $(k), 1, 1] =  $(dynmatrix[1, k, 1, 1]) : deviated = $(deviated_jobs[1, k, 1, 1])")
        for d in 2:(max_time_deviation + 1)
            deviated_jobs[1, k, 1, d] = deepcopy([ [] for i in 1:m ])
            dynmatrix[1, k, 1, d] = ff(pi, 1, k, 0, deviated_jobs[1, k, 1, d], m, n, P_bar, P_hat, w)   #  -Inf
            #println("dynmatrix[1, $(k), 1, $(d)] =  $(dynmatrix[1, k, 1, d]) : deviated = $(deviated_jobs[1, k, 1, d])")
        end
    end
    for r in 1:m
        if r == 1
            max_delayed_jobs_all_combinations = [ [] ]
        #elseif r == 2
        #    max_delayed_jobs_all_combinations = [ [ [1, 5, 2], [], [] ] ]  # [ [ [10, 1, 2], [], [] ] ]
        #elseif r == 3
        #    max_delayed_jobs_all_combinations = [ [ [1, 5, 2], [5, 10, 3], [] ] ]  # [ [ [10, 1, 2], [10, 7, 5], [] ] ]
        else
            max_delayed_jobs_all_combinations = max_delayed_jobs_all[r - 1]
            if expand_search_set
                max_delayed_jobs_all_combinations = [max_delayed_jobs_all_combinations ; delayed_jobs_all[r - 1]]
            end
        end

        max_wct[r] = 0
        max_delayed_jobs[r] = deepcopy( [ [] ] )
        max_delayed_jobs_all[r] = [ ]

        for max_delayed_jobs_prev_machine in max_delayed_jobs_all_combinations
            if r >= 2
                # calculate processing time of job in position 1
                jobk = pi[1]
                ptk = P_hat[r, jobk]
                deviated_jobs[r, 1, 1, Int64(ceil(ptk)) + 1] = deepcopy(max_delayed_jobs_prev_machine)  # deepcopy([ [] for i in 1:m ])
                dynmatrix[r, 1, 1, Int64(ceil(ptk)) + 1] = ff(pi, r, 1, 0, deviated_jobs[r, 1, 1, Int64(ceil(ptk)) + 1], m, n, P_bar, P_hat, w)
                for g in 2:(gamma + 1)
                    deviated_jobs[r, 1, g, Int64(ceil(ptk)) + 1] = deepcopy(max_delayed_jobs_prev_machine)  #  [ [] for i in 1:m ]
                    deviated_jobs[r, 1, g, Int64(ceil(ptk)) + 1][r] = [ jobk ]
                    dynmatrix[r, 1, g, Int64(ceil(ptk)) + 1] = ff(pi, r, 1, ptk, deviated_jobs[r, 1, g, Int64(ceil(ptk)) + 1], m, n, P_bar, P_hat, w)
                    #println("dynmatrix[r, 1, $(g), $(Int64(ceil(ptk)) + 1)] =  $(dynmatrix[r, 1, g, Int64(ceil(ptk)) + 1]) : deviated = $(deviated_jobs[r, 1, g, Int64(ceil(ptk)) + 1])")
                end

                deviated_jobs[r, 1, 1, 1] = deepcopy(max_delayed_jobs_prev_machine)  # deviated_jobs[r - 1, k - 1, 1, 1])  #  [ [] for i in 1:m ]
                dynmatrix[r, 1, 1, 1] = ff(pi, r, 1, 0, deviated_jobs[r, 1, 1, 1], m, n, P_bar, P_hat, w)
                for k in 2:n
                    deviated_jobs[r, k, 1, 1] = deepcopy(deviated_jobs[r, k - 1, 1, 1])
                    dynmatrix[r, k, 1, 1] = ff(pi, r, k, 0, deviated_jobs[r, k, 1, 1], m, n, P_bar, P_hat, w) # + dynmatrix[r, k - 1, 1, 1]
                    for d in 2:(max_time_deviation + 1)
                        deviated_jobs[r, k, 1, d] = deepcopy(max_delayed_jobs_prev_machine)   #  deepcopy(deviated_jobs[r - 1, k, 1, d])
                        dynmatrix[r, k, 1, d] = ff(pi, r, k, 0, deviated_jobs[r, k, 1, d], m, n, P_bar, P_hat, w) # + dynmatrix[r, k - 1, 1, 1]  # dynmatrix[r - 1, k, 1, d]
                    end
                end
            end
            for k in 2:n
                #println("*** k = $(k)")
                for g in 2:(gamma + 1)
                    #println("g = $(g)")
                    for d in 1:(max_time_deviation+1)
                        # calculate processing time of job in position k
                        jobk = pi[k]
                        ptk = Int64(ceil(P_hat[r, jobk]))
                        if d - ptk - 1 < 0
                            deviated_jobs[r, k, g, d] = deepcopy(deviated_jobs[r, k - 1, g, d])
                            dynmatrix[r, k, g, d] = ff(pi, r, k, d - 1, deviated_jobs[r, k, g, d], m, n, P_bar, P_hat, w) #+ dynmatrix[r, k - 1, g, d]
                        else
                            #println("jobk = $(jobk) ; deviated_jobs[r, k - 1, g - 1, d - ptk] = $(deviated_jobs[r, k - 1, g - 1, d - ptk])")
                            #if (dynmatrix[r, k - 1, g, d] < dynmatrix[r, k - 1, g - 1, d - ptk]) #&& !(jobk in deviated_jobs[r, k - 1, g - 1, d - ptk][r])
                                #println("Deviating job $(jobk)")
                                # deviated_jobs[r, k, g, d] = [ deviated_jobs[r, k - 1, g - 1, d - ptk] ; jobk ]
                                deviated_jobs[r, k, g, d] = deepcopy(deviated_jobs[r, k - 1, g - 1, d - ptk])
                                #println("Before change (k = $(k - 1); g = $(g - 1); d = $(d - ptk)): $(deviated_jobs[r, k - 1, g - 1, d - ptk])")
                                push!(deviated_jobs[r, k, g, d][r], jobk)
                                #println("After change (k = $(k); g = $(g); d = $(d)): $(deviated_jobs[r, k, g, d])")
                                dynmatrix[r, k, g, d] = ff(pi, r, k, d, deviated_jobs[r, k, g, d], m, n, P_bar, P_hat, w) #+ dynmatrix[r, k - 1, g - 1, d - ptk]
                                #println("CHANGE : r = $(r); k = $(k), g = $(g), d = $(d) : $(deviated_jobs[r, k, g, d]) v = $(dynmatrix[r, k, g, d])")
                            if (dynmatrix[r, k - 1, g, d] > dynmatrix[r, k, g, d]) #&& !(jobk in deviated_jobs[r, k - 1, g - 1, d - ptk][r])
                            #else
                                deviated_jobs[r, k, g, d] = deepcopy(deviated_jobs[r, k - 1, g, d])
                                dynmatrix[r, k, g, d] = ff(pi, r, k, d - 1, deviated_jobs[r, k, g, d], m, n, P_bar, P_hat, w) #+ dynmatrix[r, k - 1, g, d]
                            end
                        end
                    end
                end
            end
            for d in reverse(1:max_time_deviation+1)
                #println("dynmatrix[$(r), $(n), $(Int64(ceil(Γ[r])) + 1), $(d)] = $(dynmatrix[r, n, Int64(ceil(Γ[r])) + 1, d]) ; deviated_jobs = $(deviated_jobs[r, n, Int64(ceil(Γ[r])) + 1, d])")
                if max_wct[r] < dynmatrix[r, n, Int64(ceil(Γ[r])) + 1, d]
                    max_wct[r] = dynmatrix[r, n, Int64(ceil(Γ[r])) + 1, d]
                    max_delayed_jobs[r] = deepcopy(deviated_jobs[r, n, Int64(ceil(Γ[r])) + 1, d])
                    max_delayed_jobs_all[r] = [ deepcopy(deviated_jobs[r, n, Int64(ceil(Γ[r])) + 1, d]) ]
                elseif max_wct[r] <= dynmatrix[r, n, Int64(ceil(Γ[r])) + 1, d]  # tie
                    if !(deviated_jobs[r, n, Int64(ceil(Γ[r])) + 1, d] in max_delayed_jobs_all[r])
                        push!(max_delayed_jobs_all[r], deepcopy(deviated_jobs[r, n, Int64(ceil(Γ[r])) + 1, d]))
                    end
                elseif dynmatrix[r, n, Int64(ceil(Γ[r])) + 1, d] > 0 && expand_search_set
                    push!(delayed_jobs_all[r], deepcopy(deviated_jobs[r, n, Int64(ceil(Γ[r])) + 1, d]))
                end
            end
            #println("max_wct[$(r)] = $(max_wct[r]) ; max_delayed_jobs_all[$(r)] = $(max_delayed_jobs_all[r])")
            #println("delayed_jobs_all[$(r)] = $(delayed_jobs_all[r])")
        end
    end
    Γ_scenario_list = []
    for worst_Γ_scenario in max_delayed_jobs_all[m]
        #wct = calculate_wct_given_scenario(pi, m, n, P_bar, P_hat, w, worst_Γ_scenario, false )     # true)
        fits = true
        for i in 1:m
            if size(worst_Γ_scenario[i], 1) != Γ[i]
                fits = false
            end
        end
        #println("Given scenario = $(worst_Γ_scenario), the WCT is $(wct) ; fits = $(fits).")
        if fits
            push!(Γ_scenario_list, worst_Γ_scenario)
        end
    end
    #for r in 1:m
    #    println("max_wct[$(r)] = $(max_wct[r])")
    #end
    println("[DP] The m-machine worst-case WCT is $(max_wct[m]) ; Γ_scenario_list = $(Γ_scenario_list)\n")
    return max_wct[m], Γ_scenario_list
end

# ff(1, ptk, pi, 1, n, P_bar, P_hat, w)
# pi                : the sequence / permutation pi
# r                 : the machine index r in 1:m
# sequence_pos      : until which sequence position should the objective function be calculated ?
# current_deviation : current processing time deviation for machine r
# deviated_jobs     : the list of deviated jobs in each machine i in 1:m
# n                 : the number of jobs n
# P_har             : the nominal processing times matrix
# P_hat             : the processing time deviations matrix
# w                 : the weight array for the weighted completion time calculation
function ff(pi, r, sequence_pos, current_deviation, deviated_jobs, m, n, P_bar, P_hat, w)
    #sum = current_deviation
    # fixme
    #println("Deviated jobs : $(deviated_jobs)")
    r = m
    C = zeros(Float64, (m, n))
    k = pi[1]
    C[1, 1] = P_bar[1, k]
    if k in deviated_jobs[1]
        C[1, 1] += P_hat[1, k]
    end
    for i in 2:r
        C[i, 1] = C[i - 1, 1] + P_bar[i, k]
        if k in deviated_jobs[i]
            C[i, 1] += P_hat[i, k]
        end
    end
    for j in 2:n
        k = pi[j]
        C[1, j] = C[1, j - 1] + P_bar[1, k]
        if k in deviated_jobs[1]
            C[1, j] += P_hat[1, k]
        end
        for i in 2:r
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
    end
    if sum_wct < 0
        sum_wct = 0
    end
    return sum_wct
end     #end function ff


function test_two_machine_budget_worst_case_dp_v2(m, n, P_bar, P_hat, time_limit = 1800)
    Jobs = 1:n
    N = Jobs
    Machines = 1:m
    M = Machines
    # φ_star is the optimal Johnson solution obtained over the nominal processing time array (P_bar)
    #φ_star, pi, x = solve_deterministic_pfsp(m, n, P_bar)
    #println("Based on the nominal processing times, the optimal makespan φ_star is $(φ_star).")
    pi = [x for x in 1:n]
    println("pi: $(pi)")
    #   (1) Value of the Beta function
    #β = calculate_β(Γ1, Γ2, φ_star, x)
    println("Calculating worst-case WCT for sequence pi...")
    error_count = 0
    w = [1.0 for x in 1:n]  # For now, assume all job weights are equal to 1
    for Γ1 in 0:n
        for Γ2 in 0:n
            wct_brute, best_Γ = worst_case_wct_brute_force_2_machines(pi, n, Γ1, Γ2, P_bar, P_hat, w, true)
            Γ = [Γ1, Γ2]
            wct_dp2, best_Γ = worst_case_wct_dp_m_machines_2(pi, 2, n, Γ, P_bar, P_hat, w)
            if abs(wct_dp2 - wct_brute) > EPS
                println("ERROR : wct_dp incorrect : wct_dp2 = $(wct_dp2) x wct_brute = $(wct_brute)")
                error_count += 1
                break
            end
            println("****** Errors so far : $(error_count)")
        end
        if error_count > 0
            break
        end
    end
    if error_count > 0
        println("Found $(error_count) calculation errors!")
    else
        println("No errors found.")
    end
    println("Done.")
    return error_count
end


function test_m_machine_budget_worst_case_dp(m, n, P_bar, P_hat, pi, time_limit = 1800)
    Jobs = 1:n
    N = Jobs
    Machines = 1:m
    M = Machines
    #   (1) Value of the Beta function
    #β = calculate_β(Γ1, Γ2, φ_star, x)
    println("Calculating worst-case WCT for sequence pi...")
    error_count = 0
    w = [1.0 for x in 1:n]  # For now, assume all job weights are equal to 1
    for Γ1 in 1:2
        for Γ2 in 1:n
            for Γ3 in 1:n
                Γ = [Γ1, Γ2, Γ3]
                wct_brute, best_Γ, best_j = worst_case_wct_brute_force_3_machines(pi, 3, n, Γ, P_bar, P_hat, w)
                wct_dp, best_Γ_2 = worst_case_wct_dp_m_machines_2(pi, 3, n, Γ, P_bar, P_hat, w, false)
                if abs(wct_dp - wct_brute) > EPS
                    println("ERROR : wct_dp incorrect : wct_brute = $(wct_brute) x wct_dp = $(wct_dp)")
                    error_count += 1
                end
            end
        end
    end
    if error_count > 0
        println("Found $(error_count) calculation errors!")
    else
        println("No errors found.")
    end
    println("Done.")
    return error_count
end

function test_m_machine_budget_worst_case_dp_10jobs_2machines()
    searchdir(path,key) = filter(x->occursin(key,x), readdir(path))
    basedir = "../../../instances/robust/ying/data/10 jobs"
    error_count = 0
    file_count = 0
    filename_list = []
    error_count_list = []
    for perc in ["10%", "20%", "30%", "40%", "50%"]
        dir = "$(basedir)/$(perc)"
        for filename in searchdir(dir, ".txt")
            inputfile = "$(dir)/$(filename)"
            m, n, P_bar, P_hat = read_robust_ying_input_file(inputfile)
            file_error_count = test_two_machine_budget_worst_case_dp_v2(m, n, P_bar, P_hat)
            error_count += file_error_count
            file_count += 1
            println("****** Errors so far : $(error_count)")
            push!(filename_list, inputfile)
            push!(error_count_list, file_error_count)
        end
    end
    println("RESULT : After processing $(file_count) files, the total worst-wct calculation errors found is : $(error_count)")
    for i in 1:size(error_count_list, 1)
        println("$(filename_list[i]) : $(error_count_list[i])")
    end
end

function test_m_machine_budget_worst_case_dp_toy_3machines()
    searchdir(path,key) = filter(x->occursin(key,x), readdir(path))
    basedir = "/home/ubuntu/RPFS_Budget_TWCT/instances/robust/ying/data"
    for filename in searchdir(basedir, ".txt")
        if filename == "RB0104010_20perc_3m.txt"  # "test_toy_3m.txt"
            inputfile = "$(basedir)/$(filename)"
            m, n, P_bar, P_hat = read_robust_ying_input_file(inputfile)
            #count = test_two_machine_budget_worst_case_dp_v2(2, n, P_bar, P_hat)
            # φ_star is the optimal Johnson solution obtained over the nominal processing time array (P_bar)
            #φ_star, pi, x = solve_deterministic_pfsp(m, n, P_bar)
            #println("Based on the nominal processing times, the optimal makespan φ_star is $(φ_star).")
            pi = [x for x in 1:n]
            println("pi: $(pi)")
            count = test_m_machine_budget_worst_case_dp(m, n, P_bar, P_hat, pi)
            println("RESULT : Total worst-wct calculation errors found: $(count)")
        end
    end
end

function test_m_machine_budget_worst_case_dp_2machines_specific()
    searchdir(path,key) = filter(x->occursin(key,x), readdir(path))
    basedir = "/home/ubuntu/RPFS_Budget_TWCT/instances/robust/ying/data/10jobs/30%"
    for filename in searchdir(basedir, ".txt")
        if filename == "RB0103001.txt"
            inputfile = "$(basedir)/$(filename)"
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
            for Γ1 in 1:n
                for Γ2 in 1:n
                    println("*** Γ1 = $(Γ1); Γ2 = $(Γ2)")
                    wct_brute, best_Γ = worst_case_wct_brute_force_2_machines(pi, n, Γ1, Γ2, P_bar, P_hat, w)
                    wct_dp, best_Γ_2 = worst_case_wct_dp_m_machines_2(pi, m, n, [Γ1, Γ2], P_bar, P_hat, w, false)
                    if abs(wct_dp - wct_brute) > EPS
                        println("ERROR : wct_dp incorrect : wct_brute = $(wct_brute) x wct_dp = $(wct_dp)")
                        error_count += 1
                    end
                end
            end
            println("The total error count is $(error_count).")
            break
        end
    end
end

function test_m_machine_budget_worst_case_dp_3machines()
    searchdir(path,key) = filter(x->occursin(key,x), readdir(path))
    basedir = "/home/ubuntu/RPFS_Budget_TWCT/instances/robust/ying/data"
    for filename in searchdir(basedir, ".txt")
        if filename == "RB0104010_20perc_3m.txt"
            inputfile = "$(basedir)/$(filename)"
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
            for Γ1 in 3:n
                for Γ2 in 1:n
                    for Γ3 in 1:n
                        println("*** p_Γ1 = $(Γ1); p_Γ2 = $(Γ2); p_Γ3 = $(Γ3)")
                        #Γ1 = (p_Γ1 * n) / 100.0
                        #Γ2 = (p_Γ2 * n) / 100.0
                        wct_brute, worst_Γ_scenario_brute, worst_j = worst_case_wct_brute_force_3_machines(pi, m, n, [Γ1, Γ2, Γ3], P_bar, P_hat, w, false)
                        wct_dp, worst_Γ_scenario_dp = worst_case_wct_dp_m_machines_2(pi, m, n, [Γ1, Γ2, Γ3], P_bar, P_hat, w)
                        #println("worst_Γ_scenario_brute = $(worst_Γ_scenario_brute)")
                        if abs(wct_dp - wct_brute) > EPS
                            println("ERROR : wct_dp incorrect : wct_brute = $(wct_brute) x wct_dp = $(wct_dp)")
                            error_count += 1
                            break
                        end
                    end
                    if error_count > 0
                        break
                    end
                end
                if error_count > 0
                    break
                end
            end
        end
    end
end

#@btime test_m_machine_budget_worst_case_dp_10jobs_2machines()
# @btime test_m_machine_budget_worst_case_dp_toy_3machines()
@btime test_m_machine_budget_worst_case_dp_2machines_specific()  # 29 errors, 12,8 s, 2,2GB mem
