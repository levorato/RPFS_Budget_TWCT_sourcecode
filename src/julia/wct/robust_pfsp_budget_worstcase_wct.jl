# ===================================================================================================
# robust_pfsp_budget_worstcase_wct.jl (budgeted uncertainty)
# Worst-case determination procedure for the Robust Permutation Flow Shop Problem to Minimize the
# Weighted Sum of Job Completion Times (WCT), based on the budgeted uncertainty set (Bertsimas et al.).
# Given a permutation pi, finds the restricted worst-case scenario (i.e. for each machine r, which
# jobs will have their processing times deviated to their worst-case values), so that the
# objective function (weighted completion time - wct) reaches the worst-case (maximum) value.
# *** Dynamic Programming approach for calculating the worst-case scenario.
# ===================================================================================================

include("../config.jl")
# Include file reader and util functions
include("../pfsp_file_reader.jl")
# Deterministic PFSP functions
#include("deterministic_pfsp.jl")
include("../robust_pfsp_budget_worstcase_m_mach.jl")
include("robust_pfsp_budget_worstcase_wct_brute.jl")

using Combinatorics
using OffsetArrays  # Fortran-like arrays with arbitrary, zero or negative starting indices

EPS = 1e-5

function calculate_completion_array(m, n, pi, P_bar, P_hat, deviated_items, job_pos)
    C = zeros(Float64, (m, n))
    dev = zeros(Int64, (m, n))
    for k in 1:job_pos  # For each job position k
        for r in deviated_items[k]  # For each machine r in deviated_items[k]
            dev[r, pi[k]] = 1
        end
    end
    k = pi[1]
    C[1, 1] = P_bar[1, k] + P_hat[1, k] * dev[1, k]
    for i in 2:m
        C[i, 1] = C[i - 1, 1] + P_bar[i, k] + P_hat[i, k] * dev[i, k]
    end
    for j in 2:job_pos  # n
        k = pi[j]
        C[1, j] = C[1, j - 1] + P_bar[1, k] + P_hat[1, k] * dev[1, k]
        for i in 2:m
            C[i, j] = max(C[i, j - 1], C[i - 1, j]) + P_bar[i, k] + P_hat[i, k] * dev[i, k]
        end
    end
    println("C = $(C)")
    return C
end

# Calculates the worst-case weighted completion time given a sequence / permutation pi,
# budget parameters Γ[j] (for j in 1:m) and a weight vector w.
# Uses the dynamic programming method.
# Note: The calculation of the worst weighted completion time of each job j=pi(k)
#     depends only on the same job on the previous machine M[q-1]
#     and the previous job (pi(k-1)) on the same machine M[q].
#     So the complexity will be O( m * n^3 ).
function worst_case_wct_dp_m_machines(pi, m, n, Γ, P_bar, P_hat, w, expand_search_set = true, verbose = false)
    #Solve Separation problem using dynamic programming
    max_wct = zeros(n)
    # Create dynamic programming matrix to store results
    # First index : which machine r in 1:m
    # Second index : which position k in 1:n
    # Third index : number of deviated processing times g in 1:max_gamma+1, for a given machine, considering all jobs
    # Forth index : UNUSED (accumulated processing time deviation (discretized) d in 1:max_time_deviation+1)
    Ca = OffsetArray{Float64}(undef, 0:m, 0:(n+1), 0:(n+1), 0:n)
    Za = OffsetArray{Array{Any, 1}}(undef, 0:m, 0:(n+1), 0:(n+1))
    fill!(Ca, 0)
    fill!(Za, [])
    for y1 in 0:n  # Γ1
        for px in 1:n
            Ca[1, 0, y1, px] = 0
        end
    end
    # DP rules for machine q == 1 : we only need to look at the previous job pi(k-1) on the same machine M_1
    for y2 in 0:(ceil_T1+1)  # y2 is the budget of the machine M_1
        Ca[1, k, y2, 0] = max(Ca[1, k, y2, 0], Ca[1, k - 1, y2, 0] + P_bar[1, pi[k]])
        if Ca[1, k, y2, 0] == Ca[1, k - 1, y2, 0] + P_bar[1, pi[k]]
            Za[1, k, y2] = deepcopy(Za[1, k - 1, y2])
        end
        if y2 <= m
            Ca[1, k, y2 + 1, 0] = max(Ca[1, k, y2 + 1, 0], Ca[1, k - 1, y2, 0] + P_hat[1, pi[k]] + P_bar[1, pi[k]])
            if abs(Ca[1, k, y2 + 1, 0] - (Ca[1, k - 1, y2, 0] + P_hat[1, pi[k]] + P_bar[1, pi[k]])) < EPS
                Za[1, k, y2 + 1] = deepcopy(Za[1, k, y2])  # deepcopy(Za[k - 1, 1, y1])
                push!(Za[1, k, y2 + 1], 1)
            else
                Za[1, k, y2 + 1] = deepcopy(Za[1, k, y2])  # deepcopy(Za[k - 1, 1, y1])
            end
        end
    end
    for q in 2:m
        y1 = Int64(ceil(Γ[q - 1]))  # y1 is the budget of the previous machine M_{q-1}
        delta_Tq = ceil(Γ[q]) - floor(Γ[q])
        ceil_Tq = Int64(ceil(Γ[q]))
        # DP rules for k == 1 : we only have to look at the completion time of job pi(1) on the previous machine (q-1)
        for y2 in 0:(ceil_Tq+1)  # y2 is the budget of the current machine M_q
            Ca[q, 1, y2, 0] = max(Ca[q, 1, y2, 0], Ca[q - 1, 1, y1, 0] + P_bar[q, pi[1]])
            if Ca[q, 1, y2, 0] == Ca[q - 1, 1, y1, 0] + P_bar[q, pi[1]]
                Za[q, 1, y2] = deepcopy(Za[q - 1, 1, y1])
            end
            if y1 <= m
                Ca[q, 1, y2 + 1, 0] = max(Ca[q, 1, y2 + 1, 0], Ca[q - 1, 1, y1, 0] + P_hat[q, pi[1]] + P_bar[q, pi[1]])
                if abs(Ca[q, 1, y2 + 1, 0] - (Ca[q - 1, 1, y1, 0] + P_hat[q, pi[1]] + P_bar[q, pi[1]])) < EPS
                    Za[q, 1, y1 + 1] = deepcopy(Za[q - 1, 1, y1])
                    push!(Za[q, 1, y2 + 1], q)
                else
                    Za[q, 1, y2 + 1] = deepcopy(Za[q - 1, 1, y1])  # deepcopy(Za[k - 1, 1, y1])
                end
            end
        end
        for k in 2:n
            for y2 in 0:(ceil_Tq+1)  # y2 is the budget of the current machine M_q
                previous_value = Ca[q, k, y2, 0]
                Ca[q, k, y2, 0] = max(Ca[q, k, y2, 0], max(Ca[q, k - 1, y2, 0], Ca[q - 1, k, y1, 0]) + P_bar[q, pi[k]])
                if abs(Ca[q, k, y2, 0] - previous_value) > EPS  # Ca[q, k, y2, 0] != previous_value
                    if abs(Ca[q, k, y2, 0] - (Ca[q - 1, k, y1, 0] + P_bar[q, pi[k]])) < EPS
                        Za[q, k, y2] = deepcopy(Za[q - 1, k, y1])
                    elseif (Ca[q, k, y2, 0] - (Ca[q, k - 1, y2, 0] + P_bar[q, pi[k]])) < EPS
                        Za[q, k, y2] = deepcopy(Za[q, k - 1, y2])  # deepcopy(Za[k - 1, q, y1])  ## ??
                    end
                end
                if y2 <= m
                    previous_value = Ca[q, k, y2 + 1, 0]
                    Ca[q, k, y2 + 1, 0] = max(Ca[q, k, y2 + 1, 0], max(Ca[q - 1, k, y1, 0], Ca[q, k - 1, y2, 0]) + P_bar[q, pi[k]] + P_hat[q, pi[k]])
                    if abs(Ca[q, k, y2 + 1, 0] - previous_value) > EPS  # Ca[q, k, y2 + 1, 0] != previous_value
                        if abs(Ca[q, k, y2 + 1, 0] - (Ca[q - 1, k, y1, 0] + P_bar[q, pi[k]] + P_hat[q, pi[k]])) < EPS
                            Za[q, k, y2 + 1] = deepcopy(Za[q - 1, k, y1])
                            push!(Za[q, k, y2 + 1], q)
                        elseif abs(Ca[q, k, y2 + 1, 0] - (Ca[q, k - 1, y2, 0] + P_bar[q, pi[k]] + P_hat[q, pi[k]])) < EPS
                            Za[q, k, y2 + 1] = deepcopy(Za[q, k - 1, y2])  # deepcopy(Za[k - 1, q, y1])  # deepcopy(Za[k, q - 1, y2])
                            push!(Za[q, k, y2 + 1], q)
                        else
                            println("ERROR!")
                            ###Za[q, k, y2 + 1] = deepcopy(Za[q - 1, k, y2 + 1])  ## ??
                        end
                    end
                end
            end
        end
        println("Deviated jobs for q = 1: $(Za[q, n, ceil_Tq])")
        # Calculate the completion time matrix for jobs pos <= k, given deviated_items
        best_dev = []
        best_value = 0
        for r in 1:m
            for gamma in ceil_Tk:ceil_Tk
                println("Ca = $(Ca[k, r, gamma, 0])")
                print("Za[$(k), $(r), $(gamma)] = $(Za[k, r, gamma]) => ")
                temp_deviated_items = deepcopy(deviated_items)
                push!(temp_deviated_items, deepcopy(Za[k, r, gamma]))
                temp_ca = calculate_completion_array(m, n, pi, P_bar, P_hat, temp_deviated_items, k)
                wct = sum(w[pi[k]] * temp_ca[m, k] for k in 1:n)
                println(" $(temp_ca[m, k]) ; wct = $(wct)")
                if wct > best_value
                    best_dev = deepcopy(Za[k, r, gamma])
                    print(" [Improved r = $(r), gamma = $(gamma)] ")
                end
            end
        end
        push!(deviated_items, deepcopy(best_dev))  # Za[k, m, ceil_Tk])
        println("\nDeviated machines for k = $(k): $(deviated_items[k])")
        bestCa = calculate_completion_array(m, n, pi, P_bar, P_hat, deviated_items, k)
    end
    worst_wct = sum(w[pi[k]] * Ca[k, m, Int64(ceil(Γ[pi[k]])), 0] for k in 1:n)  # Was : ... for q in 1:(m-2)   x and x     for q in (m-1):m
    println("[DP] The m-machine worst-case Cmax is $(worst_wct).")
    worst_Γ_scenario = calculate_worst_Γ_scenario_wct(pi, m, n, P_bar, P_hat, w, Γ, Ca, Za, worst_wct)  # [ [[] for q in 2:m] ]
    if isempty(worst_Γ_scenario)
        worst_Γ_scenario = [Int64[] for x in 1:m]
    end
    println("[DP] The m-machine worst-case Cmax is $(worst_wct), with worst_Γ_scenario = $(worst_Γ_scenario).")
    worst_Γ_dev = calculate_deviation_matrix(m, n, pi, Γ, worst_Γ_scenario, [], "job")
    wct_validation = calculate_wct_given_scenario(pi, m, n, P_bar, P_hat, w, worst_Γ_dev)
    @assert abs(wct_validation - max_wct[m]) < 1e-5
    return worst_wct, worst_Γ_scenario
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

function test_wct()
    searchdir(path,key) = filter(x->occursin(key,x), readdir(path))
    basedir = joinpath(robust_ying_instances_folder, "data")
    for filename in searchdir(basedir, ".txt")
        if filename == "RB0104010_20perc_3m.txt"
            inputfile = "$(basedir)/$(filename)"
            m, n, P_bar, P_hat = read_robust_ying_input_file(inputfile)
            pi = reverse([x for x in 1:n])
            w = [1.0 for x in 1:n]
            deviated_jobs = [ [] for i in 1:m ]
            Γ = [0.4 * n for x in 1:m]
            worst_wct, worst_Γ_scenario = worst_case_wct_dp_m_machines(pi, m, n, Γ, P_bar, P_hat, w, true, true)
            println("worst_wct = $(worst_wct) ; worst_Γ_scenario = $(worst_Γ_scenario)")
        end
    end
end

#test_wct()
