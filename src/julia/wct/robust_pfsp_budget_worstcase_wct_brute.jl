# ===================================================================================================
# robust_pfsp_budget_worstcase_wct_brute.jl
# Worst-case determination procedure for the Robust Permutation Flow Shop Problem to Minimize the
# Weighted Sum of Job Completion Times (WCT), based on the budgeted uncertainty set (Bertsimas et al.).
# Given a permutation pi, finds the restricted worst-case scenario (i.e. on each machine r,
# which jobs will have their processing times deviated to their worst-case values), so that the
# objective function (weighted completion time - wct) reaches the worst-case (maximum) value.
# *** Brute force approach for calculating the worst-case scenario.
# ===================================================================================================
#include("../robust_pfsp_budget_worstcase_m_mach.jl")

using Combinatorics
using OffsetArrays  # Fortran-like arrays with arbitrary, zero or negative starting indices


# The criterion under consideration is sum_i (w_i * C_i), which is the
# the sum of the weighted job completion times.
function calculate_wct(m, n, p, w, pi, verbose = false)
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
    if verbose
        println("[WCT] C = $(C)")
    end
    return wct
end

function calculate_wct_given_scenario(pi, m, n, P_bar, P_hat, w, scenario, verbose = false)
    # The scenario parameter contains the amount of deviation of each job j on machine i
    dev = scenario
    for i in 1:m
        for j in 1:n
            P_bar[i, j] += dev[i, j] * P_hat[i, j]
        end
    end
    wct = calculate_wct(m, n, P_bar, w, pi, verbose)
    for i in 1:m
        for j in 1:n
            P_bar[i, j] -= dev[i, j] * P_hat[i, j]
        end
    end
    if verbose
        println("Given scenario $(scenario), the wct is $(wct).")
    end
    return wct
end

function calculate_deviation_matrix(m, n, pi, Γ, dev_full, dev_frac, budget_type)
    dev = zeros(m, n)
    if budget_type == "machine"
        i = 1
        for m_dev in dev_full  # for each machine i
            for j in m_dev  # for each job j
                dev[i, j] = 1
            end
            i += 1
        end
        i = 1
        for m_dev in dev_frac  # for each machine i
            for j in m_dev  # for each job j
                if ceil(Γ[i]) - floor(Γ[i]) > 0
                    dev[i, j] = ceil(Γ[i]) - floor(Γ[i])
                end
            end
            i += 1
        end
    else    # budget_type == "job"
        k = 1
        for k_dev in dev_full  # for each job in position k
            for r in k_dev
                dev[r, pi[k]] = 1
            end
            k += 1
        end
        k = 1
        for k_dev in dev_frac
            if ceil(Γ[k]) - floor(Γ[k]) > 0
                for r in k_dev
                    dev[r, pi[k]] = ceil(Γ[k]) - floor(Γ[k])
                end
            end
            k += 1
        end
    end
    return dev
end

# Calculates the worst-case weighted completion time, given a sequence / permutation pi,
# budget parameters Γ1 and Γ2 and a weight vector w. Uses brute force algorithm.
function worst_case_wct_brute_force_2_machines(pi, n, p_Γ1, p_Γ2, P_bar, P_hat, w, verbose = true)
    m = 2
    N_array = [x for x in 1:n]
    N_set = Set([x for x in 1:n])
    # e in N_minus_N1_tilde <==> !in(e, N1_tilde)
    Γ = [p_Γ1, p_Γ2]
    Γ1 = ceil(Int64, p_Γ1)
    Γ2 = ceil(Int64, p_Γ2)
    all_Γ1_combinations = combinations(N_array, Γ1)
    all_Γ2_combinations = combinations(N_array, Γ2)
    max_fo = 0
    best_Γ = []
    worst_dev = []
    all_cmax = Float64[]
    println("[Brute Force] Calculating worst-case WCT given a permutation pi = $(pi), m = $(m), n = $(n), w = $(w) and budget parameters Γ = $((Γ1, Γ2))...")
    p = P_bar
    worst_p_matrix = []
    for N1_tilde_arr in all_Γ1_combinations
        for j in N1_tilde_arr
            p[1, j] += P_hat[1, j]
        end
        for N2_tilde_arr in all_Γ2_combinations
            for j in N2_tilde_arr
                p[2, j] += P_hat[2, j]
            end
            fo = calculate_wct(m, n, p, w, pi)
            #println("*** Processing Sets : N1_tilde = $(N1_tilde_arr); N2_tilde = $(N2_tilde_arr): cmax = $(fo)")
            push!(all_cmax, fo)
            if fo > max_fo
                max_fo = fo
                best_Γ = [[N1_tilde_arr, N2_tilde_arr]]
                worst_dev = [calculate_deviation_matrix(2, n, pi, Γ, [N1_tilde_arr, N2_tilde_arr], [], "machine")]
                if verbose  worst_p_matrix = deepcopy(p)  end
            elseif abs(fo - max_fo) < EPS
                push!(best_Γ, [N1_tilde_arr, N2_tilde_arr])
                push!(worst_dev, calculate_deviation_matrix(2, n, pi, Γ, [N1_tilde_arr, N2_tilde_arr], [], "machine"))
            end
            for j in N2_tilde_arr  # undo
                p[2, j] -= P_hat[2, j]
            end
        end
        for j in N1_tilde_arr  # undo
            p[1, j] -= P_hat[1, j]
        end
    end
    #println("Found the following makespans: $(all_cmax).")
    if verbose
        println("[Brute Force] The worst-case WCT is $(max_fo) with permutation $(pi) and budget (s_Γ1, s_Γ2) = $(best_Γ).")
        #calculate_wct(m, n, worst_p_matrix, w, pi, false)
        wct = calculate_wct_given_scenario(pi, m, n, P_bar, P_hat, w, worst_dev[1], false )     # true)
        println("[Validation] Given scenario = $(best_Γ[1]), the WCT is $(wct).")
    end
    println("[Brute Force] The worst-case WCT is $(max_fo) ; permutation $(pi)")
    return max_fo, best_Γ, worst_dev
end

# Calculates the worst-case WCT given a sequence / permutation pi,
# budget parameters Γ[1], Γ[2] and Γ[3] and a weight vector w.
# Uses brute force algorithm.
function worst_case_wct_brute_force_3_machines(pi, m, n, Γ, P_bar, P_hat, w, verbose = false)
    N_array = [x for x in 1:n]
    N_set = Set([x for x in 1:n])
    # e in N_minus_N1_tilde <==> !in(e, N1_tilde)
    #Γ = [ceil(Int64, x) for x in p_Γ]
    delta_T1 = ceil(Γ[1]) - floor(Γ[1])
    delta_T2 = ceil(Γ[2]) - floor(Γ[2])
    delta_T3 = ceil(Γ[3]) - floor(Γ[3])
    ceil_T1 = Int64(ceil(Γ[1]))
    ceil_T2 = Int64(ceil(Γ[2]))
    ceil_T3 = Int64(ceil(Γ[3]))
    all_Γ1_combinations = combinations(N_array, ceil_T1)
    all_Γ2_combinations = combinations(N_array, ceil_T2)
    all_Γ3_combinations = combinations(N_array, ceil_T3)
    max_fo = 0
    worst_Γ_scenario = []
    worst_j = []
    worst_dev = []
    println("[Brute Force] Calculating worst-case WCT given a permutation pi = $(pi), m = $(m), n = $(n), w = $(w) and budget parameters Γ = $(Γ)...")
    p = P_bar
    worst_p_matrix = []
    for N1_tilde_arr in all_Γ1_combinations
        if isempty(N1_tilde_arr)  # Ignores the case when T[1] == 0
            N1_tilde_arr = [-1]
        end
        for j1 in N1_tilde_arr  # *** job j1 in Machine M_1 will only vary a % to its worst-case processing time
            if ceil_T1 > 0  # Ignores the case when T[1] == 0
                for j in N1_tilde_arr
                    if j != j1 || delta_T1 == 0
                        p[1, j] += P_hat[1, j]
                    else
                        p[1, j] += (delta_T1 * P_hat[1, j] )
                    end
                end
            end
            for N2_tilde_arr in all_Γ2_combinations
                if isempty(N2_tilde_arr)  # Treats the case when T[2] == 0
                    N2_tilde_arr = [-1]
                end
                for j2 in N2_tilde_arr  # *** job j2 in Machine M_2 will only vary a % to its worst-case processing time
                    if ceil_T2 > 0  # Ignores the case when T[2] == 0
                        for j in N2_tilde_arr
                            if j != j2  || delta_T2 == 0
                                p[2, j] += P_hat[2, j]
                            else
                                p[2, j] += (delta_T2 * P_hat[2, j] )
                            end
                        end
                    end
                    for N3_tilde_arr in all_Γ3_combinations
                        if isempty(N3_tilde_arr)  # Treats the case when T[3] == 0
                            N3_tilde_arr = [-1]
                        end
                        for j3 in N3_tilde_arr  # *** job j3 in Machine M_3 will only vary a % to its worst-case processing time
                            if ceil_T3 > 0  # Ignores the case when T[3] == 0
                                for j in N3_tilde_arr
                                    if j != j3  || delta_T3 == 0
                                        p[3, j] += P_hat[3, j]
                                    else
                                        p[3, j] += (delta_T3 * P_hat[3, j] )
                                    end
                                end
                            end
                            fo = calculate_wct(m, n, p, w, pi)
                            #println("*** Processing Sets : N1_tilde = $(N1_tilde_arr); N2_tilde = $(N2_tilde_arr): cmax = $(fo)")
                            if fo > max_fo
                                max_fo = fo
                                worst_Γ_scenario = [[N1_tilde_arr, N2_tilde_arr, N3_tilde_arr]]
                                worst_j = [[j1, j2, j3]]
                                worst_dev = [calculate_deviation_matrix(3, n, Γ, pi, [N1_tilde_arr, N2_tilde_arr, N3_tilde_arr], [j1, j2, j3], "machine")]
                            elseif abs(fo - max_fo) < EPS
                                push!(worst_Γ_scenario, [N1_tilde_arr, N2_tilde_arr, N3_tilde_arr])
                                push!(worst_j, [j1, j2, j3])
                                push!(worst_dev, calculate_deviation_matrix(3, n, pi, Γ, [N1_tilde_arr, N2_tilde_arr, N3_tilde_arr], [j1, j2, j3], "machine"))
                            end
                            if ceil_T3 > 0  # Ignores the case when T[3] == 0
                                for j in N3_tilde_arr  # undo
                                    if j != j3  || delta_T3 == 0
                                        p[3, j] -= P_hat[3, j]
                                    else
                                        p[3, j] -= (delta_T3 * P_hat[3, j] )
                                    end
                                end
                            end
                        end
                    end
                    if ceil_T2 > 0  # Ignores the case when T[2] == 0
                        for j in N2_tilde_arr  # undo
                            if j != j2 || delta_T2 == 0
                                p[2, j] -= P_hat[2, j]
                            else
                                p[2, j] -= (delta_T2 * P_hat[2, j] )
                            end
                        end
                    end
                end
            end
            if ceil_T1 > 0  # Ignores the case when T[1] == 0
                for j in N1_tilde_arr  # undo
                    if j != j1 || delta_T1 == 0
                        p[1, j] -= P_hat[1, j]
                    else
                        p[1, j] -= (delta_T1 * P_hat[1, j] )
                    end
                end
            end
        end
    end
    #println("Found the following makespans: $(all_cmax).")
    if verbose
        println("[Brute Force] The worst-case WCT is $(max_fo) ; permutation $(pi) and scenario (s_Γ) = $(worst_Γ_scenario); (s_j) = $(worst_j).")
        #calculate_wct(m, n, worst_p_matrix, w, pi, false)
        wct = calculate_wct_given_scenario(pi, m, n, P_bar, P_hat, w, worst_dev[1], false )
        println("[Validation] Given scenario = $(worst_Γ_scenario[1]), the WCT is $(wct).")
    else
        println("[Brute Force] The worst-case WCT is $(max_fo) ; permutation $(pi)")
    end
    return max_fo, worst_Γ_scenario, worst_j, worst_dev
end


# Calculates the worst-case WCT given a sequence / permutation pi,
# job budget parameters Γ[1], Γ[2], Γ[3] and Γ[4] (for jobs 1..4) and a weight vector w.
# Uses brute force algorithm.
function worst_case_wct_brute_force_4_jobs(pi, m, n, Γ, P_bar, P_hat, w, verbose = false)
    M_array = [x for x in 1:m]
    M_set = Set([x for x in 1:m])
    # e in N_minus_N1_tilde <==> !in(e, N1_tilde)
    #Γ = [ceil(Int64, x) for x in p_Γ]
    delta_T1 = ceil(Γ[pi[1]]) - floor(Γ[pi[1]])  # budget for the 1st job in permutation pi : pi(1)
    delta_T2 = ceil(Γ[pi[2]]) - floor(Γ[pi[2]])  # budget for the 2nd job in permutation pi : pi(2)
    delta_T3 = ceil(Γ[pi[3]]) - floor(Γ[pi[3]])  # budget for the 3rd job in permutation pi : pi(3)
    delta_T4 = ceil(Γ[pi[4]]) - floor(Γ[pi[4]])  # budget for the 4th job in permutation pi : pi(4)
    ceil_T1 = Int64(ceil(Γ[pi[1]]))
    ceil_T2 = Int64(ceil(Γ[pi[2]]))
    ceil_T3 = Int64(ceil(Γ[pi[3]]))
    ceil_T4 = Int64(ceil(Γ[pi[4]]))
    all_Γ1_combinations = combinations(M_array, ceil_T1)
    all_Γ2_combinations = combinations(M_array, ceil_T2)
    all_Γ3_combinations = combinations(M_array, ceil_T3)
    all_Γ4_combinations = combinations(M_array, ceil_T4)
    max_fo = 0
    worst_Γ_scenario = []
    worst_r_frac = []
    worst_dev = []
    println("[Brute Force] Calculating worst-case WCT given a permutation pi = $(pi), m = $(m), n = $(n), w = $(w) and budget parameters Γ = $(Γ)...")
    p = P_bar
    worst_p_matrix = []
    for M1_tilde_arr in all_Γ1_combinations  # All combinations of machines where job pi(1) will deviate its processing time
        if isempty(M1_tilde_arr)  # Ignores the case when T[1] == 0
            M1_tilde_arr = [-1]
        end
        for r1 in M1_tilde_arr  # *** job (r1, pi(1)) in Machine M_r will only vary a fraction (%) to its worst-case processing time
            if ceil_T1 > 0  # Ignores the case when T[1] == 0
                for r in M1_tilde_arr
                    if r != r1 || delta_T1 == 0
                        p[r, pi[1]] += P_hat[r, pi[1]]
                    else
                        p[r, pi[1]] += (delta_T1 * P_hat[r, pi[1]])
                    end
                end
            end
            for M2_tilde_arr in all_Γ2_combinations  # All combinations of machines where job pi(2) will deviate its processing time
                if isempty(M2_tilde_arr)  # Treats the case when T[2] == 0
                    M2_tilde_arr = [-1]
                end
                for r2 in M2_tilde_arr  # *** job (r2, pi(2)) in Machine M_r will only vary a fraction (%) to its worst-case processing time
                    if ceil_T2 > 0  # Ignores the case when T[2] == 0
                        for r in M2_tilde_arr
                            if r != r2 || delta_T2 == 0
                                p[r, pi[2]] += P_hat[r, pi[2]]
                            else
                                p[r, pi[2]] += (delta_T2 * P_hat[r, pi[2]])
                            end
                        end
                    end
                    for M3_tilde_arr in all_Γ3_combinations  # All combinations of machines where job pi(3) will deviate its processing time
                        if isempty(M3_tilde_arr)  # Treats the case when T[3] == 0
                            M3_tilde_arr = [-1]
                        end
                        for r3 in M3_tilde_arr  # *** job (r3, pi(3)) in Machine M_r will only vary a fraction (%) to its worst-case processing time
                            if ceil_T3 > 0  # Ignores the case when T[3] == 0
                                for r in M3_tilde_arr
                                    if r != r3 || delta_T3 == 0
                                        p[r, pi[3]] += P_hat[r, pi[3]]
                                    else
                                        p[r, pi[3]] += (delta_T3 * P_hat[r, pi[3]])
                                    end
                                end
                            end
                            for M4_tilde_arr in all_Γ4_combinations  # All combinations of machines where job pi(4) will deviate its processing time
                                if isempty(M4_tilde_arr)  # Treats the case when T[4] == 0
                                    M4_tilde_arr = [-1]
                                end
                                for r4 in M4_tilde_arr  # *** job (r4, pi(4)) in Machine M_r will only vary a fraction (%) to its worst-case processing time
                                    if ceil_T4 > 0  # Ignores the case when T[4] == 0
                                        for r in M4_tilde_arr
                                            if r != r4 || delta_T4 == 0
                                                p[r, pi[4]] += P_hat[r, pi[4]]
                                            else
                                                p[r, pi[4]] += (delta_T3 * P_hat[r, pi[4]])
                                            end
                                        end
                                    end
                                    fo = calculate_wct(m, n, p, w, pi)
                                    #println("*** Processing Sets : N1_tilde = $(N1_tilde_arr); N2_tilde = $(N2_tilde_arr): cmax = $(fo)")
                                    if fo > max_fo
                                        max_fo = fo
                                        worst_Γ_scenario = [[M1_tilde_arr, M2_tilde_arr, M3_tilde_arr, M4_tilde_arr]]
                                        worst_r_frac = [[r1, r2, r3, r4]]
                                        worst_dev = [calculate_deviation_matrix(m, n, pi, Γ, worst_Γ_scenario[1], worst_r_frac[1], "job")]
                                    elseif abs(fo - max_fo) < EPS
                                        push!(worst_Γ_scenario, [M1_tilde_arr, M2_tilde_arr, M3_tilde_arr, M4_tilde_arr])
                                        push!(worst_r_frac, [r1, r2, r3, r4])
                                        push!(worst_dev, calculate_deviation_matrix(m, n, pi, Γ, [M1_tilde_arr, M2_tilde_arr, M3_tilde_arr, M4_tilde_arr], [r1, r2, r3, r4], "job"))
                                    end
                                    if ceil_T4 > 0  # Ignores the case when T[4] == 0
                                        for r in M4_tilde_arr
                                            if r != r4 || delta_T4 == 0
                                                p[r, pi[4]] -= P_hat[r, pi[4]]
                                            else
                                                p[r, pi[4]] -= (delta_T3 * P_hat[r, pi[4]])
                                            end
                                        end
                                    end
                                end
                            end
                            if ceil_T3 > 0  # Ignores the case when T[3] == 0
                                for r in M3_tilde_arr  # undo
                                    if r != r3 || delta_T3 == 0
                                        p[r, pi[3]] -= P_hat[r, pi[3]]
                                    else
                                        p[r, pi[3]] -= (delta_T3 * P_hat[r, pi[3]])
                                    end
                                end
                            end
                        end
                    end
                    if ceil_T2 > 0  # Ignores the case when T[2] == 0
                        for r in M2_tilde_arr  # undo
                            if r != r2 || delta_T2 == 0
                                p[r, pi[2]] -= P_hat[r, pi[2]]
                            else
                                p[r, pi[2]] -= (delta_T2 * P_hat[r, pi[2]])
                            end
                        end
                    end
                end
            end
            if ceil_T1 > 0  # Ignores the case when T[1] == 0
                for r in M1_tilde_arr  # undo
                    if r != r1 || delta_T1 == 0
                        p[r, pi[1]] -= P_hat[r, pi[1]]
                    else
                        p[r, pi[1]] -= (delta_T1 * P_hat[r, pi[1]])
                    end
                end
            end
        end
    end
    if verbose
        println("[Brute Force] The worst-case WCT is $(max_fo) ; permutation $(pi) and scenario (s_Γ) = $(worst_Γ_scenario); (s_r) = $(worst_r_frac).")
        #calculate_wct(m, n, worst_p_matrix, w, pi, false)
        wct = calculate_wct_given_scenario(pi, m, n, P_bar, P_hat, w, worst_dev[1], false)
        println("[Validation] Given scenario = $(worst_Γ_scenario[1]), the WCT is $(wct).")
    else
        println("[Brute Force] The worst-case WCT is $(max_fo) ; permutation $(pi)")
    end
    return max_fo, worst_Γ_scenario, worst_r_frac, worst_dev
end


# Calculates the worst-case weighted completion time, given a sequence / permutation pi,
# global budget parameter Γ, and a weight vector w. Uses brute force algorithm.
function worst_case_wct_brute_force_global_budget(permutation, m, n, Γ, P_bar, P_hat, w, verbose = true; max_solutions = 5)
    N_array = []  # Set of all operations (r, i) concerning machine r and job i
    for r in 1:m
        for i in 1:n
            push!(N_array, (r, i))
        end
    end
    ceil_T = Int64(ceil(Γ))
    # Obtain the list of all possible scenarios with exactly Gamma oscillated operations, for a total of (N operations choose Gamma) scenarios
    all_Γ_scenarios = combinations(N_array, ceil_T)
    max_fo = 0.0
    worst_Γ_scenario = []
    println("[Brute Force] Calculating worst-case WCT given a permutation pi = $(permutation), m = $(m), n = $(n), w = $(w) and budget parameter Γ = $(Γ)...")
    p = deepcopy(P_bar)
    count = 1
    n_ops = size(N_array, 1)
    println("Total number of operations: $(n_ops)")
    # Calculate the number of combinations of n_ops, choose ceil_T elements == n_ops! / (ceil_T! (n_ops - ceil_T)! )
    total = binomial(n_ops, ceil_T)
    println("Total number of combinations: $(total)")
    prev_perc = 0
    C = zeros(Float64, (m, n))
    for scenario in all_Γ_scenarios
        for op in scenario
            r = op[1]; i = op[2]
            p[r, i] += P_hat[r, i]
        end
        # Calculate the weighted completion time for permutation
        C[1, 1] = p[1, permutation[1]]
        for i in 2:m
            C[i, 1] = C[i - 1, 1] + p[i, permutation[1]]
        end
        for j in 2:n
            C[1, j] = C[1, j - 1] + p[1, permutation[j]]
            for i in 2:m
                C[i, j] = max(C[i, j - 1], C[i - 1, j]) + p[i, permutation[j]]
            end
        end
        wct = 0.0
        for j in 1:n  # C[m, j] is the final completion time of job j
            wct += C[m, j] * w[permutation[j]]
        end
        if wct > max_fo  # Improved worst-case scenario
            max_fo = wct
            worst_Γ_scenario = [scenario]
        elseif abs(wct - max_fo) < 1e-5  # Additional worst-case scenario
			if !(scenario in worst_Γ_scenario) && (size(worst_Γ_scenario, 1) < max_solutions)
            	push!(worst_Γ_scenario, scenario)
			end
        end
        for op in scenario  # undo
            r = op[1]; i = op[2]
            p[r, i] -= P_hat[r, i]
        end
        count += 1
        # Print percentage of processed scenarios (status)
        if (100 * count / total >= prev_perc + 5)
            prev_perc = Int64(floor(100 * count / total))
            print("$(prev_perc)%    ")
        end
    end
    println("\n[Brute Force] The worst-case WCT is $(max_fo) with permutation $(permutation) and scenario (s_Γ) = $(worst_Γ_scenario).")
    return max_fo, worst_Γ_scenario
end

function calculate_wct_given_global_scenario(pi, m, n, P_bar, P_hat, w, scenario, verbose = false)
    # The scenario parameter contains the list of operations (r, i) that will deviate to
    #     the worst-case value (P_bar[r, i] + P_hat[r, i])
    for op in scenario
        r = op[1]; i = op[2]
        P_bar[r, i] += P_hat[r, i]
    end
    wct = calculate_wct(m, n, P_bar, w, pi, verbose)
    for op in scenario
        r = op[1]; i = op[2]
        P_bar[r, i] -= P_hat[r, i]
    end
    return wct
end

# Convert a job deviation list in the format [(r, i)], r in Machines, i in Jobs where job i deviates on Machine r,
# To a job deviation matrix where dev[r, i] == 1 if job i deviates on Machine r.
function convert_scenario_job_deviation_list_to_matrix(m, n, job_dev_list)
    dev = zeros(m, n)
    for (r, i) in job_dev_list
        dev[r, i] = 1
    end
    return dev
end

# Returns true if it is worth using the brute case approach when calculating
# the worst-case value (budgeted uncertainty), given an instance (m, n),
# a permutation pi and a budget Gamma.
function worst_case_brute_force_is_worth(m, n, Gamma)
    try
        n_ops = m * n
        println("Total number of operations: $(n_ops)")
        # Calculate the number of combinations of n_ops, choose Gamma elements == n_ops! / (Gamma! (n_ops - Gamma)! )
        total = binomial(n_ops, Gamma)
        println("Total number of combinations: $(total)")
        limit = 2^25  # limit for smaller values of Gamma (Gamma < (m * n) / 2)
        if (Gamma >= (m * n) / 2)
            limit = 1e6
        end
        println("Brute force limit: $(limit)")
        if (total <= limit) && ((Gamma <= 4 && n_ops <= 50) || (Gamma <= 5 && n_ops > 50))   # || (Gamma >= m * n - 4)) #(m * n > 50)
            println("Brute force approach for worst-case is worth it.")
            return true
        else
            println("Brute force approach for worst-case is NOT worth it.")
            return false
        end
    catch e
		if isa(e, OverflowError)
			println("Brute force approach for worst-case is NOT worth it.")
		else
			println("WARN: Unknown error. Returning brute force approach for worst-case is NOT worth it.")
		end
		return false
	end
end
