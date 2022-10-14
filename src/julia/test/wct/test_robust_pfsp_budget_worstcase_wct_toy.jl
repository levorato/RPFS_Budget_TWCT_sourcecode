# TEST robust_pfsp_budget_worstcase_wct.jl
# Worst-case calculation for the robust budget problem
# Weighed Completion Time Objective Function

include("../../wct/robust_pfsp_budget_worstcase_wct.jl")
include("../../wct/robust_pfsp_wilson_worstcase_wct_mip.jl")
include("../../robust_pfsp_budget_worstcase_m_mach.jl")
include("../../wct/robust_pfsp_alt_worstcase_wct_mip.jl")
include("../../config.jl")

using Combinatorics

function test_m_machine_budget_worst_case_dp_toy(desired_filename)
    basedir = robust_ying_instances_folder
    for filename in searchdir(basedir, ".txt")
        if filename == desired_filename
            inputfile = "$(basedir)/$(filename)"
            m, n, P_bar, P_hat = read_robust_ying_input_file(inputfile)
            Jobs = 1:n
            N = Jobs
            Machines = 1:m
            M = Machines
            println("Instance has n = $(n) and m = $(m).")
            w = [1.0 for x in 1:n]  # For now, assume all job weights are equal to 1
            error_count = 0
            N_array = [x for x in 1:n]
            all_pi_permutations = permutations(N_array, n)
            Γ = [2 for x in 1:m]
            for pi in all_pi_permutations
                println("Calculating worst-case makespan for sequence pi = $(pi) and Γ = $(Γ)...")
                wct_dp, worst_Γ_scenario_dp = worst_case_wct_dp_m_machines_2(pi, m, n, Γ, P_bar, P_hat, w, false)
                calculate_wct_given_scenario(pi, m, n, P_bar, P_hat, w, worst_Γ_scenario_dp[1], true)
                if m == 3
                    wct_brute, worst_Γ_scenario_brute, worst_j = worst_case_wct_brute_force_3_machines(pi, m, n, Γ, P_bar, P_hat, w, true)
                    calculate_wct_given_scenario(pi, m, n, P_bar, P_hat, w, last(worst_Γ_scenario_brute), true)
                elseif m == 2
                    wct_brute, worst_Γ_scenario_brute = worst_case_wct_brute_force_2_machines(pi, m, n, Γ, P_bar, P_hat, w, true)
                    calculate_wct_given_scenario(pi, m, n, P_bar, P_hat, w, last(worst_Γ_scenario_brute), true)
                end
                #wct_mip, worst_Γ_scenario_mip = solve_robust_model_worst_case_wilson(pi, m, n, Γ, P_bar, P_hat, w)
                #calculate_wct_given_scenario(pi, m, n, P_bar, P_hat, w, worst_Γ_scenario_mip, true)
                #cmax_dp, worst_Γ_scenario_cmax = worst_case_cmax_dp_m_machines(pi, m, n, Γ, P_bar, P_hat)
                #calculate_wct_given_scenario(pi, m, n, P_bar, P_hat, w, worst_Γ_scenario_cmax, true)

                if abs(wct_dp - wct_brute) > EPS
                    println("ERROR : wct_dp incorrect : wct_brute = $(wct_brute) x wct_dp = $(wct_dp)")
                    error_count += 1
                    break
                end
                println("wct_brute = $(wct_brute) ; wct_dp = $(wct_dp)")
            end
        end
    end
end

function test_m_machine_budget_worst_case_dp_toy_job_budget(desired_filename)
    basedir = robust_ying_instances_folder
    for filename in searchdir(basedir, ".txt")
        if filename == desired_filename
            inputfile = "$(basedir)/$(filename)"
            m, n, P_bar, P_hat = read_robust_ying_input_file(inputfile)
            Jobs = 1:n
            N = Jobs
            Machines = 1:m
            M = Machines
            println("Instance has n = $(n) and m = $(m).")
            w = [1.0 for x in 1:n]  # For now, assume all job weights are equal to 1
            error_count = 0
            N_array = [x for x in 1:n]
            all_pi_permutations = permutations(N_array, n)
            all_pi_permutations = [[2, 1, 3, 4], [4, 3, 2, 1]]
            #all_pi_permutations = [[4, 3, 2, 1]]
            Γ = [2 for x in 1:n]
            for pi in all_pi_permutations
                println("Calculating worst-case makespan for sequence pi = $(pi) and Γ = $(Γ)...")
                wct_mip_wilson, worst_Γ_scenario_mip_wilson, elapsed_time = solve_robust_model_worst_case_wilson(pi, m, n, Γ, P_bar, P_hat, w)
                calculate_wct_given_scenario(pi, m, n, P_bar, P_hat, w, worst_Γ_scenario_mip_wilson[1], true)
                wct_mip_alt, worst_Γ_scenario_mip_alt, elapsed_time = solve_robust_model_worst_case_alt(pi, m, n, Γ, P_bar, P_hat, w)
                calculate_wct_given_scenario(pi, m, n, P_bar, P_hat, w, worst_Γ_scenario_mip_alt[1], true)
                wct_brute, worst_Γ_scenario_brute, worst_r_frac_brute, worst_dev = worst_case_wct_brute_force_3_machines(pi, m, n, Γ, P_bar, P_hat, w, true)
                calculate_wct_given_scenario(pi, m, n, P_bar, P_hat, w, last(worst_dev), true)
                ###wct_dp, worst_Γ_scenario_dp = worst_case_wct_dp_m_machines(pi, m, n, Γ, P_bar, P_hat, w)
                ###calculate_wct_given_scenario(pi, m, n, P_bar, P_hat, w, worst_Γ_scenario_dp[1], true)

                #wct_mip, worst_Γ_scenario_mip = solve_robust_model_worst_case_wilson(pi, m, n, Γ, P_bar, P_hat, w)
                #calculate_wct_given_scenario(pi, m, n, P_bar, P_hat, w, worst_Γ_scenario_mip, true)
                #cmax_dp, worst_Γ_scenario_cmax = worst_case_cmax_dp_m_machines(pi, m, n, Γ, P_bar, P_hat)
                #calculate_wct_given_scenario(pi, m, n, P_bar, P_hat, w, worst_Γ_scenario_cmax, true)

                #if abs(wct_dp - wct_brute) > EPS
                #    println("ERROR : wct_dp incorrect : wct_brute = $(wct_brute) x wct_dp = $(wct_dp)")
                #    error_count += 1
                #    break
                #end
                #println("wct_brute = $(wct_brute) ; wct_dp = $(wct_dp)")
            end
        end
    end
end

#test_m_machine_budget_worst_case_dp_toy("RB0104010_20perc_3m.txt")  # "test_toy_3m.txt"
test_m_machine_budget_worst_case_dp_toy_job_budget("test_toy_3m.txt")
