# ===================================================================================================
# TEST robust_pfsp_budget_wct_separation.jl
# Solve the robust problem via separation approach
# Robust Permutation Flow Shop Problem to Minimize the Weighted Sum of Job Completion Times (WCT)
# ===================================================================================================

include("../../wct/robust_pfsp_alt_worstcase_wct_mip.jl")
include("../../wct/robust_pfsp_budget_wct_separation_manne.jl")
include("../../wct/robust_pfsp_budget_wct_separation_liao_you.jl")
include("../../wct/robust_pfsp_budget_wct_separation_hybrid.jl")
include("../../config.jl")

EXPERIMENT_NAME = "run_ccg_rpfs_wct_global"
# If possible, do not modify the lines below
EXPERIMENT_OUTPUT_FOLDER = abspath(create_full_dir(normpath(output_folder),
	[EXPERIMENT_NAME]))

function test_pfsp_budget_separation_2machines_specific()
    basedir = joinpath(robust_ying_instances_folder, "rob-pfsp-wct", "10x3")
    time_limit = 7200
    for filename in searchdir(basedir, ".txt")
        if filename == "RB0101001_010_003_10_wct_inputs.txt"
            inputfile = joinpath("$(basedir)", "$(filename)")
            m, n, P_bar, P_hat, w = read_robust_input_file(inputfile)
            println("Calculating robust PFSP budget wct solution...")
            error_count = 0
            solution_values = []
            for p_Γ in [10]  # , 20, 30]
                Γ = (p_Γ * n * m) / 100.0
                println("*** p_Γ = $(p_Γ); Γ = $(Γ)")
                #cost, pi, x, run_time, is_opt, iter, gap, lb = solve_robust_model_flowshop_wct_budget_separation_hybrid(m, n, P_bar, P_hat, w, Γ, time_limit;
                #    budget_type="global", sp_bigM_type = 7, max_cores = 2, instance_filename=inputfile, use_brute_force=true)
				cost, pi, x, run_time, is_opt, iter, gap, lb = solve_robust_model_flowshop_wct_budget_separation_liao_you(m, n, P_bar, P_hat, w, Γ, time_limit;
                    budget_type="global", sp_bigM_type = 7, max_cores = 2, instance_filename=inputfile, use_brute_force=true)
                if !is_opt
                    push!(solution_values, (Γ, cost))
                    continue
                end
                wct, Γ_scenario, elapsed_time = solve_robust_model_worst_case_alt(pi, m, n, Γ, P_bar, P_hat, w, 1;
							bigM_type = 3, max_cores = 2, solver = "CPLEX", budget_type = "global", solver_time_limit = 7200,
							single_cut_per_iteration = true, fractional_solution = false, initial_scenario = nothing)
                if abs(wct - cost) > 1e-5
                    println("ERROR : wct incorrect : cost = $(cost) x wct = $(wct)")
                    error_count += 1
                    break
                end
                push!(solution_values, (Γ, wct))
            end
            println("The solution values are: $(solution_values)")
        end
    end
end


test_pfsp_budget_separation_2machines_specific()
