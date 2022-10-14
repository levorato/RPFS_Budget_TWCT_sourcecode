# =============================================================================
# calculate_twct_of_rpfs_solution_given_nominal_P_bar.jl
# -----------------------------------------------------------------------------
# Calculate the total weighted completion time (TWCT) of a RPFS-TWCT solution
# pi_rob under the scenario of nominal processing times (P = P_bar).
# =============================================================================

# Include file reader and util functions
include("../../../pfsp_file_reader.jl")
# Deterministic PFSP functions
include("../../../deterministic_pfsp.jl")
include("../../../config.jl")

using CSV, DataFrames

CONSOLIDATED_RESULT_DIR = abspath(joinpath(home_prefix, "pfsp_experiments", "run_ccg_rpfs_wct_global"))


function read_rpfs_results_from_csv(result_filename)
    println("Reading result file: $(result_filename)")
    df = DataFrame(CSV.File(result_filename, delim=';'))
    return df
end

function calculate_twct_of_rpfs_solution_given_nominal_P_bar()
    df = read_rpfs_results_from_csv(joinpath(CONSOLIDATED_RESULT_DIR, "RPFS_TWCT_all_results.csv"))
    output_path = joinpath(CONSOLIDATED_RESULT_DIR, "RPFS_TWCT_all_results_extended.csv")
    basedir = joinpath(instances_folder, "petro", "v3")
    # Filter the dataframe to process only the instances of the case study:
    # "instance_9_4_R100_wct_inputs.txt" and "instance_15_5_R100_wct_inputs.txt"
    df = df[(((df[!, "model"] .== "hybrid-wilson") .& (df[!, "instance_name"] .== "instance_15_5_R100_wct_inputs.txt"))
            .| ((df[!, "model"] .== "hybrid-manne") .& (df[!, "instance_name"] .== "instance_9_4_R100_wct_inputs.txt"))), :]
    nrows, ncols = size(df)
    df[!, :twct_nominal] .= 0.0
    for row in 1:nrows
        instance = df[row, :instance_name]
        #alpha = df[row, :alpha]
        #n = df[row, :n]
        filepath = joinpath("$(basedir)", instance)
        m, n, P_bar, P_hat, w = read_robust_input_file(filepath)
        permutation = df[row, :permutation]
        permutation = [parse(Int, ss) for ss in split(permutation)]
        twct_nominal = calculate_wct(m, n, P_bar, w, permutation)
        println("TWCT (based on nominal processing time matrix): $(twct_nominal)")
        df[row, :twct_nominal] = twct_nominal
    end
    CSV.write(output_path, df, delim=';')
    println("Calculation successfully saved at \"$(output_path)\".")
end

calculate_twct_of_rpfs_solution_given_nominal_P_bar()
