# =============================================================================
# run_bb_rpfs_adrs_interval.jl
# =============================================================================
# Julia script to run the Branch and Bound procedure for the Robust PFS Problem
# (ADRS Cmax objective,interval processing time),based on Kouvelis et al.(2000).
# =============================================================================

using Random
using UUIDs
using Dates
include("../../src/julia/config.jl")
include("../../src/julia/pfsp_file_reader.jl")

EXPERIMENT_NAME = "run_bb_rpfs_adrs_interval"

# If possible, do not modify the lines below
EXPERIMENT_OUTPUT_FOLDER = abspath(create_full_dir(normpath(output_folder),
	[EXPERIMENT_NAME]))

function run_bb_pfsp_adrs(file_list, time_limit, dominance_rules)
    rng = MersenneTwister(1234)
    uid = uuid4(rng)
    println("$(uid)")
    batchId = Dates.format(Dates.now(), "yyyy_mm_dd-HH_MM_SS") * "-$(uid)"
    output_folder2 = abspath(create_full_dir(EXPERIMENT_OUTPUT_FOLDER, [batchId]))
	output_folder_stdout = abspath(create_full_dir(output_folder2, ["stdout"]))
    for filepath in file_list
		filename = filepath[findlast('/', filepath)+1:length(filepath)]
        println("=========================================================")
        println("      Solving instance $(filename)...")
        println("=========================================================")
		fullfilepath = joinpath(instances_folder, filepath)
        # Create an output text file to save the stdout and sterr process output
		stdout_file_path = joinpath(output_folder_stdout, "$(filename)_stdout.txt")
	    stdout_file = open(stdout_file_path, "a+")
		exit_code = -500
		try
			execpath = abspath(joinpath(project_folder, "bin", "pfsp"))
            cmd = `$(execpath) --input-file="$(fullfilepath)"
                --output-folder $(output_folder2)
                --model-version=101
				--time-limit=$(time_limit)
				--dominance=$(dominance_rules)`
			println("cmd: $(cmd)")
			err = Pipe()
		    out = Pipe()
		    proc = run(pipeline(ignorestatus(cmd),stdout=out,stderr=err))
		    close(err.in)
		    close(out.in)
			stdout_str = String(read(out))
		    println(stdout_file, "*** stdout: ",stdout_str, "\n")
		    println(stdout_file, "*** stderr: ",String(read(err)), "\n")
		    println(stdout_file, "*** exit status: ", proc.exitcode, "\n")
			flush(stdout_file)
			exit_code = proc.exitcode
        catch e
			bt = catch_backtrace()
			msg = sprint(showerror, e, bt)
			println(msg)
			println(stdout_file, "*** EXCEPTION OCCURED: $(msg)")
			flush(stdout_file)
        end
		close(stdout_file)
    end
    println("DONE.")
end

function generate_instance_list(instance_set, size_list)
    instance_folder_list = []
	if instance_set == "ying"
	    for folder in ["robust/ying/data/$(num_jobs)jobs" for num_jobs in size_list]
	        instance_folder_list = [ instance_folder_list ; [folder * "/$(alpha)%" for alpha in [10, 20, 30, 40, 50]] ]
	    end
	end
	if instance_set == "kouvelis"
    	instance_folder_list = [instance_folder_list ; ["robust/kouvelis/test/interval"]]
	end
    println("Instances folders to process: $(instance_folder_list)")
    # Randomly extract 20% of instances from each instance folder
    instances_list = []
    instance_count = 0
    for folder in instance_folder_list
        instances = ["$(joinpath(folder, file))" for file in searchdir(joinpath(instances_folder, folder), ".txt")]
        instance_count += length(instances)
        instances_list = [ instances_list; instances ]
    end
    println("All instances count: $(instance_count)")
    return instances_list
end

# Ying
#file_list = generate_instance_list("ying", ["10", "20", "50"])
file_list = generate_instance_list("ying", ["100", "150"])
#file_list = generate_instance_list("ying", ["150"])
#file_list = generate_instance_list("ying", ["200"])
# Kouvelis
#file_list = generate_instance_list("kouvelis", [])

time_limit = 7200
run_bb_pfsp_adrs(file_list, time_limit, "true")
#run_bb_pfsp_adrs(file_list, time_limit, "false")
