# ===================================================================================================
# pfsp_file_reader.jl
# Deterministic/Robust Permutation Flow Shop Problem instace file reader.
# ===================================================================================================

using DelimitedFiles

MAX_JOBS = 400  # the maximum number of jobs allowed
MAX_MACHINES = 400

function read_input_file(filepath)
    file = open(filepath)
    text = ""
    count = 1
    m = 0
    n = 0
    for ln in eachline(file)
        if length(ln) > 0 && (! occursin(ln, "!"))  # ignore comments and empty lines
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
    flush(stdout)
    return m, n, T
end

function read_robust_ying_input_file(x)
    file = open(x)
    m = 0
    i = j = 0
    P_bar = zeros(Float64, (MAX_MACHINES, MAX_JOBS))
    P_hat = zeros(Float64, (MAX_MACHINES, MAX_JOBS))
    read_m_n = 0
    for ln in eachline(file)
        if length(ln) > 0  # ignore comments and empty lines
            ln = lstrip(ln)
            if length(ln) > 0
                if occursin(ln[1], "!")
                    read_m_n = 1
                    continue
                end
                # e.g. format : 13	16	1.3	1.6  ==> p_bar_1j p_bar_2j p_hat_1j p_hat_2j
                line_array = readdlm(IOBuffer(ln))
                j += 1
                if j == 1
                    println("Read header line $(j): $(line_array)")
                    # determine the number of machines
                    if read_m_n == 1  # ! n m alpha
                        n = Int64(line_array[1])
                        m = Int64(line_array[2])
                        println("The number of jobs is $(n).")
                        println("The number of machines is $(m).")
                        P_bar = zeros(Float64, (m, n))
                        P_hat = zeros(Float64, (m, n))
                        continue
                    else
                        m = Int64(floor(length(line_array) / 2))
                        println("The number of machines is $(m).")
                        P_bar = zeros(Float64, (m, MAX_JOBS))
                        P_hat = zeros(Float64, (m, MAX_JOBS))
                    end
                end
                if read_m_n == 1
                    j = 1
                    read_m_n = 2
                end
                #println("Read line $(j): $(line_array)")
                count = 1
                for i in 1:m
                    P_bar[i, j] = line_array[count]
                    count += 1
                end
                for i in 1:m
                    P_hat[i, j] = line_array[count]
                    count += 1
                end
            end
        end
    end
    n = Int64(j)
    println("The number of jobs is $(n).")
    P_bar2 = Base.ReshapedArray(P_bar, (m,n), ())
    P_hat2 = Base.ReshapedArray(P_hat, (m,n), ())
    println("Ying robust instance file read OK : $(x).")
    flush(stdout)
    return m, n, P_bar2, P_hat2
end

function read_robust_input_file(x)
    file = open(x)
    m = 0
    n = 0
    i = j = 0
    P_bar = zeros(Float64, (MAX_MACHINES, MAX_JOBS))
    P_hat = zeros(Float64, (MAX_MACHINES, MAX_JOBS))
    w = zeros(Float64, (MAX_JOBS))
    step = 0
    for ln in eachline(file)
        if length(ln) > 0  # ignore empty lines
            line = ln
            if occursin("nJobs", line)
                println("Step 1")
                step = 1
                continue
            end
            if occursin("Weights", line)
                step = 2
                j = 1
                continue
            end
            if occursin("P_bar", line)
                step = 3
                j = 1
                continue
            end
            if occursin("P_hat", line)
                step = 4
                j = 1
                continue
            end
            ln = lstrip(ln)
            line_array = readdlm(IOBuffer(ln))
            # println("Read header line $(j): $(line) ; $(line_array)")
            # determine the number of machines
            if step == 1  # n m
                n = Int64(line_array[1])
                m = Int64(line_array[2])
                println("The number of jobs is $(n).")
                println("The number of machines is $(m).")
                P_bar = zeros(Float64, (m, n))
                P_hat = zeros(Float64, (m, n))
                w = zeros(Float64, (n))
                continue
            elseif step == 2  # job weights
                w[j] = line_array[1]
                j += 1
            elseif step == 3  # P_bar
                for i in 1:m
                    P_bar[i, j] = line_array[i]
                end
                j += 1
            elseif step == 4  # P_hat
                for i in 1:m
                    P_hat[i, j] = line_array[i]
                end
                j += 1
            end
        end
    end
    # println("M = $(m), N = $(n), P_bar = $(P_bar), P_hat = $(P_hat), w = $(w).")
    println("Flow shop robust instance file read OK : $(x).")
    flush(stdout)
    if step == 0
        return 0, 0, [], [], []
    end
    return m, n, P_bar, P_hat, w
end

function create_full_dir(basepath, folders)
    fullpath = basepath
    for folder in folders
        fullpath = joinpath(fullpath, folder)
    end
    return mkpath(fullpath)
end
