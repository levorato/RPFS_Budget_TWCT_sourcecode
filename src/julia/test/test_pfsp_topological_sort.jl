# ==============================================================================================
# test_pfsp_topological_sort.jl
# Test topological sort to infer a PFSP solution (permutation) from a binary MILP variable
# that belongs to either Manne or Liao-You PFSP MILP models.
# ==============================================================================================

include("../pfsp_util.jl")
using Random

function create_partial_order_matrix_from_permutation(perm)
    n = size(perm, 1)
    d = zeros(Int64, (n, n))
    for pos_k in 1:n  # for all jobs k in position pos_k in perm
        k = perm[pos_k]
        # d[i][k] = 1, if job i is scheduled any time before job k; 0, otherwise
        # Assign d[i][k] = 1 to all jobs i that come before job k in permutation seq_1
        for pos_i in 1:(pos_k-1)  # for each position pos_i < pos_k
            #break
            # job i is in position pos_i, and comes before job k
            i = perm[pos_i]
            # i = 1..n-1 , k = i+1..n
            if ((i < n - 1) && (k > i))
                d[i, k] = 1
            end
            if ((k < n - 1) && (i > k))
                d[k, i] = 0
            end
        end
        # for each job i that comes after job k
        for pos_i in (pos_k+1):n  # for each position pos_i > pos_k
            #break
            i = perm[pos_i]
            if ((k < n - 1) && (i > k))
                d[k, i] = 1
            end
            if ((i < n - 1) && (k > i))
                d[i, k] = 0
            end
        end
    end
    return d
end

function test_topological_sort()
    n = 10
    original_permutation = shuffle(Int64[x for x in 1:n])
    println("\noriginal_permutation = $(original_permutation)")
    d = create_partial_order_matrix_from_permutation(original_permutation)
    println("d = $(d)")
    new_permutation = topological_sort_by_dfs(d, n)
    println("new_permutation = $(new_permutation)")
end

test_topological_sort()
