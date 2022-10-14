# Return a [toplogical sort](https://en.wikipedia.org/wiki/Topological_sorting) of a directed
# graph `g` as a vector of vertices in topological order.
# nv = number of vertices
function topological_sort_by_dfs(g, nv)
    vcolor = zeros(UInt8, nv)
    verts = Vector{Int64}()
    for v in 1:nv
        vcolor[v] != 0 && continue
        S = Vector{Int64}([v])
        vcolor[v] = 1
        while !isempty(S)
            u = S[end]
            w = 0
            for n in 1:nv  # out_neighbors(g, u)
                if g[u, n] == 1
                    if vcolor[n] == 1
                        error("The input graph contains at least one loop.")
                    elseif vcolor[n] == 0
                        w = n
                        break
                    end
                end
            end
            if w != 0
                vcolor[w] = 1
                push!(S, w)
            else
                vcolor[u] = 2
                push!(verts, u)
                pop!(S)
            end
        end
    end
    return reverse(verts)
end
