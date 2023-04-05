function get_Δw(p::Float64, ϵ::Float64, w1::Float64, w2::Float64, rnd::Float64)
    if rnd < p
        Δw = ϵ * w1
    else
        if w1 >= w2
            Δw = ϵ * w2
        else
            Δw = ϵ * w1
        end
    end
    Δw
end

function update_ysm!(w::Array{Float64,1}, p::Float64, ϵ::Float64, idx1::Int64, idx2::Int64, rnd::Float64)
    """
    idx1 : giver
    idx2 : receiver

    """
    w1 = w[idx1]
    w2 = w[idx2]

    Δw = get_Δw(p, ϵ, w1, w2, rnd)

    w[idx1] -= Δw
    w[idx2] += Δw

    # if rnd < p
    #     w .+= Δw/(length(w)-1)
    #     w[idx1] -= (Δw/(length(w)-1) + Δw)
    # else
    #     w[idx1] -= Δw
    #     w[idx2] += Δw
    # end

    # return Δw
end

## Complete graph (NO network)
function change_edge_to_node_CG(i, N, E)
    i_don = Int(floor(i / (N-1) - 1e-9)) + 1
    if i_don > N
        i_don = N
    end

    i_rec = i - (N-1) * (i_don-1)
    if i_rec >= i_don
        i_rec += 1
    end

    i_don, i_rec
end

function update_ysm_series!(w::Vector{Float64}, p::Float64, ϵ::Float64, N::Int64, time::Vector{Int64}, save_time::Vector{Int64})
    sigma2 = zeros(Float64, length(time))
    w_t = zeros(Float64, N, length(save_time))

    # Pick edges
    E = N * (N-1)
    inds = zeros(Int, N)
    rnds = rand(N)

    t0 = 1
    I = 1
    II = 1
    @inbounds for t = time
        @inbounds for tt = t0:t
            # rand!(inds, N_s)
            inds .= Int.(round.(rand(N).*(E-1) .+ 1))
            rnds .= rand(N)
            # ------------- t elapsed by 1 ----------------
            @inbounds for i = 1:N
                # sample!(N_s, inds, replace=false)
                i_don, i_rec = change_edge_to_node_CG(inds[i],N,E)
                update_ysm!(w, p, ϵ, i_don, i_rec, rnds[i])
            end
            # ---------------------------------------------
            if tt in save_time
                w_t[:,II] .= w
                # log_wt[:,II] .= log10.(w)
                II += 1
            end
        end

        sigma2[I] = sig2(w, 1)[1]
        # sigma2_log[I] = sig2(log10.(w), 1)[1]
        I += 1
        t0 = t + 1
    end
    sigma2, w_t
end

## Network (edge_list)
function update_ysm_series_net!(w::Array{Float64,1}, p::Float64, ϵ::Float64, N::Int64, edge_list::Array{Int64,2}, nodes::Vector{Int64}, time::Vector{Int64}, save_time::Vector{Int64})
    sigma2 = zeros(length(time))
    N_gc = length(nodes)
    E = size(edge_list, 1)
    ind_list = zeros(Int, N_gc)
    ind = collect(1:E)
    w_t = zeros(N, length(save_time))
    rnds = rand(N_gc)
    # Δw_acc = zeros(E)
    # Δw_acc_time = zeros(E, length(save_time))

    t0 = 1
    I = 1
    II = 1
    @inbounds for t = time
        @inbounds for tt = t0:t
            rand!(ind_list, ind)
            rnds .= rand(N_gc)
            # idx_list .= edge_list[ind_list]
            # --------------------- t elapsed by 1 -----------------------
            @inbounds for i = 1:N_gc

                update_ysm!(w, p, ϵ, edge_list[ind_list[i],1], edge_list[ind_list[i],2], rnds[i])
                # Δw_acc[ind_list[i]] += Δw
            end
            # ------------------------------------------------------------
            if tt in save_time
                w_t[:,II] .= w
                # Δw_acc_time[:,II] .= Δw_acc
                II += 1
            end
        end

        sigma2[I] = sig2(w[nodes], 1)[1]
        # sigma2_log[I] = sig2(log10.(w[nodes]), 1)[1]
        I += 1
        t0 = t + 1

    end
    sigma2, w_t
end

## Regular network model
function change_edge_to_node_REG(i::Int64, N::Int64, neighbor::Int64)
    i_don = Int(floor(i / (2 * neighbor) - 1e-9)) + 1
    if i_don > N
        i_don = N
    end

    i_row = i - (2 * neighbor) * (i_don-1)
    if i_row <= neighbor
        i_rec = i_don + i_row
    else
        i_rec = i_don + i_row - 2 * neighbor - 1
    end

    if i_rec < 1
        i_rec += N
    elseif i_rec > N
        i_rec -= N
    end

    i_don, i_rec
end

function update_ysm_series_regular!(w::Array{Float64,1}, p::Float64, ϵ::Float64, N::Int64, time::Vector{Int64}, save_time::Vector{Int64}, neighbor::Int64)
    sigma2 = zeros(Float64, length(time))
    w_t = zeros(Float64, N, length(save_time))

    # Pick edges
    E = N * neighbor * 2
    inds = zeros(Int, N)
    rnds = rand(N)

    t0 = 1
    I = 1
    II = 1
    @inbounds for t = time
        @inbounds for tt = t0:t
            # rand!(inds, N_s)
            inds .= Int.(round.(rand(N).*(E-1) .+ 1))
            rnds .= rand(N)
            # ------------- t elapsed by 1 ----------------
            @inbounds for i = 1:N
                # sample!(N_s, inds, replace=false)
                i_don, i_rec = change_edge_to_node_REG(inds[i],N,neighbor)
                update_ysm!(w, p, ϵ, i_don, i_rec, rnds[i])
            end
            # ---------------------------------------------
            if tt in save_time
                w_t[:,II] .= w
                # log_wt[:,II] .= log10.(w)
                II += 1
            end
        end

        sigma2[I] = sig2(w, 1)[1]
        # sigma2_log[I] = sig2(log10.(w), 1)[1]
        I += 1
        t0 = t + 1
    end
    sigma2, w_t
end

## 1D (periodic b.c.)
# function update_ysm_series_1d!(w::Vector{Float64}, p::Float64, ϵ::Float64, N::Int64, time::Vector{Int64}, save_time::Vector{Int64})
#     sigma2 = zeros(length(time))
#     sigma2_log = zeros(length(time))
#     w_t = zeros(N, length(save_time))
#
#     # ind_c = zeros(Int,1)
#     # inc_c = zeros(Int,1)
#     ind_c = zeros(Int,N)
#     inc_c = zeros(Int,N)
#
#     ind = collect(1:N)
#     inc = [-1,1]
#
#     t0 = 1
#     I = 1
#     II = 1
#     for t = time
#         for tt = t0:t
#             rand!(ind_c, ind)
#             rand!(inc_c, inc)
#             # ------------- t elapsed by 1 ----------------
#             for i = 1:N
#                 ind_me = ind_c[i]
#                 ind_ne = ind_c[i] + inc_c[i]
#
#                 if ind_ne < 1
#                     ind_ne = N
#                 elseif ind_ne > N
#                     ind_ne = 1
#                 end
#
#                 # idx .= rand(2)
#                 # idx_int[1] = Int(ceil(idx[1] * N))
#                 # idx_int[2] = 2 * (idx[2] > 0.5) - 1
#                 # idx_int[2] = min(max(idx_int[1] + idx_int[2], 1), N)
#
#                 update_ysm!(w, p, ϵ, ind_me, ind_ne)
#             end
#             # ---------------------------------------------
#             if tt in save_time
#                 w_t[:,II] .= w
#                 II += 1
#             end
#         end
#         sigma2[I] = sig2(w, 1)[1]
#         sigma2_log[I] = sig2(log10.(w), 1)[1]
#         I += 1
#         t0 = t + 1
#     end
#     sigma2, sigma2_log, w_t
# end

######################################################################################################################
## Network construction
function configuration_model(N, γ, k₀)
    d = Pareto(γ-1,k₀)
    k = Int.(round.(rand(d, N)))
    a = [k[k .>= N]]
    while !isempty(a[1])
        b = Int.(round.(rand(d, length(a[1]))))
        k[k .>= N] .= b
        a[1] = k[k .>= N]
    end
    if sum(k) % 2 != 0
        k[1] += 1
    end

    k_lookup = [[i, k[i]] for i in 1:N]
    n_edges = Int(sum(k) / 2)
    # n_edges = sum(k)
    x = [0, 0]
    edge_list = []
    for i = 1:n_edges
        if length(k_lookup) > 1
            # A = unravel_seq(k_lookup)
            w = Weights([kl[2] for kl in k_lookup])
            # sample!(1:length(k_lookup), x, replace=false)
            sample!(1:length(k_lookup), w, x, replace=false)

            node1 = k_lookup[x[1]][1]
            node2 = k_lookup[x[2]][1]

            push!(edge_list, [node1,node2])
            push!(edge_list, [node2,node1])

            k1 = k_lookup[x[1]][2]
            k2 = k_lookup[x[2]][2]

            delete_list = []
            if k1 <= 1
                # deleteat!(k_lookup, x[1])
                push!(delete_list, x[1])
            else
                k_lookup[x[1]][2] -= 1
            end

            if k2 <= 1
                push!(delete_list, x[2])
            else
                k_lookup[x[2]][2] -= 1
            end

            sort!(delete_list)
            deleteat!(k_lookup, delete_list)
        end
    end
    unique(edge_list)
end

function unravel_seq(seq)
    """
    seq : [[1,n1], [2,n2], [3,n3],...]

    """
    A = []
    B = []
    for (i,(j,n)) in enumerate(seq)
        A = [A; fill(j,n)]
        B = [B; fill(i,n)]
    end
    [A B]
end

function star_model(N)
    edge_list = []
    for n = 2:N
        push!(edge_list, [1, n])
        push!(edge_list, [n, 1])
    end
    edge_list
end

function get_adjlist(edge_list)
    adj_list = [[] for i in 1:N]
    for edge in edge_list
        node1 = edge[1]
        node2 = edge[2]
        push!(adj_list[node1], node2)
        push!(adj_list[node2], node1)
        adj_list[node1] = unique(adj_list[node1])
        adj_list[node2] = unique(adj_list[node2])
    end
    adj_list
end

get_degree(adj_list) = [length(a) for a in adj_list]

function get_adjmatrix(adj_list)
    n = length(adj_list)
    adj_mat = zeros(Int, n, n)
    for i = 1:n
        a = adj_list[i]
        adj_mat[i,a] .= 1
    end
    adj_mat
end

zeta(N,α) = sum([1/n^α for n = 1:N])

function static_model(N::Int64, α::Float64, L::Int64)
    # edge_list = []
    edge_list = Vector{Vector{Int64}}(undef, 0)
    A = zeta(N,α)
    w = Weights([1/(n^α)/A for n=1:N])
    x = zeros(Int,2)
    count = 0
    samp = collect(1:N)
    while count < L
        sample!(samp, w, x, replace=false)
        edge = [x[1], x[2]]
        edge_rev = [x[2], x[1]]
        if (edge in edge_list) || (edge_rev in edge_list)

        else
            push!(edge_list, edge)
            push!(edge_list, edge_rev)
            count += 1
        end
    end
    edge_list
end

function neighbor_model(N)
    edge_list = Vector{Vector{Int64}}(undef, 0)
    for i = 1:N-1
        push!(edge_list, [i, i+1])
        push!(edge_list, [i+1, i])
    end
    push!(edge_list, [N, 1])
    push!(edge_list, [1, N])
    nodes = collect(1:N)
    return edge_list, nodes
end

function regular_model(N, dist)
    edge_list = Vector{Vector{Int64}}(undef, 0)
    for i = 1:N
        for d = 1:dist
            j = i + d
            if j > N
                j -= N
            end
            push!(edge_list, [i, j])
            push!(edge_list, [j, i])
        end
    end
    nodes = collect(1:N)
    return edge_list, nodes
end

function config_model_mean_deg(N, γ, k₀)
    a = (1/k₀^(γ-1) - 1/(N-1)^(γ-1)) / (γ-1)
    b = (1/k₀^(γ-2) - 1/(N-1)^(γ-2)) / (γ-2)
    return b / a
end

function neighbor_w(wt, adj_list)
    w = wt[:,end]
    k = get_degree(adj_list)
    i_rich = findall(w .> 1e-8)
    wt_neigh_sum = zeros(length(i_rich), size(wt,2))
    for (i,idx) in enumerate(i_rich)
        v = view(wt, adj_list[idx], :)
        wt_neigh_sum[i,:] .= sum(v, dims=1)[:]
    end
    wt_neigh_sum
end

function extract_giant_component(edge_list::Vector{Vector{Int64}})
    # Returns 'edge_list' of the giant component.
    g = SimpleGraph(N)
    for edge in edge_list
        add_edge!(g, edge[1], edge[2])
    end
    ind = connected_components(g)
    ind_N = [length(nodes) for nodes in ind]
    nodes = ind[ind_N .== maximum(ind_N)][1]
    new_edge_list = Vector{Vector{Int64}}(undef, 0)
    for edge in edge_list
        if (edge[1] in nodes) || (edge[2] in nodes)
            push!(new_edge_list, edge)
        end
    end

    # N_gc = length(nodes)    # Number of nodes in GC

    return new_edge_list, nodes
end

function to_matrix(edge_list::Vector{Vector{Int64}})
    return reduce(vcat,transpose.(edge_list))
end

function build_static_model(N::Int, α::Float64, L::Int)
    edge_list, nodes = extract_giant_component(static_model(N, α, L))
    edge_list = to_matrix(edge_list)
    return edge_list, nodes
end

# issimilar(num1, num2, percent) = abs((num1-num2)/num1) <= percent

# function get_save_points(p, N, time, ref)
#     # For conditions of (p, N), extract times
#     time_new = []
#     if p * N < 1.0
#         for r in ref
#             x = time[issimilar.(time ./ N, r, 0.01)]
#             if length(x) > 0
#                 push!(time_new, x[1])
#             end
#         end
#     else
#         for r in ref
#             x = time[issimilar.(time .* p, r, 0.01)]
#             if length(x) > 0
#                 push!(time_new, x[1])
#             end
#         end
#     end
#     time_new
# end

# function config_model(N, γ, k₀)
#     x = true
#     while x
#         try
#             adj = configuration_model(N, γ, k₀)
#             return adj
#         catch e
#         end
#     end
# end

# function configuration_model(N, γ, k₀)
#     d = Pareto(γ-1,k₀)
#     ks = Int.(round.(rand(d, N)))
#     a = [ks[ks .>= N]]
#     while !isempty(a[1])
#         b = Int.(round.(rand(d, length(a[1]))))
#         ks[ks .>= N] .= b
#         a[1] = ks[ks .>= N]
#     end
#     if sum(ks) % 2 != 0
#         ks[1] += 1
#     end
#
#     sort!(ks, rev=true)
#
#     # Build the graph
#     adj = [[] for i = 1:N]
#     k = zeros(Int, N)
#     for i = 1:N
#         a = union(adj[i], [i], findall(k .== ks))
#         n = ks[i] - k[i]
#         # A = setdiff(collect(1:N), a)
#         if n > 0
#             idx = sample(setdiff(collect(1:N), a), n, replace=false)
#             for ii = 1:n
#                 j = idx[ii]
#                 if k[j] < ks[j]
#                     push!(adj[i], j)
#                     push!(adj[j], i)
#                     k[i] += 1
#                     k[j] += 1
#                 end
#             end
#         end
#     end
#     adj
# end

function get_file_num()
    s = split(PROGRAM_FILE, "\\")[end]
    s = split(s, ".")[1]
    n = parse(Int64, s[end])
    # print(n)
    return n
end

function file_num(str)
    s = split(str, "\\")[end]
    s = split(s, ".")[1]
    n = parse(Int64, s[end])
    # print(n)
    return n
end


######################################################################################
### OLD ###
"""
using StatsBase

function update_ysd!(w, ϵ, p)
    idx = StatsBase.samplepair(N)

    w1 = w[idx[1]]
    w2 = w[idx[2]]

    if w1 >= w2
        sorted_idx1 = idx[1]
        sorted_idx2 = idx[2]
    else
        sorted_idx1 = idx[2]
        sorted_idx2 = idx[1]
    end

    idx_rich = sorted_idx1
    idx_poor = sorted_idx2

    w_rich = w[idx_rich]
    w_poor = w[idx_poor]

    # Choose YSM/random
    is_random = false   # Flag
    if rand() < p
        is_random = true
    end

    # Choose giver/receiver
    if rand() > 0.5
        w_giver = w_rich
        a = 1   # Label for rich -> poor
    else
        w_giver = w_poor
        a = -1   # Label for poor -> rich
    end

    if is_random
        Δw = ϵ * w_giver
    else
        Δw = ϵ * w_poor
    end

    w_rich -= a * Δw
    w_poor += a * Δw

    # Update
    w[idx_rich] = w_rich
    w[idx_poor] = w_poor
end

function update_ysd_s!(w, ϵ, p, t, ts)
    # t : repetitions between saving
    # T : How many of those t
    w_t = zeros(N, ts)
    w_t[:, 1] .= w
    for tt = 2:ts
        for t = 1:t
            update_ysd!(w, ϵ, p)
        end
        w_t[:, tt] .= w
    end
    w_t
end

"""
