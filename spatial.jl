using Plots
using StatsBase
using LaTeXStrings
using DelimitedFiles

function bfs(node, dict_nodes, dict_edges, counted, retained_link, clusterID, ID_list)

    # Breadth first search. Given the network structure, return the size of connected component that contains `node`
    # ID_list[i]: node i is in the component with id `clusterID`

    # counted: 0 not counted, 1 counted
    # retained: 0 not retained, 1 retained

    N = length(dict_nodes)

    queue_edges = Int64[]
    queue_nodes = Int64[]

    push!(queue_nodes, node)

    counted[node] = 1
    ID_list[node] = clusterID
    cluster_size = 1

    while !isempty(queue_nodes)
        node = popfirst!(queue_nodes)

        for edge in dict_nodes[node]
            if (counted[edge] == 0) & (retained_link[edge-N] == 1)
                push!(queue_edges, edge)
                counted[edge] = 1
                ID_list[edge] = clusterID
            end
        end

        while !isempty(queue_edges)
            edge = popfirst!(queue_edges)
            for node in dict_edges[edge]
                if (counted[node] == 0)
                    cluster_size += 1
                    push!(queue_nodes, node)
                    counted[node] = 1
                    ID_list[node] = clusterID
                end
            end
        end
    end
    return cluster_size, counted, ID_list, clusterID
end

function gc(dict_nodes, dict_edges, retained_link)

    # calculate the size of the maximum connected cluster

    num_nodes = length(dict_nodes)
    num_edges = length(dict_edges)
    counted = zeros(Int64, num_nodes + num_edges)
    ID_list = zeros(Int64, num_nodes + num_edges)
    max_cluster_size = 0
    max_clusterID = 0
    clusterID = 0
    for i = 1:num_nodes
        if (counted[i] == 0)
            clusterID += 1
            cluster_size, counted, ID_list, clusterID = bfs(i, dict_nodes, dict_edges, counted, retained_link, clusterID, ID_list)
            if cluster_size > max_cluster_size
                max_cluster_size = cluster_size
                max_clusterID = clusterID
            end
        else
            continue
        end
    end
    return max_cluster_size, max_clusterID, ID_list
end

function spatial_spatial_regulation(N, node_density, d0, drp, drn, c, cp, cn)
    # node_density: the number of nodes in the unit cube.
    # c, cp, cn are coefficient of the probability distribution (may not be normalized)

    # dict_nodes_stru: {nodeID -> [neighbor edgeID]}
    # dict_edges_stru: {edgeID -> [structural neighbor nodeID]} In the case of simple network, each edge always has two structural neighbors.
    # dict_edges_pos: {edgeID -> [positive regulatory neighbor nodeID]}
    # dict_edges_neg: {edgeID -> [negative regulatory neighbor nodeID]}

    # Here we consider the general case where the typical scale for positive (drp) and negative (dnp) regulations can be different


    dict_nodes_stru = Dict{Int64,Vector{Int64}}()
    dict_edges_stru = Dict{Int64,Vector{Int64}}()
    dict_edges_pos = Dict{Int64,Vector{Int64}}()
    dict_edges_neg = Dict{Int64,Vector{Int64}}()

    L_cube = sqrt(N / node_density)
    x = rand(N) * L_cube        # Coordinate of nodes
    y = rand(N) * L_cube

    # Estimate number of edges, as space need to be allocated for dict_edges_stru.

    counter = 0
    for i = 1:N
        for j = i+1:N
            dist = sqrt((x[i] - x[j])^2 + (y[i] - y[j])^2)
            prob = c * exp(-dist / d0)
            p = rand()
            if p < prob
                counter += 1
            end
        end
    end


    upper_bound = counter * 1.5   # Set an upper bound for the number of edges          

    # Allocate space for dictionary.

    for i = 1:N
        dict_nodes_stru[i] = Int64[]
    end

    for j = N+1:N+upper_bound
        dict_edges_stru[j] = Int64[]
    end

    edgeID = N + 1


    dist_avg = []

    # Construct the structural network
    for i = 1:N
        for j = i+1:N

            # Periodic boundary condition

            dx = abs(L_cube / 2 - abs(L_cube / 2 - abs(x[i] - x[j])))
            dy = abs(L_cube / 2 - abs(L_cube / 2 - abs(y[i] - y[j])))
            dist = sqrt(dx^2 + dy^2)
            push!(dist_avg, dist)

            # Caution: do not have c>1 otherwise the probability might be greater than 1.

            prob = c * exp(-dist / d0)
            p = rand()
            if p < prob
                push!(dict_nodes_stru[i], edgeID)
                push!(dict_nodes_stru[j], edgeID)
                push!(dict_edges_stru[edgeID], i)
                push!(dict_edges_stru[edgeID], j)
                edgeID += 1
            end
        end
    end
    @show mean(dist_avg)

    # delete the over-allocated space.
    max_edgeID = edgeID - 1
    for i = max_edgeID+1:N+upper_bound
        delete!(dict_edges_stru, i)
    end

    M = max_edgeID - N  # Number of edges

    for i = N+1:N+M
        dict_edges_pos[i] = Int64[]
        dict_edges_neg[i] = Int64[]
    end

    xl = zeros(M)
    yl = zeros(M)

    # Coordinates of edges
    for i = 1:M
        edgeID = N + i
        node1 = dict_edges_stru[edgeID][1]
        node2 = dict_edges_stru[edgeID][2]
        xl[i] = (x[node1] + x[node2]) / 2
        yl[i] = (y[node1] + y[node2]) / 2
    end

    # Construct the regulatory network
    for j = 1:M
        for i = 1:N
            edgeID = N + j
 
            dxl = abs(L_cube / 2 - abs(L_cube / 2 - abs(x[i] - xl[j])))
            dyl = abs(L_cube / 2 - abs(L_cube / 2 - abs(y[i] - yl[j])))
            d = sqrt(dxl^2 + dyl^2)

            # Caution: here the probabilities are not normalized. Do not let cp+cn>1.
            # Always check whether the regulations are generated as expected.

            prob_low = cp * exp(-d / drp)
            prob_high = cp * exp(-d / drp) + cn * exp(-d / drn)
      
            p = rand()
            if p < prob_low
                push!(dict_edges_pos[edgeID], i)
            elseif (p > prob_low) & (p < prob_high)
                push!(dict_edges_neg[edgeID], i)
            end
        end
    end

    @show mean([length(dict_nodes_stru[i]) for i = 1:N])
    @show mean([length(dict_edges_pos[i]) for i = N+1:length(dict_edges_pos)])
    @show mean([length(dict_edges_neg[i]) for i = N+1:length(dict_edges_neg)])

    # here return of the x and y coordiate of nodes, for plotting later
    return dict_nodes_stru, dict_edges_stru, dict_edges_pos, dict_edges_neg, x, y
end

function iteration(dict_nodes_stru, dict_edges_stru, dict_edges_pos, dict_edges_neg, p, retained_link, Tmax)

    # retained_link as a input, in order to control the initial condition of each iteration

    N = length(dict_nodes_stru)         # Numer of nodes
    M = length(keys(dict_edges_stru))   # Number of links

    active_node = ones(Int64, N)
    cluster_size_list = zeros(Tmax)
    T_list = zeros(Tmax)
    state = zeros(Int64, Tmax, N)       # state[t, i] indicates the state of node i at time t (0 for off and 1 for on)
    for t = 1:Tmax

        # Find the giant component and damage others. 
        # Only nodes in the giant component are active.

        max_cluster_size, max_clusterID, ID_list = gc(dict_nodes_stru, dict_edges_stru, retained_link)

        for i = 1:N
            if ID_list[i] == max_clusterID && (max_cluster_size >= 1)
                active_node[i] = 1
            else
                active_node[i] = 0
            end
        end

        cluster_size_list[t] = (max_cluster_size / N)
        T_list[t] = t

        state[t, :] = active_node

        for edge in keys(dict_edges_stru)
            pos_regulation = dict_edges_pos[edge]
            neg_regulation = dict_edges_neg[edge]
            sum_pos = 0
            sum_neg = 0
            for pos_node in pos_regulation
                sum_pos += active_node[pos_node]
            end

            for neg_node in neg_regulation
                sum_neg += active_node[neg_node]
            end

            if (sum_pos > 0) & (sum_neg == 0) & (rand() < p)
                retained_link[edge-N] = 1
            else
                retained_link[edge-N] = 0
            end
        end

    end
    return cluster_size_list, T_list, retained_link, state
end


function orbit_diagram()

    # Calculate the orbit diagram with respect to p.

    N = 1000               # total number of nodes
    node_density = 100      # number of nodes in a unit square
    p = 1                   # probability of retaining a link
    d0 = 0.25               # typical length of structural links
    drp = 0.25              # typical length of positive regulatory links
    drn = 0.25              # typical length of negative regulatory links
    c = 0.6                 # control parameter of average degree of structural network (it is NOT the average degree, see the function above)
    cp = 0.2                # control parameter of average positive regulatory degree
    cn = 0.2                # control parameter of average negative regulatory degree

    Tmax = 500              # Maximum time for running the dynamics
    collect = 100           # We take the last `collect` elements in the time series. This value need to be smaller than Tmax

    p_collect = zeros(collect * 101)
    R_collect = zeros(collect * 101)

    dict_nodes_stru, dict_edges_stru, dict_edges_pos, dict_edges_neg, x, y = spatial_spatial_regulation(N, node_density, d0, drp, drn, c, cp, cn)

    Threads.@threads for p_iter = 0:1:100    # Calculate for 101 different p value. 

        # Remove the macro to run sequentially.

        p = p_iter / 100
        @show p
        retained_link = rand(length(dict_edges_stru)) .< 1           # Initial condition could affect the result significantly.
        cluster_size_list, T_list, retained_link, state = iteration(dict_nodes_stru, dict_edges_stru, dict_edges_pos, dict_edges_neg, p, retained_link, Tmax)


        p_collect[p_iter*collect+1:(1+p_iter)*collect] = ones(collect) * p
        R_collect[p_iter*collect+1:(1+p_iter)*collect] = cluster_size_list[end-collect+1:end]

    end
    return p_collect, R_collect
end

function time_trajectory()


    N = 10000               # total number of nodes
    node_density = 100      # number of nodes in a unit square
    p = 1                   # probability of retaining a link (after the regulation)
    d0 = 0.25               # typical length of structural links
    drp = 0.25              # typical length of positive regulatory links
    drn = 0.25              # typical length of negative regulatory links
    c = 0.6                 # control parameter of average degree of structural network (it is NOT the average degree, see the function above)
    cp = 0.2                # control parameter of average positive regulatory degree
    cn = 0.2                # control parameter of average negative regulatory degree

    # Produce a GIF animation of the time trajectory of the spatial pattern

    pL0 = 1                 # initial condition

    Tmax = 500              # maximum time for running the dynamics
    box_size = sqrt(N / node_density)   # box size for plotting

    dict_nodes_stru, dict_edges_stru, dict_edges_pos, dict_edges_neg, x, y = spatial_spatial_regulation(N, node_density, d0, drp, drn, c, cp, cn)
    retained_link = rand(length(dict_edges_stru)) .< pL0            # Initial condition could affect the result significantly.
    @show sum(retained_link)
    cluster_size_list, T_list, retained_link, state = iteration(dict_nodes_stru, dict_edges_stru, dict_edges_pos, dict_edges_neg, p, retained_link, Tmax)

    anim = @animate for t = Tmax-100:Tmax
        state_t = state[t, :]
        active_x = filter!(z -> z != 0, state_t .* x)       # Only keep the coordinates of active nodes
        active_y = filter!(z -> z != 0, state_t .* y)

        scatter(active_x, active_y, title="t=$t", legend=false, markersize=5, framestyle=:box,
            xlabel="\$x\$", ylabel="\$y\$", xlim=(0, box_size), ylim=(0, box_size),
            guidefont=font(15), xtickfont=font(14), ytickfont=font(14), titlefont=font(14))
    end
    gif(anim, "time_trajectory.gif", fps=5)        # Generate a GIF figure of the snapshot. 

    return cluster_size_list, T_list, x, y, state
end

cluster_size_list, T_list, x, y, state = time_trajectory()
figure1 = plot(T_list, cluster_size_list)
xlabel!(L"t")
ylabel!(L"R")
display(figure1)

# p_collect, R_collect = orbit_diagram()
# figure2 = scatter(p_collect, R_collect)
# xlabel!(L"p")
# ylabel!(L"R")
# display(figure2)
