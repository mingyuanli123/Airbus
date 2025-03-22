function iterative_closest_point(map_points, pointcloud, R, t; max_iters = 10, visualize=false)
    if visualize
        f = Figure()
        ax = Axis(f[1,1], aspect=DataAspect())

        lims = ( - 5, 45)
        xlims = lims
        ylims = lims
        xlims!(ax, xlims)
        ylims!(ax, ylims)
        Mpx = [map_points[i][1] for i = 1:length(map_points)]
        Mpy = [map_points[i][2] for i = 1:length(map_points)]
        
        scatter!(ax, Mpx, Mpy, color=:blue)
        scene = display(f)

        Qt = [R*pointcloud[i]+t for i =1:length(pointcloud)]
        Qtx = [Qt[i][1] for i = 1:length(Qt)]
        Qty = [Qt[i][2] for i = 1:length(Qt)]
        sc = scatter!(ax, Qtx, Qty, color=:red)
    end

    @infiltrate 
    @assert length(pointcloud) ≤ length(map_points) 
    N = length(pointcloud)

    point_associations = Vector{Int}(1:N)

    not_converged = true
    iters = 0
    err = 0
    while not_converged
        num_changes = update_point_associations!(point_associations, pointcloud, map_points, R, t)
        err = update_point_transform!(point_associations, pointcloud, map_points, R, t)
        iters += 1
        not_converged = (num_changes > 0) && iters < max_iters
        if visualize
            Qt = [R*pointcloud[i]+t for i =1:length(pointcloud)]
            Qtx = [Qt[i][1] for i = 1:length(Qt)]
            Qty = [Qt[i][2] for i = 1:length(Qt)]
            delete!(ax, sc)
            sc = scatter!(ax, Qtx, Qty, color=:red)
            @infiltrate
        end
    end
    return (; R, t, err)
end

"""
point_associations[i] = j means that pointcloud[i] is paired with map_points[j]

Here point_associations is updated in-place.

Here we WILL allow multiple points to be associated to the same map_point. Some implementations will
not allow this, but it is fine for us.

This function returns num_changes, which is how many elements of point_associations are changed.

"""
function update_point_associations!(point_associations, pointcloud, map_points, R, t)
    num_changes = 0
    for i = 1:length(pointcloud)
        min = Inf
        index = 1
        for j = 1:length(map_points)
            dis = (map_points[j][1] - (R*pointcloud[i]+t)[1])^2 + (map_points[j][2] - (R*pointcloud[i]+t)[2])^2
            if dis <= min
                min = dis
                index = j
            end
        end
        if point_associations[i] != index
            point_associations[i] = index    
            num_changes = num_changes + 1   
        end
    end
    return num_changes
end

"""
This function updates R and t in place, to minimize 

∑ ||pᵢ - Rqᵢ - t||²2

where pᵢ is the point in map_points indicated by point_associations[i]
and qᵢ is the i-th point in pointcloud.

You will need to derive the closed-form solution to this problem. You can do this from scratch 
(I recommend you try before looking up resources), but the full approach is mostly given in 
the following reference. Note that some of the slides have some typos, but the hand-written
proof is correct. You will need to keep track of the fact that some notational differences
exist between what we call variables, etc. and what that author calls things.
See https://cs.gmu.edu/~kosecka/cs685/cs685-icp.pdf, slides 1-8.

Return err = 1/N ∑ ||pᵢ - R*qᵢ-t ||^2
"""
function update_point_transform!(point_associations, pointcloud, map_points, R, t)
    mu_p = sum(map_points[point_associations[i]] for i in 1:length(point_associations)) / length(point_associations)
    mu_q = sum(pointcloud) / length(pointcloud)

    centered_p = [map_points[point_associations[i]] - mu_p for i in 1:length(point_associations)]
    centered_q = [pointcloud[i] - mu_q for i in 1:length(pointcloud)]

    W = reduce(+, [(centered_q[i] * transpose(centered_p[i])) for i in 1:length(pointcloud)])

    F = svd(W)
    R_new = F.V * transpose(F.U)

    if abs(det(R_new) + 1) < 0.01
        R_new = F.V * [1 0; 0 -1] * transpose(F.U)
    end

    R .= R_new
    t .= mu_p - R * mu_q

    err = mean([norm(map_points[point_associations[i]] - (R * pointcloud[i] + t))^2 for i in 1:length(pointcloud)])

    return err
end


function test_ICP(; Nb = 70, visualization=false)
    map_points = readdlm("hd_map_unlabeled.csv", ',', Float64)
    map_points = eachrow(map_points) |> collect # transform into list of 2d points
    # assume pointcloud is also a list of 2d points
    
    ground_truth_map = example_map()
    dθ = -π .+ 2*π*(1:Nb)./Nb
    for i = 1:10
        x = feasible_point(ground_truth_map)
        θ = rand()*2*π-π
        Rinv = [cos(θ) sin(θ); -sin(θ) cos(θ)]
        pointcloud = []
        for d in dθ
            pt, α = compute_range(x, θ+d, ground_truth_map.segments; α_max=7.0)
            if isinf(α)
                continue
            else
                push!(pointcloud, Rinv*(pt .- x)) # transform to EGO frame
            end
        end
        θ̂ = θ + randn()*0.2
        x̂ = x + randn(2)*1.5

        R0 = [cos(θ̂) -sin(θ̂); sin(θ̂) cos(θ̂)]
        t0 = Vector{Float64}(x̂)
        (; R, t, err) = iterative_closest_point(map_points, pointcloud, R0, t0; visualize=visualization)
        x_estimated = t
        θ_estimated = atan(R[2,1], R[1,1])
        println("original error: $(norm(x-x̂)), $(θ-θ̂)")
        println("output error: $(norm(x-x_estimated)), $(θ-θ_estimated)")
        #println("x true: $x, x est: $x_estimated")
        #println("θ true: $θ, θ est: $θ_estimated")
    end
end
