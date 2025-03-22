using Infiltrator

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



"""
Unicycle model
"""
function f(x, u, ω, Δ)
    v = x[3]+0.5*Δ*(u[1]+ω[1])
    θ = x[4]+0.5*Δ*(u[2]+ω[2])
    x + Δ * [v*cos(θ), v*sin(θ), u[1]+ω[1], u[2]+ω[2]]
end

"""
Jacobian of f with respect to x, evaluated at x,u,ω,Δ.
"""
function jac_fx(x, u, ω, Δ)
    v = x[3] + 0.5*Δ*(u[1] + ω[1])
    θ = x[4] + 0.5*Δ*(u[2] + ω[2])
    return [ 1  0  Δ*cos(θ)   -Δ*v*sin(θ);
             0  1  Δ*sin(θ)    Δ*v*cos(θ);
             0  0  1           0;
             0  0  0           1 ]
end


"""
Jacobian of f with respect to u, evaluated at x,u,ω,Δ.
"""
function jac_fu(x, u, ω, Δ)
    v = x[3] + 0.5*Δ*(u[1] + ω[1])
    θ = x[4] + 0.5*Δ*(u[2] + ω[2])
    return [ 0.5*Δ^2*cos(θ)    -0.5*Δ^2*v*sin(θ);
             0.5*Δ^2*sin(θ)     0.5*Δ^2*v*cos(θ);
             Δ                 0;
             0                 Δ ]
end

"""
Jacobian of f with respect to ω, evaluated at x,u,ω,Δ.
"""
function jac_fω(x, u, ω, Δ)
    v = x[3] + 0.5*Δ*(u[1] + ω[1])
    θ = x[4] + 0.5*Δ*(u[2] + ω[2])
    return [ 0.5*Δ^2*cos(θ)    -0.5*Δ^2*v*sin(θ);
             0.5*Δ^2*sin(θ)     0.5*Δ^2*v*cos(θ);
             Δ                 0;
             0                 Δ ]
end

"""
Non-standard measurement model. Can we extract state estimate from these measurements?
"""
function h(x)
    [atan(x[2], x[1]), 
     -cos(x[4])*x[3]*(x[2]-3*x[1])]
end

"""
Jacobian of h with respect to x, evaluated at x.
"""
function jac_hx(x)
    x1, x2, v, θ = x
    r2 = x1^2 + x2^2
    # For h₁ = atan(x₂,x₁):
    dh1_dx1 = -x2/r2
    dh1_dx2 =  x1/r2
    # h₁ does not depend on v or θ.
    # For h₂ = -cos(θ)*v*(x₂-3*x₁):
    dh2_dx1 = 3*v*cos(θ)         # derivative of (x₂-3*x₁) wrt x1 gives -3, multiplied by -cosθ*v gives 3*v*cosθ.
    dh2_dx2 = -v*cos(θ)          # derivative wrt x2.
    dh2_dx3 = -cos(θ)*(x2-3*x1)
    dh2_dx4 = sin(θ)*v*(x2-3*x1)
    return [ dh1_dx1  dh1_dx2  0     0;
             dh2_dx1  dh2_dx2  dh2_dx3  dh2_dx4 ]
end


"""
Extended kalman filter implementation.

Assume that the 'true' physical update in the world is given by 

xₖ = f(xₖ₋₁, uₖ, ωₖ, Δ), where Δ is the time difference between times k and k-1.

Here, uₖ is the 'true' controls applied to the system. These controls can be assumed to be a random variable,
with probability distribution given by 𝒩 (mₖ, proc_cov) where mₖ is some IMU-like measurement, and proc_cov is a constant covariance matrix.

ωₖ is assumed to be some random disturbance which affects the system. This could be something like wind. This variable is also presumed to be random,
with probability distribution given by 𝒩 (0, dist_cov).

The process model distribution is then approximated as:

P(xₖ | xₖ₋₁, uₖ) ≈ 𝒩 ( Axₖ₋₁ + Buₖ + L*0 + c, Σ̂ )

where 
A = ∇ₓf(μₖ₋₁, mₖ, 0, Δ),
B = ∇ᵤf(μₖ₋₁, mₖ, 0, Δ),
L = ∇ω f(μₖ₋₁, mₖ, 0, Δ),
c = f(μₖ₋₁, mₖ, 0, Δ) - Aμₖ₋₁ - Bmₖ - L*0

μ̂ = Aμₖ₋₁ + Bmₖ + L*0 + c
  = f(μₖ₋₁, mₖ, 0, Δ)
Σ̂ = A Σₖ₋₁ A' + B proc_cov B' + L dist_cov L'

Further, assume that the 'true' measurement generation in the world is given by

zₖ = h(xₖ) + wₖ,

where wₖ is some additive gaussian noise with probability density function given by

𝒩 (0, meas_var).

The measurement model is then approximated as 

P(zₖ | xₖ) ≈ 𝒩 ( C xₖ + d , meas_var )

where 
C = ∇ₓ h(μ̂), 
d = h(μ̂) - Cμ̂

The extended Kalman filter update equations can be implemented as the following:

Σₖ = (Σ̂⁻¹ + C' (meas_var)⁻¹ C)⁻¹
μₖ = Σₖ ( Σ̂⁻¹ μ̂ + C' (meas_var)⁻¹ (zₖ - d) )

"""
function filter(; μ=zeros(4), Σ=Diagonal([5,5,3,1.0]), x0=zeros(4), num_steps=25, meas_freq=0.5, meas_jitter=0.025, meas_var=Diagonal([0.25,0.25]), proc_cov = Diagonal([0.2, 0.1]), dist_cov=Diagonal([0.3,0.3]), rng=MersenneTwister(5), output=true)
    gt_states = [x0,] # ground truth states that we will try to estimate
    timesteps = []
    u_constant = randn(rng) * [5.0, 0.2]
    μs = [μ,]
    Σs = Matrix{Float64}[Σ,]
    zs = Vector{Float64}[]

    u_prev = zeros(2)
    x_prev = x0

    for k = 1:num_steps
        uₖ = u_constant
        mₖ = uₖ + sqrt(proc_cov) * randn(rng, 2) # Noisy IMU measurement.
        Δ = meas_freq + meas_jitter * (2*rand(rng) - 1)
        ω_true = sqrt(dist_cov) * randn(rng, 2)
        xₖ = f(x_prev, uₖ, ω_true, Δ)
        x_prev = xₖ
        u_prev = uₖ
        zₖ = h(xₖ) + sqrt(meas_var) * randn(rng, 2)

        # TODO : perform update on Σ, μ
        # μ = ...
        # Σ = ...

        A = jac_fx(μ, mₖ, zeros(2), Δ)
        B = jac_fu(μ, mₖ, zeros(2), Δ)
        L = jac_fω(μ, mₖ, zeros(2), Δ)
        μ_hat = f(μ, mₖ, zeros(2), Δ)
        Σ_hat = A * Σ * A' + B * proc_cov * B' + L * dist_cov * L'

        C = jac_hx(μ_hat)
        d = h(μ_hat) - C*μ_hat
        Σ = inv(inv(Σ_hat) + C' * inv(meas_var) * C)
        μ = Σ * (inv(Σ_hat) * μ_hat + C' * inv(meas_var) * (zₖ - d))
        
        push!(μs, μ)
        push!(Σs, Σ)
        push!(zs, zₖ)
        push!(gt_states, xₖ)
        push!(timesteps, Δ)
        if output
            println("Ttimestep ", k, ":")
            println("   Ground truth (x,y): ", xₖ[1:2])
            println("   Estimated (x,y): ", μ[1:2])
            println("   Ground truth v: ", xₖ[3])
            println("   estimated v: ", μ[3])
            println("   Ground truth θ: ", xₖ[4])
            println("   estimated θ: ", μ[4])
            println("   measurement received: ", zₖ)
            println("   Uncertainty measure (det(cov)): ", det(Σ))
        end
    end

    (; μs, Σs)
end