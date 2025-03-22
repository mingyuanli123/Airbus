using JuMP
using Ipopt

"""
Physical properties of a vehicle.
"""
struct VehicleProperties
    half_length::Float64
    half_width::Float64
end

"""
Information regarding trajectory of a vehicle.
"""
struct Trajectory
    vehicle_properties::VehicleProperties
    x₀::Vector{Float64}
    X::Vector{Vector{Float64}}
end

"""
Defines halfspace { x: a'x ≥ b }
"""
struct HalfSpace
    a::Vector{Float64} # normal vector
    b::Float64         # offset
end
"""
Assumes x = [p₁, p₂, v, θ]
        u = [a, ω]

Return ẋ = [cosθ⋅v, sinθ⋅v, a, ω] 
"""
function f_simple_car(x, u)
    [cos(x[4])*x[3], sin(x[4])*x[3], u[1], u[2]]
end


"""
U = [u₁, ..., uₖ]
uᵢ ∈ ℝ², i ∈ {1,...,k}
x₀ ∈ ℝ⁴

Return [x₁, ..., xₖ]
where xᵢ₊₁ = xᵢ + Δ⋅f_simple_car(xᵢ, uᵢ₊₁)
"""
function generate_trajectory(x₀, U, vehicle_properties; Δ=0.25)
    X = Vector{Vector}()
    #TODO fill in X
    cur_x = x₀
    for u in U
        next_x = cur_x + Δ * f_simple_car(cur_x, u)
        push!(X, next_x)
        cur_x = next_x
    end
    
    Trajectory(vehicle_properties, x₀, X)
end


function generate_collision_free_controls()
    U = [[[0.9,0.1] for _ in 1:15]; [[0.1, -0.25] for _ in 16:25]; [[0.0,0] for _ in 26:35]]
end


"""
return true if car with state x₁ and body
properties b₁ is in collision with car x₂ with
body b₂, otherwise false.

Collision means that the two rectangular bodies are intersecting.

Assume x is of form [p₁, p₂, v, θ] and
       b are vehicle properties as defined above.
"""

function local_to_global(state::Vector{Float64}, local_point::Vector{Float64})::Vector{Float64}
    """
    Transform a point from the vehicle's local frame to the global frame.

    # Arguments
    - `state::Vector{Float64}`: [x_pos, y_pos, velocity, theta]
    - `local_point::Vector{Float64}`: [a, b] in local coordinates

    # Returns
    - `global_point::Vector{Float64}`: [x_global, y_global]
    """
    x, y, _, theta = state
    a, b = local_point

    # Rotation matrix
    R = [cos(theta) -sin(theta);
         sin(theta)  cos(theta)]

    # Global coordinates
    global_point = R * [a, b] .+ [x, y]
    return global_point
end


function check_collision(x₁::Vector{Float64},
                         b₁::VehicleProperties,
                         x₂::Vector{Float64},
                         b₂::VehicleProperties)

    ϵ = 0.0001 
    model = Model(Ipopt.Optimizer)
    set_silent(model)
    
    # Extract pos and theta for x1, x2
    x1_pos, y1_pos, _, theta1 = x₁
    x2_pos, y2_pos, _, theta2 = x₂
    
    # Define variables for points within Vehicle 1 in local coordinates
    @variable(model, -b₁.half_length <= a1 <= b₁.half_length)
    @variable(model, -b₁.half_width <= b1_local <= b₁.half_width)
    
    @variable(model, -b₂.half_length <= a2 <= b₂.half_length)
    @variable(model, -b₂.half_width <= b2_local <= b₂.half_width)
    
    cos_theta1 = cos(theta1)
    sin_theta1 = sin(theta1)
    cos_theta2 = cos(theta2)
    sin_theta2 = sin(theta2)
    
    p1_x = x1_pos + cos_theta1 * a1 - sin_theta1 * b1_local
    p1_y = y1_pos + sin_theta1 * a1 + cos_theta1 * b1_local
    p2_x = x2_pos + cos_theta2 * a2 - sin_theta2 * b2_local
    p2_y = y2_pos + sin_theta2 * a2 + cos_theta2 * b2_local
    
    # Define the squared Euclidean distance as the objective
    @objective(model, Min, (p1_x - p2_x)^2 + (p1_y - p2_y)^2)
    
    # Optimize the model
    optimize!(model)
    
    # Check if optimization was successful
    status = termination_status(model)
    if status != MOI.OPTIMAL && status != MOI.LOCALLY_SOLVED
        @warn "Optimization did not converge to an optimal solution. Status: $status"
        return false  # Assume no collision if optimization fails
    end
    
    min_dist_sq = objective_value(model)
    return min_dist_sq <= ϵ
end

"""
helper functions
"""
function check_pt_in_HalfSpace(pt, hs::HalfSpace)
    return dot(hs.a, pt) >= hs.b
end

"""
Return false if vehicle defined by state x
and vehicle properties b lies entirely within
each of the halfspaces listed in L, true otherwise.
"""
function check_lane_violations(x::Vector{Float64}, b::VehicleProperties, L::Vector{HalfSpace})
    # Extract position and orientation
    x_pos, y_pos, _, theta = x

    # Define the four corners of the vehicle in the local frame
    local_corners = [
        [-b.half_length,  b.half_width],  
        [-b.half_length, -b.half_width],  
        [ b.half_length, -b.half_width],  
        [ b.half_length,  b.half_width]  
    ]

    # Rotation matrix
    R = [cos(theta) -sin(theta);
    sin(theta)  cos(theta)]

    # Rotate and translate the corners to global coordinates
    global_corners = [R * corner .+ [x_pos, y_pos] for corner in local_corners]

    for corner in global_corners
        for hs in L
            if !check_pt_in_HalfSpace(corner, hs)
                return true
            end
        end
    end

    # If all corners are within all halfspaces, return false (no violation)
    return false
end


"""
True if trajectories are collision free and satisfy lane constraints, false otherwise.
"""
function check_trajectories(trajectories::Vector{Trajectory}, lanes::Vector{HalfSpace}) 
    states = [t.X for t in trajectories]
    for xs in zip(states...)
        for (e1,x1) in enumerate(xs)
            if check_lane_violations(x1, trajectories[e1].vehicle_properties, lanes)
                return false
            end
            for (e2,x2) in enumerate(xs)
                if e2 ≤ e1
                    continue
                else
                    if check_collision(x1, trajectories[e1].vehicle_properties,
                                       x2, trajectories[e2].vehicle_properties)
                        return false
                    end
                end
            end
        end
    end
    return true
end
