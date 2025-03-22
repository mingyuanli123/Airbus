struct DetectedObject
    id::Int
    position::SVector{2, Float64}
    velocity::SVector{2, Float64}
end

struct MyPerceptionType
    objects::Vector{DetectedObject}
end

function perception(cam_meas_channel, localization_state_channel, perception_state_channel)
    # set up stuff
    tracked_objects = Dict{Int, DetectedObject}()

    while true
        fresh_cam_meas = []
        while isready(cam_meas_channel)
            meas = take!(cam_meas_channel)
            push!(fresh_cam_meas, meas)
        end

        latest_localization_state = fetch(localization_state_channel)
        
        # process bounding boxes / run ekf / do what you think is good
        detections = DetectedObject[]

        perception_state = MyPerceptionType(0,0.0)
        if isready(perception_state_channel)
            take!(perception_state_channel)
        end
        put!(perception_state_channel, perception_state)
    end
end

