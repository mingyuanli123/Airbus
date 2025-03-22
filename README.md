# Airbus

Group member: Elsie Li, Tommy Li, Mingyuan Li


## Technical Plan
### i. Routing & Motion Planning
#### Objective
Plan a route from the vehicle's current segment to a specified target segment within a road network. Generate smooth and safe trajectories to follow the route, accounting for road constraints and dynamic objects.

#### Implementation Strategy:
1. Use polynomial trajectory generation or spline interpolation for path smoothing.

2. Consider the bicycle model for an ego vehicle, i.e. the vehicle has 4 states (position (2D), angular displacement, velocity).

3. Integrate velocity profiling for acceleration/deceleration, considering stop signs and speed limits.

4. Formulate motion planning into Model Predictive Control (MPC) optimization problem, where constraints can be thought of as stay in lane, stop at stop signs, and avoid other vehicles.

#### Key Ideas and Tools
Polyline lane segments, state dynamic model, MPC controller, collision avoidance checker.

### ii. Localization
#### Objective
Determine the vehicle’s position and orientation using GPS, IMU, and wheel encoder data.

#### Implementation Strategy:
1. Use an Extended Kalman Filter (EKF) to fuse GPS and IMU measurements.

2. Incorporate a motion model for prediction and sensor updates for correction.

3. Use map matching to improve localization accuracy in areas with poor GPS reception.

#### Key Tools
julia-based EKF implementation, sensor simulation modules.

### iii. Perception
#### Objective
Use data of ego vehicle and camera measurements to develop a perception module that estimates states of other vehicles.

#### Implementation Strategy:
1. Study the simulation code provided to us and understand what types of data are involved and how these can be used.

2. Leverage camera data to keep track of nearby vehicles, and use camera projection knowledge to estimate vehicle positions and locality differences to calculate velocity vectors (derivatives of location).
#### Key Tools
Camera Projections, Bounding Box of Vehicle, Gradients

## Testing Plan
### Unit Testing of Components
#### Routing
Simulate pathfinding in various synthetic and real map scenarios to validate shortest path accuracy and re-routing robustness.
#### Motion Planning
Use straight-line and curved tracks in VehicleSim to test basic and dynamic behaviors (e.g., stop, pull out).
#### Localization
Replay recorded sensor data and compare estimated vs. ground-truth trajectories; test under GPS drop-out conditions.
#### Perception
Use labeled datasets to evaluate detection accuracy and speed; simulate  multiple vehicle interactions and test object detection + tracking robustness.

### Integration Testing
Full stack will be deployed in a simulated environment.

End-to-end tests will include a fixed start and goal location with: Route generation, Real-time perception updates, State estimation via localization, Adaptive motion planning to reach the goal

## Timeline(Due April 19th)
Routing & Motion Control Part (Due April 2nd)

Localization & Perception Part (Due April 5th)

Modular Tests(Due April 5th)

Integrate Tests (Due April 9th)

## Division of Labor
Part 1: Routing & Motion Planning
Name: Tommy, Elsie

Part 2: Localization & Perception
Name: Mingyuan, Elsie

Part 3: Testing
Name: Mingyuan, Tommy
