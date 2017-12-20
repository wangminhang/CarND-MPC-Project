## The model

The kinematic model includes the vehicle's x and y coordinates, orientation angle and velocity, as well as the cross-track error(CTE) and psi error (epsi). Actuator outputs are acceleration and delta (steering angle). The model calculate the current state from the previous state.

## Timestep Length and Elapsed Duration 

After tried a lot of combinations, I choose them equal 15 and 0.05 respectively. 

## Polynomial Fitting and MPC Preprocessing

Transforming the waypoints to the vehicle's perspective and fit a polynomial to the them.

## Model Predictive Control with Latency

It is important to measure and incorporate that latency in the model. As for this project, the simulator have a latency of 100 ms, hence that latency is also incorporated before the actuators values are arrived and sent to vehicle.

---