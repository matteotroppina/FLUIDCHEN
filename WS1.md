# Worksheet 1 - Lid Driven Cavity Flow

Build instructions.
```shell
mkdir build && cd build
cmake ..
make
```

Run the code with the Lid Driven Cavity case.
```shell
./fluidchen ../example_cases/LidDrivenCavity/LidDrivenCavity.dat
```

## Visualizations
The following images created with ParaView depict the results of the simulation with the parameters in example_cases. The velocity field shows that the velocity of the fluid decreases with increasing distance from the moving wall on the top. The fluid directly under the moving wall has the same velocity as the moving lid. The streamlines follow the expected paths, creating a circulation in the cavity and an accumulation point at the left bottom corner. The right upper corner has the highest pressure, in contrast to the light upper corner with the lowest pressure. 

### Velocity Field
![Velocity Field](imgs/velocity_field.png)
### Direction Field
![Direction Field](imgs/direction_field.png)
### Streamlines
![Streamlines](imgs/streamlines.png)
### Pressure Field
![Pressure Field](imgs/pressure_field.png)

## SOR's behavior depending on relaxation parameter omega
The following image shows the behavior of SOR algorithm depending on various values for the relaxation parameter. It shows that with small omega such as 0.2 the convergence is very slow, more timesteps are necessary in order to decrease the iteration numbers. The value 1.7 used in the simulation for the example file converges faster than the other chosen time steps. If we run the simulation with an omega of 2.0, it leads to an unstable simulation. 

![SOR vs w](imgs/sorbehavior.jpg)



## Algorithm's behavior depending on time step dt (fixed timestep)

According to the stability conditions for the timesteps dt, the fixed dt should not exceed 0.01. If the fixed timestep is larger than 0.01, it would diverge and would not lead to a stable simulation.

## Influence of gridsize on convergence (fixed timestep)
The fixed dt should fulfill all three stability conditions mentioned in equation (12). The first stability condition would be guaranteed until imax = 64, but the second/third conditions (CFL conditions) are already not fulfilled with imax = 32. With a fixed timestep of dt = 0.05 and nu = 0.001, convergence should be only observed with imax = jmax = 16. This expectation matches with the observations made with the simulation, it diverges for all cases expect with imax = jmax = 16.

## Kinematic viscosity (with adaptive timestep)

## Difficulties


