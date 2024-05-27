# Worksheet 2 - Arbitrary Geometries and Energy Transport for the Navier-Stokes Equations

Build instructions.

```
mkdir build && cd build
cmake ..
make
```

Run different cases:

- **Plane Shear Flow**
  ```
  ./fluidchen ../example_cases/ShearFlow/ShearFlow.dat
  ```
- **Karman Vortex Street**
  ```
  ./fluidchen ../example_cases/ChannelWithObstacle/ChannelWithObstacle.dat
  ```
- **Flow over a step**
  ```
  ./fluidchen ../example_cases/ChannelWithBFS/ChannelWithBFS.dat
  ```
- **Natural Convection**
  - Case (a) - high $\nu$
    ```
    ./fluidchen ../example_cases/NaturalConvection/NaturalConvection_a.dat
    ```
  - Case (b) - low $\nu$
    ```
    ./fluidchen ../example_cases/NaturalConvection/NaturalConvection_b.dat
    ```
- **Fluid Trap**
  ```
  ./fluidchen ../example_cases/FluidTrap/FluidTrap.dat
  ```
- **Rayleigh Benard Convection**
  ```
  ./fluidchen ../example_cases/RayleighBenard/RayleighBenard.dat
  ```

## Plane Shear Flow

#### Velocity Field

![Shearflow Velocity](imgs/shearflow_velocity.png)

#### Analytical Solution Comparison (x = 5.0)

![Analytical Solution](imgs/shearflow_plot.png)

#### Pressure Field

![Velocity Field](imgs/shearflow_pressure.png)

## Karman Vortex Street

#### Velocity Field

![Shearflow Velocity](imgs/karman_velocity.png)

<!-- ### Streamlines
![Velocity Field](imgs/karman_streamlines.png) -->

#### Pressure Field

![Velocity Field](imgs/karman_pressure.png)

## Flow over a step

#### Velocity Field

![Shearflow Velocity](imgs/BFS_velocity.png)

#### Pressure Field

![Velocity Field](imgs/BFS_pressure.png)

## Natural Convection

### Case (a) - $\nu=0.001$

#### Velocity Field

![Shearflow Velocity](imgs/NC_a_velocity.png)

#### Temperature Field

![Velocity Field](imgs/NC_a_temp.png)

### Case (b) - $\nu = 0.0002$

#### Velocity Field

![Shearflow Velocity](imgs/NC_b_velocity.png)

#### Pressure Field

![Velocity Field](imgs/NC_b_temp.png)

## Fluid Trap

#### Velocity Field

![Shearflow Velocity](imgs/fluidtrap_temp.png)

## Rayleigh Benard Convection
