# Worksheet 2 - Arbitrary Geometries and Energy Transport for the Navier-Stokes Equations

Build instructions.

```shell
mkdir build && cd build
cmake ..
make
```

## Plane Shear Flow

Run the code with the "Plane Shear Flow" case:

```shell
./fluidchen ../example_cases/ShearFlow/ShearFlow.dat
```

### Velocity Field

![Shearflow Velocity](imgs/shearflow_velocity.png)

#### Analytical Solution Comparison
At x = 5.0 
![Analytical Solution](imgs/shearflow_plot.png)

### Pressure Field

![Velocity Field](imgs/shearflow_pressure.png)

## Karman Vortex Street

Run the code with the "Karman Vortex Street" case:

```shell
./fluidchen ../example_cases/ChannelWithObstacle/ChannelWithObstacle.dat
```

### Velocity Field

![Shearflow Velocity](imgs/karman_velocity.png)

<!-- ### Streamlines
![Velocity Field](imgs/karman_streamlines.png) -->

### Pressure Field

![Velocity Field](imgs/karman_pressure.png)

## Flow over a step

Run the code with the "Flow over a step" case:

```shell
./fluidchen ../example_cases/ChannelWithBFS/ChannelWithBFS.dat
```

### Velocity Field

![Shearflow Velocity](imgs/BFS_velocity.png)

### Pressure Field

![Velocity Field](imgs/BFS_pressure.png)

## Natural Convection

Run the code with the "Natural Convection" case;

```shell
./fluidchen ../example_cases/NaturalConvection/NaturalConvection.dat
```

## Fluid Trap

Run the code with the "Fluid Trap" case:

```shell
./fluidchen ../example_cases/FluidTrap/FluidTrap.dat
```

### Velocity Field

![Shearflow Velocity](imgs/fluidtrap_temp.png)

## Rayleigh Benard Convection

Run the code with the "Rayleigh Benard Convection" case:

```shell
./fluidchen ../example_cases/RayleighBenard/RayleighBenard.dat
```
