# Worksheet 2 - Arbitrary Geometries and Energy Transport for the Navier-Stokes Equations

Build instructions.

```
mkdir build && cd build
cmake ..
make
```

## Plane Shear Flow

Run the code with the "Plane Shear Flow" case:

```
./fluidchen ../example_cases/ShearFlow/ShearFlow.dat
```

#### Velocity Field

![Shearflow Velocity](imgs/shearflow_velocity.png)

#### Pressure Field

![Velocity Field](imgs/shearflow_pressure.png)

## Karman Vortex Street

Run the code with the "Karman Vortex Street" case:

```
./fluidchen ../example_cases/ChannelWithObstacle/ChannelWithObstacle.dat
```

#### Velocity Field

![Shearflow Velocity](imgs/karman_velocity.png)

<!-- ### Streamlines
![Velocity Field](imgs/karman_streamlines.png) -->

#### Pressure Field

![Velocity Field](imgs/karman_pressure.png)

## Flow over a step

Run the code with the "Flow over a step" case:

```
./fluidchen ../example_cases/ChannelWithBFS/ChannelWithBFS.dat
```

#### Velocity Field

![Shearflow Velocity](imgs/BFS_velocity.png)

#### Pressure Field

![Velocity Field](imgs/BFS_pressure.png)

## Natural Convection

### Case (a) - high $\nu$

Run the code with the "Natural Convection" case with high Kinematic Viscosisty $\nu = 0.001$ ;

```
./fluidchen ../example_cases/NaturalConvection/NaturalConvection_a.dat
```

#### Velocity Field

![Shearflow Velocity](imgs/NC_a_velocity.png)

#### Temperature Field

![Velocity Field](imgs/NC_a_temp.png)

### Case (b) - low $\nu$

Run the code with the "Natural Convection" case with low Kinematic Viscosisty $\nu = 0.0002$;

```
./fluidchen ../example_cases/NaturalConvection/NaturalConvection_b.dat
```

#### Velocity Field

![Shearflow Velocity](imgs/NC_b_velocity.png)

#### Pressure Field

![Velocity Field](imgs/NC_b_temp.png)

## Fluid Trap

Run the code with the "Fluid Trap" case:

```
./fluidchen ../example_cases/FluidTrap/FluidTrap.dat
```

#### Velocity Field

![Shearflow Velocity](imgs/fluidtrap_temp.png)

## Rayleigh Benard Convection

Run the code with the "Rayleigh Benard Convection" case:

```
./fluidchen ../example_cases/RayleighBenard/RayleighBenard.dat
```
