# Worksheet 3 - Parallelization

Build instructions:

```shell
mkdir build && cd build
cmake ..
```

Run different cases:

- **Lid-driven cavity**
  ```shell
  ../mpi_liddriven.sh
  ```
- **Fluid Trap**
  ```shell
  ../mpi_fluidtrap.sh
  ```
  **Plane Shear Flow**
  ```shell
  ../mpi_shearflow.sh
  ```
- **Karman Vortex Street**
  ```shell
  ../mpi_obstacle.sh
  ```
- **Flow over a step**
  ```shell
  ../mpi_BFS.sh
  ```
- **Natural Convection**
  - Case (a) - high $\nu$
    ```shell
    ../mpi_convectionA.sh
    ```
  - Case (b) - low $\nu$
    ```shell
    ../mpi_convectionB.sh
    ```
- **Rayleigh Benard Convection**
  ```shell
  ../mpi_rayleighBenard.sh
  ```
  If you encounter permission issues when running the `.sh` files, run the following command:

```
chmode +x [path to file]
```

then run again the file as above.

The number of processors involved and the consequent domain decomposition can be modified in the `.sh` according to the following format:

```
mpirun --oversubscribe -np [number of all processors] [path to file] [iproc] [jproc]
```

where `iproc` and `jproc` stand for the number of processors along `x` and `y` direction respectively.

`--oversubscribe` allows the oversubscription of processor cores. This means you can start more MPI processes than the number of physical CPU cores available on the machine.

## Results Validation

### Lid-driven cavity

### Fluidtrap

Simulation (1,1) Runtime: 19s

Simulation (2,3) Runtime: 4s

Speedup S(p) = 4.75

Efficiency E(p) = 0.8

Simulation (3,2) Runtime: 5s

Speedup S(p) = 3.8

Efficiency E(p) = 0.63

## Performance Analysis

### Strong Scaling with Rayleigh Benard case

maximum speedup, parallel efficiency

### Weak Scaling with Lid-driven cavity case
