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

For sequential execution:

```
mpirun --oversubscribe -np 1 [path to file] 1 1
```

## Results Validation

### Lid-driven cavity

When running the simulation with a domain size of imax = 300 and jmax = 300, the time step size becomes too small, causing the SOR method to struggle with convergence. This results in the maximum number of iterations being exceeded at each time step.

To address this issue, we carried out the simulation with a reduced domain size of `imax` = 150 and `jmax` = 150.

| Configuration        | Runtime  | Speedup S(p) | Efficiency E(p) |
| -------------------- | -------- | ------------ | --------------- |
| **Sequential (1,1)** | 647 s    | -            | -               |
| **(2,2)**            | 301 s    | 2.15         | 0.54            |
| **(1,4)**            | (broken) |              |                 |
| **(3,2)**            | 255 s    | 2.54         | 0.42            |

### Fluidtrap

| Configuration        | Runtime | Speedup S(p) | Efficiency E(p) |
| -------------------- | ------- | ------------ | --------------- |
| **Sequential (1,1)** | 19s     | -            | -               |
| **(2,3)**            | 4s      | 4.75         | 0.8             |
| **(3,2)**            | 5s      | 3.8          | 0.63            |

## Performance Analysis

### Strong Scaling with Rayleigh Benard case

### Weak Scaling with Lid-driven cavity case
