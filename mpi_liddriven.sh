cd build
make
mpirun -np 2 ./fluidchen ../example_cases/LidDrivenCavity/LidDrivenCavity.dat 2 1