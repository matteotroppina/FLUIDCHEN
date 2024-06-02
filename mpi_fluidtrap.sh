clear # clear terminal
cd build
rm ../example_cases/FluidTrap/FluidTrap_Output/*
make && mpirun --oversubscribe -np 4 ./fluidchen ../example_cases/Flui/LidDrivenCavity.dat 2 2
# mpirun -np 2 ./fluidchen ../example_cases/ShearFlow/ShearFlow.dat 2 1