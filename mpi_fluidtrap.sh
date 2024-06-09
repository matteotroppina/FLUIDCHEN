clear # clear terminal
cd build
rm ../example_cases/FluidTrap/FluidTrap_Output/*
make && mpirun --oversubscribe -np 6 ./fluidchen ../example_cases/FluidTrap/FluidTrap.dat 2 3