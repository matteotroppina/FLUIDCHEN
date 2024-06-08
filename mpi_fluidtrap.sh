clear # clear terminal
cd build
rm ../example_cases/FluidTrap/FluidTrap_Output/*
make && mpirun --oversubscribe -np 6 ./fluidchen ../example_cases/FluidTrap/FluidTrap.dat 3 2 
# mpirun -np 2 ./fluidchen ../example_cases/ShearFlow/ShearFlow.dat 2 1