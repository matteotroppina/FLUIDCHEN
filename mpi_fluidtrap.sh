clear # clear terminal
rm ../example_cases/FluidTrap/FluidTrap_Output/* #remove old output files
make && mpirun --oversubscribe -np 6 ./fluidchen ../example_cases/FluidTrap/FluidTrap.dat 3 2 