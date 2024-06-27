clear # clear terminal
rm -rf FluidTrap_Output/* #remove old output files
make && mpirun --oversubscribe -np 6 ./fluidchen FluidTrap.dat 3 2   