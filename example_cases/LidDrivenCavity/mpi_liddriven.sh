//clear # clear terminal
rm -rf LidDrivenCavity_Output/* #remove old output files
make && mpirun -np 1 ./fluidchen LidDrivenCavity.dat 1 1