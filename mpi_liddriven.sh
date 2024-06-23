//clear # clear terminal
rm -rf ../example_cases/LidDrivenCavity/LidDrivenCavity_Output/* #remove old output files
make && mpirun -np 1 ./fluidchen ../example_cases/LidDrivenCavity/LidDrivenCavity.dat 1 1