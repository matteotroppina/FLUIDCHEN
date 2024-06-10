clear # clear terminal
rm -rf ../example_cases/LidDrivenCavity/LidDrivenCavity_Output/* #remove old output files
make && mpirun --oversubscribe -np 6 ./fluidchen ../example_cases/LidDrivenCavity/LidDrivenCavity.dat 3 2