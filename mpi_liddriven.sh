clear
cd build
rm -rf ../example_cases/LidDrivenCavity/LidDrivenCavity_Output/* #remove old output files
make && mpirun --oversubscribe -np 9 ./fluidchen ../example_cases/LidDrivenCavity/LidDrivenCavity.dat 3 3
# mpirun -np 2 ./fluidchen ../example_cases/ShearFlow/ShearFlow.dat 2 1