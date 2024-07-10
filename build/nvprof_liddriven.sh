clear # clear terminal
rm -rf ../example_cases/LidDrivenCavity/LidDrivenCavity_Output/* #remove old output files
nvprof --profile-child-processes mpirun -np 1 fluidchen ../example_cases/LidDrivenCavity/LidDrivenCavity.dat