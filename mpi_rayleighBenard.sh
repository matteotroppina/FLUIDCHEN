clear # clear terminal
rm ../example_cases/RayleighBenard/RayleighBenard_Output/* #remove old output files
make && mpirun -np 5 ./fluidchen ../example_cases/RayleighBenard/RayleighBenard.dat 5 1