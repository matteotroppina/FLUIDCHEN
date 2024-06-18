clear  # clear terminal
rm ../example_cases/ShearFlow/ShearFlow_Output/* #remove old output files
make && mpirun -np 5 ./fluidchen ../example_cases/ShearFlow/ShearFlow.dat 5 1