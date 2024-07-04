clear  # clear terminal
rm ../example_cases/ShearFlow/ShearFlow_Output/* #remove old output files
make && mpirun -np 1 ./fluidchen ../example_cases/ShearFlow/ShearFlow.dat 1 1