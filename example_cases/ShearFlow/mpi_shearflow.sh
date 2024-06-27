clear  # clear terminal
rm -rf ShearFlow_Output/* #remove old output files
make && mpirun -np 5 ./fluidchen ShearFlow.dat 5 1