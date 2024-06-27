clear #clear terminal
rm -rf NaturalConvection_a_Output/* #remove old output files
make && mpirun -np 5 ./fluidchen NaturalConvection_a.dat 5 1