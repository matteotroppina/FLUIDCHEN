clear #clear terminal
rm -rf NaturalConvection_b_Output/* #remove old output files
make && mpirun -np 5 ./fluidchen NaturalConvection_b.dat 5 1