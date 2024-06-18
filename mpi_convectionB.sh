clear #clear terminal
rm ../example_cases/NaturalConvection/NaturalConvection_b_Output/* #remove old output files
make && mpirun -np 5 ./fluidchen ../example_cases/NaturalConvection/NaturalConvection_b.dat 5 1