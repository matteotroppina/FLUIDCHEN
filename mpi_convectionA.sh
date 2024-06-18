clear #clear terminal
rm ../example_cases/NaturalConvection/NaturalConvection_a_Output/* #remove old output files
make && mpirun -np 5 ./fluidchen ../example_cases/NaturalConvection/NaturalConvection_a.dat 5 1