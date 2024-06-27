clear # clear terminal
rm -rf RayleighBenard_Output/* #remove old output files
make && mpirun -np 16 ./fluidhcen RayleighBenard.dat 4 4