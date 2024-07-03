clear # clear terminal
rm ../example_cases/ChannelWithObstacle/Cylinder_Output/* #remove old output files
make && mpirun -np 2 ./fluidchen ../example_cases/ChannelWithObstacle/Cylinder.dat 2 1