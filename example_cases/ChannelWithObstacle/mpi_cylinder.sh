clear # clear terminal
rm ../example_cases/ChannelWithObstacle/Cylinder_Output/* #remove old output files
make && mpirun -np 1 ./fluidchen ../example_cases/ChannelWithObstacle/Cylinder.dat 1 1