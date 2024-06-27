clear # clear terminal
rm -rf ChannelWithObstacle_Output/* #remove old output files
make && mpirun -np 5 ./fluidchen ChannelWithObstacle.dat 5 1