clear # clear terminal
rm ../example_cases/ChannelWithObstacle/UpscaledChannelWithObstacle_Output/* #remove old output files
cmake -DCMAKE_BUILD_TYPE=GPU ..
make && mpirun -np 1 fluidchen ../example_cases/ChannelWithObstacle/UpscaledChannelWithObstacle.dat
