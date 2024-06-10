clear # clear terminal
rm ../example_cases/ChannelWithBFS/ChannelWithBFS_Output/* #remove old output files
make && mpirun --oversubscribe -np 6 ./fluidchen ../example_cases/ChannelWithBFS/ChannelWithBFS.dat 3 2