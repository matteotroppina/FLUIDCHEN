clear # clear terminal
rm ../example_cases/ChannelWithBFS/ChannelWithBFS_Output/* #remove old output files
make && mpirun --oversubscribe -np 1 ./fluidchen ChannelWithBFS.dat 1 1