clear # clear terminal
rm -rf ChannelWithBFS_Output/* #remove old output files
make && mpirun --oversubscribe -np 1 ./fluidchen ChannelWithBFS.dat 1 1