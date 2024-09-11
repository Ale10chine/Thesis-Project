#!/bin/bash
#SBATCH --job-name=job3
#SBATCH --partition=razor
#SBATCH --ntasks=1
#SBATCH --mem=14GB
# warm up processors
sudo cpupower frequency-set -g performance
sleep 0.1
stress-ng -c 4 --cpu-ops=100
# set limits
ulimit -v 16777216


# Command to execute
g++ -I/nfsd/rop/sw/ibm/cos221/concert/include -I/nfsd/rop/sw/ibm/cos221/cplex/include -Wno-deprecated-declarations -std=c++17 -c Utility.cpp -o Utility.o
g++ -I/nfsd/rop/sw/ibm/cos221/concert/include -I/nfsd/rop/sw/ibm/cos221/cplex/include -Wno-deprecated-declarations -std=c++17 -c prova.cpp -o prova.o

g++ Utility.o prova.o  -o prova -L/nfsd/rop/sw/ibm/cos221/cplex/lib/x86-64_linux/static_pic/ -L/nfsd/rop/sw/ibm/cos221/concert/lib/x86-64_linux/static_pic/  -lilocplex -lconcert -lcplex -lm -lpthread -Wno-deprecated-declarations

./prova 12345 65 

# back to power saving mode
sudo cpupower frequency-set -g powersave