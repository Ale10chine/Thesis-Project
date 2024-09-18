#!/bin/bash
#SBATCH --job-name=job9
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
./prova9 1 65 

# back to power saving mode
sudo cpupower frequency-set -g powersave