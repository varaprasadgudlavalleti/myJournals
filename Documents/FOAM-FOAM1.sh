#!/bin/bash

#SBATCH --job-name="FOAM-FOAM1"
#SBATCH --partition="q64"
#SBATCH --ntasks=64
#SBATCH --mincpus=4
#SBATCH --ntasks-per-core=1
#SBATCH --mem-per-cpu=4GB
#SBATCH --time="24:00:00"
#SBATCH --mail-user="au779625@uni.au.dk"
#SBATCH --mail-type=all

# Source the bashrc to load environment and modules
source ~/.bashrc
of2206

# Print key environment variables for verification
echo "LD_LIBRARY_PATH: $LD_LIBRARY_PATH"
echo "WM_MPLIB: $WM_MPLIB"

# Clean previous dynamic code builds to force a fresh compile
rm -rf /home/varaprasadgudlavalleti/FOAM-FOAM1/dynamicCode/*
echo "After cleaning, dynamic library directory:"
ls -la /home/varaprasadgudlavalleti/FOAM-FOAM1/dynamicCode/platforms/linux64GccDPInt32Opt/lib

# Optional: Pause briefly to allow filesystem sync
sleep 5

# Execute Allrun if needed
./Allrun

# Decompose the case for parallel run
decomposePar > log.decomposePar

# Run the simulation and capture output, including dynamic library build logs
srun -n 64 hisDriftFluxFoam -parallel > log.hisDriftFluxFoam

# Reconstruct the case
reconstructPar > log.reconstructPar

# List the dynamic library directory after the simulation
echo "After simulation, dynamic library directory:"
ls -la /home/varaprasadgudlavalleti/FOAM-FOAM1/dynamicCode/platforms/linux64GccDPInt32Opt/lib
