#!/bin/bash
#SBATCH --account=nn9887k
#SBATCH --job-name=filter3d
#SBATCH --time=0-4:00:00
#SBATCH --nodes=16
#SBATCH --ntasks-per-node=64

## Safety settings, I don't know why this is recommended in the documentation but fine
set -o errexit
set -o nounset

## Software modules
module restore system
module load foss/2021a
module load HDF5/1.10.7-gompi-2021a

# Program to run and input script
executable=program3d.Linux.64.mpic++.gfortran.OPTHIGH.MPI.ex

# Input files. Must be copied to where we run the simulation. 
input_file=filter3d.inputs
chemistry_file=air_chemistry3d.json
transport_file=bolsig_air.dat
gas_law_file=ENMSIS_Atmosphere.dat
tanh_file=tanh.dat

# Create directory, copy executable and input script
workdir=$USERWORK/${SLURM_JOB_NAME}
if [ ! -d $workdir ]; then mkdir -p $workdir; fi

# Only pout.0
export CH_OUTPUT_INTERVAL=9999999

# Copy executables to work directory
cp -f $executable $workdir/$executable
cp -f $input_file $workdir/$input_file
cp -f $chemistry_file $workdir/$chemistry_file
cp -f $transport_file $workdir/$transport_file
cp -f $gas_law_file $workdir/$gas_law_file
cp -f $tanh_file $workdir/$tanh_file

# Lustre striping
lfs setstripe --stripe-count 8 --stripe-size 64M $workdir

# Navigate to directory and run
cd $workdir
mpirun ./$executable $input_file Driver.restart=0 FieldSolverMultigrid.filter = 1
