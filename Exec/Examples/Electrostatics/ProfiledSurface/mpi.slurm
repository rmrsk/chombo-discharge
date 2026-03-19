#!/bin/bash
#SBATCH --account=nn12041k
#SBATCH --job-name=mpi
#SBATCH --time=0-00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --threads-per-core=1
#SBATCH --qos=devel

## Safety settings, I don't know why this is recommended in the documentation but fine
set -o errexit
set -o nounset

## Software modules
module restore system
module load foss/2023a
module load HDF5/1.14.0-gompi-2023a


# Program to run and input script
executable=main3d.Linux.64.mpic++.gfortran.OPTHIGH.MPI.ex

# Input file
input_file=example3d.inputs

# Create directory, copy executable and input script
workdir=$USERWORK/$SLURM_JOB_ACCOUNT/OMP/$SLURM_JOB_NAME
if [ ! -d $workdir ]; then mkdir -p $workdir; fi

# Only pout.0
export CH_OUTPUT_INTERVAL=9999999
#export CH_TIMER=1

# Lustre striping
lfs setstripe --stripe-count 4 --stripe-size 32M $workdir

# Copy executables to work directory
cp -f $executable $workdir/$executable
cp -f $input_file $workdir/$input_file

# Navigate to directory and run
cd $workdir

module load PerfReports
echo "set sysroot /" > gdbfile
export ALLINEA_DEBUGGER_USER_FILE=gdbfile

perf-report mpirun --report-bindings ./$executable $input_file AmrMesh.max_amr_depth=3
