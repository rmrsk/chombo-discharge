#!/bin/bash
#SBATCH --account=nn12041k
#SBATCH --job-name=hybrid
#SBATCH --time=0-00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=16
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
executable=main3d.Linux.64.mpic++.gfortran.OPTHIGH.MPI.OPENMPCC.ex

# Input file
input_file=example3d.inputs

# Create directory, copy executable and input script
workdir=$USERWORK/$SLURM_JOB_ACCOUNT/OMP/$SLURM_JOB_NAME
if [ ! -d $workdir ]; then mkdir -p $workdir; fi

# Only pout.0
export CH_TIMER=1

# Lustre striping
lfs setstripe --stripe-count 4 --stripe-size 32M $workdir

# Copy executables to work directory
cp -f $executable $workdir/$executable
cp -f $input_file $workdir/$input_file

# Navigate to directory and run
cd $workdir

export OMP_SCHEDULE="dynamic,1"
export OMP_NUM_THREADS=16
export OMP_PROC_BIND=true
export OMP_PLACES=cores

# we added these two lines:
module load PerfReports/23.1.2
echo "set sysroot /" > gdbfile
export ALLINEA_DEBUGGER_USER_FILE=gdbfile

perf-report mpirun --report-bindings --rank-by core --map-by numa:PE=16 --bind-to core ./$executable $input_file AmrMesh.max_amr_depth=3

  # perf-report srun --report-bindings bash -c '
  #   RANK=$SLURM_PROCID
  #   START_CPU=$((RANK * 16))
  #   END_CPU=$((START_CPU + 15))
    
  #   # Create comma-separated CPU list: 0-15, 16-31, etc.
  #   export OMP_PLACES="{$START_CPU:16:1}"
    
  #   ./main3d.Linux.64.mpic++.gfortran.OPTHIGH.MPI.OPENMPCC.ex example3d.inputs AmrMesh.max_amr_depth=3
  # '

