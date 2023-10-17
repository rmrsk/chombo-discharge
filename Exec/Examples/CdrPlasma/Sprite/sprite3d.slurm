#!/bin/bash
#SBATCH --account=nn9887k
#SBATCH --job-name=sprite3d_no_filter_coarsen
#SBATCH --time=0-0:30:00
#SBATCH --nodes=300
#SBATCH --ntasks-per-node=128

## Safety settings, I don't know why this is recommended in the documentation but fine
set -o errexit
set -o nounset

## Software modules
module restore system
module load foss/2020b
module load HDF5/1.10.7-gompi-2020b

# Program to run and input script
executable=program3d.Linux.64.mpic++.gfortran.OPTHIGH.MPI.ex

# Input files. Must be copied to where we run the simulation. 
input_file=sprite3d.inputs
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
lfs setstripe --stripe-count 64 --stripe-size 16M $workdir

# Navigate to directory and run
cd $workdir

#mpirun ./$executable $input_file CdrPlasmaJSON.chemistry_file=$chemistry_file  Driver.max_plot_depth=7 Driver.restart=1200 Driver.initial_regrids=1 Driver.max_steps=1201 Driver.plot_interval=-1 Driver.checkpoint_interval=-1 Driver.profile=true
#mpirun ./$executable $input_file CdrPlasmaJSON.chemistry_file=$chemistry_file Driver.max_steps=2000 Driver.restart=1500 Driver.initial_regrids=1 CdrPlasmaStreamerTagger.buffer = 0 FieldSolverMultigrid.filter=1 AmrMesh.max_amr_depth = 10

# Parse options into file and run.
echo "" >> $input_file
echo "# Simulations-specific options" >> $input_file
echo Driver.max_steps=2500 >> $input_file
echo Driver.restart=2250 >> $input_file
echo Driver.initial_regrids=1 >> $input_file
echo CdrPlasmaGodunovStepper.filter_rho=0 >> $input_file
echo CdrPlasmaStreamerTagger.max_coarsen_lvl=5 >> $input_file
echo McPhoto.num_sampling_packets=2 >> $input_file
echo CdrPlasmaStreamerTagger.buffer=0 >> $input_file
echo AmrMesh.max_amr_depth=12 >> $input_file
mpirun ./$executable $input_file 

