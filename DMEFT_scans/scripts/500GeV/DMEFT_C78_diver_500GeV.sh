#!/bin/bash
#
#SBATCH -N 20                              # Number of nodes
#SBATCH --ntasks-per-node 8                # Number of MPI processes per node.
#SBATCH -c 3                               # Number of cores per MPI process.
#SBATCH -t 3-00:00                         # Runtime in D-HH:MM
#SBATCH -p plgrid                          # Partition to submit to
#SBATCH -A gambitv5                        # Account to charge
#SBATCH --signal=SIGUSR1@300                                           # Signal and time before walltime to send
#SBATCH -o runs/DMEFT_C78_diver_500GeV/DMEFT_C78_diver_500GeV.out      # File to which STDOUT will be written
#SBATCH -e runs/DMEFT_C78_diver_500GeV/DMEFT_C78_diver_500GeV.err      # File to which STDERR will be written
#SBATCH --mail-type=ALL                                                # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=sanjay.bloor12@imperial.ac.uk                      # Email to which notifications will be sent

module load plgrid/tools/cmake/3.12.3 plgrid/tools/intel/18.0.0 plgrid/tools/impi/2018 plgrid/libs/gsl/2.4 plgrid/libs/hdf5/1.8.19-serial plgrid/libs/eigen/3.2.7 plgrid/libs/boost/1.66.0 plgrid/tools/python/2.7.14
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
mpiexec -outfile-pattern=runs/DMEFT_C78_diver_500GeV/DMEFT_C78_diver_500GeV.out-%r -errfile-pattern=runs/DMEFT_C78_diver_500GeV/DMEFT_C78_diver_500GeV.err-%r ./gambit -rf DMEFT_scans/yaml_files/500GeV/DMEFT_C78_diver_500GeV.yaml