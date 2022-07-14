#!/bin/bash

#SBATCH --nodes=12
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 20
#SBATCH --time=60:00:00


module load ipd/3.6-2019.3.075

# Generate secret for Parallel Python authentication between server and client
secret=`cat /dev/urandom | tr -cd 'a-f0-9' | head -c 32`

# Pick IP port for Parallel Python communication (IANA: use 49152 to 65535 for dynamic/private ports)
port=`python -c 'import random; print (random.randint(49152,65535))'`

# Add additional directories to the Python search path.
# This must be done before starting ppserver.py!
#export PYTHONPATH="$PYTHONPATH:/users/modellers/hpo/code/mizer2/mizer/python"

# define initial part of filename prodcued after gathering AMM7 data 
name=archerfull_gather_forcing

cd ~/code/fabm-mizer/python

export MKL_NUM_THREADS=1

# run the servers
srun ~/.local/bin/ppserver.py -w $SLURM_CPUS_ON_NODE -t 3600 -s $secret -p $port -r &

# launch the parallel processing
python3 -u run_offline_amm7.py ~/projects/COMFORT_final/${name}.nc ~/projects/COMFORT_final/${name}-offline_mizer_w80_F08 --ppservers $SLURM_JOB_NODELIST:$port --ncpus 0 --secret=$secret --debug --parameters=/users/modellers/hpo/projects/COMFORT_final/mizer-amm7.yaml 

