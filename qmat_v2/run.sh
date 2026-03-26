#!/bin/bash
#SBATCH -o job.%j.out
#SBATCH --partition=C032M0128G
#SBATCH --qos=high
#SBATCH -J qmat
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32

module load intel/2018.1
module add openmpi/3.0.0-intel2018.0

DIR="/gpfs/share/home/{你的目录}"
cd "$DIR"

HOST_FILE="./slurm.hosts"
srun hostname -s | sort -n > "$HOST_FILE"

mpirun -n 32 -machinefile "$HOST_FILE" ./qmat

