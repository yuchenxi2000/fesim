#!/bin/bash
#SBATCH -o job.%j.out
#SBATCH --partition=C032M0128G
#SBATCH --qos=high
#SBATCH -J BTO
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32

module load intel/2018.1

DIR="/gpfs/share/home/{你的目录}"
cd "$DIR"

SYSTEM_N=10
T=350
STEPS=10000

run_sim() {
    mkdir "./$1"
    case $1 in
    1)
        cat > "./sim.ini" << EOF
[sys]
Nx = $SYSTEM_N
Ny = $SYSTEM_N
Nz = $SYSTEM_N

[init]
u = random
v = random
eta_h0 = 0.012
eta_h1 = 0.012
eta_h2 = 0.012
eta_h3 = 0.000
eta_h4 = 0.000
eta_h5 = 0.000
q = ./cache_q_10

[sim]
T = $T
steps = $STEPS
pressure = -5.0e+9

[monitor]
eta_h = ./$1/monitor_eta_h
; e = ./$1/monitor_e
eta_h_step = 1
e_step = 0

[out]
u = ./$1/out_u.bin
v = ./$1/out_v.bin
eta_h = ./$1/out_eta_h.bin

EOF
        ./main
        cp "./$1/out_u.bin" "./out_u.bin"
        cp "./$1/out_v.bin" "./out_v.bin"
        cp "./$1/out_eta_h.bin" "./out_eta_h.bin"
        ;;
    *)
        cat > "./sim.ini" << EOF
[sys]
Nx = $SYSTEM_N
Ny = $SYSTEM_N
Nz = $SYSTEM_N

[init]
u = ./out_u.bin
v = ./out_v.bin
eta_h = ./out_eta_h.bin
q = ./cache_q_10

[sim]
T = $T
steps = $STEPS
pressure = -5.0e+9

[monitor]
eta_h = ./$1/monitor_eta_h
; e = ./$1/monitor_e
eta_h_step = 1
e_step = 0

[out]
u = ./$1/out_u.bin
v = ./$1/out_v.bin
eta_h = ./$1/out_eta_h.bin

EOF
        ./main
        cp "./$1/out_u.bin" "./out_u.bin"
        cp "./$1/out_v.bin" "./out_v.bin"
        cp "./$1/out_eta_h.bin" "./out_eta_h.bin"
        ;;
    esac
}


for i in {1..60}; do
    echo "temperature step $i: T = $T"
    run_sim $i
    T=$(expr $T - 5)
done
