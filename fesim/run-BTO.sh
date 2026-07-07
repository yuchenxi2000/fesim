#!/bin/bash
# Cooling-curve simulation: 350 K down to 55 K, 5 K per step.
# Usage: fill in the variables below, then submit with sbatch or run directly.

# ===== user settings =====
CALC_DIR="/path/to/your/workdir"        # working directory (output goes here)
MAIN_BIN="/path/to/build/fesim/main"    # FESim binary
CACHE_Q="/path/to/cache_q_10"           # dipole matrix from QMat
# =========================

# Load environment (uncomment and adjust for your cluster)
# module load intel openmpi

SYSTEM_N=10
T=350
STEPS=10000

mkdir -p "$CALC_DIR"
cd "$CALC_DIR" || exit 1

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
q = $CACHE_Q

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
        "$MAIN_BIN" ./sim.ini
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
q = $CACHE_Q

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
        "$MAIN_BIN" ./sim.ini
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
