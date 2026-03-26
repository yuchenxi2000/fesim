# FESim

BaTiO3 ferroelectric transition simulation based on effective hamiltonian theory

[中文](https://github.com/yuchenxi2000/fesim/blob/main/README.md)

## Usage

The main program named FESim is under the `fesim` directory, implementing the Monte Carlo simulation of BaTiO3 supercell. The hamiltonian of the system is the effective hamiltonian in the paper [Zhang et al., 1995][1]. CBLAS is used for matrix related calculations, and OpenMP is used for calculating the hamiltonian in parallel. Due to the difference in the CBLAS implementations and considering that the program is simple, the Makefile is not given. However, there is a script `make.sh` which can be used as a reference (though it only works on the super computer I used).

The configuration file for FESim is written in ini format. `fesim/run-BTO.sh` is an example of configuration file, which implements the simulation of ferroelectric phase transition under cooling. The detailed implementation of cooling is that we divide the simulation steps into several stages. For each stage we write a configuration file, and run the simulation under fixed temperature.

> `fesim/run-BTO.sh` will cost about one day (on a 64 core machine). You can decrease STEPS to make it faster, however, STEPS should not be lesser than 500 to get meaningful result.

The program under `qmat` directory named QMat calculates the dipole interaction matrix. For each Monte Carlo step, FESim needs to calculate the system energy which requires the dipole interaction matrix given by QMat. In effective hamiltonian theory, the dipole interaction matrix is approximated to be constant, therefore, we calculated it before simulation to save computational effort. This program depends on CBLAS and OpenMPI library. Similar with FESim, I did not give the Makefile due to the differences in the CBLAS and OpenMPI implementations. However, the scripts under directory can be used as references.

The configuration file for QMat is written in ini format. `qmat_v2/dipole.ini` is an example of configuration file. It should be noted that to calculate the dipole interation matrix used for FESim simulation, the lattice constants ax, ay, az must be set to 1.0. This is because FESim automatically handles lattice constants properly. Also, the system size Nx, Ny, Nz should be the same as that in the configuration file of FESim. 

> [inih][2], a single C header library, is used to handle ini format, which was included in the project.

## Theory

Here I give a brief introduction to the theory behind. Firstly we expand the hamiltonian into polynomials of local phonon modes and strain tensor. The coefficients in the expansion can be calculated by first principle calculation using density functional theory. The local phonon modes and strain tensor can be seen as a representation of atom displacements inside the crystal. The reason to choose phonon modes and strain instead of atom displacements is two-fold. First, the phonon modes have relatively high symmetry which can simplify the expanded hamiltonian. Second, based on soft mode theory of ferroelectricity, the soft phonon mode governs the ferroelectric transition, therefore, we can neglect high-order phonons which reduces computational cost.

After calculating the hamiltonian, we can use Monte Carlo methods (used in the paper [Zhang et al., 1995][1]) or molecular dynamics to simulate the phase transition.

## Results

![](https://github.com/yuchenxi2000/fesim/blob/main/fesim/fig.png)

The strain (which reflects lattice constants) of a 10x10x10 supercell as the function of temperature during cooling process. This result is obtained by `fesim/run-BTO.sh` with plotting script `fesim/plot.py`. It shows the phase transition from C, T, O to R phase when cooling.

![](https://github.com/yuchenxi2000/fesim/blob/main/fesim/fig_prb.png)

As a comparison, the above figure is the result of [Zhang et al., 1995][1]. By the way, the order of strain components in Voigt notation seems wrong in the figure. In section A of the paper, the six strain components are defined to be exx, eyy, ezz, eyz, ezx, exy. However, in the figure they seem to be exx, eyy, ezz, exy, eyz, ezx.

[1]: https://journals.aps.org/prb/abstract/10.1103/PhysRevB.52.6301 "Effective hamiltonian theory"
[2]: https://github.com/benhoyt/inih "Ini library"

