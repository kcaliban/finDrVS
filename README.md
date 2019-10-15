# finDrVS

finDrVS is a distributed computing application to perform a virtual screening on several computers.

finDrVS is part of finDr, a toolset developed to find peptide binders
for reflect, the iGEM Team of Freiburg 2019.

## Description

finDrVS simply runs a rigid AutoDock Vina docking for every file in the
library folder specified in `config.ini`, distributing equally amongst
all computing nodes and running one docking per processor per node
at a time.

## Dependencies

* An implementation of the MPI standard, for distribution on cluster nodes, but also required to run on a single computer. We used [MPICH](https://www.mpich.org/)
* [AutoDock Vina](http://vina.scripps.edu/), for binding affinity calculations
* [MGLTools](http://mgltools.scripps.edu/), for generation of files required by AutoDock Vina
* [PyMOL (Open Source)](https://sourceforge.net/projects/pymol/), for generation of PDB files


## Installation

finDrVS is written for Linux. Make sure you have all dependencies installed, then
move on to the following steps:

Go to any directory you like and clone this repository, change into its directory
and compile using make
```bash
cd ilikethisdirectory
git clone https://github.com/kcaliban/finDrVS.git
cd finDrVS
make
```

Before you can use finDrVS you have to configure it. Take a look at `config.ini`
and change the settings accordingly, making sure all directories you
specify exist.

## Usage

### Single computer

Dvelopr is written for computer clusters, it can however be executed on a single
computer.

Change to the directory you created in the installation step and run
the following command:
```bash
mpirun -np 1 ./finDrVS : -np 1 ./finDrVS -w
```

### Computer cluster

For computation on a computing cluster, you have to specify how many
individual computing nodes (not threads!) you can use.

Change to the directory you created in the installation step and run
the following command:
```bash
mpirun -np 1 ./finDrVS : -np NUMNODES ./finDrVS -w
```

## License

See LICENSE file

Used libraries:
* inih is written by Ben Hoyt, see src/inih/LICENSE
* cxxopts is written by Jarryd Beck, see src/cxxopts/LICENSE
