# Dockr

Dockr is a distributed computing application to perform a (mirror image)
phage display in silico.

## Description

Dockr simply runs a rigid AutoDock Vina docking for every file in the
library folder specified in `config.ini`, distributing equally amongst
all computing nodes and running one docking per processor per node
at a time.

## Dependencies

* An implementation of the MPI standard, for distribution on cluster nodes, but also required to run on a single computer. We used [MPICH](https://www.mpich.org/)
* [AutoDock Vina](http://vina.scripps.edu/), for binding affinity calculations
* [MGLTools](http://mgltools.scripps.edu/), for generation of files required by AutoDock Vina
* [PyMOL (Open Source)](https://sourceforge.net/projects/pymol/), for generation of PDB files


## Installation

Dockr is written for Linux. Make sure you have all dependencies installed, then
move on to the following steps:

Go to any directory you like and clone this repository, change into its directory
and compile using make
```bash
cd ilikethisdirectory
git clone https://github.com/kcaliban/Dockr.git
cd Dockr
make
```

Before you can use Dockr you have to configure it. Take a look at `config.ini`
and change the settings accordingly, making sure all directories you
specify exist.

## Usage

### Single computer

Dvelopr is written for computer clusters, it can however be executed on a single
computer.

```bash
mpirun -np 1 ./Dockr : -np 1 ./Dockr -w
```

### Computer cluster

For computation on a computing cluster, you have to specify how many
individual computing nodes (not threads!) you can use:

```bash
mpirun -np 1 ./Dockr : -np NUMNODES ./Dockr -w
```

## License

See LICENSE file

Used libraries:
* inih is written by Ben Hoyt, see src/inih/LICENSE
* cxxopts is written by Jarryd Beck, see src/cxxopts/LICENSE
