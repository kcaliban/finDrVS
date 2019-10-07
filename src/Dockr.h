/* Copyright 2019 Fabian Krause
 *
 * Dockr Main class
 *
 * Performs a mirror image phage display in a distributed computing cluster
 * (or one computer if one is very very patient)
*/
#ifndef SRC_DOCKR_H_
#define SRC_DOCKR_H_
#include <omp.h>
#include <mpi.h>
#include <sys/stat.h>
#include <fstream>
#include "inih/INIReader.h"
#include "cxxopts/cxxopts.hpp"
#include "VinaInstance/VinaInstance.h"
std::string mgltoolstilitiesPath;
std::string pythonShPath;
std::string workDir;
std::string receptors;
std::string library;
std::string vina;
int exhaustiveness;
Info * info;
bool receptorsprep;

#endif  // SRC_DOCKR_H_
