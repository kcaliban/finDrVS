/* Copyright 2019 Fabian Krause
 *
 * finDrVS Main class
 *
 * Performs a mirror image phage display in a distributed computing cluster
 * (or one computer if one is very very patient)
*/
#ifndef SRC_FINDRVS_H_
#define SRC_FINDRVS_H_
#include <math.h>
#include <omp.h>
#include <mpi.h>
#include <sys/stat.h>
#include <fstream>
#include "inih/INIReader.h"
#include "cxxopts/cxxopts.hpp"
#include "VinaInstance/VinaInstance.h"
#include "Serialization/Serialization.h"
#include "Communication.h"
std::string mgltoolstilitiesPath;
std::string pythonShPath;
std::string workDir;
std::string receptorsdir;
std::vector<std::string> receptors;
std::string library;
std::string vina;
int exhaustiveness;
Info * info;
bool receptorsprep;

#endif  // SRC_FINDRVS_H_
