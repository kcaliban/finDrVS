// Copyright 2019 iGEM Team Freiburg 2019
#include "finDrVS.h"

void mImage(const std::string pdb) {
  // Read receptor file
  std::ifstream t(pdb);
  t.seekg(0, std::ios::end);
  size_t size = t.tellg();
  std::string buffer(size, ' ');
  t.seekg(0);
  t.read(&buffer[0], size);
  t.close();

  std::string output;
  // Read line per line, set min and max accordingly
  std::string line;
  std::stringstream receptorStream(buffer);
  while (std::getline(receptorStream, line, '\n')) {
    std::string newline = line;
    if (line.substr(0, 4) == "ATOM" || line.substr(0, 6) == "HETATM") {
      float x = - stof(line.substr(30, 8));
      char invertedX[9];  // Null terminating char at the end => 9
      sprintf(invertedX, "%8.3f", x);
      for (unsigned int i = 0; i < 8; i++) {
        newline[30 + i] = invertedX[i];
      }
    }
    output.append(newline + "\n");
  }

  std::ofstream outf(pdb, std::ofstream::out | std::ofstream::trunc);
  outf << output;
  outf.close();
}

std::string PDBtoFASTA(std::string filename) {
  std::ifstream file(filename);
  std::unordered_map<std::string, std::string> AA({
    {"ALA", "A"}, {"ARG", "R"}, {"ASN", "N"}, {"ASP", "D"}, {"CYS", "C"},
    {"GLU", "E"}, {"GLN", "Q"}, {"GLY", "G"}, {"HIS", "H"}, {"ILE", "I"},
    {"LEU", "L"}, {"LYS", "K"}, {"MET", "M"}, {"PHE", "F"}, {"PRO", "P"},
    {"SER", "S"}, {"THR", "T"}, {"TRP", "W"}, {"TYR", "Y"}, {"VAL", "V"}
  });

  std::string FASTA;
  std::string line;
  int previd = -1;
  while (getline(file, line)) {
    if (line.substr(0, 4) == "ATOM" || line.substr(0, 6) == "HETATM") {
      std::string aa = line.substr(17, 3);
      int id = stoi(line.substr(22, 4));
      if (id != previd) {
        FASTA.append(AA[aa]);
      }
      previd = id;
    }
  }
  return FASTA;
}


// Get ligand filenames
std::vector<std::string> getLibrary(std::string dir) {
  // Read all pdb files in a directory
  std::string command;
  command.append("ls ");
  command.append(dir);
  command.append(" | cat | grep -v \"\\.pdbqt\"");
  std::string output;
  FILE * lsOutputStream = popen(command.c_str(), "r");
  char buf[1024];
  while (fgets(buf, 1024, lsOutputStream)) {
    output += buf;
  }
  pclose(lsOutputStream);
  // Return vector of every filename
  std::vector<std::string> result;
  std::string line;
  std::stringstream outputStream(output);
  while (std::getline(outputStream, line, '\n')) {
    result.push_back(line);
  }
  return result;
}

// Get receptor filenames
std::vector<std::string> getReceptorsM(std::string dir, bool prep = false) {
  // Read all pdb files in a directory
  std::string command;
  command.append("ls ");
  command.append(dir);
  command.append(" | cat ");
  if (!prep) {
    command.append("| grep -v .pdbqt | grep -v conf");
  } else {
    command.append("| grep -v conf");
  }
  std::string output;
  FILE * lsOutputStream = popen(command.c_str(), "r");
  char buf[1024];
  while (fgets(buf, 1024, lsOutputStream)) {
    output += buf;
  }
  pclose(lsOutputStream);
  // Output of receptors
  // Return vector of every filename
  std::vector<std::string> result;
  std::string line;
  std::stringstream outputStream(output);
  while (std::getline(outputStream, line, '\n')) {
    result.push_back(dir + "/" + line);
  }
  return result;
}

void prepareConfig(std::string receptor) {
  // Generate conf file for docking

  // Read receptor file
  std::ifstream t(receptor);
  t.seekg(0, std::ios::end);
  size_t size = t.tellg();
  std::string buffer(size, ' ');
  t.seekg(0);
  t.read(&buffer[0], size);
  t.close();

  // std::cout << "Read receptor file!" << std::endl;

  // For max and min no sorting is required
  float xmin = std::numeric_limits<float>::infinity();
  float ymin = std::numeric_limits<float>::infinity();
  float zmin = std::numeric_limits<float>::infinity();
  float xmax = - std::numeric_limits<float>::infinity();
  float ymax = - std::numeric_limits<float>::infinity();
  float zmax = - std::numeric_limits<float>::infinity();

  // Read line per line, set min and max accordingly
  std::string line;
  std::stringstream receptorStream(buffer);
  while (std::getline(receptorStream, line, '\n')) {
    if (line.substr(0, 4) == "ATOM" || line.substr(0, 6) == "HETATM") {
      float x = stof(line.substr(31, 8));
      xmin = (x < xmin) ? x : xmin;
      xmax = (x > xmax) ? x : xmax;
      float y = stof(line.substr(39, 8));
      ymin = (y < ymin) ? y : ymin;
      ymax = (y > ymax) ? y : ymax;
      float z = stof(line.substr(47, 8));
      zmin = (z < zmin) ? z : zmin;
      zmax = (z > zmax) ? z : zmax;
    }
  }

  float sizex = xmax - xmin;
  float sizey = ymax - ymin;
  float sizez = zmax - zmin;

  // float offsetx = xmin + sizex / 2;
  float offsetx = xmin + sizex / 2.0;
  float offsety = ymin + sizey / 2.0;
  float offsetz = zmin + sizez / 2.0;

  std::string outfile = receptor + "_conf";

  std::ofstream confFile;
  confFile.open(outfile.c_str(), std::ios::trunc);
  if (!confFile) {
    std::cout << "Could not open config file!" << std::endl;
    exit(-1);
  }
  confFile << "center_x = " << offsetx << std::endl;
  confFile << "center_y = " << offsety << std::endl;
  confFile << "center_z = " << offsetz << std::endl;
  confFile << std::endl;
  confFile << "size_x = " << sizex + 30 << std::endl;
  confFile << "size_y = " << sizey + 30 << std::endl;
  confFile << "size_z = " << sizez + 30 << std::endl;
  confFile.close();
}

void prepareReceptor(std::string receptor) {
  // Generate a PDBQT
  std::string command;
  command.append(pythonShPath);
  command.append(" ");
  command.append(mgltoolstilitiesPath);
  command.append("/prepare_receptor4.py -r ");
  command.append(receptor);
  command.append(" -A bonds_hydrogens -U nphs -o ");
  command.append(receptor);
  command.append("qt");
  int success = system(command.c_str());
  if (success != 0) {
    throw VinaException("Could not generate pdbqt file for receptor", receptor);
  }
}

void prepareLigand(std::string ligand) {
  // Generate a PDBQT
  std::string command;
  command.append(pythonShPath);
  command.append(" ");
  command.append(mgltoolstilitiesPath);
  command.append("/prepare_ligand4.py -l ");
  command.append(ligand);
  command.append(" -Z -A bonds_hydrogens -U nphs -o ");
  command.append(ligand);
  command.append("qt >/dev/null");
  int success = system(command.c_str());
  if (success != 0) {
    throw VinaException("Could not generate pdbqt file for ligand",
                        ligand,
                        "PQT");
  }
}

float genDock(std::string file) {
  float affinity = 10;
  // Prepare ligand
  try {
    prepareLigand(file);
  } catch (...) {
    throw;
  }
  // Do a docking for each receptor
  for (unsigned int i = 0; i < receptors.size(); i++) {
    VinaInstance vinaInstance(vina.c_str(), receptors.at(i).c_str(),
                              file.c_str(), info);
    float recaffinity = vinaInstance.calculateBindingAffinity(exhaustiveness,
                                                              5);
    if (recaffinity < affinity) { affinity = recaffinity; }
  }

  return affinity;
}

void check(const std::string p) {
  struct stat st;
  if (stat(p.c_str(), &st) != 0) {
    std::cout << "File or path \"" + p + "\" does not exist\n"
                   "Check your config.ini" << std::endl;
    exit(-1);
  }
}

void checkExecutable(const std::string e) {
  std::string command;
  command.append(e);
  command.append(" --help");
  command.append(" 2>/dev/null 1>&2");
  int success = system(command.c_str());
  if (success != 0) {
    std::cout << "Cannot start " + e + "!\n"
                   "Check your config.ini" << std::endl;
    exit(-1);
  }
}

std::vector<std::vector<std::string>> distribute(std::vector<std::string>& fils,
                                                 int world_size) {
  /* Equal distribution for lack of innovation at the time of coding */
  std::vector<std::vector<std::string>> distribution;
  unsigned int vecPos = 0;
  for (int i = 0; i < world_size - 1; i++) {
    std::vector<std::string> bucket;
    for (unsigned int j = 0; j < (unsigned int)
                                 ceil((float) fils.size() /
                                      (float) (world_size - 1)); j++) {
      if (vecPos > fils.size() - 1) { break; };
      bucket.push_back(fils.at(vecPos++));
    }
    distribution.push_back(bucket);
    if (vecPos > fils.size() - 1) { break; };
  }
  return distribution;
};

void work(int world_size, int world_rank) {
  /* Receive workload, spread jobs amongst OpenMP threads */
  // Get receptors
  receptors = getReceptorsM(receptorsdir, receptorsprep);
  // Wait to receive a vector of files
  std::vector<std::string> FILES;
  unsigned int FILESSize = 0;

  MPI_Request request;
  MPI_Irecv(&FILESSize, 1, MPI_INT, 0, SENDFILESSIZE,
           MPI_COMM_WORLD, &request);
  MPI_Wait(&request, MPI_STATUS_IGNORE);

  char * tmp = new char[FILESSize];
  MPI_Irecv(&tmp[0], FILESSize, MPI_BYTE, 0, SENDFILESCONT,
           MPI_COMM_WORLD, &request);
  MPI_Wait(&request, MPI_STATUS_IGNORE);
  deserialize(FILES, tmp, FILESSize);
  free(tmp);

  info->infoMsg("Worker #" + std::to_string(world_rank) + "'s workload:" +
                std::to_string(FILES.size()));
  // Do the right thing
  std::vector<std::pair<std::string, float>> results;
  #pragma omp parallel
  #pragma omp for ordered
  for (unsigned int j = 0; j < FILES.size(); j++) {
    // Send back after every 5 results
    // Makes all threads wait until here, then executed sequentially.
    // To avoid accidental removal of result
    // (e.g. A pushes back result, B clears immediately after)
    #pragma omp ordered
    {
      if (results.size() >= 5) {
        // Send back the results
        unsigned int resultsSize;
        tmp = serialize(results, &resultsSize);
        MPI_Isend(&resultsSize, 1, MPI_INT, 0, SENDAFFINSIZE, MPI_COMM_WORLD,
                  &request);
        MPI_Wait(&request, MPI_STATUS_IGNORE);
        MPI_Isend(&tmp[0], resultsSize, MPI_BYTE, 0, SENDAFFINCONT, MPI_COMM_WORLD,
                  &request);
        MPI_Wait(&request, MPI_STATUS_IGNORE);
        free(tmp);
        // Clear results
        results.clear();
      }
    }
    // Dock
    try {
      float aff = genDock(library + "/" + FILES.at(j));
      // Add result
      #pragma omp critical
      results.push_back(std::make_pair(FILES.at(j), aff));
    } catch (...) {
      info->errorMsg("Docking for " + FILES.at(j) +
                     " failed, skipping...", false);
      continue;
    }
  }
  // Send back whatever is left over
  unsigned int resultsSize;
  tmp = serialize(results, &resultsSize);
  MPI_Isend(&resultsSize, 1, MPI_INT, 0, SENDAFFINSIZE, MPI_COMM_WORLD,
            &request);
  MPI_Wait(&request, MPI_STATUS_IGNORE);
  MPI_Isend(&tmp[0], resultsSize, MPI_BYTE, 0, SENDAFFINCONT, MPI_COMM_WORLD,
            &request);
  MPI_Wait(&request, MPI_STATUS_IGNORE);
  free(tmp);
  // Clear results
  results.clear();
  info->infoMsg("Worker #" + std::to_string(world_rank) + " done.");
}

void letOthersWork(int world_size, int world_rank, bool mirrorImage) {
  /* Prepare receptors */
  receptors = getReceptorsM(receptorsdir, receptorsprep);
  if (mirrorImage || !receptorsprep) {
    for (auto s : receptors) {
      if (mirrorImage) {
        mImage(s);
      }
      prepareConfig(s);
      try {
        prepareReceptor(s);
      } catch (std::exception& e) {
        info->errorMsg(e.what(), true);
      }
    }
  }
  /**************/
  /* Get and distribute ligands */
  // Prepare workloads
  std::vector<std::string> ligands = getLibrary(library);
  info->infoMsg("Library size: " + std::to_string(ligands.size()));
  auto distr = distribute(ligands, world_size);
  int kArrSize = world_size - 1;
  // Send buckets to subprocesses
  info->infoMsg("Master is sending his work...");
  MPI_Request requests[kArrSize];
  unsigned int bucketSizes[kArrSize];
  char * bucketBin[kArrSize];
  // Send size in bytes
  for (int i = 1; i < world_size; i++) {
    bucketBin[i - 1] = serialize(distr.at(i - 1), &bucketSizes[i - 1]);
    MPI_Isend(&bucketSizes[i - 1], 1, MPI_INT, i,
              SENDFILESSIZE, MPI_COMM_WORLD, &requests[i - 1]);
  }
  MPI_Waitall(world_size - 1, requests, MPI_STATUS_IGNORE);
  // Send binary data
  for (int i = 1; i < world_size; i++) {
    MPI_Isend(&bucketBin[i - 1][0], bucketSizes[i - 1], MPI_BYTE, i,
              SENDFILESCONT, MPI_COMM_WORLD, &requests[i - 1]);
  }
  MPI_Waitall(world_size - 1, requests, MPI_STATUS_IGNORE);
  // Free memory
  for (int i = 1; i < world_size; i++) {
    free(bucketBin[i - 1]);
  }
  /* Cache results and write every X results */
  std::ofstream ofstream (workDir + "/RESULTS", std::ios::app);
  unsigned int count = 0;
  while (42) {
    if (count >= ligands.size()) {
      break;
    }
    // Get results of each bucket
    std::vector<std::vector<std::pair<std::string, float>>> bucketResults;
    // Get filesize from each
    unsigned int resSize[world_size - 1];
    for (int i = 1; i < world_size; i++) {
      MPI_Irecv(&resSize[i - 1], 1, MPI_INT, i,
                SENDAFFINSIZE, MPI_COMM_WORLD, &requests[i - 1]);
    }
    MPI_Waitall(world_size - 1, requests, MPI_STATUS_IGNORE);
    // Get contents from each
    char * resBin[world_size - 1];
    for (int i = 1; i < world_size; i++) {
      resBin[i - 1] = new char[resSize[i - 1]];
      MPI_Irecv(&resBin[i - 1][0], resSize[i - 1], MPI_BYTE, i, SENDAFFINCONT,
                MPI_COMM_WORLD, &requests[i - 1]);
    }
    MPI_Waitall(world_size - 1, requests, MPI_STATUS_IGNORE);
    // Deserialize binary messages
    for (int i = 1; i < world_size; i++) {
      std::vector<std::pair<std::string, float>> results;
      deserialize(results, resBin[i - 1], resSize[i - 1]);
      bucketResults.push_back(results);
      free(resBin[i - 1]);
    }
    for (int i = 1; i < world_size; i++) {
      std::cout << "From Worker #" << i << ":" << std::endl;
      for (auto j : bucketResults.at(i - 1)) {
        std::cout << j.first << ": " << j.second << std::endl;
      }
    }
    // Write to file
    for (auto buck : bucketResults) {
      for (auto res : buck) {
        ofstream << res.first << "\t" << PDBtoFASTA(library + "/" + res.first)
                 << "\t" << std::to_string(res.second) << "\n";
        count++;
      }
    }
    ofstream.flush();
  }
  ofstream.close();
  info->infoMsg("All results back, check " + workDir + "/RESULTS");
}

int main(int argc, char *argv[]) {
  /* Determine whether to run in Worker mode */
  cxxopts::Options options("finDrVS", "In silico phage display");
  options.add_options()
    ("w", "Launch process in worker-mode", cxxopts::value<bool>()->default_value("false"))
    ("mi",
     "(optional) Convert the target into its mirror-image (L to D or D to L)\n"
     "Target has to be unprepared (just the .pdb file, no .pdbqt and conf)\n"
     "Attention: Original target gets overwritten!"
     , cxxopts::value<bool>()->default_value("false"));
  bool worker;
  bool mirrorImage;
  try {
    auto result = options.parse(argc, argv);
    worker = result["w"].as<bool>();
    mirrorImage = result["mi"].as<bool>();

  } catch (std::exception& e) {
    std::cout << e.what() << std::endl;
    std::cout << options.help() << std::endl;
    exit(-1);
  }
  /*********************/
  /* Read and check config.ini */
  INIReader reader("config.ini");
  if (reader.ParseError() != 0) {
        std::cout << "Can't load 'config.ini'\n"
                     "Check if it exists in the same dir as finDrVS";
        return 1;
  }
  workDir = reader.Get("finDrVS", "workdir", "");
  check(workDir);
  vina = reader.Get("finDrVS", "vina", "vina");
  checkExecutable(vina);
  mgltoolstilitiesPath = reader.Get("finDrVS", "MGLToolsUtilities", "");
  check(mgltoolstilitiesPath);
  pythonShPath = reader.Get("finDrVS", "pythonsh", "");
  checkExecutable(pythonShPath);
  library = reader.Get("finDrVS", "library", "");
  check(library);
  receptorsdir = reader.Get("finDrVS", "receptors", "");
  check(receptorsdir);
  exhaustiveness = reader.GetInteger("finDrVS", "exhaustiveness", 8);
  receptorsprep = reader.GetBoolean("finDrVS", "receptorsprep", false);
  /*********************/
  /* Initialize OpenMPI */
  int world_size, world_rank;
  // Initialize the MPI environment
  MPI_Init(&argc , &argv);
  // Get the number of processes
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  // Get the rank of the process
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  /*********************/
  info = new Info(true, true, workDir + "/LOG");
  if (worker) {
    work(world_size, world_rank);
  } else {
    letOthersWork(world_size, world_rank, mirrorImage);
  }
  /*********************/
  MPI_Finalize();
  return 0;
}

