// Copyright 2019 Fabian Krause
#include "Dockr.h"

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
  command.append(" | cat ");
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
  // Return vector of every filename
  std::vector<std::string> result;
  std::string line;
  std::stringstream outputStream(output);
  while (std::getline(outputStream, line, '\n')) {
    result.push_back(line);
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
  command.append(" -h");
  command.append(" 2>/dev/null 1>&2");
  int success = system(command.c_str());
  if (success != 0) {
    std::cout << "Cannot start AutoDock Vina!\n"
                   "Check your config.ini" << std::endl;
    exit(-1);
  }
}

void work() {
  /* Prepare receptors */
  std::vector<std::string> receptorfiles = getReceptorsM(receptors,
                                                        receptorsprep);
  if (!receptorsprep) {
    for (auto s : receptorfiles) {
      prepareConfig(s);
      try {
        prepareReceptor(s);
      } catch (std::exception& e) {
        info->errorMsg(e.what(), true);
      }
    }
  }
  /**************/
  /* Get ligands */
  std::vector ligands = getLibrary(library);
}

void letOthersWork() {
  /* Receive workload, spread jobs amongst OpenMP threads */
}

int main(int argc, char *argv[]) {
  /* Determine whether to run in Worker mode */
  cxxopts::Options options("Dockr", "In silico phage display");
  options.add_options()
    ("w", "Launch process in worker-mode", cxxopts::value<bool>()->default_value("false"));
  bool worker;
  try {
    auto result = options.parse(argc, argv);
    worker = result["w"].as<bool>();

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
                     "Check if it exists in the same dir as Dockr";
        return 1;
  }
  workDir = reader.Get("Dockr", "workdir", "");
  check(workDir);
  vina = reader.Get("Dockr", "vina", "vina");
  checkExecutable(vina);
  mgltoolstilitiesPath = reader.Get("Dockr", "MGLToolsUtilities", "");
  check(mgltoolstilitiesPath);
  pythonShPath = reader.Get("Dockr", "pythonsh", "");
  checkExecutable(pythonShPath);
  library = reader.Get("Dockr", "library", "");
  check(library);
  receptors = reader.Get("Dockr", "receptors", "");
  check(receptors);
  exhaustiveness = reader.GetInteger("Dockr", "exhaustiveness", 8);
  receptorsprep = reader.GetBoolean("Dockr", "receptorsprep", false);
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
    work();
  } else {
    letOthersWork();
  }
  /*********************/
  MPI_Finalize();
  return 0;
}

