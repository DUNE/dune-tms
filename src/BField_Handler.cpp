#include "BField_Handler.h"

void BField_Handler::SetFile(std::string filename) {
  std::cout << "BField handler loading file " << filename << "..." << std::endl;
  InputFileName = filename;

  std::string line;
  std::ifstream Input(InputFileName.c_str());

  if (Input.is_open()) {
    while (std::getline(Input, line)) {
      // Check for a comment, or lines that are just a space
      if (line[0] == '#' || line.empty()) continue;
      // Chop the line up and get the x,y,z,B values
      double x, y, z = -999.999;
      // Each entry separated by a space
      int nentries = 0;
      while (line.find_first_of(" ") != std::string::npos) {
        int newpos = line.find_first_of(" ");
        std::string substring = line.substr(0, newpos);
        line = line.substr(newpos+1, line.size());
        // First entry is x, second y, third z, etc
        if      (nentries == 0) x = std::atof(substring.c_str());
        else if (nentries == 1) y = std::atof(substring.c_str());
        else if (nentries == 2) z = std::atof(substring.c_str());
        nentries++;
      }
      // Final entry is the bfield
      double bfield = std::atof(line.c_str());

      xvec.push_back(x);
      yvec.push_back(y);
      zvec.push_back(z);
      bvec.push_back(bfield);

      // Check if there are four entries (i.e. three spaces)
      if (nentries != 3) {
        std::cerr << "Found more than four entries (three spaces) in B field map " << InputFileName << std::endl;
        std::cerr << "Exiting..." << std::endl;
        throw;
      }
    }
  } else {
    std::cout << "Unable to open B-field file " << filename << std::endl;
    throw;
  }

}


