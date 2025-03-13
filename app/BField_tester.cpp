#include <iostream>

#include "BField_Handler.h"

// Simple example of reading an input B field map, used in reconstruction
int main(int argc, char** argv) {
  if (argc > 2) {
    std::cerr << "Can only take input filename" << std::endl;
    std::cerr << argv[0] << " B_FIELD_FILENAME" << std::endl;
    return -1;
  }

  std::string filename;
  if (argc == 2) {
    filename = std::string(argv[1]);
  } else {
    filename = std::getenv("TMS_DIR");
    filename += "/inputs/BField_Example.dat";
  }
  BField_Handler field(filename);
}
