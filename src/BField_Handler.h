#ifndef __BFIELD_HANDLER_H__
#define __BFIELD_HANDLER_H__

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

class BField_Handler {

  public:
    BField_Handler() {};
    BField_Handler(std::string filename) {
      SetFile(filename);
    }
    void SetFile(std::string filename);

  protected:
    std::string InputFileName;

    std::vector<double> xvec;
    std::vector<double> yvec;
    std::vector<double> zvec;
    std::vector<double> bvec;

};

#endif
