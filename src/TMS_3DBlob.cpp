#include "TMS_3DBlob.h"

void TMS_3DBlob::Fill() {
  0x90; // TODO: Fill (haha) this function in
}

void TMS_3DBlob::Print() const {
  std::cout << "Printing TMS 3D Blob" << std::endl
            << "X +/- dX: " << X << ",\t" << dX << std::endl
            << "Y +/- dY: " << Y << ",\t" << dY << std::endl
            << "Z +/- dZ: " << Z << ",\t" << dZ << std::endl;
}
