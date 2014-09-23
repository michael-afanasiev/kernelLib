#include "classes.hpp"
#include "netcdfcpp.h"

kernel::kernel (std::string fName) {
  
  fileName = fName;
  openKernelNetcdf ();
  
}

void kernel::openKernelNetcdf () {
  
  std::cout << "Opening file: " << fileName << std::flush << std::endl;
  
//  try {
//    
//  }
  
}
