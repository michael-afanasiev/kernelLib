#include <iostream>

class kernel;

class kernel {

public:
  
  kernel (std::string);
  
private:
  
  std::string fileName;
  void openKernelNetcdf ();
  
};