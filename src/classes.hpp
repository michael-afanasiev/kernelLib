#include <iostream>
#include <vector>
#include <cmath>
#include "mpi.h"

// Colour codes for pretty writes.
static const char *rst = "\x1b[0m";
static const char *red = "\x1b[31m";
static const char *grn = "\x1b[32m";
static const char *mgn = "\x1b[35m";
static const char *blu = "\x1b[36m";
// End colour codes.

// Utility functions.
float rad2deg        (float);
float deg2rad        (float);
void radThetaPhi2xyz (float, float, float, float &, float &, float &);
void xyz2RadThetaPhi (float &, float &, float &, float, float, float);

void singlePrint (std::string);
// End utility functions.

// Class definitions.
class kernel;

class kernel {

public:
  
  kernel  (std::string);
  ~kernel ();
  
private:
  
  // Netcdf filename.
  std::string fileName;
  
  // Netcdf details.
  int numWroteProcs;
  int numGLLPoints;
  
  // Kernel values.
  float *rawKernel;
  
  // Kernel dimensions.
  float *radius;
  float *theta;
  float *phi;

  // For rotations.
  float radCenter;
  float phiCenter;
  float thetaCenter;
    
  // Internal functions.
  void openCoordNetcdf      ();
  void openKernelNetcdf     ();
  void findChunkDimensions  ();
  void rotateZaxis          ();
  void rotateYaxis          ();
  void rotateXaxis          ();
    
};