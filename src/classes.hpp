#include <algorithm>
#include <iostream>
#include <vector>
#include <cmath>
#include <ctime>

#include "kdtree.h"
#include "mpi.h"

// Global variables.
static const float R_EARTH=6371.0;
// End global variables.

// Colour codes for pretty writes.
static const char *rst = "\x1b[0m";
static const char *red = "\x1b[31m";
static const char *grn = "\x1b[32m";
static const char *yel = "\x1b[33m";
static const char *mgn = "\x1b[35m";
static const char *blu = "\x1b[36m";
// End colour codes.

// Utility functions.
float rad2deg        (float);
float deg2rad        (float);

bool checkRadius (float, float, float);
bool checkTheta  (float, float, float);
bool checkPhi    (float, float, float);

void radThetaPhi2xyz (float &, float &, float &, float &, float &, float &);
void xyz2RadThetaPhi (float &, float &, float &, float &, float &, float &);

void singlePrint (std::string);

void writeExodus (float *regMeshArr, float *regX, float *regY, float *regZ, 
                  int &nx, int &ny, int &nz);
// End utility functions.

// Class definitions.
class kernel;
class mesh;

class kernel {
  
  friend class mesh;

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
  float *theta; // latitude.
  float *phi;   // longitude.
  
  // Spherical extremes;
  float radiusMin, thetaMin, phiMin;
  float radiusMax, thetaMax, phiMax;
  
  // Cartesian extremes.
  float xMin, yMin, zMin;
  float xMax, yMax, zMax;

  // For rotations.
  float radCenter;
  float phiCenter;
  float thetaCenter;
    
  // KDTree.
  int *KDdat;
  kdtree *tree;
    
  // Internal functions.
  void openCoordNetcdf      ();
  void openKernelNetcdf     ();
  void findChunkDimensions  ();
  void rotateZaxis          ();
  void rotateYaxis          ();
  void rotateXaxis          ();
  void createKDTree         ();
    
};

class mesh {

public:
  
  mesh (kernel &);
  
private:
  
  // Discritization.
  float dx, dy, dz;
  int nx, ny, nz;
  long gridSize;
  
  // Internal functions.
  void createMesh (kernel &);
  
  // Coordinate arrays.
  float *x;
  float *y;
  float *z;
  
  // Parameter.
  float *value;

};