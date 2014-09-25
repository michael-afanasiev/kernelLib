#include "classes.hpp"

void singlePrint (std::string string) {
  
  // For use with MPI -- prints from a single core.
  
  if (MPI::COMM_WORLD.Get_rank () == 0)
    std::cout << string << std::flush << std::endl;
  
}

float rad2deg (float angleInRad) {
  
  // Converts radians to degrees.
  
  float angleInDeg = angleInRad * 180. / M_PI;
  return angleInDeg;
  
}

float deg2rad (float angleInDeg) {

  // Converts degrees to radians.
  
  float angleInRad = angleInDeg * M_PI / 180.;
  return angleInRad;
  
}

void radThetaPhi2xyz (float &rad, float &theta, float &phi,
                      float &x,   float &y,     float &z) {
                         
  // Converts a point in spherical coordinates to cartesian coordinates.
                         
  x = rad * cos (phi) * sin (theta);
  y = rad * sin (phi) * sin (theta);
  z = rad * cos (theta);                     
                         
}

void xyz2RadThetaPhi (float &rad,  float &theta,  float &phi,
                      float &x,    float &y,      float &z) {
                        
  // Converts a point in cartesian coordiantes to spherical coordinates.
                         
  rad = sqrt (x*x + y*y + z*z);
  theta = acos (z / rad);
  phi = atan2 (y, x);

}

bool checkRadius (float min, float max, float val) {
  
  bool within=false;
  if ((val <= max) && (val >= min))
    within = true;
  
  return within;
  
}

bool checkTheta (float min, float max, float val) {
  
  bool within=false;
  if ((val <= max) && (val >= min))
    within = true;
  
  return within;
  
}

bool checkPhi (float min, float max, float val) {
  
  bool within=false;
  if ((val <= max) && (val >= min))
    within = true;
  
  return within;
  
}