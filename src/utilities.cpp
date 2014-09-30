#include "classes.hpp"

void singlePrint (std::string string) {
  
  // For use with MPI -- prints from a single core.
  
  MPI::COMM_WORLD.Barrier ();  
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
                        
  if (theta < 1)
    theta = deg2rad (1);
                         
  rad = sqrt (x*x + y*y + z*z);
  theta = acos (z / rad);
  phi = atan2 (y, x);

}

bool checkRadius (float min, float max, float val) {
  
  // Tells us if we're between radMin and radMax.
  
  bool within=false;
  if ((val <= max) && (val >= min))
    within = true;
  
  return within;
  
}

bool checkTheta (float min, float max, float val) {

  // Tells us if we're between thetaMin and thetaMax.
  
  bool within=false;
  if ((val <= max) && (val >= min))
    within = true;
  
  return within;
  
}

bool checkPhi (float min, float max, float val) {
  
  // Tells us if we're between phiMin and phiMax.
  
  bool within=false;
  if ((val <= max) && (val >= min))
    within = true;
  
  
  return within;
  
}

float distFromPoint (float &xPoint, float &yPoint, float &zPoint,
                     float &xTest,  float &yTest,  float &zTest) {
                       
  // Determines the euclidian distance from an arbitrary point.
                       
  float x = xPoint - xTest;
  float y = yPoint - yTest;
  float z = zPoint - zTest;
  
  float arg      = x*x + y*y + z*z;
  float distance = sqrt (arg);
  
  return distance;
                                                                                         
} 

std::vector<float> getNormalVector (std::vector<float> &A,
                                    std::vector<float> &B,
                                    std::vector<float> &C) {
                                      
   // Gets the normal vector to 3 points in 3-dimensions. Used to determine the equation of a plane.
                                      
  std::vector<float> AB;
  std::vector<float> AC;
  std::vector<float> n;
  
  AB.resize (3);
  AC.resize (3);
  n.resize  (3);
  
  AB[0] = B[0] - A[0];
  AB[1] = B[1] - A[1];
  AB[2] = B[2] - A[2];
  
  AC[0] = C[0] - A[0];
  AC[1] = C[1] - A[1];
  AC[2] = C[2] - A[2];
  
  n[0] = AB[1] * AC[2] - AB[2] * AC[1];
  n[1] = AB[0] * AC[2] - AB[2] * AC[0] * (-1);
  n[2] = AB[0] * AC[1] - AB[1] * AC[0];
  
  float magnitude = sqrt (n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
  n[0] = n[0] / magnitude;
  n[1] = n[1] / magnitude;
  n[2] = n[2] / magnitude;
  
  return n;
  
}

float projWonV_Dist (std::vector<float> &x, std::vector<float> &v, std::vector<float> &x0) {
  
  // Projects a vector x - x0 onto the plane v.
  
  float dotVW = v[0] * (x[0]-x0[0]) + v[1] * (x[1]-x0[1]) + v[2] * (x[2] - x0[2]);
  float magV  = sqrt (v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
  
  return dotVW / magV;
  
}
                                      
                                      

float oneDimDist (float &xPoint, float &yPoint, float &zPoint,
                  float &xTest,  float &yTest,  float &zTest) {
                    
  // TODO INCOMPLETE.
  
  // float d[3];
  float d1 = xPoint - xTest;
  float d2 = yPoint - yTest;
  float d3 = zPoint - zTest;
  
  // return *std::min_element (d, d+3);
  
}