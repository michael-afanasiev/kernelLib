#include "classes.hpp"

using namespace std;

mesh::mesh (kernel &kern) {
  
  singlePrint ("\nGenerating regular mesh.");
  createMesh (kern);
  
}  

void mesh::createMesh (kernel &kern) {
  
  int myRank = MPI::COMM_WORLD.Get_rank ();
  
  dx = 50 / R_EARTH;
  dy = 50 / R_EARTH;
  dz = 10 / R_EARTH;
  
  nx = int (((kern.xMax + dx - kern.xMin) / dx) + 1);
  ny = int (((kern.yMax + dx - kern.yMin) / dy) + 1);
  nz = int (((kern.zMax + dx - kern.zMin) / dz) + 1);
    
  gridSize = nx * ny * nz;
  
  if (myRank == 0 )
    std::cout << mgn << "\tTotal size:\t" << gridSize 
      << "\n\tnx:\t\t" << nx 
      << "\n\tny:\t\t" << ny 
      << "\n\tnz:\t\t" << nz << rst << std::endl;
      
  value = new float [gridSize]();
  x     = new float [nx]();
  y     = new float [ny]();
  z     = new float [nz]();
  
  singlePrint ("\n\x1b[33mInterpolating.\x1b[0m");
  clock_t begin = std::clock ();
  for (size_t i=0; i<nx; i++) {
    for (size_t j=0; j<ny; j++) {
      for (size_t k=0; k<nz; k++) {
        
        float xLoc = kern.xMin + i * dx;
        float yLoc = kern.yMin + j * dy;
        float zLoc = kern.zMin + k * dz;

        float rad, theta, phi;
        xyz2RadThetaPhi (rad, theta, phi, xLoc, yLoc, zLoc);

        bool checkR = checkRadius (kern.radiusMin, kern.radiusMax, rad);
        bool checkT = checkTheta  (kern.thetaMin, kern.thetaMax, theta);
        bool checkP = checkPhi    (kern.phiMin, kern.phiMax, phi);
        
        float data = 0.;
        if (checkR && checkT && checkP) {
          
          kdres *set = kd_nearest3 (kern.tree, xLoc, yLoc, zLoc);
          void  *ind = kd_res_item_data (set);
          int pnt    = * (int *) ind;
          data       = kern.rawKernel[pnt];

          kd_res_free (set);
          
        }
        
        int index    = k + j * nx + i * (nz * ny);
        value[index] = data;
        x[i]         = xLoc;
        y[j]         = yLoc;
        z[k]         = zLoc;
        
      }            
    }    
  }
    
  MPI::COMM_WORLD.Barrier ();      
  clock_t end = std::clock();
  double elapsed = double (end - begin) / CLOCKS_PER_SEC;  
  // cout << "Loop took: " << elapsed << " seconds for processor " << MPI::COMM_WORLD.Get_rank ()
    // << std::endl;

  if (myRank == 0)
    std::cout << "\x1b[32mDone.\x1b[0m (" << elapsed << " seconds)\n";
  
  if (myRank == 0)
    writeExodus (value, x, y, z, nx, ny, nz);
  
}  
