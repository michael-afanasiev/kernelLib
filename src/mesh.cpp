#include "classes.hpp"

using namespace std;

mesh::mesh (kernel &kern) {
  
  singlePrint ("\nGenerating regular mesh.");
  createMesh (kern);
  
}  

void mesh::createMesh (kernel &kern) {
  
  // MPI variables.
  int myRank = MPI::COMM_WORLD.Get_rank ();
  
  // Determine the discritization in each direction. TODO make variables.
  dx = 500 / R_EARTH;
  dy = 500 / R_EARTH;
  dz = 500 / R_EARTH;
  
  // Determines the number of points in each direction.
  nx = int (((kern.xMax + dx - kern.xMin) / dx) + 1);
  ny = int (((kern.yMax + dx - kern.yMin) / dy) + 1);
  nz = int (((kern.zMax + dx - kern.zMin) / dz) + 1);
    
  // Total gridlock!
  gridSize = nx * ny * nz;
  
  // Report whats up in pretty colours.
  if (myRank == 0 )
    std::cout << mgn << "\tTotal size:\t" << gridSize 
      << "\n\tnx:\t\t" << nx 
      << "\n\tny:\t\t" << ny 
      << "\n\tnz:\t\t" << nz << rst << std::endl;
      
    
  // Set up interpolation variables. TODO Chunk might be unused....
  value = new float [gridSize]();
  chunk = new bool  [gridSize];
  x     = new float [nx]();
  y     = new float [ny]();
  z     = new float [nz]();
  
  // Begin interpolation.
  singlePrint ("\n\x1b[33mInterpolating.\x1b[0m");
  clock_t begin = std::clock ();
  for (size_t i=0; i<nx; i++) {
    for (size_t j=0; j<ny; j++) {
      for (size_t k=0; k<nz; k++) {
        
        // Position in regular mesh.
        float xLoc = kern.xMin + i * dx;
        float yLoc = kern.yMin + j * dy;
        float zLoc = kern.zMin + k * dz;

        // Spherical position in regular mesh.
        float rad, theta, phi;
        xyz2RadThetaPhi (rad, theta, phi, xLoc, yLoc, zLoc);

        // Check whether we're within the kernel bounds.
        bool checkR = checkRadius (kern.radiusMin, kern.radiusMax, rad);
        bool checkT = checkTheta  (kern.thetaMin, kern.thetaMax, theta);
        bool checkP = checkPhi    (kern.phiMin, kern.phiMax, phi);
        
        // Initialize interpolation values. 0 if we're outside (should already be 0).
        // Calculate index stride as well. TODO within may be unused...
        bool within = false;
        float data  = 0.;
        int index   = k + j * nz + i * (nz * ny);
        if (checkR && checkT && checkP) {
          
          // Extract the point from the KDtree.
          kdres *set = kd_nearest3 (kern.tree, xLoc, yLoc, zLoc);
          void  *ind = kd_res_item_data (set);
          int pnt    = * (int *) ind;
          data       = kern.rawKernel[pnt];
          within     = true;

          kd_res_free (set);
                    
        }
              
        // Save the point kernel in value, the physical location in xyz.
        value[index] = data;        
        x[i]         = xLoc;
        y[j]         = yLoc;
        z[k]         = zLoc;
        
      }            
    }    
  }   
   
  // Report aww yeah. 
  MPI::COMM_WORLD.Barrier ();      
  clock_t end = std::clock();
  double elapsed = double (end - begin) / CLOCKS_PER_SEC;  
  if (myRank == 0)
    std::cout << "\x1b[32mDone.\x1b[0m (" << elapsed << " seconds)\n";
  
  if (myRank == 0)
    writeExodus (value, x, y, z, nx, ny, nz);
  
}