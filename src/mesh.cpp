#include "classes.hpp"

using namespace std;

mesh::mesh (kernel &kern) {
  
  singlePrint ("\nGenerating regular mesh.");
  createMesh (kern);
  smoothMesh (kern);
  
}  

void mesh::createMesh (kernel &kern) {
  
  // MPI variables.
  int myRank = MPI::COMM_WORLD.Get_rank ();
  
  // Determine the discritization in each direction. TODO make variables.
  dx = 100 / R_EARTH;
  dy = 100 / R_EARTH;
  dz = 10 / R_EARTH;
  
  // Determines the number of points in each direction.
  nx = int (((kern.xMax + dx - (kern.xMin - dx)) / dx) + 1);
  ny = int (((kern.yMax + dx - (kern.yMin - dx)) / dy) + 1);
  nz = int (((kern.zMax + dx - (kern.zMin - dx)) / dz) + 1);
    
  // Total gridlock!
  gridSize = nx * ny * nz;
  
  // Report whats up in pretty colours.
  if (myRank == 0 )
    std::cout << mgn << "Total size:\t" << gridSize 
      << "\nnx:\t\t" << nx 
      << "\nny:\t\t" << ny 
      << "\nnz:\t\t" << nz << rst << std::endl;
      
    
  // Set up interpolation variables. TODO Chunk might be unused....
  value = new float [gridSize]();
  reg   = new bool  [gridSize]();
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
        // Calculate index stride as well.
        bool within = false;
        float data  = 0.;
        int index   = k + j * nz + i * (nz * ny);
        if (checkR && checkT && checkP) {
          
          // Extract the point from the KDtree.
          kdres *set = kd_nearest3 (kern.tree, xLoc, yLoc, zLoc);
          void  *ind = kd_res_item_data (set);
          int pnt    = * (int *) ind;
          data       = kern.rawKernel[pnt];
          within     = kern.inReg[pnt];

          kd_res_free (set);
                    
        }
              
        // Save the point kernel in value, the physical location in xyz.
        value[index] = data;  
        reg[index]   = within;      
        x[i]         = xLoc;
        y[j]         = yLoc;
        z[k]         = zLoc;
        
      }            
    }    
  }   

  // Here we can play with the intmerpolated kernels (for debugging).
  for (size_t i=0; i<gridSize; i++) {
      value[i] = 0.;
    }

  for (size_t i=0; i<gridSize; i++) {
    if (reg[i]) {
      value[i] = 1.;
      break;
    }
  }
  //
  // Report aww yeah. 
  MPI::COMM_WORLD.Barrier ();      
  clock_t end = std::clock();
  double elapsed = double (end - begin) / CLOCKS_PER_SEC;  
  if (myRank == 0)
    std::cout << "\x1b[32mDone.\x1b[0m (" << elapsed << " seconds)\n";
  
  if (myRank == 0)
    writeExodus (value, x, y, z, nx, ny, nz, "unSmoothed.ex2");
  
}

void mesh::smoothMesh (kernel &kern) {
  
  int myRank = MPI::COMM_WORLD.Get_rank ();
  
  float xVarSquare = 100 / R_EARTH;
  float yVarSquare = 100 / R_EARTH;
  float zVarSquare = 100 / R_EARTH;
  
  // MPI variables.
  int rank      = MPI::COMM_WORLD.Get_rank ();
  int worldSize = MPI::COMM_WORLD.Get_size ();
  
  // Array which will hold smoothed image. Initialized to zero.
  smoothValue    = new float [gridSize]();
  float *xSmooth = new float [gridSize]();
  float *ySmooth = new float [gridSize]();
  
  // Initialize normalization factors.
  float *normFactors = new float [gridSize]();
  
  // Time.
  clock_t begin = std::clock();
  
  singlePrint ("\n\x1b[33mSmoothing.\x1b[035m\nPass 1.");  
  for (size_t i=0; i<nx; i++) {
    for (size_t j=0; j<ny; j++) {
      for (size_t k=0; k<nz; k++) {
        
        // Mesh index.
        size_t index = k + nz * (j + i * ny);

        // Only go in here is we're within the bounds of the original kernel.
        if (reg[index]) {
        
          // Loop over the x-axis.
          for (size_t ipass=0; ipass<nx; ipass++ ) {
          
            // Parameters for 1D gaussian.
            float dist = x[i] - x[ipass];
            float dist2 = dist*dist;
            float arg = (-1) * dist2 / (2*xVarSquare);
            float shape = exp (arg);
          
            // Convolution along x-axis.
            size_t axIndex  = k + nz * (j + ipass * ny);
            xSmooth[index] += shape * value[axIndex];
          
            // Store normalization factor.
            normFactors[index] += shape;                            
          
          }                  
        }
      }
    }
  }
  
  singlePrint ("Pass 2.");  
  // Smooth the smoothed x-array in the y-direction.
  for (size_t i=0; i<nx; i++) {
    for (size_t j=0; j<ny; j++) {
      for (size_t k=0; k<nz; k++) {
        
        // Mesh index.
        size_t index = k + nz * (j + i * ny);     
                        
        // Only go in here is we're within the bounds of the original kernel.
        if (reg[index]) {

          // Loop over the y-axis.
          for (size_t jpass=0; jpass<ny; jpass++ ) {
          
            // Parameters for 1D gaussian.
            float dist = y[j] - y[jpass];
            float dist2 = dist*dist;
            float arg = (-1) * dist2 / (2*yVarSquare);
            float shape = exp (arg);
          
            // Convolution along y-axis.
            size_t axIndex = k + nz * (jpass + i * ny);
            ySmooth[index] += shape * xSmooth[axIndex];
            
            // Store normalization factor.
            normFactors[index] += shape;      
          
          }                                
        }        
      }
    }
  }
  
  singlePrint ("Pass 3.");  
  // Smooth the smooth x-y along the z-axis.
  for (size_t i=0; i<nx; i++) {
    for (size_t j=0; j<ny; j++) {
      for (size_t k=0; k<nz; k++) {
        
        // Mesh index.
        size_t index = k + nz * (j + i * ny);

        // Only go in here is we're within the bounds of the original kernel.
        if (reg[index]) {
        
          // Loop over the z-axis.
          for (size_t kpass=0; kpass<nz; kpass++ ) {
          
            // Parameters for 1D gaussian.
            float dist = z[k] - z[kpass];
            float dist2 = dist*dist;
            float arg = (-1) * dist2 / (2*zVarSquare);
            float shape = exp (arg);
          
            // Convolution along y-axis.
            size_t axIndex      = kpass + nz * (j + i * ny);
            smoothValue[index] +=  shape * ySmooth[axIndex];
          
            // Store normalization factor.
            normFactors[index] += shape;
          
          }
                    
          // smoothValue[index] = smoothValue[index] / normFactors[index];
        }        
      }
    }
  }
  
  MPI::COMM_WORLD.Barrier ();      
  clock_t end = std::clock();
  double elapsed = double (end - begin) / CLOCKS_PER_SEC;  
  
  if (myRank == 0)
    std::cout << "\x1b[32mDone.\x1b[0m (" << elapsed << " seconds)\n";
  
  float valueSum = 0.;
  for (size_t i=0; i<gridSize; i++) {
    if (reg[i]) {
    valueSum += smoothValue[i];
    // smoothValue[i] = normFactors[i];  
  }
  }
  
  if (myRank==0) {
    cout << "value sum: " << valueSum << endl;
    cout << "max value: " << *std::max_element (smoothValue, smoothValue+gridSize) << endl;
  }
  
  if (myRank == 0)
    writeExodus (smoothValue, x, y, z, nx, ny, nz, "smooth.ex2");
  

  
}