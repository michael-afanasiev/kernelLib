#include "classes.hpp"
#include <netcdf>

using namespace std;

kernel::kernel (std::string fName) {
  
  // Initialize and split kernels.
  
  fileName = fName;
  openCoordNetcdf      ();
  openKernelNetcdf     ();
  
  singlePrint ("Rotating chunks to z axis.");
  findChunkDimensions ();
  rotateZaxis         ();
  findChunkDimensions ();
  rotateYaxis         ();
  findChunkDimensions ();
    
  singlePrint ("Creating KD-tree.");
  createKDTree ();
  
  MPI::COMM_WORLD.Barrier ();
  
}

kernel::~kernel () {
  
  delete [] rawKernel;
  delete [] radius;
  delete [] theta;
  delete [] phi;
  
}

void kernel::createKDTree () {
  
  // Initialize tree.
  tree = kd_create (3);
  
  // Initialize the data array.
  KDdat = new int [numGLLPoints];
  
  // Local xyz values. May not be appropriate.
  float x, y, z;
  for (size_t i=0; i<numGLLPoints; i++) {
    
    radThetaPhi2xyz (radius[i], theta[i], phi[i], x, y, z);
    KDdat[i] = i;
    kd_insert3 (tree, x, y, z, &KDdat[i]);
      
  }
  
}

void kernel::rotateXaxis () {
  
  // Rotates individual chunks around the x-axis.
  
  float xAng = thetaCenter * (1);
  
  // Define rotation matrix about axis.
  float rx11 = 1;
  float rx12 = 0;
  float rx13 = 0;
  float rx21 = 0;
  float rx22 = cos (xAng);
  float rx23 = sin (xAng) * (-1);
  float rx31 = 0; 
  float rx32 = sin (xAng);
  float rx33 = cos (xAng);

  for ( size_t i=0; i<numGLLPoints; i++ ) {
    
    // Local values for rotations.
    float x, y, z;
    
    // Convert (saved) spherical coordinate points to cartesian.
    radThetaPhi2xyz (radius[i], theta[i], phi[i], x, y, z);

    float xNew = rx11 * x + rx12 * y + rx13 * z;
    float yNew = rx21 * x + rx22 * y + rx23 * z;
    float zNew = rx31 * x + rx32 * y + rx33 * z;
    
    // Convert back, overwriting initial values.
    xyz2RadThetaPhi (radius[i], theta[i], phi[i], xNew, yNew, zNew);
    
  }
  
}

void kernel::rotateYaxis () {
  
  // Rotates individual chunks around the y-axis.
  
  float yAng = thetaCenter * (-1);
  
  // Define rotation matrix about axis.
  float ry11 = cos (yAng);
  float ry12 = 0;
  float ry13 = sin (yAng);
  float ry21 = 0;
  float ry22 = 1;
  float ry23 = 0;
  float ry31 = sin (yAng) * (-1);
  float ry32 = 0; 
  float ry33 = cos (yAng);

  for ( size_t i=0; i<numGLLPoints; i++ ) {
    
    // Local values for rotations.
    float x, y, z;
    
    // Convert (saved) spherical coordinate points to cartesian.
    radThetaPhi2xyz (radius[i], theta[i], phi[i], x, y, z);

    float xNew = ry11 * x + ry12 * y + ry13 * z;
    float yNew = ry21 * x + ry22 * y + ry23 * z;
    float zNew = ry31 * x + ry32 * y + ry33 * z;
    
    // Convert back, overwriting initial values.    
    xyz2RadThetaPhi (radius[i], theta[i], phi[i], xNew, yNew, zNew);
    
  }
  
}

void kernel::rotateZaxis () {
  
  // Rotates individual chunks around the z-axis.
  
  float zAng = phiCenter * (-1);
  
  // Define rotation matrix about axis.
  float rz11 = cos (zAng);
  float rz12 = sin (zAng) * (-1);
  float rz13 = 0; 
  float rz21 = sin (zAng);
  float rz22 = cos (zAng);
  float rz23 = 0;
  float rz31 = 0;
  float rz32 = 0;
  float rz33 = 1;

  for ( size_t i=0; i<numGLLPoints; i++ ) {
    
    // Local values for rotations.
    float x, y, z;

    // Convert (saved) spherical coordinate points to cartesian.
    radThetaPhi2xyz (radius[i], theta[i], phi[i], x, y, z);

    // Rotate.
    float xNew = rz11 * x + rz12 * y + rz13 * z;
    float yNew = rz21 * x + rz22 * y + rz23 * z;
    float zNew = rz31 * x + rz32 * y + rz33 * z;
    
    // Convert back, overwriting initial values.
    xyz2RadThetaPhi (radius[i], theta[i], phi[i], xNew, yNew, zNew);
    
  }
  
}

void kernel::findChunkDimensions () {
  
  // Find the middle of a chunk. This is done by averaging all cartesian vectors in the chunk,
  // normalizing the result, and transforming this normalized point to spherical coordinates.
  
  // Local xyz for coordinate transform.  
  float x, y, z;
  
  // Initizlize component averages.
  float xSum=0;
  float ySum=0;
  float zSum=0;
  
  // Initialize rad averages.
  float rSum=0;
  
  // Initialize xyzValues.  
  radThetaPhi2xyz (radius[0], theta[0], phi[0], xMin, yMin, zMin);
  radThetaPhi2xyz (radius[0], theta[0], phi[0], xMax, yMax, zMax);
  
  for (size_t i=0; i<numGLLPoints; i++) {
    
    radThetaPhi2xyz (radius[i], theta[i], phi[i], x, y, z);
    
    // Cartesian box extremes.
    if (x < xMin)
      xMin = x;
    
    if (y < yMin)
      yMin = y;
    
    if (z < zMin)
      zMin = z;
    
    if (x > xMax)
      xMax = x;
    
    if (y > yMax)
      yMax = y;
    
    if (z > zMax)
      zMax = z;        
    
    // Sum radius for average rad.
    rSum += radius[i];
    
    // Sum components.
    xSum += x;
    ySum += y;
    zSum += z;
    
  }
    
  // Calculate magnitude and normalize.
  float magnitude = xSum*xSum + ySum*ySum + zSum*zSum;
  float xCenter   = xSum / magnitude;
  float yCenter   = ySum / magnitude;
  float zCenter   = zSum / magnitude;
  
  // Return center point.
  xyz2RadThetaPhi (radCenter, thetaCenter, phiCenter, xCenter, yCenter, zCenter);
  
  // Get average radius.
  radCenter = rSum / numGLLPoints;
  
  // Get extreme spherical coordinates.
  radiusMin = *std::min_element (radius, radius+numGLLPoints);
  thetaMin  = *std::min_element (theta, theta+numGLLPoints);
  phiMin    = *std::min_element (phi, phi+numGLLPoints);

  radiusMax = *std::max_element (radius, radius+numGLLPoints);
  thetaMax  = *std::max_element (theta, theta+numGLLPoints);
  phiMax    = *std::max_element (phi, phi+numGLLPoints);
  
}

void kernel::openCoordNetcdf () {
  
  // Open the NetCDF coordinate file output from the solver.
  
  using namespace netCDF;
  using namespace netCDF::exceptions;

  // Local mpi variables.
  int myRank = MPI::COMM_WORLD.Get_rank ();
  int worldSize = MPI::COMM_WORLD.Get_size ();
  
  if (myRank == 0)
    std::cout << "Opening coordinate file: " << blu << "./krn/xyzCrustMantle.nc"  << rst 
      << " with " << worldSize << " processors." << std::flush << std::endl;
  
  try {

    std::string coordName = "./krn/xyzCrustMantle.nc";

    // Open the file.
    NcFile dataFile (coordName, NcFile::read);
    
    // Get variable.
    NcVar NcRadius = dataFile.getVar ("radius");
    NcVar NcTheta  = dataFile.getVar ("theta");
    NcVar NcPhi    = dataFile.getVar ("phi");
    
    // Get array sizes.
    NcDim procDim  = NcRadius.getDim (0);
    NcDim coordDim = NcRadius.getDim (1);  
    numWroteProcs  = procDim.getSize ();
    numGLLPoints   = coordDim.getSize ();
  
    if (myRank == 0)
      std::cout << mgn << "\tNumber of solver processers:\t " << numWroteProcs
        << "\n\tNumber of GLL points:\t\t " << numGLLPoints << rst << "\n" << std::endl;
    
    // Current error handling. Only can have as many cores as we did for the simulation.
    if (worldSize > numWroteProcs) {
      
      if (myRank == 0)
        std::cout << red << "Currently, you can only use as many cores for processing as were used " 
          << "in the forward run. Exiting.\n" << rst << std::flush << std::endl;
      
      MPI::COMM_WORLD.Abort (EXIT_FAILURE);
      exit (EXIT_FAILURE);
      
    }
    
    // Set up the MPI read chunk array.
    std::vector<size_t> start;
    std::vector<size_t> count;
    start.resize(2);
    count.resize(2);
    
    // Row major MPI read. Start at [myRank, 0]
    start[0] = myRank;
    start[1] = 0;
    
    // Read until end of line [myrank, numGLLPoints]
    count[0] = 1;
    count[1] = numGLLPoints;
    
    // Of course only read in with the number of processors used to create the file.
    if (myRank < numWroteProcs) {
      radius = new float [numGLLPoints];
      theta  = new float [numGLLPoints];
      phi    = new float [numGLLPoints];
      
      NcRadius.getVar (start, count, radius);
      NcTheta.getVar  (start, count, theta);
      NcPhi.getVar    (start, count, phi);
    }
    
    // Destructor will close file.
        
  } catch (NcException &error) {
    
    std::cout << error.what() << std::endl;
    std::cout << red << "Failure reading: " << fileName << std::endl;
    std::exit (EXIT_FAILURE);
    
  }
  
}

void kernel::openKernelNetcdf () {
  
  // Open the kernel file output from the solver.
  
  using namespace netCDF;
  using namespace netCDF::exceptions;

  // Local mpi variables.
  int myRank = MPI::COMM_WORLD.Get_rank ();
  int worldSize = MPI::COMM_WORLD.Get_size ();
  
  if (myRank == 0)
    std::cout << "Opening kernel file: " << blu << fileName << rst << " with " << worldSize 
      << " processors." << std::flush << std::endl;
  
  try {

    // Open the file.
    NcFile dataFile (fileName, NcFile::read);
    
    // Get variable.
    NcVar NcKernel = dataFile.getVar ("rawKernel");
    
    // Get array sizes.
    NcDim procDim = NcKernel.getDim (0);
    NcDim kernDim = NcKernel.getDim (1);  
    numWroteProcs = procDim.getSize ();
    numGLLPoints  = kernDim.getSize ();
  
    if (myRank == 0)
      std::cout << mgn << "\tNumber of solver processers:\t " << numWroteProcs
        << "\n\tNumber of GLL points:\t\t " << numGLLPoints << rst << "\n" << std::endl;
    
    // Set up the MPI read chunk array.
    std::vector<size_t> start;
    std::vector<size_t> count;
    start.resize(2);
    count.resize(2);
    
    // Row major MPI read. Start at [myRank, 0]
    start[0] = myRank;
    start[1] = 0;
    
    // Read until end of line [myrank, numGLLPoints]
    count[0] = 1;
    count[1] = numGLLPoints;
    
    // Of course only read in with the number of processors used to create the file.
    if (myRank < numWroteProcs) {
      rawKernel = new float [numGLLPoints];
      NcKernel.getVar (start, count, rawKernel);
    }
        
    // Destructor will close file.
    
  } catch (NcException &error) {
    
    std::cout << error.what() << std::endl;
    std::cout << red << "Failure reading: " << fileName << std::endl;
    std::exit (EXIT_FAILURE);
    
  }
    
}