#include "classes.hpp"
#include <netcdf>

using namespace std;

kernel::kernel (std::string fName) {
  
  // Initialize and split kernels.
  
  fileName = fName;
  openCoordNetcdf      ();
  openKernelNetcdf     ();
  
  // Ugly, but need to rotate to  z-axis to get haze boundaries.
  findChunkDimensions ();  
  rotateZaxis         ();
  findChunkDimensions ();
  rotateYaxis         ();
  findChunkDimensions ();
  findSideSets        ();
  
  singlePrint ("Rotating chunks to z axis.");
  rotateZaxis         ();
  findChunkDimensions ();
  rotateYaxis         ();
  findChunkDimensions ();
    
  singlePrint  ("Creating KD-tree.");
  createKDTree ();
  
  MPI::COMM_WORLD.Barrier ();
  
}

kernel::~kernel () {
  
  delete [] rawKernel;
  delete [] radius;
  delete [] theta;
  delete [] phi;
  
}

void kernel::findSideSets () {
  
  // Find the sidesets of a mesh chunk -- that is the nodes close enough to an edge that might 
  // require communicating.
  
  // Allocate sideSet array.
  sideSet = new bool [numGLLPoints]();
  
  // Re-usable vectors for edges of chunk (for face plane).  
  std::vector<float> A, B, C;
  A.resize (3);
  B.resize (3);
  C.resize (3);

  // Normal face plane vectors.
  std::vector <float> n1, n2, n3, n4;  
  n1.resize (3);
  n2.resize (3);
  n3.resize (3);
  n4.resize (3);
  
  // Save edges of chunks for use to define plane later.
  std::vector<float> p1, p2, p3, p4;
  p1.resize (3);
  p2.resize (3);
  p3.resize (3);
  p4.resize (3);
  
  // Face 1.
  A[0] = xMax; A[1] = yMin; A[2] = zMax;
  B[0] = xMin; B[1] = yMin; B[2] = zMax;
  C[0] = xMax; C[1] = yMin; C[2] = zMin;
  n1 = getNormalVector (A, B, C);
  p1 = A;

  // Face 2.
  A[0] = xMax; A[1] = yMax; A[2] = zMax;
  B[0] = xMax; B[1] = yMin; B[2] = zMax;
  C[0] = xMax; C[1] = yMax; C[2] = zMin;
  n2 = getNormalVector (A, B, C);
  p2 = A;

  // Face 3.
  A[0] = xMin; A[1] = yMax; A[2] = zMax;
  B[0] = xMax; B[1] = yMax; B[2] = zMax;
  C[0] = xMin; C[1] = yMax; C[2] = zMin;  
  n3 = getNormalVector (A, B, C);
  p3 = A;

  // Face 4.
  A[0] = xMin; A[1] = yMin; A[2] = zMax;
  B[0] = xMin; B[1] = yMax; B[2] = zMax;
  C[0] = xMin; C[1] = yMin; C[2] = zMin;
  n4 = getNormalVector (A, B, C);
  p4 = A;
  
  // Setup vector to test all GLL points -- face distance.
  std::vector<float> xTest;
  xTest.resize (3);
  
  // Loop over all points.
  for (size_t i=0; i<numGLLPoints; i++) {
    
    xTest[0] = xStore[i];
    xTest[1] = yStore[i];
    xTest[2] = zStore[i];
    
    // Take dot product with face plane to get distance. Expand to physical dimensions.
    float dFace1 = abs (projWonV_Dist (xTest, n1, p1) * R_EARTH);
    float dFace2 = abs (projWonV_Dist (xTest, n2, p2) * R_EARTH);
    float dFace3 = abs (projWonV_Dist (xTest, n3, p3) * R_EARTH);
    float dFace4 = abs (projWonV_Dist (xTest, n4, p4) * R_EARTH);    
    
    // Add point to sideSet.
    if (dFace1 < dHaze || dFace2 < dHaze || dFace3 < dHaze || dFace4 < dHaze)
      sideSet[i] = true;
  
  }
}

void kernel::exploreGaussianHaze () {
  
  // Uses the side sets found in getSideSets and copies to each processor some of the neighboring
  // chunk that might be needed in the gaussian smoother.
  
  // MPI variables.
  int myRank    = MPI::COMM_WORLD.Get_rank ();
  int worldSize = MPI::COMM_WORLD.Get_size ();
  
  // Buffer to transfer the side sets. TODO might also transfer ibool haze array, might not.
  float *xStoreBuffer = new float [numGLLPoints];
  float *yStoreBuffer = new float [numGLLPoints];
  float *zStoreBuffer = new float [numGLLPoints];
  
  singlePrint ("\x1b[33mExploring Gaussian Haze.\x1b[0m");
  
  // Loop over all processors. This could be avoided if some way was invented to tell which of those
  // were just too far away. Will have to see how fast this is.
  clock_t begin = std::clock ();
  for (size_t i=0; i<worldSize; i++) {
        
    if (i == myRank) {
      
      // Throw the current processor's coordStores to all others.
      xStoreBuffer = xStore;    
      yStoreBuffer = yStore;
      zStoreBuffer = zStore;
      
    }
                    
    // coordinate store broadcast.
    MPI::COMM_WORLD.Bcast (&xStoreBuffer[0], numGLLPoints, MPI_FLOAT, i);
    MPI::COMM_WORLD.Bcast (&yStoreBuffer[0], numGLLPoints, MPI_FLOAT, i);
    MPI::COMM_WORLD.Bcast (&zStoreBuffer[0], numGLLPoints, MPI_FLOAT, i);
    
    // Once these are broadcast, work through all interesting points and see if they should be
    // included in the gaussian haze. LOCAL ARRAY. TODO.
    for (size_t j=0; j<numGLLPoints; j++) {

      float xLoc = xStore[j];
      float yLoc = yStore[j];
      float zLoc = zStore[j];
            
      // Work through all interesting points of REMOTE ARRAY. TODO might only need to broadcast
      // the hazes in the first place. Like or with a bool array that only picks out of the hazes
      // of both chunks and compares the distances. Yes i think that should work.
      for (size_t k=0; k<numGLLPoints; k++) {
          
        float xRem = xStoreBuffer[k];
        float yRem = yStoreBuffer[k];
        float zRem = zStoreBuffer[k];

        float distance = oneDimDist (xLoc, yLoc, zLoc, xRem, yRem, zRem);
        
      }
            
    }

        
  }
  
  clock_t end = std::clock();
  double elapsed = double (end - begin) / CLOCKS_PER_SEC;  
  
  if (myRank == 0)
    std::cout << "\x1b[32mDone.\x1b[0m (" << elapsed << " seconds)\n";
    
}

void kernel::createKDTree () {
  
  // Create kdTree of the kernel.
  
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
    
    // Store the cartesian coordinates.
    xStore[i] = x;
    yStore[i] = y;
    zStore[i] = z;
    
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
    
    // Preallocate cartesian arrays.
    xStore = new float [numGLLPoints];
    yStore = new float [numGLLPoints];
    zStore = new float [numGLLPoints];
    
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