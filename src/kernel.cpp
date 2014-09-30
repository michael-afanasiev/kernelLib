#include "classes.hpp"
#include <netcdf>

using namespace std;

kernel::kernel (std::string fName) {
  
  // Initialize and split kernels.
  
  fileName = fName;
  openCoordNetcdf      ();
  openKernelNetcdf     ();
  
  findChunkDimensions ();
  findNeighbours  ();
  constructMaster ();
  
  singlePrint ("Rotating chunks to z axis.");
  findChunkDimensions ();
  rotateZaxis         ();
  findChunkDimensions ();
  rotateYaxis         ();
  findChunkDimensions ();
    
  createKDTree ();
  
  MPI::COMM_WORLD.Barrier ();
  
}

kernel::~kernel () {
  
  delete [] rawKernel;
  delete [] radius;
  delete [] theta;
  delete [] phi;
  
}


void kernel::constructMaster () {
  
  // This function assembles a master chunk from the surrounding ones, and resets the kernel 
  // definitions.
  
  singlePrint ("Constructing master chunk.");
  
  // MPI variables.
  int myRank     = MPI::COMM_WORLD.Get_rank ();
  int RAD_TAG    = 0;
  int THETA_TAG  = 1;
  int PHI_TAG    = 2;
  int KERNEL_TAG = 3;
  
  // New number of GLL points is the number of one chunk times the number of neighbouring chunks.
  int newNumGLLPoints = numGLLPoints * (neighbours.size()+1);
  
  // These are scratch arrays that will be copied to.
  float *scratchRadius    = new float [newNumGLLPoints];
  float *scratchTheta     = new float [newNumGLLPoints];
  float *scratchPhi       = new float [newNumGLLPoints];
  float *scratchRawKernel = new float [newNumGLLPoints];
      
  // Initialize region tag.
  bool *scratchInReg = new bool [newNumGLLPoints]();
  
  // These are the MPI recv buffer arrays that will be recieved from each of the nieghbours.
  float *recvBufRadius    = new float [numGLLPoints];
  float *recvBufTheta     = new float [numGLLPoints];
  float *recvBufPhi       = new float [numGLLPoints];
  float *recvBufRawKernel = new float [numGLLPoints];  
    
  // Send the neighbouring arrays. Non-blocking -- LARGE BUFFERS. Send to the processor stored in
  // neighbours[i]. Tag the file with the current rank.
  for (size_t i=0; i<neighbours.size(); i++) {
    
    MPI::COMM_WORLD.Isend (&radius[0],    numGLLPoints, MPI::FLOAT, neighbours[i], RAD_TAG);
    MPI::COMM_WORLD.Isend (&theta[0],     numGLLPoints, MPI::FLOAT, neighbours[i], THETA_TAG);
    MPI::COMM_WORLD.Isend (&phi[0],       numGLLPoints, MPI::FLOAT, neighbours[i], PHI_TAG);
    MPI::COMM_WORLD.Isend (&rawKernel[0], numGLLPoints, MPI::FLOAT, neighbours[i], KERNEL_TAG);    
  
  // Recieve the neighbouring arrays from all interesting processes.
    MPI::COMM_WORLD.Recv (recvBufRadius,    numGLLPoints, MPI::FLOAT, neighbours[i], RAD_TAG);
    MPI::COMM_WORLD.Recv (recvBufTheta,     numGLLPoints, MPI::FLOAT, neighbours[i], THETA_TAG);
    MPI::COMM_WORLD.Recv (recvBufPhi,       numGLLPoints, MPI::FLOAT, neighbours[i], PHI_TAG);
    MPI::COMM_WORLD.Recv (recvBufRawKernel, numGLLPoints, MPI::FLOAT, neighbours[i], KERNEL_TAG);
    
    // Here we fill up the scratch arrays. k loops over the receive buffer, and copies the receive
    // buffer into the scratch array, which is dimensioned by its position in num_neighbours.    
    size_t k=0;
    for (size_t j=(i)*numGLLPoints; j<(i+1)*numGLLPoints; j++) {
      
      scratchRadius[j]    = recvBufRadius[k];
      scratchTheta[j]     = recvBufTheta[k];
      scratchPhi[j]       = recvBufPhi[k];
      scratchRawKernel[j] = recvBufRawKernel[k];
      k++;
      
    }
                
  }
  
  // Still need to copy one last value -- the original values belonging to that processor.
  size_t k=0;
  for (size_t i=neighbours.size()*numGLLPoints; i<(neighbours.size()+1)*numGLLPoints; i++) {    
    
    scratchRadius[i]    = radius[k];
    scratchTheta[i]     = theta[k];
    scratchPhi[i]       = phi[k];
    scratchRawKernel[i] = rawKernel[k];
    scratchInReg[i]     = true;
    k++;
        
  }
  
  // Phew. Free a whole bunch of now-useless memory.
  delete [] phi;
  delete [] theta;
  delete [] radius;
  delete [] xStore;
  delete [] yStore;
  delete [] zStore;
  delete [] rawKernel;
  delete [] recvBufPhi;
  delete [] recvBufTheta;
  delete [] recvBufRadius;
  delete [] recvBufRawKernel;
  
  // re-allocate the arrays that we're replacing with the new, expanded values.
  phi       = new float [newNumGLLPoints];
  inReg     = new bool  [newNumGLLPoints];
  theta     = new float [newNumGLLPoints];
  radius    = new float [newNumGLLPoints];
  xStore    = new float [newNumGLLPoints];
  yStore    = new float [newNumGLLPoints];
  zStore    = new float [newNumGLLPoints];
  rawKernel = new float [newNumGLLPoints];
  
  // Copy the scratch arrays into the original default arrays.
  for (size_t i=0; i<newNumGLLPoints; i++) {
    
    radius[i]    = scratchRadius[i];
    theta[i]     = scratchTheta[i];
    phi[i]       = scratchPhi[i];
    rawKernel[i] = scratchRawKernel[i];   
    inReg[i]     = scratchInReg[i]; 
    
  }
  
  // Reset the number of GLL points. TODO see if anything else needs resetting.
  numGLLPoints = newNumGLLPoints;
  
  // Free the scratch memory.
  delete [] scratchPhi;
  delete [] scratchInReg;
  delete [] scratchTheta;
  delete [] scratchRadius;
  delete [] scratchRawKernel;
      
}

void kernel::resetRotations () {
    
  // Resets coordinates to their original values. This is to save floating point errors when
  // rotating like a madman.
  
  // Copy original arrays.
  for (size_t i=0; i<numGLLPoints; i++) {
    
    radius[i] = radiusOrig[i];
    theta[i]  = thetaOrig[i];
    phi[i]    = phiOrig[i];
    
  }
  
  findChunkDimensions ();
    
}

void kernel::findNeighbours () {
  
  // MPI variables.
  int myRank    = MPI::COMM_WORLD.Get_rank ();
  int worldSize = MPI::COMM_WORLD.Get_size ();
  
  // These vectors are used to determine neighboring chunks by distances to center.
  std::vector<float> centerDistances;
  std::vector<float> centerDistancesIndex;
    
  singlePrint ("Finding neighbours.");  
  float xCenterBuf;
  float yCenterBuf;
  float zCenterBuf;
  
  // Loop over all processors. This could be avoided if some way was invented to tell which of those
  // were just too far away. Will have to see how fast this is.
  for (size_t i=0; i<worldSize; i++) {
        
    if (i == myRank) {
      
      // Throw the current processor's centers to all others.
      xCenterBuf = xCenter;
      yCenterBuf = yCenter;
      zCenterBuf = zCenter;
      
    }
    
    // Broadcast centers.
    MPI::COMM_WORLD.Bcast (&xCenterBuf, 1, MPI::FLOAT, i);
    MPI::COMM_WORLD.Bcast (&yCenterBuf, 1, MPI::FLOAT, i);
    MPI::COMM_WORLD.Bcast (&zCenterBuf, 1, MPI::FLOAT, i);
    
    // Determine distance.
    float distFromCenter = distFromPoint (xCenter, yCenter, zCenter, xCenterBuf, yCenterBuf, 
      zCenterBuf) * R_EARTH;
      
    // Save distances. There are two arrays as one will be used an an index, and the other will 
    // be sorted.
    centerDistances.push_back      (distFromCenter);
    centerDistancesIndex.push_back (distFromCenter);      
                              
  }  
  
  // Sort the distance array.
  std::sort (centerDistances.begin(), centerDistances.end());
  
  // Loop through sorted distance array, taking n closest centers as the neighbouring chunks
  // (n=8 usually).
  for (size_t i=0; i<centerDistances.size(); i++) {    
    for (size_t j=0; j<centerDistancesIndex.size(); j++) {
  
      if ((centerDistances[i] == centerDistancesIndex[j]) && (centerDistances[i] != 0))
        neighbours.push_back (j);      
  
    }
  }
  
  // Sort by processor number (helps future broadcast).
  std::sort (neighbours.begin(), neighbours.end());
    
}

void kernel::createKDTree () {
  
  // Create kdTree of the kernel.
  
  int myRank = MPI::COMM_WORLD.Get_rank ();
    
  singlePrint ("\n\x1b[33mCreating KD-tree. ");
  double begin = std::clock ();  
  
  // Initialize tree.
  tree = kd_create (3);
  
  // Initialize the data array.
  KDdat = new int [numGLLPoints];  
  
  // Local xyz values. May not be appropriate.
  float x, y, z;
  for (size_t i=0; i<numGLLPoints; i++) {

    x = xStore[i];
    y = yStore[i];
    z = zStore[i];
    // radThetaPhi2xyz (radius[i], theta[i], phi[i], x, y, z);
    KDdat[i] = i;
    kd_insert3 (tree, x, y, z, &KDdat[i]);
    
  }

  MPI::COMM_WORLD.Barrier ();      
  clock_t end = std::clock();
  double elapsed = double (end - begin) / CLOCKS_PER_SEC;  
  if (myRank == 0)
    std::cout << "\x1b[32mDone.\x1b[0m (" << elapsed << " seconds)\n";
  
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
    
    if (z < zMin)
      zMin = z;
    
    // Sum radius for average rad.
    rSum += radius[i];
    
    // Sum components.
    xSum += x;
    ySum += y;
    zSum += z;
    
    
  }
    
  // Calculate magnitude and normalize.
  xCenterPhys  = xSum / numGLLPoints;
  yCenterPhys  = ySum / numGLLPoints;
  zCenterPhys  = zSum / numGLLPoints;
    
  float magnitude = sqrt (xCenterPhys*xCenterPhys + yCenterPhys*yCenterPhys + 
    zCenterPhys*zCenterPhys);
  
  xCenter         = xCenterPhys / magnitude;
  yCenter         = yCenterPhys / magnitude;
  zCenter         = zCenterPhys / magnitude;
  
      
  // Return center point.
  xyz2RadThetaPhi (radCenter, thetaCenter, phiCenter, xCenter, yCenter, zCenter);
  
  // Get average radius.
  radCenter = rSum / numGLLPoints;

  // Get extreme spherical coordinates (after possible rotation).
  radiusMin     = *std::min_element (radius, radius+numGLLPoints);
  thetaMin      = *std::min_element (theta, theta+numGLLPoints);
  phiMin        = *std::min_element (phi, phi+numGLLPoints);

  radiusMax     = *std::max_element (radius, radius+numGLLPoints);
  thetaMax      = *std::max_element (theta, theta+numGLLPoints);
  phiMax        = *std::max_element (phi, phi+numGLLPoints);
  
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
      std::cout << mgn << "Number of solver processers:\t " << numWroteProcs
        << "\nNumber of GLL points:\t\t " << numGLLPoints << rst << "\n" << std::endl;
    
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
    
    // Save original arrays.
    radiusOrig = new float [numGLLPoints];
    thetaOrig  = new float [numGLLPoints];
    phiOrig    = new float [numGLLPoints];    
    for (size_t i=0; i<numGLLPoints; i++) {
      
      radiusOrig[i] = radius[i];
      thetaOrig[i]  = theta[i];
      phiOrig[i]    = phi[i];
      
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
      std::cout << mgn << "Number of solver processers:\t " << numWroteProcs
        << "\nNumber of GLL points:\t\t " << numGLLPoints << rst << "\n" << std::endl;
    
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

void kernel::quickSortCenter (int i1st, int i2nd) {
  
  int pivotElement;
  
  float d1st = distFromPoint (xStore[i1st], yStore[i1st], zStore[i1st], xCenterPhys, yCenterPhys, 
    zCenterPhys);
    
  float d2nd = distFromPoint (xStore[i2nd], yStore[i2nd], zStore[i2nd], xCenterPhys, yCenterPhys, 
    zCenterPhys);
  
  if (i1st < i2nd) {
    pivotElement = pivot (i1st, i2nd, d1st, d2nd);
    quickSortCenter (i1st, pivotElement-1);
    quickSortCenter (pivotElement+1, i2nd);
  }
  
}

int kernel::pivot (int &i1st, int &i2nd, float &d1st, float &d2nd) {
              
  int p              = i1st;
  float pivotElement = d1st;
  
  for (int i=i1st+1; i<=i2nd; i++) {
    
    float dTest = distFromPoint (xStore[i], yStore[i], zStore[i], 
      xCenterPhys, yCenterPhys, zCenterPhys);
      
    if (dTest <= pivotElement) {
      
      p++;
      std::swap (inReg[i],     inReg[p]);
      std::swap (xStore[i],    xStore[p]);
      std::swap (yStore[i],    yStore[p]);
      std::swap (zStore[i],    zStore[p]);
      std::swap (rawKernel[i], rawKernel[p]);
      
    }
  }
  
  std::swap (inReg[p],     inReg[i1st]);
  std::swap (xStore[p],    xStore[i1st]);
  std::swap (yStore[p],    yStore[i1st]);
  std::swap (zStore[p],    zStore[i1st]);
  std::swap (rawKernel[p], rawKernel[i1st]);
  
  return p;
  
}