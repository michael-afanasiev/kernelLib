#include "classes.hpp"
#include "exodusII.h"

void exodusErrorCheck (int, std::string);

void writeExodus (float *regMeshArr, float *regX, float *regY, float *regZ, 
                  int &nx, int &ny, int &nz, std::string fName)
{

  // Exodus parameters.
  int  comp_ws      = sizeof(float);
  int  io_ws        = 0;
  int  nVars        = 1;
  int  nDim         = 3;
  int  nBlock       = 1;
  int  nNodeSet     = 0;
  int  nSideSet     = 0;
  int  nNodePerElem = 8;
  int  numNodes     = nx*ny*nz;
  int  numElem      = (nx-1)*(ny-1)*(nz-1);

  char *varNames[1];
  varNames[0] = (char*) "Sensitivity";
  
  int *nodeNumArr = new int [numNodes];
  float *nodeCorZ = new float [numNodes];
  float *nodeCorY = new float [numNodes];
  float *nodeCorX = new float [numNodes];

  // Unpack coordinate arrays.
  int it = 0;
  for ( int i=0; i<nx; i++ ) {
    for ( int j=0; j<ny; j++ ) {
      for ( int k=0; k<nz; k++ ) {
    
        nodeNumArr[it] = it+1;
        nodeCorZ[it]   = regZ[k];
        nodeCorY[it]   = regY[j];
        nodeCorX[it]   = regX[i];
        it++;

      }
    }
  }

  int *connect = new int [numElem*nNodePerElem];

  // Construct connectivity array.
  int count=0;
  for ( int i=0; i<nx-1; i++ ) {
    for ( int j=0; j<ny-1; j++ ) {
      for ( int k=0; k<nz-1; k++ ) {

        connect[count]   = nodeNumArr[k+(nz)*(j+i*ny)];
        connect[count+1] = nodeNumArr[k+(nz)*(j+i*ny)+nz];
        connect[count+2] = nodeNumArr[k+(nz)*(j+i*ny)+ny*nz+nz];
        connect[count+3] = nodeNumArr[k+(nz)*(j+i*ny)+ny*nz];
        connect[count+4] = nodeNumArr[(k+1)+(nz)*(j+i*ny)];
        connect[count+5] = nodeNumArr[(k+1)+(nz)*(j+i*ny)+nz];
        connect[count+6] = nodeNumArr[(k+1)+(nz)*(j+i*ny)+ny*nz+nz];
        connect[count+7] = nodeNumArr[(k+1)+(nz)*(j+i*ny)+ny*nz];
        count           += 8;

      }
    }
  }

  // std::cout << "Writing exodus file." << std::flush << std::endl;

  // Interpolated array.
  int idexo = ex_create ( fName.c_str(), EX_CLOBBER, &comp_ws, &io_ws );
  exodusErrorCheck ( ex_put_init( idexo, "Kernel", nDim, numNodes, numElem, nBlock, nNodeSet, 
    nSideSet ), "ex_put_init" );
  exodusErrorCheck ( ex_put_coord ( idexo, nodeCorX, nodeCorY, nodeCorZ ), "ex_put_coord" );
  exodusErrorCheck ( ex_put_elem_block ( idexo, 1, "HEX", numElem, nNodePerElem, 0 ), 
    "ex_put_elem_block" );
  exodusErrorCheck ( ex_put_node_num_map ( idexo, nodeNumArr ), "ex_put_node_num_map" );
  exodusErrorCheck ( ex_put_elem_conn  ( idexo, 1, connect ), "ex_put_elem_conn" );
  exodusErrorCheck ( ex_put_var_param ( idexo, "n", nVars ), "ex_put_var_param" );
  exodusErrorCheck ( ex_put_var_names ( idexo, "n", nVars, varNames ), "ex_put_var_names" );
  exodusErrorCheck ( ex_put_nodal_var ( idexo, 1, 1, numNodes, regMeshArr ), "ex_put_nodal_var" ); 
  exodusErrorCheck ( ex_close ( idexo ), "ex_close" );

}

void exodusErrorCheck ( int ier, std::string function )
{

  if ( ier != 0 )
  {
    std::cout << "ERROR in " << red << function << rst << std::flush << std::endl;
    exit (EXIT_FAILURE);
  }

}