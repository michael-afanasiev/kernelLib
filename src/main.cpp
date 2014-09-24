#include "classes.hpp"

int main () {
  
  MPI::Init ();
  
  if (MPI::COMM_WORLD.Get_rank () == 0)
    std::cout << rst << std::endl;
  
  kernel kern ("./krn/betahKernelCrustMantle.nc");
  
  MPI::Finalize ();

}