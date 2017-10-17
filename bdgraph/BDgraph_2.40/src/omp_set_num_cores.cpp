#include <sstream>
#include <string>        // std::string, std::to_string
#include <vector>        // for using vector
#include <limits>        // for std::numeric_limits<double>::max()
#include <Rmath.h>
#include <Rconfig.h>
#include <omp.h>
#include <R.h>

extern "C" {
void omp_set_num_cores( int *cores ) 
{
   // we removed the portion which tests for the 
   // SUPPORT_OPENMP definition, since this seemed not to work
   // for us despite specifiying -fopenm for the R compiler options
   omp_set_num_threads( *cores );
}
}
