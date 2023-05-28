#ifndef IFSI2D_VARIABLES_H_
#define IFSI2D_VARIABLES_H_

#include <iostream>
#include <iomanip>

/* ------------------------------------ */
using size_type  = int;
using value_type = double;
/* ------------------------------------ */

const size_type nnodl = 4 ; /* number of bi-linear interpolation functions    */
const size_type nnodp = 9 ; /* number of bi-quadratic interpolation functions */
const size_type ngaus = 3 ; /* number of Gaussian integration rule            */ 
const size_type nevab = 27; /* the size of each local fluid matrix (u,v,p)    */

#endif  // IFSI2D_VARIABLES_H_