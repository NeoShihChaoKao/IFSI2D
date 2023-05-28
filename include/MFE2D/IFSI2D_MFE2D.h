#ifndef INCLUDE_IFSI2D_INPUT_IFSI2D_INPUT_H_
#define INCLUDE_IFSI2D_INPUT_IFSI2D_INPUT_H_

#include <IFSI2D_ConstVariables.h>
#include <inlakit.h>

#include <iostream>
#include <iomanip>
#include <vector>


namespace ifsi2d {
namespace mfe2d  {    

//===============================================

template<typename T>
class MFE2D {

  //----------------------------------------

  template<typename S> friend class DenseVector;
  template<typename S> friend class DenseMatrix;
   
  //----------------------------------------

  using size_type  = int;
  using value_type = double;

  //----------------------------------------
  
  public:
    MFE2D(); // default constructor
    virtual ~MFE2D(); // destructor
  private:
    
    size_type               melem; // number of total elements of domain
    size_type               mpoin; // number of total nodes    of domain
    size_type               mtotv; // number of total variables 
    
    std::vector<size_type>  lboud; // label of boundary condition
    std::vector<value_type> boudv; // value of boundary condition
    
}

//===============================================

}  // namespace mfe2d
}  // namespace ifsi2d

#endif  // INCLUDE_IFSI2D_INPUT_H_
