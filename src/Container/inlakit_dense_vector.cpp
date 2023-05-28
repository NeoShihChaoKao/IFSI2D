#ifndef SRC_CONTAINER_INLAKIT_DENSE_VECTOR_CPP_
#define SRC_CONTAINER_INLAKIT_DENSE_VECTOR_CPP_

#include <inlakit_container_dense_vector.h>
#include <vector_traits.h>
#include <inlakit_container_sparse_matrix.h>
#include <inlakit_container_sparse_matrix_csr.h>
#include <inlakit_container_sparse_matrix_coo.h>

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <vector>

/* =================================================================================== */

namespace inlakit   {
namespace container {

// #####################################################################################

template<typename T>
DenseVector<T>::DenseVector() {
 /*
  *   @Function : DenseVector
  *   @Author   : NeoKao
  *   @Brief    : Default constructor of DenseVector class
  */
  m_nnu = 0;
}

// #####################################################################################

template<typename T>
DenseVector<T>::DenseVector(const size_type &nnu) {
 /*
  *   @Function : DenseVector
  *   @Author   : NeoKao
  *   @Brief    : Constructor of DenseVector class
  */
  if ( nnu >= 0 ) m_nnu = nnu;
  m_val.resize(m_nnu);
}

// #####################################################################################

template<typename T>
DenseVector<T>::DenseVector(const self_type &DenVec) {
 /*
  *   @Function : DenseVector
  *   @Author   : NeoKao
  *   @Brief    : Copy constructor of DenseVector class
  */
  // check if their vector size are identical
  if ( m_nnu == DenVec.m_nnu ) {
    m_nnu = DenVec.m_nnu;
    m_val.resize(m_nnu);
    for ( size_type i = 0 ; i < m_nnu ; i++ ) m_val[i] = DenVec.m_val[i];
  }
}

// #####################################################################################

template<typename T>
DenseVector<T>::~DenseVector() {
 /*
  *   @Function : ~DenseVector
  *   @Author   : NeoKao
  *   @Brief    : Destructor
  */
  m_nnu = 0;
  m_val.clear();
  std::vector<T>().swap(m_val);  // free the memory
}

// #####################################################################################

template<typename T> inline
void DenseVector<T>::create(const size_type &nnu_) {
  if ( nnu_ > 0 ) m_nnu = nnu_;
  m_val.resize(m_nnu);
}

// #####################################################################################

template<typename T> inline
void DenseVector<T>::AxpBy(const value_type &alpha, const self_type &x,
                           const value_type &beta, const self_type &y) {
  for ( size_type i = 0 ; i < m_nnu ; i++ ) {
    m_val[i] = alpha * x.m_val[i] + beta * y.m_val[i];
  }
}

// #####################################################################################

template<typename T> inline
void DenseVector<T>::fill_with_zero() {
  if ( m_nnu > 0 ) {
    for ( size_type i = 0 ; i < m_nnu ; i++ ) m_val[i] = 0 ;
  }
}

// #####################################################################################

template<typename T> inline
void DenseVector<T>::fill_with_one() {
  if ( m_nnu > 0 ) {
    for ( size_type i = 0 ; i < m_nnu ; i++ ) m_val[i] = 1 ;
  }
}

// #####################################################################################

template<typename T> inline
void DenseVector<T>::fill_with_value(T& value) {
  if ( m_nnu > 0 ) {
    for ( size_type i = 0 ; i < m_nnu ; i++ ) m_val[i] = value ;
  }
}

// #####################################################################################

template<typename T> inline
void DenseVector<T>::copy(const self_type &DenseVector_) {
  if ( DenseVector_.m_nnu > 0 ) {
    m_nnu = DenseVector_.m_nnu;
    for ( size_type i = 0 ; i < m_nnu ; i++ ) {
      m_val[i] = DenseVector_.m_val[i];
    }
  }
}

// #####################################################################################

template<typename T> inline
void DenseVector<T>::swap(self_type &DenseVector_) {
 /*
  *   @Function : swap
  *   @Author   : NeoKao <neokao@moldex3d.com>
  *   @Brief    : Exchange the DesneVector object with DenseVector
  */
  if (m_nnu == DenseVector_.m_nnu) {
    m_val.swap(DenseVector_.m_val);
  } else  {
    std::cout << " Error happens at [swap] " << std::endl;
    std::cout << " The size is not consistent !! "      << std::endl;
  }
}

// #####################################################################################

template<typename T> inline
void DenseVector<T>::swap(const size_type &i, const size_type &j) {
 /*
  *   Function : swap
  *   Author   : NeoKao <neokao@moldex3d.com>
  *   Brief    : swap the value DenseVector[i] and DenseVector[j]
  */   
  T tmpValue;
  tmpValue = m_val[i];
  m_val[i] = m_val[j];
  m_val[j] = tmpValue;
}

// #####################################################################################

template<typename T> inline
void DenseVector<T>::clear() {
 /*
  *   @Function : clear
  *   @Author   : NeoKao <neokao@moldex3d.com>
  *   @Brief    : Clean the dense vector object
  */
  m_val.clear();
  std::vector<T> ().swap(m_val);  // clear and recycle the space
}

// #####################################################################################

template<typename T> inline
void DenseVector<T>::show() {
 /*
  *   @Function : show
  *   @Author   : NeoKao <neokao@moldex3d.com>
  *   @Brief    : Output the DenseVector profile
  */
  std::cout.precision(5);
  std::cout.setf(std::ios_base::fixed);
  std::cout << std::endl;
  for (size_type i = 0; i < m_nnu; i++) std::cout << "  " << m_val[i];
  std::cout << std::endl;
}

// #####################################################################################

template<typename T> inline
T* DenseVector<T>::getValues() {
  return m_val.data();
}

// #####################################################################################

// template<typename T> inline
// size_type DenseVector<T>::getSize() {
//   return m_nnu;
// }

// #####################################################################################

template<typename T> inline
T DenseVector<T>::calcInnerProduct() {
  T InnerProduct = 0;
  for ( size_type i = 0 ; i < m_nnu ; i++ ) {
    InnerProduct += m_val[i] * m_val[i];
  }
  return InnerProduct;
}

// #####################################################################################

template<typename T> inline
T DenseVector<T>::calcInnerProduct(const self_type &DenseVector_) {
  T InnerProduct = 0;
  if ( m_nnu == DenseVector_.m_nnu ) {
    for ( size_type i = 0 ; i < m_nnu ; i++ ) {
      InnerProduct += m_val[i] * DenseVector_.m_val[i];
    }
  }
  return InnerProduct;
}

// #####################################################################################

template<typename T> inline
T DenseVector<T>::calcNorm1() {
  T Norm1 = 0;
  for ( size_type i = 0 ; i < m_nnu ; i++ ) {
    Norm1 += m_val[i];
  }
  Norm1 /= (m_nnu * 1.0);
  return Norm1;
}

// #####################################################################################

template<typename T> inline
T DenseVector<T>::calcNorm2() {
  T Norm2 = 0;
  for ( size_type i = 0 ; i < m_nnu ; i++ ) {
    Norm2 += m_val[i] * m_val[i];
  }
  Norm2 /= (m_nnu * 1.0);
  Norm2 = std::sqrt(Norm2);
  return Norm2;
}

// #####################################################################################

template<typename T> inline
T DenseVector<T>::calcNormInf() {
  T max = *std::max_element(m_val.begin(), m_val.end());
  return max;
}

// #####################################################################################

template<typename T> inline
T DenseVector<T>::getMaxValue() {
  T max = *std::max_element(m_val.begin(), m_val.end());
  return max;
}

// #####################################################################################

template<typename T> inline
T DenseVector<T>::getMinValue() {
  T min = *std::min_element(m_val.begin(), m_val.end());
  return min;
}

// #####################################################################################

template<typename T> inline
T DenseVector<T>::calcMeanValue() {
  value_type MeanValue = 0.0;
  if ( m_nnu > 0 ) {
    for (size_type i = 0 ; i < m_nnu ; i++) {
      MeanValue += m_val[i];
    }
  }
  MeanValue /= m_nnu;
  return MeanValue;
}

// #####################################################################################

template<typename T> inline
void DenseVector<T>::resetSize(const size_type &nnu) {
  if (m_nnu > nnu || m_nnu < nnu) {
    m_nnu = nnu;
    m_val.clear();
    std::vector<T> ().swap(m_val);  // clear and recycle the space
    m_val.resize(m_nnu);
  }
}

// #####################################################################################

template<typename T> inline
void DenseVector<T>::setValue(const size_type &index, const T &value) {
  if ( index >= 0 ) m_val[index] = value;
}

// #####################################################################################

// template<typename T> inline
// T DenseVector<T>::getValue(const size_type &index) {
//   return m_val.data()[index];
// }

// #####################################################################################

template<typename T> inline
T DenseVector<T>::calcSum() {
  T Sum = 0;
  if ( m_nnu > 0 ) {
    for ( size_type i = 0 ; i < m_nnu ; i++ ) {
      Sum += m_val[i];
    }
  }
  return Sum;
}

// #####################################################################################

template<typename T> inline
T DenseVector<T>::calcAbsoluteSum() {
  T AbsoluteSum = 0;
  if ( m_nnu > 0 ) {
    for ( size_type i = 0 ; i < m_nnu ; i++ ) {
      AbsoluteSum += std::abs(m_val[i]);
    }
  }
  return AbsoluteSum;
}

// #####################################################################################

template<typename T> inline
void DenseVector<T>::addValue(const T &value) {
if (m_capacity < m_val.capacity()) {
    m_val[m_capacity] = value;
    m_capacity++;
  } else {
    m_val.push_back(value);
    m_nnu = m_val.size();
  }
}

// #####################################################################################

template<typename T> inline
DenseVector<T> DenseVector<T>::operator + (const DenseVector<T> &DenVec) {
 /*
  *   @Function : operator +
  *   @Author   : NeoKao <neokao@moldex3d.com>
  *   @Brief    : Overload the operator +
  */
  inlakit::container::DenseVector<T> tmp_vec;
  if ( m_nnu != DenVec.m_nnu ) {
  } else {
    //------------------------------------
    tmp_vec.create(m_nnu);
    //------------------------------------
    for ( size_type i = 0 ; i < m_nnu ; i++ ) {
      tmp_vec.m_val[i] = m_val[i] + DenVec.m_val[i];
    }
  }
  return tmp_vec;
}

// #####################################################################################

template<typename T> inline
DenseVector<T> DenseVector<T>::operator - (const DenseVector<T> &DenVec) {
 /* 
  *   @Function : operator -
  *   @Author   : NeoKao <neokao@moldex3d.com>
  *   @Brief    : Overload the operator -
  */
  inlakit::container::DenseVector<T> tmp_vec;
  if ( m_nnu != DenVec.m_nnu ) {
  } else {
    //------------------------------------
    tmp_vec.create(m_nnu);
    //------------------------------------
    for ( size_type i = 0 ; i < m_nnu ; i++ ) {
      tmp_vec.m_val[i] = m_val[i] - DenVec.m_val[i];
    }
  }
  return tmp_vec;
}

// #####################################################################################

template<typename T> inline
DenseVector<T> & DenseVector<T>::operator = (const DenseVector<T> &DenseVector_) {
 /*
  *   @Function : operator =
  *   @Author   : NeoKao <neokao@moldex3d.com>
  *   @Brief    : Overload the operator =
  */   
  if ( this != &DenseVector_ ) {
    this -> m_nnu = DenseVector_.m_nnu;
      m_val.resize(m_nnu);
      for ( size_type i = 0 ; i < m_nnu ; i++ ) m_val[i] = DenseVector_.m_val[i];
  }
  return (*this);
}

// #####################################################################################

template<typename T>
DenseVector<T> &DenseVector<T>::operator *= (const T scaling) {
  for (size_type i = 0; i < m_nnu; i++) {
    m_val[i] *= scaling;
  }
  return *this;
}

// #####################################################################################

// template<typename T>
// const DenseVector<T> &operator * (const T &scale_factor, const DenseVector<T> &DenVec) {
//  /*
//   *   @Function : oeprator *
//   *   @Author   : NeoKao <neokao@moldex3d.com>
//   *   @Brief    : Overload the operator *
//   */
//   return DenseVector<T>(DenVec) *= scale_factor;
// }

// #####################################################################################

// template<typename T>
// const DenseVector<T> &operator * (const DenseVector<T> &DenVec, const T &scale_factor) {
//  /*
//   *   @Function : oeprator *
//   *   @Author   : NeoKao <neokao@moldex3d.com>
//   *   @Brief    : Overload the operator *
//   */
//   return DenseVector<T>(DenVec) *= scale_factor;
// }

// #####################################################################################

template<typename T> inline
T DenseVector<T>::operator()(size_type index) {
 /*
  *   @Function : operator ()
  *   @Author   : NeoKao <neokao@moldex3d.com>
  *   @Brief    : Overload the operator ()
  */ 
  return m_val[index];
}

// #####################################################################################

template<typename T> inline
void DenseVector<T>::BubbleSorting() {
 /*
  *   @Function : BubbleSorting
  *   @Author   : NeoKao <neokao@moldex3d.com>
  *   @Brief    : Sort the dense vector using bubble sorting algorithm
  */ 
  T tmp_value = 0.0;
  // ----------------------------------------------------
  for (size_type stage = 0; stage < m_nnu-1; stage++) {
    for (size_type j = 0; j < m_nnu - stage-1; j++) {
      if (m_val[j] > m_val[j+1]) {
        tmp_value  = m_val[j];
        m_val[j]   = m_val[j+1];
        m_val[j+1] = tmp_value;
      }
    }  // end j-loop
  }  // end stage-loop
  // ----------------------------------------------------
}

// #####################################################################################

template<typename T> inline
void DenseVector<T>::BubbleSorting_d() {
 /*
  *   @Function : BubbleSorting
  *   @Author   : NeoKao <neokao@moldex3d.com>
  *   @Brief    : Sort the dense vector using bubble sorting algorithm
  */ 
  T tmp_value = 0.0;
  // ----------------------------------------------------
  for (size_type stage = 0; stage < m_nnu - 1; stage++) {
    for (size_type j = 0; j < m_nnu - stage - 1; j++) {
      if (m_val[j] < m_val[j+1]) {
        tmp_value  = m_val[j];
        m_val[j]   = m_val[j+1];
        m_val[j+1] = tmp_value;
      }
    }  // end j-loop
  }  // end stage-loop
  // ----------------------------------------------------
}

// #####################################################################################

template<typename T> inline
void DenseVector<T>::SelectionSorting() {
 /*
  *   @Function : SelectionSorting()
  *   @Author   : NeoKao <neokao@moldex3d.com>
  *   @Brief    : Sorting the dense vector using selection sorting algorithm
  */ 
  // to be added ....
}

// #####################################################################################

template<typename T> inline
void DenseVector<T>::InsertionSorting() {
 /*
  *   @Function : InsertionSorting()
  *   @Author   : NeoKao <neokao@moldex3d.com>
  *   @Brief    : Sorting the dense vector using insertion sorting algorithm
  */ 
  // to be added ....
}

// #####################################################################################

template<typename T> inline
void DenseVector<T>::QuickSorting() {
 /*
  *   @Function : operator ()
  *   @Author   : NeoKao <neokao@moldex3d.com>
  *   @Brief    : Sorting the dense vector using quick sorting algorithm
  */ 
}

// ---------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------
// Instance (To export the dll library)
# pragma warning(disable:4661) /* avoid the 4661 warning message */
template class DenseVector<double>;
// ---------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------

}  // namespace container
}  // namespace inlakit

#endif  // SRC_CONTAINER_INLAKIT_DENSE_VECTOR_CPP_