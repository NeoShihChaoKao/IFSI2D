#ifndef SRC_CONTAINER_INLAKIT_DENSE_MATRIX_CPP_
#define SRC_CONTAINER_INLAKIT_DENSE_MATRIX_CPP_

/* --------------- header files ---------------- */
#include <inlakit_container_dense_vector.h>
#include <inlakit_container_dense_matrix.h>
#include <vector_traits.h>
#include <matrix_traits.h>
/* --------------------------------------------- */

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <vector>

/* =================================================================================== */

namespace inlakit   {
namespace container {

// #####################################################################################

template<typename T>
DenseMatrix<T>::DenseMatrix() {
 /*
  *   @Function : DenseMatrix
  *   @Author   : NeoKao <neokao@moldex3d.com>
  *   @Brief    : Default constructor of DenseMatrix class
  */ 
  m_rows = 0;
  m_cols = 0;
}

// #####################################################################################

template<typename T>
DenseMatrix<T>::DenseMatrix(const size_type &rowSize, const size_type &colSize) {
 /*
  *   @Function : DenseMatrix
  *   @Author   : NeoKao <neokao@moldex3d.com>
  *   @Brief    : Constructor of DenseMatrix class
  */ 
  if (rowSize > 0 && colSize > 0) {
    //-----------------------------------------------------------
    this -> m_rows = rowSize;
    this -> m_cols = colSize;
    //-----------------------------------------------------------
    m_val.resize(m_rows);
    typename std::vector< std::vector<T> >::iterator iter;
    for (iter = m_val.begin(); iter != m_val.end(); ++iter) {
      iter -> resize(m_cols);
    }  // end iter-loop
    //-----------------------------------------------------------
  }
  //-----------------------------------------
  this -> fill_with_zero();
  //-----------------------------------------
}

// #####################################################################################

template<typename T>
DenseMatrix<T>::DenseMatrix(const self_type &DenMat) {
 /* 
  *   @Function : DenseMatrix
  *   @Author   : NeoKao <neokao@moldex3d.com>
  *   @Brief    : Copy constructor of DenseMatrix class
  */ 
  if (DenMat.m_rows > 0 && DenMat.m_cols > 0) {
    this -> m_rows = DenMat.m_rows;
    this -> m_cols = DenMat.m_cols;
    this -> create(m_rows, m_cols);
  }
  //-----------------------------------------
  for (size_type i = 0; i < m_rows; i++) {
    for (size_type j = 0; j < m_cols; j++) {
      m_val[i][j] = DenMat.m_val[i][j];
    }  // end j-loop
  }  // end i-loop
  //-----------------------------------------
}

// #####################################################################################

template<typename T>
DenseMatrix<T>::~DenseMatrix() {
 /* 
  *   @Function : DenseMatrix
  *   @Author   : NeoKao <neokao@moldex3d.com>
  *   @Brief    : Destructor of DenseMatrix class
  */ 
  m_rows = 0;
  m_cols = 0;
}

// #####################################################################################

template<typename T> inline
bool DenseMatrix<T>::checkIndex(const size_type &rowIndex,
                                const size_type &colIndex) {
 /*
  *   @Function : checkIndex
  *   @Author   : NeoKao <neokao@moldex3d.com>
  *   @Brief    : check the (row,col) index
  */  
  if (rowIndex >= 0 && colIndex >= 0) {
    return true;
  } else {
    std::cout << " Error happens in function [checkDimension] " << std::endl;
    std::cout << " check the parameter rowIndex or colIndex "   << std::endl;
    exit(1);
  }
}

// #####################################################################################

template<typename T> inline
void DenseMatrix<T>::create(const size_type& rowSize,
                            const size_type& colSize) {
 /*
  *   @Function : create
  *   @Author   : NeoKao <neokao@moldex3d.com>
  *   @Brief    : Generate the matrix object in Matrix class
  */  
  if (rowSize > 0 && colSize > 0) {
    this -> m_rows = rowSize;
    this -> m_cols = colSize;
    //-----------------------------------------------------------
    m_val.resize(m_rows);
    typename std::vector< std::vector<T> >::iterator iter;
    for (iter = m_val.begin(); iter != m_val.end(); ++iter) {
      iter -> resize(m_cols);
    }
  }
  // ------------------------------------------------------------
  //   for (size_type i = 0; i < m_rows; i++) {
  //     for (size_type j = 0; j < m_cols; j++) {
  //       m_val[i][j] = 0;
  //     }  // end j-loop
  //   }  // end i-loop
  this -> fill_with_zero();
  // ------------------------------------------------------------
}

// #####################################################################################

template<typename T> inline
void DenseMatrix<T>::setValue(const size_type  &rowIndex,
                              const size_type  &colIndex,
                              const value_type &value) {
 /* 
  *   Function : setValue
  *   Author   : NeoKao <neokao@moldex3d.com>
  *   Brief    : Set the DenseMatrix data at the specific location (rowIndex,colIndex)
  */  
  if (checkIndex(rowIndex , colIndex)) {
    m_val[rowIndex][colIndex] = value;
  }
}

// #####################################################################################

template<typename T> inline
void DenseMatrix<T>::copy(const self_type &DenMat) {
 /*   
  *   Function : copy
  *   Author   : NeoKao <neokao@moldex3d.com>
  *   Brief    : Copy the DenseMatrix object
  */  
  if (m_rows != DenMat.m_rows && m_cols != DenMat.m_cols) {
    std::cout << " Error happens in MatrixDataCopy" << std::endl;
    std::cout << " The matrix size is not matched"  << std::endl;
    exit(1);
  } else {
    // ------------------------------------------------------------
    for (size_type i = 0; i < m_rows; i++) {
      for (size_type j = 0; j < m_cols; j++) {
        m_val[i][j] = DenMat.m_val[i][j];
      }  // end j-loop
    }  // end i-loop
    // ------------------------------------------------------------
  }
}

// #####################################################################################

template<typename T> inline
void DenseMatrix<T>::transpose() {
  T **tmpDenseMatrix = new T*[m_rows];
  //------------------------------------------
  for ( size_type i = 0 ; i < m_rows ; i++ ) {
    tmpDenseMatrix[i] = new T[m_cols];
  }
  //------------------------------------------
  for ( size_type i = 0 ; i < m_rows ; i++ ) {
    for ( size_type j = 0 ; j < m_cols ; j++ ) {
      tmpDenseMatrix[i][j] = m_val[i][j];
    }
  }
  //------------------------------------------
  for ( size_type i = 0 ; i < m_rows ; i++ ) {
    for ( size_type j = 0 ; j < m_cols ; j++ ) {
      m_val[i][j] = tmpDenseMatrix[j][i];
    }
  }
  //----- free the matrix --------------------
  for (size_type i =  0 ; i < m_rows ; i++) delete[] tmpDenseMatrix[i];
  delete[] tmpDenseMatrix;
}

// #####################################################################################

template<typename T> inline
void DenseMatrix<T>::fill_with_zero() {
 /*   
  *   Function : fill_with_zero()
  *   Author   : NeoKao <neokao@moldex3d.com>
  *   Brief    : Fill the matrix data with entry 0
  */  
  if (m_rows > 0 && m_cols > 0) {
    // -------------------------------------------
    for (size_type i = 0; i < m_rows; i++) {
      for (size_type j = 0; j < m_cols; j++) {
        m_val[i][j] = 0;
      }  // end j-loop
    }  // end i-loop
    // -------------------------------------------
  }
}

// #####################################################################################

template<typename T> inline
void DenseMatrix<T>::fill_with_one() {
 /*   
  *   Function : fill_with_zero()
  *   Author   : NeoKao <neokao@moldex3d.com>
  *   Brief    : Fill the matrix data with entry 0
  */  
  if (m_rows > 0 && m_cols > 0) {
    // -------------------------------------------
    for (size_type i = 0; i < m_rows; i++) {
      for (size_type j = 0; j < m_cols; j++) {
        m_val[i][j] = 1;
      }  // end j-loop
    }  // end i-loop
    // -------------------------------------------
  }
}

// #####################################################################################

template<typename T> inline
void DenseMatrix<T>::fill_with_value(const value_type& value) {
 /*   
  *   Function : fill_with_value()
  *   Author   : NeoKao <neokao@moldex3d.com>
  *   Brief    : Fill the matrix data with entry 1
  */  
  if (m_rows > 0 && m_cols > 0) {
    // -------------------------------------------
    for (size_type i = 0; i < m_rows; i++) {
      for (size_type j = 0; j < m_cols; j++) {
        m_val[i][j] = value;
      }  // end j-loop
    }  // end i-loop
    // -------------------------------------------
  }
}

// #####################################################################################

template<typename T> inline
void DenseMatrix<T>::createIdentity() {
 /*   
  *   Function : createIdentity()
  *   Author   : NeoKao <neokao@moldex3d.com>
  *   Brief    : Create the identity matrix
  */  
  this -> fill_with_zero();
  for (size_type i = 0; i < m_rows; i++) m_val[i][i] = 1;
}

template<typename T> inline
void DenseMatrix<T>::scale(const value_type& scale_factor) {
 /*   
  *   Function : scale
  *   Author   : NeoKao <neokao@moldex3d.com>
  *   Brief    : Scaling the DenseMatrix entries by the factor scaling_
  */  
  for (size_type i = 0; i < m_rows; i++) {
    for (size_type j = 0; j < m_cols; j++) {
      m_val[i][j] *= scale_factor;
    }
  }
}

// #####################################################################################

template<typename T> inline
void DenseMatrix<T>::scale(const size_type  &rowIndex,
                           const size_type  &colIndex,
                           const value_type &scale_factor) {
 /*   
  *   Function : scale
  *   Author   : NeoKao <neokao@moldex3d.com>
  *   Brief    : scale the specific matrix entry by the factor coeff
  */  
  m_val[rowIndex][colIndex] *= scale_factor;
}

// #####################################################################################

template<typename T> inline
void DenseMatrix<T>::scale_row(const size_type &rowIndex, const value_type &scale_factor) {
 /*   
  *   Function : scale_row
  *   Author   : NeoKao <neokao@moldex3d.com>
  *   Brief    : scale the matrix entry at the specific row index
  */  
  for (size_type jcol = 0; jcol < m_cols; jcol++) {
    m_val[rowIndex][jcol] *= scale_factor;
  }
}

// #####################################################################################

template<typename T> inline
void DenseMatrix<T>::scale_col(const size_type &colIndex, const value_type &scale_factor) {
/*   
  *   Function : scale
  *   Author   : NeoKao <neokao@moldex3d.com>
  *   Brief    : scale the matrix entry at the specific column index
  */  
  for (size_type irow = 0; irow < m_rows; irow++) {
    m_val[irow][colIndex] *= scale_factor;
  }
}

// #####################################################################################

template<typename T> inline
T DenseMatrix<T>::calcRowSumAtRowIndex(const size_type &rowIndex) {
  T RowSum = 0;
  for (size_type index = 0 ; index < m_cols ; index++) {
    RowSum += m_val[rowIndex][index];
  }
  return RowSum;
}

// #####################################################################################

template<typename T> inline
T DenseMatrix<T>::calcColSumAtColIndex(const size_type &colIndex) {
  T ColSum = 0;
  for (size_type index = 0 ; index < m_rows ; index++) {
    ColSum += m_val[index][colIndex];
  }
  return ColSum;
}

// #####################################################################################

template<typename T> inline
T DenseMatrix<T>::calcTrace() {
  T Trace = 0;
  for ( size_type index = 0; index < m_rows; index++ ) Trace += m_val[index][index];
  return Trace;
}

// #####################################################################################

template<typename T> inline
void DenseMatrix<T>::show() {
 /*   
  *   Function : show
  *   Author   : NeoKao <neokao@moldex3d.com>
  *   Brief    : output the profile of DenseMatrix
  */
  std::cout.precision(6);
  std::cout.setf(std::ios_base::fixed);
  //---------------------------------------------------
  for (size_type i = 0; i < m_rows; i++) {
    std::cout << std::endl;
    std::cout << " |";
    for (size_type j = 0; j < m_cols; j++) {
      std::cout << "  " << std::showpos << m_val[i][j];
    }
    std::cout << "  | ";
  }
  std::cout << std::endl;
  std::cout << std::endl;
}

// #####################################################################################

template<typename T> inline
std::vector<T>& DenseMatrix<T>::operator[](size_type id) {
 /*   
  *   Function : operator[]
  *   Author   : NeoKao <neokao@moldex3d.com>
  *   Brief    : overload the operator []
  */
  return m_val[id];
}

// #####################################################################################

template<typename T> inline
const std::vector<T>& DenseMatrix<T>::operator[](size_type id) const {
 /*   
  *   Function : operator[]
  *   Author   : NeoKao <neokao@moldex3d.com>
  *   Brief    : overload the operator []
  */
  return m_val[id];
}

// #####################################################################################

template<typename T> inline
void DenseMatrix<T>::DenseVectorMultiplyTransposeDenseVector(DenseVectorType& DenVec1,
                                                             DenseVectorType& DenVec2) {
 /*                 
  *   Function : vector_multiply_trans_vector
  *   Author   : NeoKao <neokao@moldex3d.com>
  *   Brief    : get the matrix profile by multiplying a vector 
  *              with a transpose vector
  */  
  for (size_type irow = 0; irow < m_rows; irow++) {
    for (size_type jcol = 0; jcol < m_cols; jcol++) {
      m_val[irow][jcol] = DenVec1[irow] * DenVec2[jcol];
    }  // end jcol-loop
  }  // end irow-loops
}

// #####################################################################################

template<typename T> inline
void DenseMatrix<T>::TransposeDenseMatrixMultiplyDenseVector(DenseVectorType &v, DenseVectorType &TranAv) {
 /*   
  *   Function : TransposeDenseMatrixMultiplyDenseVector
  *   Author   : NeoKao <neokao@moldex3d.com>
  *   Brief    : Perform the transpose matrix-vector multiplication
  */  
  T tmpValue = 0;
  TranAv.fill_with_zero();
  // check the consistence
  if (m_rows == TranAv.getSize()) {
    for (size_type jcol = 0; jcol < m_cols; jcol++) {
      tmpValue = 0;
      for (size_type irow = 0; irow < m_rows; irow++) {
        tmpValue += m_val[irow][jcol] * v.getValue(irow);
        TranAv.setValue(jcol, tmpValue);
      }  // end irow-loop
    }  // end jcol-loop
  }
}

// #####################################################################################

template<typename T>
inline void DenseMatrix<T>::DenseMatrixMultiplyDenseMatrix(self_type& DenMatA,
                                                           self_type& DenMatB) {
 /*
  *   Function : DenseMatrixMultiplyDenseMatrix
  *   Author   : NeoKao
  *   Brief    : Perform the dense matrix multiply dense matrix operation
  */
  this -> fill_with_zero();
  //-------------------------------------------
  for (size_type l = 0; l < m_cols; l++) {
    for (size_type j = 0; j < m_cols; j++) {
      for (size_type k = 0; k < m_cols; k++) {
        m_val[l][j] += DenMatA.m_val[l][k] * DenMatB.m_val[k][j];
      }  // end k-loop
    }  // end j-loop
  }  // end i-loop
}

// #####################################################################################

template<typename T>
inline void DenseMatrix<T>::TransposeDenseMatrixMultiplyDenseMatrix(self_type& DenMatA,
                                                                    self_type& DenMatB) {
 /*
  *   Function : TransposeDenseMatrixMultiplyDenseMatrix
  *   Author   : NeoKao
  *   Brief    : Perform the dense matrix multiply dense matrix operation
  */
  this -> fill_with_zero();
  //-------------------------------------------
  for (size_type l = 0; l < m_cols; l++) {
    for (size_type j = 0; j < m_cols; j++) {
      for (size_type k = 0; k < m_cols; k++) {
        m_val[l][j] += DenMatA.m_val[k][l] * DenMatB.m_val[k][j];
      }  // end k-loop
    }  // end j-loop
  }  // end i-loop
}

// #####################################################################################

template<typename T> inline
void DenseMatrix<T>::DenseMatrixMultiplyDenseVector(DenseVectorType& v, DenseVectorType& Av) {
 /*   
  *   Function : MultiplyVector
  *   Author   : NeoKao <neokao@moldex3d.com>
  *   Brief    : Perform the matrix-vector multiplication
  */  
  T tmpValue;
  Av.fill_with_zero();
  // check the consistence
  if (m_cols == Av.getSize()) {
    for (size_type irow = 0; irow < m_rows; irow++) {
      tmpValue = 0;
      for (size_type jcol = 0; jcol < m_cols; jcol++) {
        tmpValue += m_val[irow][jcol] * v.getValue(jcol);
        Av.setValue(irow, tmpValue);
      }  // end jcol-loop
    }  // end irow-loop
  }
}

// #####################################################################################

template<typename T> inline
DenseMatrix<T> DenseMatrix<T>::operator + (const DenseMatrix<T>& DenMat) {
 /*   
  *   Function : operator +
  *   Author   : NeoKao <neokao@moldex3d.com>
  *   Brief    : overload the operator "+"
  */
  inlakit::container::DenseMatrix<T> tmp_val;
  if (DenMat.m_rows > 0 && DenMat.m_cols > 0) {
    tmp_val.create(DenMat.m_rows , DenMat.m_cols);
    for (size_type i = 0; i < DenMat.m_rows; i++) {
      for (size_type j = 0; j < DenMat.m_cols; j++) {
        tmp_val[i][j] = m_val[i][j] + DenMat.m_val[i][j];
      }  // end j-loop
    }  // end i-loop
  }
  return tmp_val;
}

// #####################################################################################

template<typename T> inline
DenseMatrix<T> DenseMatrix<T>::operator - (const DenseMatrix<T>& DenMat) {
 /*   
  *   Function : operator -
  *   Author   : NeoKao <neokao@moldex3d.com>
  *   Brief    : overload the operator "-"
  */
  inlakit::container::DenseMatrix<T> tmp_val;
  if (DenMat.m_rows > 0 && DenMat.m_cols > 0) {
    tmp_val.create(DenMat.m_rows , DenMat.m_cols);
    for (size_type i = 0; i < DenMat.m_rows; i++) {
      for (size_type j = 0; j < DenMat.m_cols; j++) {
        tmp_val[i][j] = m_val[i][j] - DenMat.m_val[i][j];
      }  // end j-loop
    }  // end i-loop
  }
  return tmp_val;
}

// #####################################################################################

template<typename T> inline
DenseMatrix<T>& DenseMatrix<T>::operator = (const DenseMatrix<T>& DenMat) {
 /*   
  *   Function : operator =
  *   Author   : NeoKao <neokao@moldex3d.com>
  *   Brief    : overload the operator "="
  */
  if (this != &DenMat) {
    this -> m_rows = DenMat.m_rows;
    this -> m_cols = DenMat.m_cols;
    ///----------------------------------------------------------
    m_val.resize(m_rows);  // resize the m_val
    typename std::vector< std::vector<T> >::iterator iter;
    for (iter = m_val.begin(); iter != m_val.end(); ++iter) {
      iter -> resize(m_cols);
    }  // end iter-loop
  }
  for (size_type i = 0; i < m_rows; i++) {
    for (size_type j = 0; j < m_cols; j++) {
      m_val[i][j] = DenMat.m_val[i][j];
    }  // end j-loop
  }  // end i-loop
  return (*this);
}

// #####################################################################################

template<typename T> inline
DenseMatrix<T>& DenseMatrix<T>::operator *= (const T& scale_factor) {
  for (size_type i = 0; i < m_rows; i++) {
    for (size_type j = 0; j < m_cols; j++) {
      m_val[i][j] *= scale_factor;
    }  // end j-loop
  }  // end i-loop
  return *this;
}

// #####################################################################################

// template<typename T> inline
// DenseMatrix<T> &DenseMatrix<T>::operator*(const T scale_factor, const DenseMatrix<T> &DenMat) {
//  /*
//   *   Function : oeprator *
//   *   Author   : NeoKao <neokao@moldex3d.com>
//   *   Brief    : Overload the operator *
//   */
//   return DenseMatrix<T>(DenMat) *= scale_factor;
// }

// // #####################################################################################

// template<typename T> inline
// DenseMatrix<T> &DenseMatrix<T>::operator*(const DenseMatrix<T> &DenMat , const T scale_factor) {
//  /*
//   *   Function : oeprator *
//   *   Author   : NeoKao <neokao@moldex3d.com>
//   *   Brief    : Overload the operator *
//   */
//   return DenseMatrix<T>(DenMat) *= scale_factor;
// }

// template<typename T> inline
// DenseMatrix<T>& DenseMatrix<T>::operator * (const T& scale_factor, const DenseMatrix<T>& DenMat) {
//  /*
//   *   Function : oeprator *
//   *   Author   : NeoKao <neokao@moldex3d.com>
//   *   Brief    : Overload the operator *
//   */
//   return DenseMatrix<T>(DenMat) *= scale_factor;
// }

// // #####################################################################################

// template<typename T> inline
// DenseMatrix<T>& DenseMatrix<T>::operator * (const DenseMatrix<T>& DenMat , const T& scale_factor) {
//  /*
//   *   Function : oeprator *
//   *   Author   : NeoKao <neokao@moldex3d.com>
//   *   Brief    : Overload the operator *
//   */
//   return DenseMatrix<T>(DenMat) *= scale_factor;
// }

// #####################################################################################

template<typename T> inline
T DenseMatrix<T>::operator()(size_type rowIndex, size_type colIndex) {
 /*
  *   Function : operator ()
  *   Author   : NeoKao <neokao@moldex3d.com>
  *   Brief    : Overload the operator ()
  */ 
  return m_val[rowIndex][colIndex];
}

// ---------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------

// Instance (To export the dll library)
# pragma warning(disable:4661) /* avoid the 4661 warning message */
template class DenseMatrix<double>;

// ---------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------

}  // namespace container
}  // namespace inlakit

#endif  // SRC_CONTAINER_INLAKIT_DENSE_MATRIX_CPP_