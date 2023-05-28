/************************************************************************************/
//
//  FILENAME  : inlakit_dense_matrix.h
//  PURPOSE   : Definition of "DENSE" matrix class in INLAKit
//  PROJECT   : INLAKIT
//  COPYRIGHT : (C) 2023, CoreTech System CO.,LTD.(http://www.moldex3d.com)
//              Intellegent property of CoreTech System, restricted to developers
//              of RD team of CoreTech specified people. Use, reproduction,
//              copying, transfer and disclosure to others of the codes are
//              prohibited except by written permission of CoreTech System
//
//  AUTHOR   : NeoKao (neokao@moldex3d.com)
//
//  HISTORY  : VERSION   NOTE                                      DATE
//             1.0.0     Release the 1st version                   2022/7/15
//  ------------------------------------------------------------------------
//
/************************************************************************************/

#ifndef  SRC_INC_INLAKIT_CONTAINER_DENSE_MATRIX_H_
#define  SRC_INC_INLAKIT_CONTAINER_DENSE_MATRIX_H_

/* ---------- header files ---------- */
#include <inlakit_container_dense_vector.h>

#include <inlakit_container_sparse_matrix.h>
#include <inlakit_container_sparse_matrix_coo.h>
#include <inlakit_container_sparse_matrix_csr.h>

#include <matrix_traits.h>
#include <vector_traits.h>
/* ---------------------------------- */

#include <inlakit_algorithm_transformation.h>
#include <inlakit_calculator_dense_matrix.h>

/* ------------------------------------------ */

#include <inlakit_direct_solver.h>

/* ------------------------------------------ */

#include <iostream>
#include <algorithm>
#include <vector>

namespace inlakit   {
namespace container {

template<typename T>
class DenseMatrix {

  // -------------------------------------------------------

  template<typename S> friend class DenseVector;
  template<typename S> friend class DenseMatrix;
  template<typename S> friend class SparseMatrixCOO;
  template<typename S> friend class SparseMatrixCSR;

  // -------------------------------------------------------

  template<typename S> friend class MatrixTransformation;

  // -------------------------------------------------------

  template<typename S> friend class DenseMatrixCalc;
  template<typename S> friend class SparseMatrixCalc;

  // -------------------------------------------------------

  template<typename S> friend class DirectSolver;

  // -------------------------------------------------------

  using self_type = DenseMatrix<T>;
  using base_type = std::vector<std::vector<T>>;

 public :

  using size_type  = int;
  using value_type = T;

  // -------------------------------------------------------

  using DenseVectorType   = inlakit::container::template DenseVector<T>;
  using vector_value_type = typename inlakit::traits::VectorTraits<DenseVectorType>::vector_value_type;
  using matrix_size_type  = std::vector<std::vector<size_type>>;
  using matrix_value_type = std::vector<std::vector<value_type>>;
  using data_type         = base_type;

  // --------------------------------------------------------------------------------------------------

  /* Default constructor */
  DenseMatrix();

  /* Constructor */
  DenseMatrix(const size_type &rowSize, const size_type &colSize);

  /* Copy constructor */
  DenseMatrix(const self_type &DenMat);  // NOLINT

  /* Destructor */
  virtual ~DenseMatrix();

  /* Check the validity of matrix (row,col) index >>> (unit testing is passed) */
  bool checkIndex(const size_type &rowIndex, const size_type &colIndex);

  void create(const size_type &rowSize, const size_type &colSize);

  /* Fill all the dense matrix compoents with value 0 >>> (unit testing is passed) */
  void fill_with_zero();

  /* Fill all the dense matrix components with value 1 >>> (unit testing is passed) */
  void fill_with_one();

  /* Fill all the dense matrix components with the specific value >>> (unit testing is passed) */
  void fill_with_value(const value_type &value);

  /* Create the identity matrix >>> (unit testing is passed) */
  void createIdentity();

  /* Cooy the dense matrix object to the self dense matrix object >>> (unit testing is passed) */
  void copy(const self_type &DenMat);

  /* Scale all the dense matrix components with DenseMatrixScaling >>> (unit testing is passed) */
  void scale(const value_type &scale_factor);

  /* Scale the specific dense matrix (RowID,ColID) components with DenseMatrixScaling */
  void scale(const size_type &rowIndex, const size_type &colIndex, const value_type &scale_factor);

  /* Scale the matrix for the specific row index */
  void scale_row(const size_type &rowIndex, const value_type &scale_factor);

  /* Scale the matrix for the specific column index */
  void scale_col(const size_type &colIndex, const value_type &scale_factor);

  // To do list .....
  //  value_type getDenseMatrixNorm1()    ; /* get the 1-norm     of dense matrix */
  //  value_type getDenseMatrixNorm2()    ; /* get the 2-norm     of dense matrix */
  //  value_type getDenseMatrixNormInfty(); /* get the infty-norm of dense matrix */
  //  value_type getDenseMatrixConditionNumber();

  /* Set the value of dense matrix at the specific location >>> (unit testing is passed) */
  void setValue(const size_type &rowIndex, const size_type &colIndex, const value_type &value);

  /* get the row size of dense matrix >>> (unit testing is passed) */
  size_type getRowSize() { return m_rows; }

  /* get the column size of dense matrix >>> (unit testing is passed) */
  size_type getColSize() { return m_cols; }

  /* get the transpose of the dense matrix */
  void transpose();

  /* get the value of dense matrix at the specific location >>> (unit testing is passed) */
  T getValue(const size_type &rowIndex, const size_type &colIndex) { return m_val[rowIndex][colIndex]; }

  /* get the row sum of dense matrix along the specific row ID >>> (unit testing is passed) */
  T calcRowSumAtRowIndex(const size_type &rowIndex);

  /* get the column sum of dense matrix along the specific column ID >>> (unit testing is passed) */
  T calcColSumAtColIndex(const size_type &colIndex);

  /* get the diagonal entries of dense matrix >>> (unit testing is passed) */
  void getDiagonalValues(DenseVectorType &diagonalValues);

  /* get the trace of dense matrix (the sum of diagonal entries) */
  T calcTrace();

  /* get the diagonal entry of dense matrix at the specific location >>> (unit testing is passed) */
  T getDiagonalValue(const size_type &Index) { return m_val[Index][Index]; }

  //  void getDenseMatrixMaxNonDiagonalAbsoluteValue(size_type  &rowIndex,
  //                                                         size_type  &colIndex,
  //                                                         value_type &value);

  /* show dense matrix profile */
  void show();

  /* calculate the matrix-vector multuplication >>> (unit testing is passed) */
  void DenseMatrixMultiplyDenseVector(DenseVectorType &v, DenseVectorType &Av);

  /* calculate the trapsose matrix-vector multiplication >>> (unit testing is passed) */
  void TransposeDenseMatrixMultiplyDenseVector(DenseVectorType &v , DenseVectorType &TranAv);

  /* calculate the dense matrix multiple dense matrix >>> (unit testing is passed) */
  void DenseMatrixMultiplyDenseMatrix(self_type &A , self_type &B);

  /* calculate the transpose dense matrix multiple dense matrix >>> (unit testing is passed) */
  void TransposeDenseMatrixMultiplyDenseMatrix(self_type &A , self_type &B);

  /* get the cross-product using two vectors >>> (unit testing is passed) */
  void DenseVectorMultiplyTransposeDenseVector(DenseVectorType &DenVec1 , DenseVectorType &DenVec2);

  ///  ----------------------------
  ///      overload operator
  ///  ----------------------------

  DenseMatrix<T>  operator+(const DenseMatrix<T>& DenMat);
  DenseMatrix<T>  operator-(const DenseMatrix<T>& DenMat);
  DenseMatrix<T> &operator=(const DenseMatrix<T>& DenMat);

  DenseMatrix<T> &operator*=(const T& scaling);

  friend DenseMatrix<T> operator*(DenseMatrix<T> DenMat, double scale_factor) {DenMat *= scale_factor; return DenMat;} // NOLINT
  friend DenseMatrix<T> operator*(double scale_factor, DenseMatrix<T> DenMat) {DenMat *= scale_factor; return DenMat;} // NOLINT

  template<typename S>
  friend DenseMatrix<T>&operator*(const S coeff , const DenseMatrix<T>& mat);

  template<typename S>
  friend DenseMatrix<T>&operator*(const DenseMatrix<T>& mat , const S coeff);

        std::vector<T>& operator[](size_type id);
  const std::vector<T>& operator[](size_type id) const;

  T operator()(size_type rowIndex, size_type colIndex);

 private:
  size_type         m_rows;  /* row size */
  size_type         m_cols;  /* col size */
  matrix_value_type m_val;   /* dense matrix data */
};

}  // namespace container
}  // namespace inlakit

#endif  // SRC_INC_INLAKIT_CONTAINER_DENSE_MATRIX_H_