#ifndef SRC_CONTAINER_INLAKIT_SPARSE_MATRIX_COO_CPP_
#define SRC_CONTAINER_INLAKIT_SPARSE_MATRIX_COO_CPP_

//------------------------------------------
#include <inlakit_container_sparse_matrix.h>
#include <inlakit_container_sparse_matrix_coo.h>
#include <inlakit_container_sparse_matrix_type.h>
//------------------------------------------

#include <iostream>
#include <iomanip>
#include <vector>

namespace inlakit   {
namespace container {

// #####################################################################################

template<typename T>
SparseMatrixCOO<T>::SparseMatrixCOO() {
  m_nnu   = 0;
  m_nna   = 0;
  m_incre = 0;
}

// #####################################################################################

template<typename T>
SparseMatrixCOO<T>::SparseMatrixCOO(const size_type &nnu) {
  if ( nnu > 0 ) this -> m_nnu = nnu;
  m_ia.reserve(m_nnu+1);

  //**************************************************************

  matrix_column_index = NULL;
  matrix_column_index = new std::vector<std::set<size_type>>(m_nnu);

  //**************************************************************

  m_incre = 0;

  // Generate the L2G mapping

  m_L2G.reserve(m_nnu);
  m_L2G.resize(m_nnu);

  for (size_type i = 0; i < m_nnu; i++) m_L2G[i] = i;
}

// #####################################################################################

template<typename T>
SparseMatrixCOO<T>::SparseMatrixCOO(const self_type &SpaMatCOO) {
 /*
  *   Function : SparseMatrixCOO
  *   Author   : NeoKao <neokao@moldex3d.com>
  *   Brief    : Copy constructor of Sparse matrix class
  */
  if (SpaMatCOO.m_nnu > 0 && SpaMatCOO.m_nna > 0) {
    m_nnu = SpaMatCOO.m_nnu;
    m_nna = SpaMatCOO.m_nna;
    //--------------------------------------
    m_jr.resize(m_nna);
    m_jc.resize(m_nna);
    m_aa.resize(m_nna);
    //--------------------------------------
    for (size_type i = 0; i < m_nna; i++) {
      m_jr[i] = SpaMatCOO.m_jr[i];
      m_jc[i] = SpaMatCOO.m_jc[i];
      m_aa[i] = SpaMatCOO.m_aa[i];
    }  // end i-loop
  }
}

// #####################################################################################

template<typename T>
SparseMatrixCOO<T>::~SparseMatrixCOO() {
 /*
  *   Function : ~SparseMatrixCOO
  *   Author   : NeoKao <neokao@moldex3d.com>
  *   Brief    : Destructor of Sparse matrix class
  */
  m_nnu = 0;
  m_nna = 0;
  m_jr.clear();
  m_jc.clear();
  m_aa.clear();
}

// #####################################################################################

template<typename T> inline
void SparseMatrixCOO<T>::addIndex(const size_type& rowIndex , const size_type& colIndex) {
 /*
  *  Function : addIndex
  *  Author   : NeoKao <neokao@moldex3d.com>
  *  Brief    : add the sparse matrix space
  */
  if (rowIndex >= 0 && colIndex >= 0) {
    /* Add index (row_id/col_id) */
    matrix_column_index -> at(rowIndex).insert(colIndex);
    /* set not compress */
    IsCompressed = false;
  }
}

// #####################################################################################

template<typename T> inline
bool SparseMatrixCOO<T>::compressIndex() {
 /*
  *  Function : compressIndex
  *  Author   : NeoKao <neokao@moldex3d.com>
  *  Brief    : compress the matrix into the sparse matrix
  */
  this -> compressIndexCOO();
  return true;
}

// #####################################################################################

template<typename T> inline
bool SparseMatrixCOO<T>::compressIndexCOO() {
 /*
  *  Function : compressIndexCOO
  *  Author   : NeoKao <neokao@moldex3d.com>
  *  Brief    : compress the input matrix into the COO format
  */
  if (!m_jr.empty()) { m_jr.clear(); }
  if (!m_jc.empty()) { m_jc.clear(); }
  if (!m_aa.empty()) { m_aa.clear(); }

  if (!m_ia.empty()) { m_ia.clear(); }

  /* ----- construct m_ia ----- */

  m_ia.resize(m_nnu+1);

  int AccumulateCount = 0;

  m_ia[0] = AccumulateCount;

  for (size_type irow = 0; irow < m_nnu; irow++) {
    AccumulateCount += int(matrix_column_index -> at(irow).size());  //  NOLINT
    m_ia[irow+1] = AccumulateCount;
  }

  /* ----- construct m_coo_row , m_coo_col ----- */

  m_nna = 0;

  for (int irow = 0; irow < m_nnu; irow++) {
    /* --- for each row, add the number of non-zero entries --- */
    m_nna += int(matrix_column_index -> at(irow).size());  //  NOLINT
  }

  /* ---------- clean aa ---------- */

  m_jr.resize(m_nna);
  m_jc.resize(m_nna);
  m_aa.resize(m_nna);

  int icr = 0;

  for (int irow = 0; irow < m_nnu; irow++) {
    /* ------------------------------------------ */
    for (std::set<size_type>::iterator sit = matrix_column_index -> at(irow).begin();
      sit != matrix_column_index -> at(irow).end(); ++sit) {
      m_jr[icr] = irow;
      m_jc[icr] = (*sit);
      icr++;
    }
  }

  //m_nnu_global = m_nnu;
  //m_nnu_local  = m_nnu_global;

  //m_nna_global = m_nna;
  //m_nna_local  = m_nna_global;

  // ----- clean m_aa -----

  //  m_aa.clear();

  // ----- set compress -----

  IsCompressed = true;
  return IsCompressed;
}

template<typename T> inline
void SparseMatrixCOO<T>::spmv(DenseVectorType &v, DenseVectorType &Av) {
 /*
  *  Function : spmv
  *  Author   : NeoKao
  *  Brief    : Perform the sparse matrix-vector multiplication (spmv)
  *             The sparse matrix is stored in COO format
  */
  Av.fill_with_zero();
  for (size_type inna = 0; inna < m_nna; inna++) {
    Av.m_val[m_jr[inna]] += m_aa[inna] * v.m_val[m_jc[inna]];
  }
}

// ---------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------

// Instance (To export the dll library)
# pragma warning(disable:4661) /* avoid the 4661 warning message */
template class SparseMatrixCOO<double>;

// ---------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------

}  // namespace container
}  // namespace inlakit

#endif  // SRC_CONTAINER_INLAKIT_SPARSE_MATRIX_COO_CPP_