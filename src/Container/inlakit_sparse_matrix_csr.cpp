#ifndef SRC_CONTAINER_INLAKIT_SPARSE_MATRIX_CSR_CPP_
#define SRC_CONTAINER_INLAKIT_SPARSE_MATRIX_CSR_CPP_

//----------------------------------------------------------
#include <inlakit_container_dense_vector.h>
#include <inlakit_container_sparse_matrix.h>
#include <inlakit_container_sparse_matrix_csr.h>
#include <inlakit_container_sparse_matrix_type.h>
//----------------------------------------------------------

#include <iostream>
#include <iomanip>
#include <vector>

namespace inlakit   {
namespace container {

// #####################################################################################

template<typename T>
SparseMatrixCSR<T>::SparseMatrixCSR() {
  m_nnu = 0;
  m_nna = 0;
}

// #####################################################################################

template<typename T>
SparseMatrixCSR<T>::SparseMatrixCSR(const size_type &nnu) {
  if ( nnu > 0 ) this -> m_nnu = nnu;
  m_ia.reserve(m_nnu + 1);

  //--------------------------------------------------------------

  matrix_column_index = NULL;
  matrix_column_index = new std::vector<std::set<size_type>>(m_nnu);

  // Generate the L2G mapping

  m_L2G.reserve(nnu);
  m_L2G.resize(nnu);

  for (size_type i = 0; i < nnu; i++) m_L2G[i] = i;
}

// #####################################################################################

template<typename T>
SparseMatrixCSR<T>::SparseMatrixCSR(const self_type &SpaMatCSR) {

 /*
  *   Function : SparseMatrixCSR
  *   Author   : NeoKao <neokao@moldex3d.com>
  *   Brief    : Copy constrctor of class SparseMatrixCSR
  */

  if (SpaMatCSR.m_nnu > 0 && SpaMatCSR.m_nna > 0) {
    m_nnu = SpaMatCSR.m_nnu;
    m_nna = SpaMatCSR.m_nna;
    //--------------------------------------
    m_ia.resize(m_nnu + 1);
    m_ja.resize(m_nna);
    m_aa.resize(m_nna);
    //--------------------------------------
    for (size_type i = 0; i < m_nnu+1; i++) m_ia[i] = SpaMatCSR.m_ia[i];
      for (size_type i = 0; i < m_nna  ; i++) {
        m_ja[i] = SpaMatCSR.m_ja[i];
        m_aa[i] = SpaMatCSR.m_aa[i];
    }
  }
}

// #####################################################################################

template<typename T>
SparseMatrixCSR<T>::~SparseMatrixCSR() {
 /*
  *   Function : SparseMatrixCSR
  *   Author   : NeoKao <neokao@moldex3d.com>
  *   Brief    : Destructor of class SparseMatrixCSR
  */
  m_nnu = 0;
  m_nna = 0;
  m_ia.clear();
  m_ja.clear();
  m_aa.clear();
}

// #####################################################################################

template<typename T> inline
bool SparseMatrixCSR<T>::addValue(const size_type  &rowIndex,
                                  const size_type  &colIndex,
                                  const value_type &value) {
 /*
  *   Function : addValue
  *   Author   : NeoKao <neokao@moldex3d.com>
  *   Brief    : add the (row,col) index and value to the CSR format
  */

  for (size_type index = m_ia[rowIndex]; index < m_ia[rowIndex+1]; ++index) {
    if ((colIndex) == m_ja[index]) {
      m_aa[index] += value;
       return true;
    }
  }  // end index-loop
  return true;
}

// #####################################################################################

template<typename T> inline
bool SparseMatrixCSR<T>::copy(self_type & SparseMatrix) {
  if (SparseMatrix.getSize() > 0 && SparseMatrix.getNonZeroSize() > 0) {
    m_nnu = SparseMatrix.getSize();
    m_nna = SparseMatrix.getNonZeroSize();
    //--------------------------------------
    m_ia.resize(m_nnu + 1);
    m_ja.resize(m_nna);
    m_aa.resize(m_nna);
    //--------------------------------------
    for (size_type i = 0; i < m_nnu + 1; i++) m_ia[i] = SparseMatrix.getValue_IA(i);
    for (size_type i = 0; i < m_nna    ; i++) m_ja[i] = SparseMatrix.getValue_JA(i);
    for (size_type i = 0; i < m_nna    ; i++) m_aa[i] = SparseMatrix.getValue_AA(i);
    //--------------------------------------
  }
  return true;
}

// #####################################################################################

template<typename T>
void SparseMatrixCSR<T>::setValue_IA(const size_type &id , const size_type &value) {
  if (id >= 0) m_ia[id] = value;
}

// #####################################################################################

template<typename T>
void SparseMatrixCSR<T>::setValue_JA(const size_type &id , const size_type &value) {
  if (id >= 0) m_ja[id] = value;
}

// #####################################################################################

template<typename T>
void SparseMatrixCSR<T>::setValue_AA(const size_type &id , const value_type &value) {
  if (id >= 0) m_aa[id] = value;
}

// #####################################################################################

template<typename T>
void SparseMatrixCSR<T>::setValues_IA(const  size_type *ia) {
  m_ia.assign(ia, ia + m_ia.size());
}

// #####################################################################################

template<typename T>
void SparseMatrixCSR<T>::setValues_JA(const  size_type *ja) {
  m_ja.assign(ja, ja + m_ja.size());
}

// #####################################################################################

template<typename T>
void SparseMatrixCSR<T>::setValues_AA(const value_type *aa) {
  m_aa.assign(aa, aa + m_aa.size());
}

// #####################################################################################

template<typename T> inline
void SparseMatrixCSR<T>::show() {
 /*
  *   Function : ShowSparseMatrixProfile()
  *   Author   : NeoKao <neokao@moldex3d.com>
  *   Brief    : Output the sparse matrix profile (CSR format)
  */  

  std::cout << std::endl;
  std::cout << " >>> The sparse matrix format : " << this -> getFormat() << std::endl;
  std::cout << "       # ia size = " << m_ia.size() << std::endl;
  std::cout << "       # ja size = " << m_ja.size() << std::endl;
  std::cout << "       # aa size = " << m_aa.size() << std::endl;
  // ++++++++++ output ia ++++++++++
  std::cout << std::endl;
  std::cout << " >>> ia : " << std::endl;
  std::cout << std::endl;
  std::cout << " |";
  for (size_type i = 0; i < m_nnu+1; i++) std::cout << " " << m_ia[i];
  std::cout << " |";
  // ++++++++++ output ja ++++++++++
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << " >>> ja : " << std::endl;
  std::cout << std::endl;
  std::cout << " |";
  for (size_type i = 0; i < m_nna; i++) std::cout << " " << m_ja[i];
  std::cout << " |";
  // ++++++++++ output aa ++++++++++
  std::cout.precision(3);
  std::cout.setf(std::ios_base::fixed);
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << " >>> aa : " << std::endl;
  std::cout << std::endl;
  std::cout.precision(5);
  std::cout.setf(std::ios_base::fixed);
  std::cout << " | ";
  for (size_type i = 0; i < m_nna; i++) std::cout << m_aa[i] << " ";
  std::cout << "|";
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;
}

// #####################################################################################

template<typename T> inline
void SparseMatrixCSR<T>::ZeroToOneBasedIndexing() {
 /*
  *   @Function : ZeroToOneBasedIndexing
  *   @Author   : NeoKao
  *   @Brief    : change the index of CSR matrix from zero to one based 
  */
  if (m_ia[0] == 0) {
    for (size_type i = 0; i < m_nnu+1; i++) m_ia[i]++;
    for (size_type i = 0; i < m_nna  ; i++) m_ja[i]++;
  }
}

// #####################################################################################

template<typename T> inline
void SparseMatrixCSR<T>::OneToZeroBasedIndexing() {
 /*
  *   @Function : OneToZeroBasedIndexing
  *   @Author   : NeoKao
  *   @Brief    : change the index of CSR matrix from one to zero based 
  */
  if (m_ia[0] == 1) {
    for (size_type i = 0; i < m_nnu+1; i++) m_ia[i]--;
    for (size_type i = 0; i < m_nna  ; i++) m_ja[i]--;
  }
}

// #####################################################################################

template<typename T> inline
void SparseMatrixCSR<T>::swap(std::vector<size_type> &ia,
                              std::vector<size_type> &ja,
                              std::vector<value_type> &aa) {
  std::vector<size_type > ().swap(m_ia);  // clear and recycle the space
  std::vector<size_type > ().swap(m_ja);  // clean and recycle the space
  std::vector<value_type> ().swap(m_aa);   // clean and recycle the space

  m_ia.resize(ia.size());
  m_ja.resize(ja.size());
  m_aa.resize(aa.size());

  m_ia.swap(ia);
  m_ja.swap(ja);
  m_aa.swap(aa);

  m_nnu = m_ia.size()-1;
  m_nna = m_ja.size();

  /* 1-based -> 0-based */

  for (size_type i = 0; i < ia.size(); i++) m_ia[i]-- ;
  for (size_type i = 0; i < ja.size(); i++) m_ja[i]-- ;

  // ----------------------------------------------------------------

  std::vector<size_type> ().swap(ia);  // clear and recycle the space
  std::vector<size_type> ().swap(ja);  // clean and recycle the space
  std::vector<value_type> ().swap(aa);  // clean and recycle the space
}

// #####################################################################################

template<typename T> inline
void SparseMatrixCSR<T>::addIndex(const size_type &rowIndex,
                                          const size_type &colIndex) {
 /*
  *  Function : addIndex
  *  Author   : NeoKao <neokao@moldex3d.com>
  *  Brief    : add the sparse matrix space
  */
  if (rowIndex >= 0 && colIndex >= 0) {
    // Add index (row_id/col_id)
    matrix_column_index -> at(rowIndex).insert(colIndex);
    // set not compress
    IsCompressed = false;
  }
}

// #####################################################################################

template<typename T> inline
bool SparseMatrixCSR<T>::compressIndex() {
 /*
  *  Function : compressIndex
  *  Author   : NeoKao <neokao@moldex3d.com>
  *  Brief    : compress the matrix into the sparse matrix
  */
  this -> compressIndexCSR();
  return true;
}

// #####################################################################################

template<typename T> inline
bool SparseMatrixCSR<T>::compressIndexCSR() {
 /*
  *  Function : compressIndexCSR
  *  Author   : NeoKao <neokao@moldex3d.com>
  *  Brief    : compress the input matrix into the CSR format
  */

  if (!m_ia.empty()) { m_ia.clear(); }
  if (!m_ja.empty()) { m_ja.clear(); }

  // ----- construct m_ia -----

  m_ia.resize(m_nnu+1);

  size_type AccumulateCount = 0;

  m_ia[0] = AccumulateCount;

  for (size_type irow = 0; irow < m_nnu; irow++) {
    AccumulateCount += int(matrix_column_index -> at(irow).size());  // NOLINT
    m_ia[irow+1] = AccumulateCount;
  }

  // number of non-zero entries

  m_nna = AccumulateCount;

  // ----- m_ja , m_aa -----

  //  m_ja.resize(m_nna);
  //  m_aa.resize(m_nna);

  m_ja.resize(m_nna);
  m_aa.resize(m_nna);

  // -----------------------

  int tempCount = 0;

  for (size_type irow = 0; irow < m_nnu; ++irow) {
    tempCount = m_ia[irow];
    /* ----- insert m_csr_col index ----- */
    for (std::set<size_type>::iterator sit = matrix_column_index -> at(irow).begin();
      sit != matrix_column_index -> at(irow).end(); ++sit) {
      m_ja[tempCount] = (*sit);
      ++tempCount;
    }
  }

  /* ----- aa clear ----- */

  for (int i = 0; i < AccumulateCount; ++i) { m_aa[i] = 0.0; }

  /* ----- set compress ----- */

  IsCompressed = true;

  return IsCompressed;
}

template<typename T> inline
void SparseMatrixCSR<T>::spmv(DenseVectorType &v , DenseVectorType& Av) {
 /*
  *  Function : spmv
  *  Author   : NeoKao
  *  Brief    : Perform the sparse matrix-vector multiplication (spmv)
  *             The sparse matrix is stored in CSR format
  */
  T tmpValue;
  
  /* Initialization the vector */
  Av.fill_with_zero();

  for (size_type irow = 0; irow < m_nnu; irow++) {
    tmpValue = 0;
    for (size_type jcol = m_ia[irow]; jcol < m_ia[irow+1]; jcol++) {
      tmpValue += m_aa[jcol] * v.getValue(m_ja[jcol]);
    }
    Av.setValue(irow, tmpValue);
  }
}

// ---------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------
// Instance (To export the dll library)
# pragma warning(disable:4661) /* avoid the 4661 warning message */
template class SparseMatrixCSR<double>;
// ---------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------

}  // namespace container
}  // namespace inlakit

#endif  // SRC_CONTAINER_INLAKIT_SPARSE_MATRIX_CSR_CPP_