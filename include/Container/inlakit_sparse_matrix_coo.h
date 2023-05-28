#ifndef IFSI2D_CONTAINER_INLAKIT_SPARSE_MATRIX_COO_H_
#define IFSI2D_CONTAINER_INLAKIT_SPARSE_MATRIX_COO_H_

//------------------------------------------
#include <inlakit_container_dense_vector.h>
#include <inlakit_container_sparse_matrix_type.h>
#include <vector_traits.h>
#include <matrix_traits.h>

//#include <inlakit_calculator_sparse_matrix.h>
//------------------------------------------

#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <utility>
#include <set>
#include <string>

namespace inlakit   {
namespace container {

template<typename T>
class SparseMatrixCOO {

    template<typename S> friend class DenseVector;
    template<typename S> friend class DenseMatrix;
    template<typename S> friend class SparseMatrixCalc;

    using self_type = SparseMatrixCOO<T>;

 public:

    using size_type  = int;
    using value_type = T;

    // ----------------------------------------------------------------

    using DenseVectorType = inlakit::container::template DenseVector<T>;

    // ----------------------------------------------------------------

    using vector_size_type  = typename inlakit::traits::VectorTraits<DenseVectorType>::vector_size_type;
    using vector_value_type = typename inlakit::traits::VectorTraits<DenseVectorType>::vector_value_type;

    /* Default constructor */
    SparseMatrixCOO();

    /* Constructor */
    SparseMatrixCOO(const size_type &nnu);  // NOLINT

    /* Copy constructor */
    SparseMatrixCOO(const self_type &SparseMatrixCOO); // NOLINT

    /* Destructor */
    virtual ~SparseMatrixCOO();

    // ---------------------------------------------------------------

    /* get the format of sparse matrix */
    inline std::string getFormat() { return m_format; }

    /* set the global sparse matrix size >>> unit testing is passed */
    // inline void setGlobalSize(const size_type &nnu_global) { m_nnu_global = nnu_global; }

    /* set the total number of non-zero entries in sparse matrix >>> unit testing is passed */
    // inline void setGlobalNonZeroSize(const size_type &nna_global) { m_nna_global = nna_global; }

    /* set the local sparse matrix size(in parallel solver) */
    // inline void setLocalSize(const size_type &nnu_local) { m_nnu_local = nnu_local; }

    /* set the local number of non-zero entries in sparse matrix(in parallel solver) */
    // inline void setLocalNonZeroSize(const size_type &nna_local) { m_nna_local = nna_local; }

    /* set the indexing base of sparse matrix */
    inline void setIndexBase(const std::string IndexBase) { m_IndexBase = IndexBase; }

    /* get the size of sparse matrix >>> unit testing is passed */
    inline size_type getSize() { return m_nnu; }

    /* get the number of non-zero entries in sparse matrix >>> unit testing is passed */
    inline size_type getNonZeroSize() { return m_nna; }

    /* get the number of non-zero entries in upper triangular sparse matrix */
    // inline size_type getNonZeroSize_uTri() { return m_nna_uTri; }

    /* get the number of non-zero entries in lower triangular sparse matrix */
    // inline size_type getNonZeroSize_lTri() { return m_nna_lTri; }

    /* get the size of the non-zero entries at the specific row index >>> unit testing is passed */
    inline size_type getNonZeroSizeAtOneRow(const size_type &rowIndex) { return int(matrix_column_index -> at(rowIndex).size()); } // NOLINT

    /* get the global matrix size */
    // inline size_type getGlobalSize() { return m_nnu_global; }

    /* get the global number of non-zero entries */
    // inline size_type getGlobalNonZeroSize() { return m_nna_global; }

    /* get the local matrix size */
    // inline size_type getLocalSize() { return m_nnu_local; }

    /* get the local number of non-zero matrix entries */
    // inline size_type getLocalNonZeroSize() { return m_nna_local; }

    /* get the format of sparse matrix */
    inline std::string getIndexBase() { return m_IndexBase; }

    /* compress the sparse matrix space >>> unit testing is passed */
    bool compressIndex();
    bool compressIndexCOO();

    /* add the sparse matrix space >>> unit testing is passed */
    void addIndex(const size_type &rowIndex, const size_type &colIndex);

    /* add the value to the sparse matrix [COO] >>> unit testing is passed */
    bool addValue(const size_type &rowIndex, const size_type &colIndex, const value_type &value);

    /* perform the sparse matrix-vector multiplication [COO format] >>> unit testing is passed */
    void spmv(DenseVectorType &v, DenseVectorType &Av);  // spmv function

    /* exchange the sparse matrix data */
    void swapCOO(std::vector<size_type>  &jr, std::vector<size_type> &jc,
                 std::vector<value_type> &aa, std::vector<size_type> &ia);

    /* exchange the Local to global mapping table */
    void swapL2G(std::vector<size_type> &L2G);

    /* get the upper triangular profile of COO sparse matrix */
    bool GetUpperTriangularMatrix();

    /* get the lower triangular profile of COO sparse matrix */
    bool GetLowerTriangularMatrix();

    /* get the IA value at the Idx location, this IA is identical to the one in CSR format */
    inline size_type getValue_IA(const size_type &Index) { return m_ia[Index]; }

    /* get the JR value [row-location] of general COO sparse matrix at the location Idx >>> unit test is passed */
    inline size_type getValue_JR(const size_type &Index) { return m_jr[Index]; }

    /* get the JC value [column-location] of general COO sparse matrix at the location Idx >>> unit test is passed */
    inline size_type getValue_JC(const size_type &Index) { return m_jc[Index]; }

    /* get the AA value [non-zero matrix value] of general COO sparse matrix at the location Idx */
    inline value_type getValue_AA(const size_type &Index) { return m_aa[Index]; }

    /* get the JR value [row-location] of upper triangular COO sparse matrix at the location Idx */
    // inline size_type getValue_uTri_JR(const size_type &Index) { return m_jr_uTri[Index]; }

    /* get the JC value [column-location] of upper triangular COO sparse matrix at the location Idx */
    // inline size_type getValue_uTri_JC(const size_type &Index) { return m_jc_uTri[Index]; }

    /* get the AA value [non-zero matrix value] of upper triangular COO sparse matrix at the location Idx */
    // inline value_type getValue_uTri_AA(const size_type &Index) { return m_aa_uTri[Index]; }

    /* get the JR value [row-location] of upper triangular COO sparse matrix at the location Idx */
    // inline size_type getValue_lTri_JR(const size_type &Index) { return m_jr_lTri[Index]; }

    /* get the JC value [column-location] of lower triangular COO sparse matrix at the location Idx */
    // inline size_type getValue_lTri_JC(const size_type &Index) { return m_jc_lTri[Index]; }

    /* get the AA value [non-zero matrix value] of lower triangular COO sparse matrix at the location Idx */
    // inline value_type getValue_lTri_AA(const size_type &Index) { return m_aa_lTri[Index]; }

    /* get the local to global mapping at the location Idx */
    inline size_type getValue_L2G(const size_type &Index) { return m_L2G[Index]; }

    /* output the profile of general COO sparse matrix */
    void show();

    /* output the profile of upper triangular COO sparse matrix */
    // void show_u_tri();

    /* read the matrix profile with name MatrixName */
    // bool read(INLAKIT_INT &cpuID, INLAKIT_STRING MatrixName, INLAKIT_STRING IndexBase);

 private:
    std::string             m_format{"COO"};   // Sparse matrix format
    std::string             m_IndexBase{"C"};  // The base of index(C or Fortran)
    size_type               m_nnu;  // matrix size
    size_type               m_nna;  // the number of non-zero entries
    std::vector<size_type>  m_jr;   // the row index of non-zero entry
    std::vector<size_type>  m_jc;   // the column index of non-zero entry
    std::vector<value_type> m_aa;   // the value of non-zero entry

    std::vector<size_type> m_ia;
    std::vector<size_type> m_L2G;   // local to global mapping

    /// ========== Compress ID tool ==========

    std::vector< std::set<size_type> > *matrix_column_index;
    size_type m_incre;
    bool      IsCompressed = false;

};

}  // namespace container
}  // namespace inlakit


#endif  // IFSI2D_CONTAINER_INLAKIT_SPARSE_MATRIX_COO_H_