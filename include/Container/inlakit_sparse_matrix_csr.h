#ifndef IFSI2D_CONTAINER_INLAKIT_SPARSE_MATRIX_CSR_H_
#define IFSI2D_CONTAINER_INLAKIT_SPARSE_MATRIX_CSR_H_

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
#include <string>
#include <set>
#include <algorithm>
#include <utility>

namespace inlakit   {
namespace container {

template<typename T>
class SparseMatrixCSR {

    template<typename S> friend class DenseVector;
    template<typename S> friend class DenseMatrix;
    template<typename S> friend class SparseMatrixCalc;

    using self_type = SparseMatrixCSR<T>;

 public:

    using size_type  = int;
    using value_type = T;

    // ----------------------------------------------------------------

    using DenseVectorType = inlakit::container::template DenseVector<T>;
    
    // ----------------------------------------------------------------

    using vector_size_type  = typename inlakit::traits::VectorTraits<DenseVectorType>::vector_size_type;
    using vector_value_type = typename inlakit::traits::VectorTraits<DenseVectorType>::vector_value_type;

    /* Default constructor */
    SparseMatrixCSR();

    /* Constructor */
    SparseMatrixCSR(const size_type &nnu);  // NOLINT

    /* Copy constructor */
    SparseMatrixCSR(const self_type &SparseMatrixCSR); // NOLINT

    /* Destructor */
    virtual ~SparseMatrixCSR();

    // ---------------------------------------------------------------

    /* return the entire m_ia(vector type) */
    inline std::vector<size_type > getValues_ia() { return m_ia; }

    /* return the entire m_ja(vector type) */
    inline std::vector<size_type > getValues_ja() { return m_ja; }

    /* return the entire m_aa(vector type) */
    inline std::vector<value_type> getValues_aa() { return m_aa; }

    /* return the entire m_ia(array type) */
    inline size_type* getValues_ia_data() { return m_ia.data(); }

    /* return the entire m_ja(array type) */
    inline size_type* getValues_ja_data() { return m_ja.data(); }

    /* return the entire m_aa(array type) */
    inline value_type* getValues_aa_data() { return m_aa.data(); }

    /* get the size of the non-zero entries at the specific row index >>> unit testing is passed */
    inline size_type getNonZeroSizeAtOneRow(const size_type &rowIndex_) { return m_nnzPerRow[rowIndex_]; }

    /* get the sparse matrix size >>> unit testing is passed */
    inline size_type getSize() { return m_nnu; }

    /* get the number of non-zero entries of sparse matrix >>> unit testing is passed */
    inline size_type getNonZeroSize() { return m_nna; }

    /* get the value of ia of CSR sparse matrix at the location Idx */
    inline size_type  getValue_IA(const size_type &idx) { return m_ia[idx]; }

    /* get the value of ja at CSR sparse matrix at the location Idx */
    inline size_type  getValue_JA(const size_type &idx) { return m_ja[idx]; }

    /* get the value of aa at CSR sparse matrix at the location Idx */
    inline value_type getValue_AA(const size_type &idx) { return m_aa[idx]; }

    void setValue_L2G(size_type cpuID, size_type nproc, size_type &nnuG);

    /* get the local to global mapping */
    size_type getValue_L2G(const size_type &idx) { return m_L2G[idx]; }

    /* set the "ia" value of CSR sparse matrix at the location Idx */
    void setValue_IA(const size_type &idx, const size_type &value);

    /* set the "ja" value of CSR sparse matrix at the location Idx */
    void setValue_JA(const size_type &idx, const size_type &value);

    /* set the "aa" value of CSR sparse matrix at the location Idx */
    void setValue_AA(const size_type &idx, const value_type &value);

    /* set all the "ia" values of CSR sparse matrix */
    void setValues_IA(const size_type *ia);

    /* set all the "ja" values of CSR sparse matrix */
    void setValues_JA(const size_type *ja);

    /* set all the "aa" values of CSR sparse matrix */
    void setValues_AA(const value_type *aa);

    /* swap the CSR sparse matrix */
    void swap(std::vector<size_type> &ia, std::vector<size_type> &ja, std::vector<value_type> &aa);

    /* swap the local to global mapping */
    void swapL2G(std::vector<size_type> &L2G);

    /* swap value */
    void swap_value(int &a, int &b);
    void swap_value(T   &a, T   &b);

    /* sort the column index and the corresponding value */
    void SortColumnIndexAndValue();

    void quickSorting(const int left, const int right, int *ja, T* aa);

    void ZeroToOneBasedIndexing();
    void OneToZeroBasedIndexing();

    /* Output the matrix file */
    void ExportMatrixProfile();

    /* get the format name of sparse matrix */
    inline std::string getFormat() { return m_format; }

    /* compress the CSR sparse matrix >>> unit testing is passed */
    bool compressIndex();  // compress the sparse matrix
    bool compressIndexCSR();

    /* add the sparse matrix space >>> unit testing is passed */
    void addIndex(const size_type &rowIndex, const size_type &colIndex);

    /* add the value to the sparse matrix [CSR] >>> unit testing is passed */
    bool addValue(const size_type &rowIndex, const size_type &colIndex, const value_type &value);

    /* perform the sparse matrix-vector multiplication [CSR format] */
    void spmv(DenseVectorType &v, DenseVectorType &Av);  // spmv function for CSR format

    /* copy the sparse matrix data */
    bool copy(self_type &SparseMatrixCSR);

    /* output the profile of general CSR sparse matrix */
    void show();

 private:
    size_type               m_nnu;  // matrix size
    size_type               m_nna;  // the number of non-zero entries
    std::vector<size_type>  m_ia;   // the first non-zero entry point for each row
    std::vector<size_type>  m_ja;   // the column index of each non-zero entry
    std::vector<value_type> m_aa;   // the value of non-zero entry

    std::string m_format{"CSR"};  // sparse matrix format

    size_type *m_nnzPerRow;  // number of non-zero entries per row

    std::vector<size_type> m_L2G;  // The local to global mapping(To get the first row index of the distributed matrix)

    std::string m_MatrixName;  // The matrix coefficient name

    std::vector<std::set<size_type>> *matrix_column_index;

    bool IsCompressed = false;
};

}  // namespace container
}  // namespace inlakit


#endif  // IFSI2D_CONTAINER_INLAKIT_SPARSE_MATRIX_CSR_H_