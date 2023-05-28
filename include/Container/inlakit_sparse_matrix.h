#ifndef IFSI2D_CONTAINER_INLAKIT_SPARSE_MATRIX_H_
#define IFSI2D_CONTAINER_INLAKIT_SPARSE_MATRIX_H_

#include <inlakit_container_sparse_matrix.h>
#include <inlakit_container_sparse_matrix_csr.h>
#include <inlakit_container_sparse_matrix_coo.h>
#include <inlakit_container_sparse_matrix_type.h>

//#include <inlakit_calculator_sparse_matrix.h>

// ---------------------------------------------------------------

namespace inlakit   {
namespace container {

template <typename T, inlakit::container::SparseMatrixType SpaMatType>
struct _SparseMatrix;

/* ============ COO format ============ */

template <typename T>
struct _SparseMatrix < T, inlakit::container::SparseMatrixType::COO > {
  using type = inlakit::container::SparseMatrixCOO<T>;
};

/* ============ CSR format ============ */

template < typename T>
struct _SparseMatrix < T, inlakit::container::SparseMatrixType::CSR > {
  using type = inlakit::container::SparseMatrixCSR<T>;
};

/* =============== Wrapper =============== */
template < typename T, inlakit::container::SparseMatrixType SpaMatType >
using SparseMatrix = typename _SparseMatrix <T, SpaMatType>::type;
/* ======================================= */

}  // namespace container
}  // namespace inlakit

#endif  // IFSI2D_CONTAINER_INLAKIT_SPARSE_MATRIX_H_