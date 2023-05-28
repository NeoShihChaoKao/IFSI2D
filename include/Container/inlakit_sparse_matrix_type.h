#ifndef IFSI2D_CONTAINER_INLAKIT_SPARSE_MATRIX_TYPE_H_
#define IFSI2D_CONTAINER_INLAKIT_SPARSE_MATRIX_TYPE_H_

#include <string>
#include <cstdint>

namespace inlakit   {
namespace container {

  /* Define the sparse matrix type */
  enum class SparseMatrixType {
    COO, CSR, ELL
  };

}  // namespace container
}  // namespace inlakit

#endif  // IFSI2D_CONTAINER_INLAKIT_SPARSE_MATRIX_TYPE_H_