#ifndef SRC_INC_MATRIX_TRAITS_H_
#define SRC_INC_MATRIX_TRAITS_H_

namespace inlakit {
namespace traits  {

template<typename T>
class MatrixTraits {
 public :
  using size_type         = typename T::size_type;
  using value_type        = typename T::value_type;
  using matrix_size_type  = typename T::matrix_size_type;
  using matrix_value_type = typename T::matrix_value_type;
};

}  // namespace traits
}  // namespace inlakit

#endif  // SRC_INC_MATRIX_TRAITS_H_