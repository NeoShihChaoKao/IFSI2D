#ifndef SRC_INC_VECTOR_TRAITS_H_
#define SRC_INC_VECTOR_TRAITS_H_

namespace inlakit {
namespace traits  {

template<typename T>
class VectorTraits {
 public :
  using size_type         = typename T::size_type;
  using value_type        = typename T::value_type;
  using vector_size_type  = typename T::vector_size_type;
  using vector_value_type = typename T::vector_value_type;
};

}  // namespace traits
}  // namespace inlakit

#endif  // SRC_INC_VECTOR_TRAITS_H_