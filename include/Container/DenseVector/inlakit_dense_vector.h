/************************************************************************************/
//
//  FILENAME  : inlakit_dense_vector.h
//  PURPOSE   : Definition of "DENSE" vector class in INLAKit
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

#ifndef  SRC_INC_INLAKIT_CONTAINER_DENSE_VECTOR_H_
#define  SRC_INC_INLAKIT_CONTAINER_DENSE_VECTOR_H_

#include <inlakit_container_dense_matrix.h>
#include <inlakit_container_sparse_matrix.h>
#include <inlakit_container_sparse_matrix_coo.h>
#include <inlakit_container_sparse_matrix_csr.h>

// -----------------------------------------------------------

#include <vector_traits.h>
#include <matrix_traits.h>

// -----------------------------------------------------------

#include <inlakit_algorithm_transformation.h>

// -----------------------------------------------------------

#include <inlakit_calculator_dense_matrix.h>
#include <inlakit_calculator_sparse_matrix.h>

// -----------------------------------------------------------

#include <inlakit_direct_solver.h>

/************************************************************/

#include <iostream>
#include <iomanip>
#include <iterator>
#include <vector>
#include <algorithm>
#include <utility>
#include <string>

namespace inlakit   {
namespace container {

/* template class definition */
template<typename T>
class DenseVector {

    //-------------------------------------------------------

    template<typename S> friend class DenseVector;
    template<typename S> friend class DenseMatrix;
    template<typename S> friend class SparseMatrixCOO;
    template<typename S> friend class SparseMatrixCSR;

    // ------------------------------------------------------

    template<typename S> friend class MatrixTransformation;

    // ------------------------------------------------------

    template<typename S> friend class DenseMatrixCalc;
    template<typename S> friend class SparseMatrixCalc;

    //-------------------------------------------------------

    template<typename S> friend class DirectSolver;

    //-------------------------------------------------------

    using self_type = DenseVector<T>;
    using base_type = std::vector<T>;

 public:

    using size_type  = int;
    using value_type = T;

    using vector_size_type  = std::vector<size_type>;
    using vector_value_type = std::vector<value_type>;
    using data_type         = base_type;

    /* Default constructor */
    DenseVector();

    /* Constructor */
    DenseVector(const size_type &nnu);  // NOLINT

    /* Copy constructor */
    DenseVector(const self_type &DenVec);  // NOLINT

    /* Destructor */
    virtual ~DenseVector();

    /*--------------------------------------------------------------------*/
    /*                                                                    */
    /*                           Functions list                           */
    /*                                                                    */
    /*--------------------------------------------------------------------*/

    void create(const size_type &nnu);

    /* Perform the operation z = alpha * x + beta * y */
    void AxpBy(const value_type &alpha, const self_type &x, const value_type &beta, const self_type &y);

    /* Fill all the components in dense vector with 0 */
    void fill_with_zero();

    /* Fill all the componetns in dense vector with 1 */
    void fill_with_one();

    /* Fill all the components in dense vector with specific value */
    void fill_with_value(T& value);

    /* Copy the dense vector */
    void copy(const self_type &DenVec);

    /* Swap the dense vector */
    void swap(self_type &DenVec);

    /* Swap the single component in a dense vector */
    void swap(const size_type &i, const size_type &j);

    /* Clean the components in a dense vector */
    void clear();

    /* Show the components in a dense vector */
    void show();

    /* add the value "value_" to the dense vector */
    void addValue(const T &value);

    /* get all the values of dense vector */
    T* getValues();

    /* get the size of dense vector >>> (unit testing is passed) */
    size_type  getSize() { return m_nnu; }

    /* get the inner-product of the self-dense vector >>> (unit testing is passed) */
    T calcInnerProduct();

    /* get the inner-product with the dense vector DenseVector >>> (unit testing is passed) */
    T calcInnerProduct(const self_type &DenVec);

    /* get the 1-norm of dense vector >>> (unit testing is passed) */
    T calcNorm1();

    /* get the 2-norm of dense vector >>> (unit testing is passed) */
    T calcNorm2();

    /* get the infinity-norm of dense vector >>> (unit testing is passed) */
    T calcNormInf();

    /* get the maximal value of dense vector >>> (unit testing is passed) */
    T getMaxValue();

    /* get the minimal value of dense vector >>> (unit testing is passed) */
    T getMinValue();

    /* get the mean value of dense vector >>> (unit testing is passed) */
    T calcMeanValue();

    /* reset the dense vector size >>> (unit testing is passed) */
    void resetSize(const size_type &nnu);

    /* set the value to the dense vector at the location id_ >>> (unit testing is passed) */
    void setValue(const size_type &index, const T &value);

    /* get the value of dense vector at the location id_ */
    T getValue(const size_type &index) { return m_val[index] ; }

    /* get the sum of dense vector >>> (unit testing is passed) */
    T calcSum();

    /* get the sum of absolute dense vector >>> (unit testing is passed) */
    T calcAbsoluteSum();

    /* overload operator */

    DenseVector<T>  operator + (const DenseVector<T>& DenVec);
    DenseVector<T>  operator - (const DenseVector<T>& DenVec);
    DenseVector<T>& operator = (const DenseVector<T>& DenVec);

    DenseVector<T>& operator *= (const T scale_factor);

    friend DenseVector<T> operator*(DenseVector<T> DenVec, double scale_factor) {DenVec *= scale_factor; return DenVec;} // NOLINT
    friend DenseVector<T> operator*(double scale_factor, DenseVector<T> DenVec) {DenVec *= scale_factor; return DenVec;} // NOLINT

    /* sorting function : bubble sorting (ascending order) */
    void BubbleSorting();

    /* sorting function : bubble sorting (descending order) */
    void BubbleSorting_d();

    /* sorting function : selection sorting */
    void SelectionSorting();

    /* sorting function : insertion sorting */
    void InsertionSorting();

    /* sorting function : quick sorting */
    void QuickSorting();

    /* Read the file */
    //  bool read_ID0(size_type &nnu, size_type &cpuID, INLAKIT_STRING VectorName);

    /* Parallel partition */
    //  bool ParPartition(size_type &cpuID, size_type &nproc, INLAKIT_STRING VectorName);

    T operator() (size_type index);

          T& operator[] (size_type ij) {return m_val[ij];}
    const T& operator[] (size_type ij) const {return m_val[ij];}

 private:
    size_type      m_nnu = 0; /* vector size */
    size_type      m_capacity = 0; /* vector capacity */
    std::vector<T> m_val; /* vector values */
    std::string    m_VectorName; /* The vector coefficient name */
};

}  // namespace container
}  // namespace inlakit

#endif  // SRC_INC_INLAKIT_CONTAINER_DENSE_VECTOR_H_