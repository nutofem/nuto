//
// Created by Thomas Titscher on 1/10/17.
//
#pragma once
#include <iostream>
#include <eigen3/Eigen/Core>
#include "base/Exception.h" // base/Exception is wrapped for python, Exception not.

namespace NuTo
{
namespace Test
{

//! @brief test the python (numpy interface) of dense matrix types
class NumpyDense
{
public:
    //! @brief checks correctness of double rMatrix (from python) and returns a matrix to python
    //! @param rMatrix np.array([[1., 2.],[3., 4.], [5., 6.]])
    //! @return np.array([[2., 4.],[6., 8.], [10., 12.])
    static Eigen::MatrixXd CheckMatrixXd(const Eigen::MatrixXd& rMatrix)
    {
        return CheckMatrixX<double>(rMatrix);
    }

    //! @brief checks correctness of integer rMatrix (from python) and returns a matrix to python
    //! @param rMatrix np.array([[1, 2],[3, 4], [5., 6.]])
    //! @return np.array([[2, 4],[6, 8], [10, 12])
    static Eigen::MatrixXi CheckMatrixXi(const Eigen::MatrixXi& rMatrix)
    {
        return CheckMatrixX<int>(rMatrix);
    }

private:
    //! @brief checks correctness of rMatrix (from python) and returns a matrix to python
    //! @param rMatrix np.array([[1, 2],[3, 4], [5, 6]])
    //! @return np.array([[2, 4],[6, 8], [10, 12])
    template <typename T>
    static Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
    CheckMatrixX(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& rMatrix)
    {
        if (rMatrix.rows() != 3)
            throw Exception(__PRETTY_FUNCTION__, "rMatrix.rows = " + std::to_string(rMatrix.rows()) + ", should be 3");
        if (rMatrix.cols() != 2)
            throw Exception(__PRETTY_FUNCTION__, "rMatrix.cols = " + std::to_string(rMatrix.cols()) + ", should be 2");

        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> correct(3, 2);
        correct << 1, 2, 3, 4, 5, 6;

        if ((correct - rMatrix).cwiseAbs().maxCoeff() > 1.e-10)
        {
            std::cout << "got \n" << rMatrix << std::endl;
            std::cout << "expected \n" << correct << std::endl;
            throw Exception(__PRETTY_FUNCTION__, "Incorrect values");
        }
        return correct * 2;
    }
};


} // namespace Test
} // namespace NuTo
