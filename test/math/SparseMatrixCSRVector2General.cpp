#include "BoostUnitTest.h"

#include <eigen3/Eigen/Core>
#include "math/SparseMatrixCSRVector2General.h"

using namespace NuTo;

BOOST_AUTO_TEST_CASE(concatenate)
{
    auto B = SparseMatrixCSRVector2General<double>(4, 4);

    B.AddValue(1, 1, 3.0);
    B.AddValue(0, 1, 8.0);
    B.AddValue(0, 0, 2.0);
    B.AddValue(1, 3, 3.0);
    B.AddValue(1, 2, 5.0);
    B.AddValue(3, 2, 1.0);
    B.AddValue(3, 3, 7.0);
    B.AddValue(3, 0, 9.0);
    B.AddValue(3, 1, 4.0);

    auto C = SparseMatrixCSRVector2General<double>(4,1);

    C.AddValue(0, 0, 4.0);
    C.AddValue(1, 0, 2.0);
    C.AddValue(3, 0, 3.0);

    auto D = SparseMatrixCSRVector2General<double>(1,5);

    D.AddValue(0, 1, 4.0);
    D.AddValue(0, 4, 1.0);

    B.ConcatenateColumns(C);
    B.ConcatenateRows(D);

    Eigen::MatrixXd B_FullRef(5, 5);
    B_FullRef << 2.0, 8.0, 0.0, 0.0, 4.0,
                  0.0, 3.0, 5.0, 3.0, 2.0,
                  0.0, 0.0, 0.0, 0.0, 0.0,
                  9.0, 4.0, 1.0, 7.0, 3.0,
                  0.0, 4.0, 0.0, 0.0, 1.0;

    auto B_Full = B.ConvertToFullMatrix();
    BOOST_CHECK_SMALL((B_FullRef - B_Full).maxCoeff(), 1e-12);
}
