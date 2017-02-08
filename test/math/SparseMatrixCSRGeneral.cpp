#include "BoostUnitTest.h"

#include <eigen3/Eigen/Core>
#include "math/SparseMatrixCSRGeneral.h"

using namespace NuTo;

BOOST_AUTO_TEST_CASE(insertValues)
{
    Eigen::MatrixXd A_full_ref(5, 5);
    A_full_ref << 1.0, -1.0, 0.0, -3.0,  0.0,
                 -2.0,  5.0, 0.0,  0.0,  0.0,
                  0.0,  0.0, 4.0,  6.0,  4.0,
                 -4.0,  0.0, 2.0,  7.0,  0.0,
                  0.0,  8.0, 0.0,  0.0, -5.0;

    auto A_sparse = SparseMatrixCSRGeneral<double>(5, 5);
    A_sparse.AddValue(0, 0,  1.0);
    A_sparse.AddValue(1, 0, -2.0);
    A_sparse.AddValue(3, 0, -4.0);
    A_sparse.AddValue(0, 1, -1.0);
    A_sparse.AddValue(1, 1,  5.0);
    A_sparse.AddValue(4, 1,  8.0);
    A_sparse.AddValue(2, 2,  4.0);
    A_sparse.AddValue(3, 2,  2.0);
    A_sparse.AddValue(0, 3, -3.0);
    A_sparse.AddValue(2, 3,  6.0);
    A_sparse.AddValue(3, 3,  7.0);
    A_sparse.AddValue(2, 4,  4.0);
    A_sparse.AddValue(4, 4, -5.0);
    
    auto A_full = A_sparse.ConvertToFullMatrix();
    BOOST_CHECK_SMALL((A_full_ref - A_full).maxCoeff(), 1e-12);
}
