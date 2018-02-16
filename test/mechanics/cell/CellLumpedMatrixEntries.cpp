#include "BoostUnitTest.h"
#include <fakeit.hpp>
#include "mechanics/cell/CellLumpedMatrixEntries.h"

using namespace NuTo;

BOOST_AUTO_TEST_CASE(AssemblerLumpedMass)
{
    NuTo::DofType d("0", 1);

    fakeit::Mock<CellInterface> mockCell;
    Method(mockCell, DofNumbering) = Eigen::Vector3i(1, 2, 3);

    DofMatrix<double> m;
    m(d, d).resize(3, 3);
    m(d, d) << 11, 12, 13, 21, 22, 23, 31, 32, 33;
    OverloadedMethod(mockCell, Integrate, DofMatrix<double>(CellInterface::MatrixFunction)) = m;

    Group<CellInterface> g{mockCell.get()};
    CellLumpedMatrixEntries entries(g, CellInterface::MatrixFunction(), {d});


    /*
     *   11    12    13 ||
     *                  ||
     *   21    22    23 ||             active J
     *                  ||
     *   31    32    44 || 12    13
     *   ===========================
     *               21 || 22    23
     *                  ||             dependent K
     *               31 || 32    33
     *
     *     active J         dependent K
     */


    double totalMass = 11 + 12 + 13 + 21 + 22 + 23 + 31 + 32 + 33;
    Eigen::Vector3d localMassDiagonal = (Eigen::Vector3d(11, 22, 33) * (totalMass / (11. + 22. + 33.)));
    BoostUnitTest::CheckEigenMatrix(entries.Get().first[d], localMassDiagonal);
}
