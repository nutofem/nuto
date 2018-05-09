#include "BoostUnitTest.h"
#include <fakeit.hpp>
#include "nuto/mechanics/cell/SimpleAssembler.h"
#include "nuto/mechanics/cell/CellInterface.h"

fakeit::Mock<NuTo::CellInterface> MockCell(const NuTo::DofType& dof, Eigen::Vector3i numbering)
{
    fakeit::Mock<NuTo::CellInterface> cell;

    Eigen::VectorXi mockNumbering;
    NuTo::DofVector<double> mockGradient;
    NuTo::DofMatrix<double> mockHessian;

    mockNumbering = numbering;
    mockGradient[dof] = Eigen::Vector3d(11, 22, 33);
    mockHessian(dof, dof).resize(3, 3);
    mockHessian(dof, dof) << 11, 12, 13, 21, 22, 23, 31, 32, 33;

    Method(cell, DofNumbering) = mockNumbering;
    fakeit::When(OverloadedMethod(cell, Integrate, NuTo::DofVector<double>(NuTo::CellInterface::VectorFunction)))
            .AlwaysReturn(mockGradient);
    fakeit::When(OverloadedMethod(cell, Integrate, NuTo::DofMatrix<double>(NuTo::CellInterface::MatrixFunction)))
            .AlwaysReturn(mockHessian);

    return cell;
}
BOOST_AUTO_TEST_CASE(DefaultCtor)
{
    NuTo::SimpleAssembler assembler;

    NuTo::DofType d("0", 1);
    auto mockCell0 = MockCell(d, Eigen::Vector3i(0, 1, 2));
    BOOST_CHECK_THROW(assembler.BuildVector({mockCell0.get()}, {d}, NuTo::CellInterface::VectorFunction()),
                      NuTo::Exception);

    NuTo::DofInfo dofInfo;
    dofInfo.numIndependentDofs[d] = 3;
    dofInfo.numDependentDofs[d] = 0;

    assembler.SetDofInfo(dofInfo);
    BOOST_CHECK_NO_THROW(assembler.BuildVector({mockCell0.get()}, {d}, NuTo::CellInterface::VectorFunction()));

    NuTo::SimpleAssembler assembler2(dofInfo);
    BOOST_CHECK_NO_THROW(assembler2.BuildVector({mockCell0.get()}, {d}, NuTo::CellInterface::VectorFunction()));
}

NuTo::SimpleAssembler SetupAssembler(NuTo::DofType d)
{
    NuTo::DofInfo dofInfo;

    dofInfo.numIndependentDofs[d] = 3;
    dofInfo.numDependentDofs[d] = 2;

    return NuTo::SimpleAssembler(dofInfo);
}

BOOST_AUTO_TEST_CASE(AssemblerGradient)
{
    NuTo::DofType d("0", 1);
    NuTo::SimpleAssembler assembler = SetupAssembler(d);

    auto mockCell0 = MockCell(d, Eigen::Vector3i(0, 1, 2));
    auto mockCell1 = MockCell(d, Eigen::Vector3i(2, 3, 4));

    /*
     *    |  11  22  33  ||
     *  +          | 11  ||  22  33 |
     *                   ||
     *                   ||
     *       11  22  44  ||  22  33
     *      active       ||  dependent
     */

    NuTo::DofVector<double> gradient =
            assembler.BuildVector({mockCell0.get(), mockCell1.get()}, {d}, NuTo::CellInterface::VectorFunction());

    BoostUnitTest::CheckEigenMatrix(gradient[d], (Eigen::VectorXd(5)<< 11, 22, 44, 22, 33).finished());
}

BOOST_AUTO_TEST_CASE(AssemblerHessian)
{
    NuTo::DofType d("0", 1);
    NuTo::SimpleAssembler assembler = SetupAssembler(d);

    auto mockCell0 = MockCell(d, Eigen::Vector3i(0, 1, 2));
    auto mockCell1 = MockCell(d, Eigen::Vector3i(2, 3, 4));

    NuTo::DofInfo dofInfo;
    dofInfo.numDependentDofs[d] = 2;
    dofInfo.numIndependentDofs[d] = 3;
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

    auto hessian =
            assembler.BuildMatrix({mockCell0.get(), mockCell1.get()}, {d}, dofInfo, NuTo::CellInterface::MatrixFunction());

    Eigen::MatrixXd hessianE = (Eigen::MatrixXd(5,5) << 11, 12, 13,  0,  0,
                                                       21, 22, 23,  0,  0,
                                                       31, 32, 44, 12, 13,
                                                        0,  0, 21, 22, 23,
                                                        0,  0, 31, 32, 33).finished();


    BoostUnitTest::CheckEigenMatrix(Eigen::MatrixXd(hessian), hessianE);
}

BOOST_AUTO_TEST_CASE(AssemblerLumpedMass)
{
    NuTo::DofType d("0", 1);
    NuTo::SimpleAssembler assembler = SetupAssembler(d);

    auto mockCell0 = MockCell(d, Eigen::Vector3i(0, 1, 2));
    auto mockCell1 = MockCell(d, Eigen::Vector3i(2, 3, 4));


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

    NuTo::DofVector<double> massLumped = assembler.BuildDiagonallyLumpedMatrix({mockCell0.get(), mockCell1.get()}, {d},
                                                                             NuTo::CellInterface::MatrixFunction());

    double totalMass = 11 + 12 + 13 + 21 + 22 + 23 + 31 + 32 + 33;
    Eigen::Vector3d localMassDiagonal = (Eigen::Vector3d(11, 22, 33) * (totalMass / (11. + 22. + 33.)));
    Eigen::VectorXd globalMassDiagonal(5);
    globalMassDiagonal << localMassDiagonal[0], localMassDiagonal[1], localMassDiagonal[2] + localMassDiagonal[0],
            localMassDiagonal[1], localMassDiagonal[2];

    BoostUnitTest::CheckEigenMatrix(Eigen::MatrixXd(massLumped[d]), globalMassDiagonal);

}
