#include "BoostUnitTest.h"
#include <fakeit.hpp>
#include "mechanics/cell/SimpleAssember.h"
#include "mechanics/cell/CellInterface.h"

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
            .Return(mockGradient);
    fakeit::When(OverloadedMethod(cell, Integrate, NuTo::DofMatrix<double>(NuTo::CellInterface::MatrixFunction)))
            .Return(mockHessian);

    return cell;
}

NuTo::SimpleAssembler SetupAssembler(const NuTo::DofType& d)
{
    NuTo::DofContainer<int> numIndependent;
    NuTo::DofContainer<int> numDependent;

    numIndependent[d] = 3;
    numDependent[d] = 2;

    return NuTo::SimpleAssembler(numIndependent, numDependent);
}

BOOST_AUTO_TEST_CASE(AssemberGradient)
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

    NuTo::GlobalDofVector gradient =
            assembler.BuildVector({mockCell0.get(), mockCell1.get()}, {d}, NuTo::CellInterface::VectorFunction());

    BoostUnitTest::CheckEigenMatrix(gradient.J[d], Eigen::Vector3d(11, 22, 44));
    BoostUnitTest::CheckEigenMatrix(gradient.K[d], Eigen::Vector2d(22, 33));
}

BOOST_AUTO_TEST_CASE(AssemberHessian)
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

    NuTo::GlobalDofMatrixSparse hessian =
            assembler.BuildMatrix({mockCell0.get(), mockCell1.get()}, {d}, NuTo::CellInterface::MatrixFunction());
    Eigen::Matrix3d JJ = (Eigen::Matrix3d() << 11, 12, 13, 21, 22, 23, 31, 32, 44).finished();
    Eigen::Matrix2d KK = (Eigen::Matrix2d() << 22, 23, 32, 33).finished();
    Eigen::MatrixXd JK = (Eigen::MatrixXd(3, 2) << 0, 0, 0, 0, 12, 13).finished();
    Eigen::MatrixXd KJ = (Eigen::MatrixXd(2, 3) << 0, 0, 21, 0, 0, 31).finished();

    BoostUnitTest::CheckEigenMatrix(Eigen::MatrixXd(hessian.JJ(d, d)), JJ);
    BoostUnitTest::CheckEigenMatrix(Eigen::MatrixXd(hessian.JK(d, d)), JK);
    BoostUnitTest::CheckEigenMatrix(Eigen::MatrixXd(hessian.KJ(d, d)), KJ);
    BoostUnitTest::CheckEigenMatrix(Eigen::MatrixXd(hessian.KK(d, d)), KK);
}
