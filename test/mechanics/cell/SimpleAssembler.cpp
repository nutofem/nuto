#include "BoostUnitTest.h"
#include <fakeit.hpp>
#include "mechanics/cell/SimpleAssembler.h"

using namespace NuTo;

DofVector<double> MockGradient(DofType d)
{
    DofVector<double> v;
    v[d] = Eigen::Vector3d(11, 22, 33);
    return v;
}

DofMatrix<double> MockHessian(DofType d)
{
    DofMatrix<double> m;
    m(d, d).resize(3, 3);
    m(d, d) << 11, 12, 13, 21, 22, 23, 31, 32, 33;
    return m;
}

DofVector<int> MockNumbering(DofType d, Eigen::Vector3i numbering)
{
    DofVector<int> v;
    v[d] = numbering;
    return v;
}

DofInfo GetDofInfo(NuTo::DofType d)
{
    NuTo::DofInfo dofInfo;

    dofInfo.numIndependentDofs[d] = 3;
    dofInfo.numDependentDofs[d] = 2;

    return dofInfo;
}

BOOST_AUTO_TEST_CASE(AssemblerGradient)
{
    NuTo::DofType d("0", 1);
    NuTo::SimpleAssembler assembler;

    /*
     *    |  11  22  33  ||
     *  +          | 11  ||  22  33 |
     *                   ||
     *                   ||
     *       11  22  44  ||  22  33
     *      active       ||  dependent
     */
    VectorEntry entry1 = {MockGradient(d), MockNumbering(d, Eigen::Vector3i(0, 1, 2))};
    VectorEntry entry2 = {MockGradient(d), MockNumbering(d, Eigen::Vector3i(2, 3, 4))};

    fakeit::Mock<DofVectorGenerator> entries;

    fakeit::When(Method(entries, Dofs)).AlwaysReturn({d});
    fakeit::When(Method(entries, Next)).AlwaysReturn();

    // Why true, true, true, true, false? Well, this is pure trial and error.
    // boost::iterator_facade calls IsValid() _somehow/somewhen_. But at some
    // point, we have to return false to stop the generator.
    fakeit::When(Method(entries, IsValid)).Return(true, true, true, true, false);
    fakeit::When(Method(entries, Get)).Return(entry1).Return(entry2);

    // no dof numbering set
    BOOST_CHECK_THROW(assembler.BuildVector(entries.get()), Exception);

    assembler.SetDofInfo(GetDofInfo(d));
    NuTo::GlobalDofVector gradient = assembler.BuildVector(entries.get());

    BoostUnitTest::CheckEigenMatrix(gradient.J[d], Eigen::Vector3d(11, 22, 44));
    BoostUnitTest::CheckEigenMatrix(gradient.K[d], Eigen::Vector2d(22, 33));
}

BOOST_AUTO_TEST_CASE(AssemblerHessian)
{
    NuTo::DofType d("0", 1);
    NuTo::SimpleAssembler assembler(GetDofInfo(d));

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

    MatrixEntry entry1 = {MockHessian(d), MockNumbering(d, Eigen::Vector3i(0, 1, 2))};
    MatrixEntry entry2 = {MockHessian(d), MockNumbering(d, Eigen::Vector3i(2, 3, 4))};

    fakeit::Mock<DofMatrixGenerator> entries;

    fakeit::When(Method(entries, Dofs)).AlwaysReturn({d});
    fakeit::When(Method(entries, Next)).AlwaysReturn();
    fakeit::When(Method(entries, IsValid)).Return(true, true, true, true, false);
    fakeit::When(Method(entries, Get)).Return(entry1).Return(entry2);

    NuTo::GlobalDofMatrixSparse hessian = assembler.BuildMatrix(entries.get());
    Eigen::Matrix3d JJ = (Eigen::Matrix3d() << 11, 12, 13, 21, 22, 23, 31, 32, 44).finished();
    Eigen::Matrix2d KK = (Eigen::Matrix2d() << 22, 23, 32, 33).finished();
    Eigen::MatrixXd JK = (Eigen::MatrixXd(3, 2) << 0, 0, 0, 0, 12, 13).finished();
    Eigen::MatrixXd KJ = (Eigen::MatrixXd(2, 3) << 0, 0, 21, 0, 0, 31).finished();

    BoostUnitTest::CheckEigenMatrix(Eigen::MatrixXd(hessian.JJ(d, d)), JJ);
    BoostUnitTest::CheckEigenMatrix(Eigen::MatrixXd(hessian.JK(d, d)), JK);
    BoostUnitTest::CheckEigenMatrix(Eigen::MatrixXd(hessian.KJ(d, d)), KJ);
    BoostUnitTest::CheckEigenMatrix(Eigen::MatrixXd(hessian.KK(d, d)), KK);
}

BOOST_AUTO_TEST_CASE(AssemblerLumpedMass)
{
    // NuTo::DofType d("0", 1);
    // NuTo::SimpleAssembler assembler = SetupAssembler(d);

    // auto mockCell0 = MockCell(d, Eigen::Vector3i(0, 1, 2));
    // auto mockCell1 = MockCell(d, Eigen::Vector3i(2, 3, 4));


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

    // NuTo::GlobalDofVector massLumped = assembler.BuildDiagonallyLumpedMatrix({mockCell0.get(), mockCell1.get()}, {d},
    //                                                                         NuTo::CellInterface::MatrixFunction());
    //
    // double totalMass = 11 + 12 + 13 + 21 + 22 + 23 + 31 + 32 + 33;
    // Eigen::Vector3d localMassDiagonal = (Eigen::Vector3d(11, 22, 33) * (totalMass / (11. + 22. + 33.)));
    // Eigen::VectorXd globalMassDiagonal(5);
    // globalMassDiagonal << localMassDiagonal[0], localMassDiagonal[1], localMassDiagonal[2] + localMassDiagonal[0],
    //        localMassDiagonal[1], localMassDiagonal[2];
    //
    // Eigen::VectorXd J = globalMassDiagonal.head(3);
    // Eigen::VectorXd K = globalMassDiagonal.tail(2);
    //
    // BoostUnitTest::CheckEigenMatrix(Eigen::MatrixXd(massLumped.J[d]), J);
    // BoostUnitTest::CheckEigenMatrix(Eigen::MatrixXd(massLumped.K[d]), K);
}
