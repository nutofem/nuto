#include "BoostUnitTest.h"
#include "nuto/mechanics/cell/Assembler.h"

using namespace NuTo;

BOOST_AUTO_TEST_CASE(VectorAssembly)
{
    NuTo::DofType d("d", 1);
    DofVector<double> v;
    v[d] = Eigen::Vector3d(1, 2, 3);

    DofVector<int> numbering0, numbering1;
    numbering0[d] = Eigen::Vector3i(0, 1, 2);
    numbering1[d] = Eigen::Vector3i(2, 3, 4);

    DofContainer<int> size;
    size[d] = 5;
    VectorAssembler vectorAssembler(size);
    vectorAssembler.Add(v, numbering0);
    vectorAssembler.Add(v, numbering1);

    BoostUnitTest::CheckVector(vectorAssembler.Get()[d], std::vector<double>{1, 2, 3 + 1, 2, 3}, 5);

    vectorAssembler.Reset();
    BoostUnitTest::CheckVector(vectorAssembler.Get()[d], std::vector<double>{0, 0, 0, 0, 0}, 5);
}


BOOST_AUTO_TEST_CASE(MatrixAssembly)
{
    NuTo::DofType d("d", 1);
    DofMatrix<double> m;
    m(d, d).resize(3, 3);
    m(d, d) << 11, 12, 13, 21, 22, 23, 31, 32, 33;

    DofVector<int> numbering0, numbering1;
    numbering0[d] = Eigen::Vector3i(0, 1, 2);
    numbering1[d] = Eigen::Vector3i(2, 3, 4);

    DofContainer<int> size;
    size[d] = 5;
    MatrixAssembler matrixAssembler(size);
    matrixAssembler.Add(m, numbering0);
    matrixAssembler.Add(m, numbering1);

    matrixAssembler.Finish();

    Eigen::MatrixXd expected(5, 5);
    expected.setZero();
    expected.block<3, 3>(0, 0) = m(d, d);
    expected.block<3, 3>(2, 2) = m(d, d);
    expected(2, 2) = 44;

    BoostUnitTest::CheckEigenMatrix(Eigen::MatrixXd(matrixAssembler.Get()(d, d)), expected);

    matrixAssembler.Reset();
    BoostUnitTest::CheckEigenMatrix(Eigen::MatrixXd(matrixAssembler.Get()(d, d)), Eigen::MatrixXd::Zero(5, 5));
}
