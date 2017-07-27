#include "BoostUnitTest.h"

#include "base/Timer.h"
#include <iostream>

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>


#include "math/SparseMatrixCSRVector2General.h"


#include "mechanics/dofSubMatrixStorage/BlockScalar.h"
#include "mechanics/dofSubMatrixStorage/BlockFullVector.h"
#include "mechanics/dofSubMatrixStorage/BlockFullMatrix.h"
#include "mechanics/dofSubMatrixStorage/BlockSparseMatrix.h"
#include "mechanics/structures/StructureOutputBlockMatrix.h"
#include "mechanics/elements/ElementOutputBlockMatrixDouble.h"
#include "mechanics/dofSubMatrixStorage/DofStatus.h"
#include "base/Exception.h"
#include "mechanics/nodes/NodeEnum.h"


//! @brief BlockFullVectorTest [BVT]
//! @remark tests are done with int-vectors for easier comparison without epsilons...
//! @remark assertion: NuTo::FullVector calculations are correct.
BOOST_AUTO_TEST_CASE(BlockFullVector)
{
    NuTo::DofStatus s;
    s.SetDofTypes({NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::TEMPERATURE});
    s.SetActiveDofTypes({NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::TEMPERATURE});

    NuTo::BlockFullVector<int> v1(s), v2(s);
    v1[NuTo::Node::eDof::DISPLACEMENTS] = Eigen::Vector3i({1, 1, 1});
    v2[NuTo::Node::eDof::DISPLACEMENTS] = Eigen::Vector3i({2, 2, 2});

    v1[NuTo::Node::eDof::TEMPERATURE] = Eigen::Vector2i({10, 10});
    v2[NuTo::Node::eDof::TEMPERATURE] = Eigen::Vector2i({20, 20});

    // Export
    Eigen::VectorXi ev1 = v1.Export();
    Eigen::VectorXi ev2 = v2.Export();
    BOOST_CHECK_EQUAL(ev1.rows(), 5);
    BOOST_CHECK_EQUAL(ev1.sum(), 23);
    BOOST_CHECK_EQUAL(ev2.sum(), 46);

    // Addition
    BOOST_CHECK(v1 + v2 == v2 + v1);
    BOOST_CHECK((v1 + v2).Export() == ev2 + ev1);

    // Subtraction
    BOOST_CHECK(v1 - v1 == v2 - v2);
    BOOST_CHECK((v1 - v2).Export() == ev1 - ev2);

    // ScalarMultiplication
    BOOST_CHECK(v1 + v1 == v1 * 2);
    BOOST_CHECK((v1 + v1).Export() == ev1 * 2);

    // Chaining
    NuTo::BlockFullVector<int> result1(s), result2(s);
    result1 = v1 * 3 + v1 - v2 - v2 + v2 * 42;

    BOOST_TEST_MESSAGE("Chained operations done.");

    result2 = v1;
    result2 *= 3;
    result2 += v1;
    result2 -= v2;
    result2 -= v2;
    result2 += v2 * 42;

    Eigen::VectorXi eresult = ev1 * 3 + ev1 - ev2 - ev2 + ev2 * 42;
    BOOST_CHECK(result1 == result2);
    BOOST_CHECK(result1.Export() == eresult);

    // ActiveDofTypes
    s.SetActiveDofTypes({NuTo::Node::eDof::DISPLACEMENTS});
    BOOST_CHECK_EQUAL(v1.Export().rows(), 3);
    v1.Info();
    s.SetActiveDofTypes({NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::TEMPERATURE});

    // import
    // allocate BlockFullVector with the right dimensions.
    NuTo::BlockFullVector<int> imported(v1);

    // create new random vector with the global dimensions of v1
    Eigen::VectorXi toImport = Eigen::VectorXi::Random(v1.GetNumRows());

    // import and check if equal.
    imported.Import(toImport);

    BoostUnitTest::CheckEigenMatrix(toImport, imported.Export());
}


//! @brief BlockFullMatrixTest [BMT]
//! @remark just test the access operators
BOOST_AUTO_TEST_CASE(BlockFullMatrix)
{
    NuTo::DofStatus s;
    s.SetDofTypes({NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::TEMPERATURE});
    s.SetActiveDofTypes({NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::TEMPERATURE});

    NuTo::BlockFullMatrix<double> m(s);

    int numD = 3;
    int numT = 2;

    m(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::DISPLACEMENTS) = Eigen::MatrixXd::Random(numD, numD);

    m(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::TEMPERATURE).resize(numD, numT);
    m(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::TEMPERATURE) << 1, 2, 3, 1, 2, 3;

    m(NuTo::Node::eDof::TEMPERATURE, NuTo::Node::eDof::DISPLACEMENTS) =
            m(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::TEMPERATURE).transpose();

    auto& matrixTT = m(NuTo::Node::eDof::TEMPERATURE, NuTo::Node::eDof::TEMPERATURE);
    matrixTT.resize(numT, numT);

    m.Info();
    BOOST_TEST_MESSAGE("" << m);
    BOOST_CHECK_NO_THROW(m.CheckDimensions());
}

//! @brief BlockSparseMatrixTest
//! @remark tests access operations, vector*matrix operation,
BOOST_AUTO_TEST_CASE(BlockSparseMatrix)
{
    NuTo::DofStatus s;
    s.SetDofTypes({NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::TEMPERATURE});
    s.SetActiveDofTypes({NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::TEMPERATURE});
    s.SetIsSymmetric(NuTo::Node::eDof::DISPLACEMENTS, true);

    NuTo::BlockSparseMatrix m(s);

    size_t numD = 5;
    size_t numT = 2;
    double density = 1;

    m(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::DISPLACEMENTS) =
            NuTo::SparseMatrixCSRVector2Symmetric<double>::Random(numD, numD, density);
    m(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::TEMPERATURE) =
            NuTo::SparseMatrixCSRVector2General<double>::Random(numD, numT, density);
    m(NuTo::Node::eDof::TEMPERATURE, NuTo::Node::eDof::DISPLACEMENTS) =
            NuTo::SparseMatrixCSRVector2General<double>::Random(numT, numD, density);
    m(NuTo::Node::eDof::TEMPERATURE, NuTo::Node::eDof::TEMPERATURE) =
            NuTo::SparseMatrixCSRVector2General<double>::Random(numT, numT, density);
    BOOST_CHECK_NO_THROW(m.CheckDimensions());
    m.Info();


    NuTo::BlockFullVector<double> v(s);
    v[NuTo::Node::eDof::DISPLACEMENTS] = Eigen::VectorXd::Random(numD);
    v[NuTo::Node::eDof::TEMPERATURE] = Eigen::VectorXd::Random(numT);

    auto result = m * v * 4 + m * v;

    Eigen::MatrixXd exportReference = m.ExportToFullMatrix();
    Eigen::MatrixXd exportCSRVector2 = m.ExportToCSRVector2General().ConvertToFullMatrix();
    Eigen::MatrixXd exportCSR = m.ExportToCSRGeneral().ConvertToFullMatrix();

    BoostUnitTest::CheckEigenMatrix(exportCSRVector2, exportReference);
    BoostUnitTest::CheckEigenMatrix(exportCSR, exportReference);

    s.SetActiveDofTypes({NuTo::Node::eDof::DISPLACEMENTS});

    auto CSR = m.ExportToCSR();
    BOOST_CHECK(CSR->IsSymmetric());
    
    Eigen::MatrixXd exportCSRSymm = CSR->ConvertToFullMatrix();
    BoostUnitTest::CheckEigenMatrix(exportCSRSymm, m.ExportToFullMatrix());
}

//! @brief StructureOutputBlockMatrixTest
//! @remark allocates random sparse block matrices
void StructureOutputBlockMatrixTestGeneral(int rNumDAct, int rNumTAct, int rNumDDep, int rNumTDep, double rDensity)
{
    NuTo::Timer timer("StructureOutputBlockMatrixTestGeneral::DefineRandomMatrices");
    NuTo::DofStatus s;
    s.SetDofTypes({NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::TEMPERATURE});
    s.SetActiveDofTypes({NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::TEMPERATURE});


    NuTo::StructureOutputBlockMatrix BM4(s);

    BM4.JJ(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::DISPLACEMENTS) =
            NuTo::SparseMatrixCSRVector2General<double>::Random(rNumDAct, rNumDAct, rDensity);
    BM4.JJ(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::TEMPERATURE) =
            NuTo::SparseMatrixCSRVector2General<double>::Random(rNumDAct, rNumTAct, rDensity);
    BM4.JJ(NuTo::Node::eDof::TEMPERATURE, NuTo::Node::eDof::DISPLACEMENTS) =
            NuTo::SparseMatrixCSRVector2General<double>::Random(rNumTAct, rNumDAct, rDensity);
    BM4.JJ(NuTo::Node::eDof::TEMPERATURE, NuTo::Node::eDof::TEMPERATURE) =
            NuTo::SparseMatrixCSRVector2General<double>::Random(rNumTAct, rNumTAct, rDensity);

    BM4.JK(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::DISPLACEMENTS) =
            NuTo::SparseMatrixCSRVector2General<double>::Random(rNumDAct, rNumDDep, rDensity);
    BM4.JK(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::TEMPERATURE) =
            NuTo::SparseMatrixCSRVector2General<double>::Random(rNumDAct, rNumTDep, rDensity);
    BM4.JK(NuTo::Node::eDof::TEMPERATURE, NuTo::Node::eDof::DISPLACEMENTS) =
            NuTo::SparseMatrixCSRVector2General<double>::Random(rNumTAct, rNumDDep, rDensity);
    BM4.JK(NuTo::Node::eDof::TEMPERATURE, NuTo::Node::eDof::TEMPERATURE) =
            NuTo::SparseMatrixCSRVector2General<double>::Random(rNumTAct, rNumTDep, rDensity);

    BM4.KJ(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::DISPLACEMENTS) =
            NuTo::SparseMatrixCSRVector2General<double>::Random(rNumDDep, rNumDAct, rDensity);
    BM4.KJ(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::TEMPERATURE) =
            NuTo::SparseMatrixCSRVector2General<double>::Random(rNumDDep, rNumTAct, rDensity);
    BM4.KJ(NuTo::Node::eDof::TEMPERATURE, NuTo::Node::eDof::DISPLACEMENTS) =
            NuTo::SparseMatrixCSRVector2General<double>::Random(rNumTDep, rNumDAct, rDensity);
    BM4.KJ(NuTo::Node::eDof::TEMPERATURE, NuTo::Node::eDof::TEMPERATURE) =
            NuTo::SparseMatrixCSRVector2General<double>::Random(rNumTDep, rNumTAct, rDensity);

    BM4.KK(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::DISPLACEMENTS) =
            NuTo::SparseMatrixCSRVector2General<double>::Random(rNumDDep, rNumDDep, rDensity);
    BM4.KK(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::TEMPERATURE) =
            NuTo::SparseMatrixCSRVector2General<double>::Random(rNumDDep, rNumTDep, rDensity);
    BM4.KK(NuTo::Node::eDof::TEMPERATURE, NuTo::Node::eDof::DISPLACEMENTS) =
            NuTo::SparseMatrixCSRVector2General<double>::Random(rNumTDep, rNumDDep, rDensity);
    BM4.KK(NuTo::Node::eDof::TEMPERATURE, NuTo::Node::eDof::TEMPERATURE) =
            NuTo::SparseMatrixCSRVector2General<double>::Random(rNumTDep, rNumTDep, rDensity);

    BM4.CheckDimensions();


    NuTo::BlockSparseMatrix CMatrix(s, false);
    CMatrix(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::DISPLACEMENTS) =
            NuTo::SparseMatrixCSRVector2General<double>::Random(rNumDDep, rNumDAct, rDensity);
    CMatrix(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::TEMPERATURE).Resize(rNumDDep, rNumTAct);
    CMatrix(NuTo::Node::eDof::TEMPERATURE, NuTo::Node::eDof::DISPLACEMENTS).Resize(rNumTDep, rNumDAct);
    CMatrix(NuTo::Node::eDof::TEMPERATURE, NuTo::Node::eDof::TEMPERATURE) =
            NuTo::SparseMatrixCSRVector2General<double>::Random(rNumTDep, rNumTAct, rDensity);

    NuTo::BlockSparseMatrix hessian(s, true);
    hessian = BM4.JJ; // just to get every sub matrix ...
    hessian.SetZero(); // ... to the right dimensions


    auto CMatrixExp = CMatrix.ExportToCSRVector2General();
    auto JJExp = BM4.JJ.ExportToCSRVector2General();
    auto JKExp = BM4.JK.ExportToCSRVector2General();
    auto KJExp = BM4.KJ.ExportToCSRVector2General();
    auto KKExp = BM4.KK.ExportToCSRVector2General();
    auto hessianExp = JJExp;
    hessianExp.SetZeroEntries();

    timer.Reset("StructureOutputBlockMatrixTestGeneral::ApplyCMatrix_NEW");

    BM4.ApplyCMatrixScal(hessian, CMatrix, 1);

    timer.Reset("StructureOutputBlockMatrixTestGeneral::ApplyCMatrix_OLD");

    hessianExp += JJExp;
    hessianExp += (CMatrixExp.Transpose() * KKExp * CMatrixExp - CMatrixExp.Transpose() * KJExp - JKExp * CMatrixExp);

    timer.Reset("StructureOutputBlockMatrixTestGeneral::cleanup");


    auto diff = hessian.ExportToCSRVector2General() - hessianExp;
    double tolerance = 1.e-8;
    double diffMaxMin = diff.Max() - diff.Min();
    if (diffMaxMin > tolerance)
    {
        std::cout << diffMaxMin << std::endl;
        throw NuTo::Exception("[StructureOutputBlockMatrixTestGeneral] ApplyCMatrix incorrect.");
    }


    BM4.AddScal(BM4, 2.7);

    diff = (JJExp * 3.7 - BM4.JJ.ExportToCSRVector2General());
    if (diff.Max() - diff.Min() > tolerance)
        throw NuTo::Exception("[StructureOutputBlockMatrixTestGeneral] AddScal for JJ incorrect.");

    if (rNumDDep > 0)
    {
        diff = (JKExp * 3.7 - BM4.JK.ExportToCSRVector2General());
        if (diff.Max() - diff.Min() > tolerance)
            throw NuTo::Exception("[StructureOutputBlockMatrixTestGeneral] AddScal for JK incorrect.");

        diff = (KJExp * 3.7 - BM4.KJ.ExportToCSRVector2General());
        if (diff.Max() - diff.Min() > tolerance)
            throw NuTo::Exception("[StructureOutputBlockMatrixTestGeneral] AddScal for KJ incorrect.");

        diff = (KKExp * 3.7 - BM4.KK.ExportToCSRVector2General());
        if (diff.Max() - diff.Min() > tolerance)
            throw NuTo::Exception("[StructureOutputBlockMatrixTestGeneral] AddScal for KK incorrect.");
    }
}

BOOST_AUTO_TEST_CASE(SparseGeneral)
{
    StructureOutputBlockMatrixTestGeneral(10, 8, 0, 0, 1); // cmat == 0
    StructureOutputBlockMatrixTestGeneral(10, 8, 4, 2, 1);
};

//! @brief StructureOutputBlockMatrixTest
//! @remark allocates random sparse block matrices
void StructureOutputBlockMatrixTestSymmetric(int rNumDAct, int rNumDDep, double rDensity)
{
    NuTo::Timer timer("StructureOutputBlockMatrixTestSymmetric::DefineRandomMatrices");
    NuTo::DofStatus s;
    s.SetDofTypes({NuTo::Node::eDof::DISPLACEMENTS});
    s.SetActiveDofTypes({NuTo::Node::eDof::DISPLACEMENTS});
    s.SetIsSymmetric(NuTo::Node::eDof::DISPLACEMENTS, true);

    NuTo::StructureOutputBlockMatrix BM4(s);

    BM4.JJ(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::DISPLACEMENTS) =
            NuTo::SparseMatrixCSRVector2Symmetric<double>::Random(rNumDAct, rDensity);
    BM4.JK(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::DISPLACEMENTS) =
            NuTo::SparseMatrixCSRVector2General<double>::Random(rNumDAct, rNumDDep, rDensity);
    BM4.KJ(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::DISPLACEMENTS) =
            BM4.JK(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::DISPLACEMENTS)
                    .AsSparseMatrixCSRVector2General()
                    .Transpose();
    BM4.KK(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::DISPLACEMENTS) =
            NuTo::SparseMatrixCSRVector2Symmetric<double>::Random(rNumDDep, rDensity);
    BM4.CheckDimensions();

    timer.Reset("StructureOutputBlockMatrixTestSymmetric::Export");

    NuTo::BlockSparseMatrix CMatrix(s, false);
    CMatrix(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::DISPLACEMENTS) =
            NuTo::SparseMatrixCSRVector2General<double>::Random(rNumDDep, rNumDAct, rDensity);

    NuTo::BlockSparseMatrix hessian(s, true);
    hessian = BM4.JJ; // just to get every sub matrix ...
    hessian.SetZero(); // ... to the right dimensions


    auto CMatrixExp = CMatrix.ExportToCSRVector2General();
    auto JJExp = BM4.JJ.ExportToCSRVector2General();
    auto JKExp = BM4.JK.ExportToCSRVector2General();
    auto KJExp = BM4.KJ.ExportToCSRVector2General();
    auto KKExp = BM4.KK.ExportToCSRVector2General();
    auto hessianExp = JJExp;
    hessianExp.SetZeroEntries();

    timer.Reset("StructureOutputBlockMatrixTestSymmetric::ApplyCMatrix_NEW");

    BM4.ApplyCMatrixScal(hessian, CMatrix, 1);

    timer.Reset("StructureOutputBlockMatrixTestSymmetric::ApplyCMatrix_OLD");

    hessianExp += JJExp;
    hessianExp += (CMatrixExp.Transpose() * KKExp * CMatrixExp - CMatrixExp.Transpose() * KJExp - JKExp * CMatrixExp);

    timer.Reset("StructureOutputBlockMatrixTestSymmetric::cleanup");


    auto diff = hessian.ExportToCSRVector2General() - hessianExp;
    double tolerance = 1.e-8;
    double diffMaxMin = diff.Max() - diff.Min();
    if (diffMaxMin > tolerance)
    {
        std::cout << diffMaxMin << std::endl << diff.ConvertToFullMatrix() << std::endl;
        throw NuTo::Exception("[StructureOutputBlockMatrixTestSymmetric] ApplyCMatrix incorrect.");
    }
}

BOOST_AUTO_TEST_CASE(SparseSymmetric)
{
    StructureOutputBlockMatrixTestSymmetric(10, 0, 1); // cmat == 0
    StructureOutputBlockMatrixTestSymmetric(10, 2, 1);
}

BOOST_AUTO_TEST_CASE(BlockScalarTest)
{
    NuTo::DofStatus s;
    s.SetDofTypes(
            {NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::TEMPERATURE, NuTo::Node::eDof::WATERVOLUMEFRACTION});
    s.SetActiveDofTypes({NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::TEMPERATURE});

    NuTo::BlockScalar A(s);
    A.DefineDefaultValueToIninitializedDofTypes(-3.);
    NuTo::BlockScalar B(s);
    B.DefineDefaultValueToIninitializedDofTypes(2.);
    NuTo::BlockScalar C(s);
    NuTo::BlockScalar D = B * 2;
    NuTo::BlockScalar E = D / 2;
    C = B;

    A.Info();
    B.Info();
    C.Info();
    D.Info();
    E.Info();

    BOOST_CHECK(A.CheckDofWiseLessActivDofs(B));
    BOOST_CHECK(A < B);
    BOOST_CHECK(A != B);
    BOOST_CHECK(B == C);
    BOOST_CHECK(B.CheckDofWiseLessActivDofs(D));
    BOOST_CHECK(B < D);
    BOOST_CHECK(C == E);

    B[NuTo::Node::eDof::DISPLACEMENTS] = 10;
    BOOST_CHECK(!(B.CheckDofWiseLessActivDofs(D)));
}


