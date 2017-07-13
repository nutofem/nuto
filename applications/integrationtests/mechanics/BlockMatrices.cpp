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


//! @brief converts a NuTo SparseMatrix in Vector2 format to a Eigen SparseMatrix
void SparseNuToToEigen(NuTo::SparseMatrixCSRVector2General<double>& rNuTo, Eigen::SparseMatrix<double>& rEigen)
{

    const std::vector<std::vector<int>>& columns = rNuTo.GetColumns();
    const std::vector<std::vector<double>>& values = rNuTo.GetValues();

    rEigen.resize(rNuTo.GetNumRows(), rNuTo.GetNumColumns());

    std::vector<Eigen::Triplet<double>> tripletList;


    // insert every nonzero element...
    for (unsigned int row = 0; row < columns.size(); row++)
        for (unsigned int col_count = 0; col_count < columns[row].size(); col_count++)
        {
            int col = columns[row][col_count];
            double val = values[row][col_count];
            tripletList.push_back(Eigen::Triplet<double>(row, col, val));
        }

    rEigen.setFromTriplets(tripletList.begin(), tripletList.end());
    rEigen.makeCompressed();
}

//! @brief BlockFullVectorTest [BVT]
//! @remark tests are done with int-vectors for easier comparison without epsilons...
//! @remark assertion: NuTo::FullVector calculations are correct.
void BlockFullVectorTest()
{
    NuTo::Timer timer("BVT - BlockFullVectorTest:Init");

    NuTo::DofStatus s;
    s.SetDofTypes({NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::TEMPERATURE});
    s.SetActiveDofTypes({NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::TEMPERATURE});


    NuTo::BlockFullVector<int> v1(s), v2(s);
    v1[NuTo::Node::eDof::DISPLACEMENTS] = Eigen::Vector3i({1, 1, 1});
    v2[NuTo::Node::eDof::DISPLACEMENTS] = Eigen::Vector3i({2, 2, 2});

    v1[NuTo::Node::eDof::TEMPERATURE] = Eigen::Vector2i({10, 10});
    v2[NuTo::Node::eDof::TEMPERATURE] = Eigen::Vector2i({20, 20});

    /*
     * Export
     */
    timer.Reset("BVT:Export");
    Eigen::VectorXi ev1 = v1.Export();
    Eigen::VectorXi ev2 = v2.Export();

    if (ev1.rows() != 5)
        throw NuTo::Exception("[BVT:Export] Exported v1 has wrong size.");
    if (ev1.sum() != 23)
        throw NuTo::Exception("[BVT:Export] Exported v1 has wrong sum.");
    if (ev2.sum() != 46)
        throw NuTo::Exception("[BVT:Export] Exported v2 has wrong sum.");

    /*
     * Addition
     */
    timer.Reset("BVT:Addition");
    if (v1 + v2 != v2 + v1)
        throw NuTo::Exception("[BVT:Addition] v1 + v2 != v2 + v1");
    if ((v1 + v2).Export() != ev1 + ev2)
        throw NuTo::Exception("[BVT:Addition] e(v1 + v2) != ev1 + ev2");

    /*
     * Subtraction
     */
    timer.Reset("BVT:Subtraction");
    if (v1 - v1 != v2 - v2)
        throw NuTo::Exception("[BVT:Subtraction] v1 - v1 != v2 - v2");
    if ((v1 - v2).Export() != ev1 - ev2)
        throw NuTo::Exception("[BVT:Subtraction] e(v1 - v2) != ev1 - ev2");

    /*
     * ScalarMultiplication
     */
    timer.Reset("BVT:ScalarMultiplication");
    if (v1 + v1 != v1 * 2)
        throw NuTo::Exception("[BVT:ScalarMultiplication] v1 + v1 != v1 * 2");
    if ((v1 + v1).Export() != ev1 * 2)
        throw NuTo::Exception("[BVT:ScalarMultiplication] e(v1 + v1) != ev1 * 2");

    /*
     * Chaining
     */
    timer.Reset("BVT:Chaining");
    NuTo::BlockFullVector<int> result1(s), result2(s);
    result1 = v1 * 3 + v1 - v2 - v2 + v2 * 42;

    std::cout << "Chained operations done." << std::endl;

    result2 = v1;
    result2 *= 3;
    result2 += v1;
    result2 -= v2;
    result2 -= v2;
    result2 += v2 * 42;

    Eigen::VectorXi eresult = ev1 * 3 + ev1 - ev2 - ev2 + ev2 * 42;

    if (result1 != result2)
        throw NuTo::Exception("[BVT:Chaining] went wrong ... ");
    if (result1.Export() != eresult)
        throw NuTo::Exception("[BVT:Chaining] went wrong ... ");


    /*
     * ActiveDofTypes
     */
    timer.Reset("BVT:ActiveDofTypes");
    s.SetActiveDofTypes({NuTo::Node::eDof::DISPLACEMENTS});

    if (v1.Export().rows() != 3)
        throw NuTo::Exception("[BVT:ActiveDofTypes] Exported v1 has wrong size.");


    if (v1[NuTo::Node::eDof::TEMPERATURE] != (v1 + v1)[NuTo::Node::eDof::TEMPERATURE])
        throw NuTo::Exception("[BVT:ActiveDofTypes] Addition changes inactive dofs.");
    if (v1[NuTo::Node::eDof::TEMPERATURE] != (v1 - v1)[NuTo::Node::eDof::TEMPERATURE])
        throw NuTo::Exception("[BVT:ActiveDofTypes] Subtraction changes inactive dofs.");
    if (v1[NuTo::Node::eDof::TEMPERATURE] != (v1 * 2)[NuTo::Node::eDof::TEMPERATURE])
        throw NuTo::Exception("[BVT:ActiveDofTypes] Multiplication changes inactive dofs.");

    v1.Info();


    s.SetActiveDofTypes({NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::TEMPERATURE});

    /*
     * import
     */

    // allocate BlockFullVector with the right dimensions.
    NuTo::BlockFullVector<int> imported(v1);

    // create new random vector with the global dimensions of v1
    Eigen::VectorXi toImport = Eigen::VectorXi::Random(v1.GetNumRows());

    // import and check if equal.
    imported.Import(toImport);

    if (toImport != imported.Export())
    {
        std::cout << "to Import: \n " << toImport << std::endl;
        std::cout << "imported: \n ";
        std::cout << imported.Export();
        throw NuTo::Exception("[BVT:Import] failed.");
    }
}


//! @brief BlockFullMatrixTest [BMT]
//! @remark just test the access operators
void BlockFullMatrixTest()
{
    NuTo::Timer timer("BlockFullMatrixTest:Init");

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
    m.CheckDimensions();
}

//! @brief BlockSparseMatrixTest
//! @remark tests access operations, vector*matrix operation,
void BlockSparseMatrixTest()
{
    NuTo::Timer timer("BMT - BlockSparseMatrixTest:Init");

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
    m.CheckDimensions();
    m.Info();

    timer.Reset("BMT - vector*matrix");

    NuTo::BlockFullVector<double> v(s);
    v[NuTo::Node::eDof::DISPLACEMENTS] = Eigen::VectorXd::Random(numD);
    v[NuTo::Node::eDof::TEMPERATURE] = Eigen::VectorXd::Random(numT);

    auto result = m * v * 4 + m * v;

    timer.Reset("BMT - Export");

    Eigen::MatrixXd exportReference = m.ExportToFullMatrix();
    Eigen::MatrixXd exportCSRVector2 = m.ExportToCSRVector2General().ConvertToFullMatrix();
    Eigen::MatrixXd exportCSR = m.ExportToCSRGeneral().ConvertToFullMatrix();

    if ((exportCSRVector2 - exportReference).norm() > 1.e-8)
    {
        std::cout << "Reference \n" << exportReference << std::endl;
        std::cout << "Export to CSRVector2 \n" << exportCSRVector2 << std::endl;
        throw NuTo::Exception("[BlockSparseMatrixTest] Export to CSRVector2 failed.");
    }

    if ((exportCSR - exportReference).norm() > 1.e-8)
    {
        std::cout << "Reference \n" << exportReference << std::endl;
        std::cout << "Export to CSR \n" << exportCSR << std::endl;
        throw NuTo::Exception("[BlockSparseMatrixTest] Export to CSR failed.");
    }

    s.SetActiveDofTypes({NuTo::Node::eDof::DISPLACEMENTS});

    auto CSR = m.ExportToCSR();
    if (not CSR->IsSymmetric())
        throw NuTo::Exception("[BlockSparseMatrixTest] Symmetric export to CSR failed.");

    Eigen::MatrixXd exportCSRSymm = CSR->ConvertToFullMatrix();
    if ((exportCSRSymm - m.ExportToFullMatrix()).norm() > 1.e-8)
    {
        std::cout << "Reference \n" << exportReference << std::endl;
        std::cout << "Export to CSRSymm \n" << exportCSRSymm << std::endl;
        throw NuTo::Exception("[BlockSparseMatrixTest] Export to CSRSymm failed.");
    }
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

    s.SetHasInteractingConstraints(true);
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

    s.SetHasInteractingConstraints(true);
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


//! @brief CSR vs CSRVector2 ... simulates the calculation of the hessian in the time integration
//! @remark conjecture confirmed: Vector2+conversion is much (>> factor 1000 for big matrices) faster.
void RandomAddition_CSR_vs_CSRVector2()
{
    NuTo::Timer timer("CSR_vs_CSRVector2");

    size_t dim = 1e4;
    double density = 5. / dim;

    auto CSRVector2a = NuTo::SparseMatrixCSRVector2General<double>::Random(dim, dim, density, 17);
    auto CSRVector2b = NuTo::SparseMatrixCSRVector2General<double>::Random(dim, dim, density, 1337);
    auto CSRVector2c = NuTo::SparseMatrixCSRVector2General<double>::Random(dim, dim, density, 42);


    NuTo::SparseMatrixCSRGeneral<double> CSRa(CSRVector2a);
    NuTo::SparseMatrixCSRGeneral<double> CSRb(CSRVector2b);
    NuTo::SparseMatrixCSRGeneral<double> CSRc(CSRVector2c);

    timer.Reset("CSR");

    NuTo::SparseMatrixCSRGeneral<double> resultCSR = CSRa;
    resultCSR += CSRb;
    resultCSR += CSRc;


    timer.Reset("CSRVector2");
    NuTo::SparseMatrixCSRVector2General<double> tmp = CSRVector2a;
    tmp += CSRVector2b;
    tmp += CSRVector2c;
    NuTo::SparseMatrixCSRGeneral<double> resultCSRVector2(tmp);

    timer.Reset("cleanup.");
}


void SimpleTestResult(bool rResult, std::string rTestName)
{
    if (rResult)
    {
        std::cout << rTestName << " --- passed" << std::endl;
    }
    else
    {
        throw NuTo::Exception(rTestName + " --- failed.");
    }
}

//! @brief Test the block scalar class
void BlockScalarTest()
{
    NuTo::Timer timer("BlockScalarTest");
    NuTo::DofStatus s;
    s.SetDofTypes(
            {NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::TEMPERATURE, NuTo::Node::eDof::WATERVOLUMEFRACTION});
    s.SetActiveDofTypes({NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::TEMPERATURE});

    std::cout << "Active Dofs: DISPLACEMENTS / TEMPERATURE" << std::endl << std::endl;
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

    SimpleTestResult(A.CheckDofWiseLessActivDofs(B), std::string("[BlockScalarTest] A.CheckDofWiseLessActivDofs(B)"));
    SimpleTestResult(A < B, std::string("[BlockScalarTest] A < B"));
    SimpleTestResult(A != B, std::string("[BlockScalarTest] A!=B"));
    SimpleTestResult(B == C, std::string("[BlockScalarTest] B==C"));
    SimpleTestResult(B.CheckDofWiseLessActivDofs(D), std::string("[BlockScalarTest] B.CheckDofWiseLessActivDofs(D)"));
    SimpleTestResult(B < D, std::string("[BlockScalarTest] B < D"));
    SimpleTestResult(C == E, std::string("[BlockScalarTest] C==E"));

    B[NuTo::Node::eDof::DISPLACEMENTS] = 10;
    SimpleTestResult(!(B.CheckDofWiseLessActivDofs(D)),
                     std::string("[BlockScalarTest] B.CheckDofWiseLessActivDofs(D) modified"));
}


int main()
{
    BlockFullVectorTest();
    BlockFullMatrixTest();
    BlockScalarTest();
    BlockSparseMatrixTest();
    StructureOutputBlockMatrixTestGeneral(10, 8, 0, 0, 1); // cmat == 0
    StructureOutputBlockMatrixTestGeneral(10, 8, 4, 2, 1);

    StructureOutputBlockMatrixTestSymmetric(10, 0, 1); // cmat == 0
    StructureOutputBlockMatrixTestSymmetric(10, 2, 1);

    return EXIT_SUCCESS;
}
