#include "nuto/base/Timer.h"
#include <iostream>

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>


#include "nuto/math/FullMatrix.h"
#include "nuto/math/SparseMatrixCSRVector2General.h"


#include "nuto/mechanics/dofSubMatrixStorage/BlockScalar.h"
#include "nuto/mechanics/dofSubMatrixStorage/BlockFullVector.h"
#include "nuto/mechanics/dofSubMatrixStorage/BlockFullMatrix.h"
#include "nuto/mechanics/dofSubMatrixStorage/BlockSparseMatrix.h"
#include "nuto/mechanics/structures/StructureOutputBlockMatrix.h"
#include "nuto/mechanics/elements/ElementOutputBlockMatrixDouble.h"
#include "nuto/mechanics/dofSubMatrixStorage/DofStatus.h"
#include "nuto/mechanics/MechanicsException.h"


//! @brief converts a NuTo SparseMatrix in Vector2 format to a Eigen SparseMatrix
void SparseNuToToEigen(NuTo::SparseMatrixCSRVector2General<double>& rNuTo, Eigen::SparseMatrix<double>& rEigen)
{

    const std::vector<std::vector<int> >& columns = rNuTo.GetColumns();
    const std::vector<std::vector<double>>& values = rNuTo.GetValues();

    rEigen.resize(rNuTo.GetNumRows(), rNuTo.GetNumColumns());

    std::vector<Eigen::Triplet<double> > tripletList;


    // insert every nonzero element...
    for (unsigned int row=0; row < columns.size(); row++)
        for (unsigned int col_count=0; col_count < columns[row].size(); col_count++)
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
    s.SetDofTypes       ({NuTo::Node::DISPLACEMENTS, NuTo::Node::TEMPERATURES});
    s.SetActiveDofTypes ({NuTo::Node::DISPLACEMENTS, NuTo::Node::TEMPERATURES});


    NuTo::BlockFullVector<int> v1(s), v2(s);
    v1[NuTo::Node::DISPLACEMENTS] = NuTo::FullVector<int, Eigen::Dynamic>({1,1,1});
    v2[NuTo::Node::DISPLACEMENTS] = NuTo::FullVector<int, Eigen::Dynamic>({2,2,2});

    v1[NuTo::Node::TEMPERATURES] = NuTo::FullVector<int, Eigen::Dynamic>({10,10});
    v2[NuTo::Node::TEMPERATURES] = NuTo::FullVector<int, Eigen::Dynamic>({20,20});

    /*
     * Export
     */
    timer.Reset("BVT:Export");
    NuTo::FullVector<int, Eigen::Dynamic> ev1 = v1.Export();
    NuTo::FullVector<int, Eigen::Dynamic> ev2 = v2.Export();

    if (ev1.GetNumRows() != 5)              throw NuTo::MechanicsException("[BVT:Export] Exported v1 has wrong size.");
    if (ev1.Sum() != 23)                    throw NuTo::MechanicsException("[BVT:Export] Exported v1 has wrong sum.");
    if (ev2.Sum() != 46)                    throw NuTo::MechanicsException("[BVT:Export] Exported v2 has wrong sum.");

    /*
     * Addition
     */
    timer.Reset("BVT:Addition");
    if (v1 + v2 != v2 + v1)                 throw NuTo::MechanicsException("[BVT:Addition] v1 + v2 != v2 + v2");
    if ((v1 + v2).Export() != ev1 + ev2)    throw NuTo::MechanicsException("[BVT:Addition] e(v1 + v2) != ev1 + ev2");

    /*
     * Subtraction
     */
    timer.Reset("BVT:Subtraction");
    if (v1 - v1 != v2 - v2)                 throw NuTo::MechanicsException("[BVT:Subtraction] v1 - v1 != v2 - v2");
    if ((v1 - v2).Export() != ev1 - ev2)    throw NuTo::MechanicsException("[BVT:Subtraction] e(v1 - v2) != ev1 - ev2");

    /*
     * ScalarMultiplication
     */
    timer.Reset("BVT:ScalarMultiplication");
    if (v1 + v1 != v1 * 2)                  throw NuTo::MechanicsException("[BVT:ScalarMultiplication] v1 + v1 != v1 * 2");
    if ((v1 + v1).Export() != ev1 * 2)      throw NuTo::MechanicsException("[BVT:ScalarMultiplication] e(v1 + v1) != ev1 * 2");

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
    result2 += v2*42;

    NuTo::FullVector<int, Eigen::Dynamic> eresult = ev1 * 3 + ev1 - ev2 - ev2 + ev2 * 42;

    if (result1 != result2)                 throw NuTo::MechanicsException("[BVT:Chaining] went wrong ... ");
    if (result1.Export() != eresult)        throw NuTo::MechanicsException("[BVT:Chaining] went wrong ... ");


    /*
     * ActiveDofTypes
     */
    timer.Reset("BVT:ActiveDofTypes");
    s.SetActiveDofTypes({NuTo::Node::DISPLACEMENTS});

    if (v1.Export().GetNumRows() != 3)      throw NuTo::MechanicsException("[BVT:ActiveDofTypes] Exported v1 has wrong size.");


    if (v1[NuTo::Node::TEMPERATURES] != (v1+v1)[NuTo::Node::TEMPERATURES])
                                            throw NuTo::MechanicsException("[BVT:ActiveDofTypes] Addition changes inactive dofs.");
    if (v1[NuTo::Node::TEMPERATURES] != (v1-v1)[NuTo::Node::TEMPERATURES])
                                            throw NuTo::MechanicsException("[BVT:ActiveDofTypes] Subtraction changes inactive dofs.");
    if (v1[NuTo::Node::TEMPERATURES] != (v1*2)[NuTo::Node::TEMPERATURES])
                                            throw NuTo::MechanicsException("[BVT:ActiveDofTypes] Multiplication changes inactive dofs.");

    v1.Info();



    s.SetActiveDofTypes({NuTo::Node::DISPLACEMENTS, NuTo::Node::TEMPERATURES});

    /*
     * import
     */

    // allocate BlockFullVector with the right dimensions.
    NuTo::BlockFullVector<int> imported(v1);

    // create new random vector with the global dimensions of v1
    NuTo::FullVector<int, Eigen::Dynamic> toImport = NuTo::FullVector<int, Eigen::Dynamic>::Random(v1.GetNumRows());

    // import and check if equal.
    imported.Import(toImport);

    if (toImport != imported.Export())
    {
        std::cout << "to Import: \n ";
        toImport.Info();
        std::cout << "imported: \n ";
        imported.Export().Info();
        throw NuTo::MechanicsException("[BVT:Import] failed.");
    }



}


//! @brief BlockFullMatrixTest [BMT]
//! @remark just test the access operators
void BlockFullMatrixTest()
{
    NuTo::Timer timer("BlockFullMatrixTest:Init");

    NuTo::DofStatus s;
    s.SetDofTypes       ({NuTo::Node::DISPLACEMENTS, NuTo::Node::TEMPERATURES});
    s.SetActiveDofTypes ({NuTo::Node::DISPLACEMENTS, NuTo::Node::TEMPERATURES});

    NuTo::BlockFullMatrix<double> m(s);


    int numD = 3;
    int numT = 2;

    m(NuTo::Node::DISPLACEMENTS, NuTo::Node::DISPLACEMENTS) = NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>::Random(numD, numD);

    m(NuTo::Node::DISPLACEMENTS, NuTo::Node::TEMPERATURES ).Resize(numD,numT);
    m(NuTo::Node::DISPLACEMENTS, NuTo::Node::TEMPERATURES ) << 1, 2, 3, 1, 2, 3;

    m(NuTo::Node::TEMPERATURES,  NuTo::Node::DISPLACEMENTS) = m(NuTo::Node::DISPLACEMENTS, NuTo::Node::TEMPERATURES ).Trans();

    auto& matrixTT = m(NuTo::Node::TEMPERATURES,  NuTo::Node::TEMPERATURES );
    matrixTT.Resize(numT,numT);

    m.Info();
    m.CheckDimensions();
}

//! @brief BlockSparseMatrixTest
//! @remark tests access operations, vector*matrix operation,
void BlockSparseMatrixTest()
{
    NuTo::Timer timer("BMT - BlockSparseMatrixTest:Init");

    NuTo::DofStatus s;
    s.SetDofTypes       ({NuTo::Node::DISPLACEMENTS, NuTo::Node::TEMPERATURES});
    s.SetActiveDofTypes ({NuTo::Node::DISPLACEMENTS, NuTo::Node::TEMPERATURES});

    NuTo::BlockSparseMatrix m(s);

    size_t numD = 5;
    size_t numT = 2;
    double density = 1;

    m(NuTo::Node::DISPLACEMENTS, NuTo::Node::DISPLACEMENTS) = NuTo::SparseMatrixCSRVector2General<double>::Random(numD,numD,density);
    m(NuTo::Node::DISPLACEMENTS, NuTo::Node::TEMPERATURES ) = NuTo::SparseMatrixCSRVector2General<double>::Random(numD,numT,density);
    m(NuTo::Node::TEMPERATURES,  NuTo::Node::DISPLACEMENTS) = NuTo::SparseMatrixCSRVector2General<double>::Random(numT,numD,density);
    m(NuTo::Node::TEMPERATURES,  NuTo::Node::TEMPERATURES ) = NuTo::SparseMatrixCSRVector2General<double>::Random(numT,numT,density);
    m.CheckDimensions();
    m.Info();

    timer.Reset("BMT - vector*matrix");

    NuTo::BlockFullVector<double> v(s);
    v[NuTo::Node::DISPLACEMENTS] = NuTo::FullVector<double, Eigen::Dynamic>::Random(numD);
    v[NuTo::Node::TEMPERATURES ] = NuTo::FullVector<double, Eigen::Dynamic>::Random(numT);

    auto result = m*v*4 + m*v;

    timer.Reset("BMT - Export");

    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> exportReference        (m.ExportToFullMatrix());
    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> exportCSRVector2       (m.ExportToCSRVector2General());
    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> exportCSR              (m.ExportToCSRGeneral());

    if ((exportCSRVector2 - exportReference).norm() > 1.e-8)
    {
        std::cout << "Reference \n";
        exportReference.Info();
        std::cout << "Export to CSRVector2 \n";
        exportCSRVector2.Info();
        throw NuTo::MechanicsException("[BlockSparseMatrixTest] Export to CSRVector2 failed.");
    }

    if ((exportCSR - exportReference).norm() > 1.e-8)
    {
        std::cout << "Reference \n";
        exportReference.Info();
        std::cout << "Export to CSR \n";
        exportCSR.Info();
        throw NuTo::MechanicsException("[BlockSparseMatrixTest] Export to CSR failed.");
    }
}

//! @brief StructureOutputBlockMatrixTest
//! @remark allocates random sparse block matrices
void StructureOutputBlockMatrixTestGeneral(int rNumDAct, int rNumTAct, int rNumDDep, int rNumTDep, double rDensity)
{
    NuTo::Timer timer("StructureOutputBlockMatrixTestGeneral::DefineRandomMatrices");
    NuTo::DofStatus s;
    s.SetDofTypes       ({NuTo::Node::DISPLACEMENTS, NuTo::Node::TEMPERATURES});
    s.SetActiveDofTypes ({NuTo::Node::DISPLACEMENTS, NuTo::Node::TEMPERATURES});


    NuTo::StructureOutputBlockMatrix BM4(s);

    BM4.JJ(NuTo::Node::DISPLACEMENTS, NuTo::Node::DISPLACEMENTS) = NuTo::SparseMatrixCSRVector2General<double>::Random(rNumDAct, rNumDAct, rDensity);
    BM4.JJ(NuTo::Node::DISPLACEMENTS, NuTo::Node::TEMPERATURES ) = NuTo::SparseMatrixCSRVector2General<double>::Random(rNumDAct, rNumTAct, rDensity);
    BM4.JJ(NuTo::Node::TEMPERATURES,  NuTo::Node::DISPLACEMENTS) = NuTo::SparseMatrixCSRVector2General<double>::Random(rNumTAct, rNumDAct, rDensity);
    BM4.JJ(NuTo::Node::TEMPERATURES,  NuTo::Node::TEMPERATURES ) = NuTo::SparseMatrixCSRVector2General<double>::Random(rNumTAct, rNumTAct, rDensity);

    BM4.JK(NuTo::Node::DISPLACEMENTS, NuTo::Node::DISPLACEMENTS) = NuTo::SparseMatrixCSRVector2General<double>::Random(rNumDAct, rNumDDep, rDensity);
    BM4.JK(NuTo::Node::DISPLACEMENTS, NuTo::Node::TEMPERATURES ) = NuTo::SparseMatrixCSRVector2General<double>::Random(rNumDAct, rNumTDep, rDensity);
    BM4.JK(NuTo::Node::TEMPERATURES,  NuTo::Node::DISPLACEMENTS) = NuTo::SparseMatrixCSRVector2General<double>::Random(rNumTAct, rNumDDep, rDensity);
    BM4.JK(NuTo::Node::TEMPERATURES,  NuTo::Node::TEMPERATURES ) = NuTo::SparseMatrixCSRVector2General<double>::Random(rNumTAct, rNumTDep, rDensity);

    BM4.KJ(NuTo::Node::DISPLACEMENTS, NuTo::Node::DISPLACEMENTS) = NuTo::SparseMatrixCSRVector2General<double>::Random(rNumDDep, rNumDAct, rDensity);
    BM4.KJ(NuTo::Node::DISPLACEMENTS, NuTo::Node::TEMPERATURES ) = NuTo::SparseMatrixCSRVector2General<double>::Random(rNumDDep, rNumTAct, rDensity);
    BM4.KJ(NuTo::Node::TEMPERATURES,  NuTo::Node::DISPLACEMENTS) = NuTo::SparseMatrixCSRVector2General<double>::Random(rNumTDep, rNumDAct, rDensity);
    BM4.KJ(NuTo::Node::TEMPERATURES,  NuTo::Node::TEMPERATURES ) = NuTo::SparseMatrixCSRVector2General<double>::Random(rNumTDep, rNumTAct, rDensity);

    BM4.KK(NuTo::Node::DISPLACEMENTS, NuTo::Node::DISPLACEMENTS) = NuTo::SparseMatrixCSRVector2General<double>::Random(rNumDDep, rNumDDep, rDensity);
    BM4.KK(NuTo::Node::DISPLACEMENTS, NuTo::Node::TEMPERATURES ) = NuTo::SparseMatrixCSRVector2General<double>::Random(rNumDDep, rNumTDep, rDensity);
    BM4.KK(NuTo::Node::TEMPERATURES,  NuTo::Node::DISPLACEMENTS) = NuTo::SparseMatrixCSRVector2General<double>::Random(rNumTDep, rNumDDep, rDensity);
    BM4.KK(NuTo::Node::TEMPERATURES,  NuTo::Node::TEMPERATURES ) = NuTo::SparseMatrixCSRVector2General<double>::Random(rNumTDep, rNumTDep, rDensity);

    BM4.CheckDimensions();


    NuTo::BlockSparseMatrix CMatrix(s, false);
    CMatrix(NuTo::Node::DISPLACEMENTS, NuTo::Node::DISPLACEMENTS) = NuTo::SparseMatrixCSRVector2General<double>::Random(rNumDDep, rNumDAct, rDensity);
    CMatrix(NuTo::Node::DISPLACEMENTS, NuTo::Node::TEMPERATURES ).Resize(                                               rNumDDep, rNumTAct);
    CMatrix(NuTo::Node::TEMPERATURES,  NuTo::Node::DISPLACEMENTS).Resize(                                               rNumTDep, rNumDAct);
    CMatrix(NuTo::Node::TEMPERATURES,  NuTo::Node::TEMPERATURES ) = NuTo::SparseMatrixCSRVector2General<double>::Random(rNumTDep, rNumTAct, rDensity);

    NuTo::BlockSparseMatrix hessian(s, true);
    hessian = BM4.JJ;  // just to get every sub matrix ...
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
    double diffMaxMin = diff.Max()-diff.Min();
    if (diffMaxMin > tolerance)
    {
        std::cout << diffMaxMin << std::endl;
        throw NuTo::MechanicsException("[StructureOutputBlockMatrixTestGeneral] ApplyCMatrix incorrect.");
    }


    BM4.AddScal(BM4,2.7);

    diff = (JJExp*3.7 - BM4.JJ.ExportToCSRVector2General());
    if (diff.Max()-diff.Min() > tolerance)   throw NuTo::MechanicsException("[StructureOutputBlockMatrixTestGeneral] AddScal for JJ incorrect.");

    if (rNumDDep > 0)
    {
        diff = (JKExp*3.7 - BM4.JK.ExportToCSRVector2General());
        if (diff.Max()-diff.Min() > tolerance)   throw NuTo::MechanicsException("[StructureOutputBlockMatrixTestGeneral] AddScal for JK incorrect.");

        diff = (KJExp*3.7 - BM4.KJ.ExportToCSRVector2General());
        if (diff.Max()-diff.Min() > tolerance)   throw NuTo::MechanicsException("[StructureOutputBlockMatrixTestGeneral] AddScal for KJ incorrect.");

        diff = (KKExp*3.7 - BM4.KK.ExportToCSRVector2General());
        if (diff.Max()-diff.Min() > tolerance)   throw NuTo::MechanicsException("[StructureOutputBlockMatrixTestGeneral] AddScal for KK incorrect.");
    }

}

//! @brief StructureOutputBlockMatrixTest
//! @remark allocates random sparse block matrices
void StructureOutputBlockMatrixTestSymmetric(int rNumDAct, int rNumDDep, double rDensity)
{
    NuTo::Timer timer("StructureOutputBlockMatrixTestSymmetric::DefineRandomMatrices");
    NuTo::DofStatus s;
    s.SetDofTypes       ({NuTo::Node::DISPLACEMENTS});
    s.SetActiveDofTypes ({NuTo::Node::DISPLACEMENTS});
    s.SetIsSymmetric    (NuTo::Node::DISPLACEMENTS, true);

    NuTo::StructureOutputBlockMatrix BM4(s);

    BM4.JJ(NuTo::Node::DISPLACEMENTS, NuTo::Node::DISPLACEMENTS) = NuTo::SparseMatrixCSRVector2Symmetric<double>::Random(rNumDAct, rDensity);
    BM4.JK(NuTo::Node::DISPLACEMENTS, NuTo::Node::DISPLACEMENTS) = NuTo::SparseMatrixCSRVector2General<double>::Random(rNumDAct, rNumDDep, rDensity);
    BM4.KJ(NuTo::Node::DISPLACEMENTS, NuTo::Node::DISPLACEMENTS) = BM4.JK(NuTo::Node::DISPLACEMENTS, NuTo::Node::DISPLACEMENTS).AsSparseMatrixCSRVector2General().Transpose();
    BM4.KK(NuTo::Node::DISPLACEMENTS, NuTo::Node::DISPLACEMENTS) = NuTo::SparseMatrixCSRVector2Symmetric<double>::Random(rNumDDep, rDensity);
    BM4.CheckDimensions();

    timer.Reset("StructureOutputBlockMatrixTestSymmetric::Export");

    NuTo::BlockSparseMatrix CMatrix(s, false);
    CMatrix(NuTo::Node::DISPLACEMENTS, NuTo::Node::DISPLACEMENTS) = NuTo::SparseMatrixCSRVector2General<double>::Random(rNumDDep, rNumDAct, rDensity);

    NuTo::BlockSparseMatrix hessian(s, true);
    hessian = BM4.JJ;  // just to get every sub matrix ...
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
    double diffMaxMin = diff.Max()-diff.Min();
    if (diffMaxMin > tolerance)
    {
        std::cout << diffMaxMin << std::endl;
        NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> (diff).Info();
        throw NuTo::MechanicsException("[StructureOutputBlockMatrixTestSymmetric] ApplyCMatrix incorrect.");
    }



}



//! @brief CSR vs CSRVector2 ... simulates the calculation of the hessian in the time integration
//! @remark conjecture confirmed: Vector2+conversion is much (>> factor 1000 for big matrices) faster.
void RandomAddition_CSR_vs_CSRVector2()
{
    NuTo::Timer timer("CSR_vs_CSRVector2");

    size_t dim = 1e4;
    double density = 5./dim;

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
    if(rResult)
    {
        std::cout << rTestName << " --- passed" << std::endl;
    }
    else
    {
        throw NuTo::MechanicsException(rTestName +" --- failed.");
    }
}

//! @brief Test the block scalar class
void BlockScalarTest()
{
    NuTo::Timer timer("BlockScalarTest::Init");
    NuTo::DofStatus s;
    s.SetDofTypes       ({NuTo::Node::DISPLACEMENTS, NuTo::Node::TEMPERATURES, NuTo::Node::WATERVOLUMEFRACTION});
    s.SetActiveDofTypes ({NuTo::Node::DISPLACEMENTS, NuTo::Node::TEMPERATURES});

    std::cout << "Active Dofs: DISPLACEMENTS / TEMPERATURES" << std::endl << std::endl;
    NuTo::BlockScalar A(s); A.DefineDefaultValueToIninitializedDofTypes(-3.);
    NuTo::BlockScalar B(s); B.DefineDefaultValueToIninitializedDofTypes(2.);
    NuTo::BlockScalar C(s);
    NuTo::BlockScalar D = B*2;
    NuTo::BlockScalar E = D/2;
    C = B;

    A.Info();
    B.Info();
    C.Info();
    D.Info();
    E.Info();

    SimpleTestResult(A.CheckDofWiseLessActivDofs(B) ,std::string("[BlockScalarTest] A.CheckDofWiseLessActivDofs(B)"));
    SimpleTestResult(A < B                          ,std::string("[BlockScalarTest] A < B"));
    SimpleTestResult(A!=B                           ,std::string("[BlockScalarTest] A!=B"));
    SimpleTestResult(B==C                           ,std::string("[BlockScalarTest] B==C"));
    SimpleTestResult(B.CheckDofWiseLessActivDofs(D) ,std::string("[BlockScalarTest] B.CheckDofWiseLessActivDofs(D)"));
    SimpleTestResult(B < D                          ,std::string("[BlockScalarTest] B < D"));
    SimpleTestResult(C==E                           ,std::string("[BlockScalarTest] C==E"));

    B[NuTo::Node::eDof::DISPLACEMENTS] = 10;
    SimpleTestResult(!(B.CheckDofWiseLessActivDofs(D)) ,std::string("[BlockScalarTest] B.CheckDofWiseLessActivDofs(D) modified"));
}

int main()
{


#define TRY_CATCH

    try
    {
        BlockFullVectorTest();
        BlockFullMatrixTest();
        BlockScalarTest();
        BlockSparseMatrixTest();
        StructureOutputBlockMatrixTestGeneral(10, 8, 0, 0, 1); // cmat == 0
        StructureOutputBlockMatrixTestGeneral(10, 8, 4, 2, 1);

        StructureOutputBlockMatrixTestSymmetric(10, 0, 1);       // cmat == 0
        StructureOutputBlockMatrixTestSymmetric(10, 2, 1);

//        int dim = 1e5;
//        StructureOutputBlockMatrixTestGeneral(2*dim, dim, 1000, 1000, 0.0001);
//        StructureOutputBlockMatrixTestSymmetric(2*dim, 1000, 0.0001);

    }
    catch (NuTo::MechanicsException& e)
    {
        std::cout << "\n\n\n NuTo Mechanics errors occurred: \n\n\n";
        std::cout << e.what();
        return EXIT_FAILURE;
    }
    catch (NuTo::MathException& e)
    {
        std::cout << "\n\n\n NuTo Math errors occurred: \n\n\n";
        std::cout << e.what();
        return EXIT_FAILURE;
    }
    catch (...)
    {
        std::cout << "\n\n\n non-NuTo errors occurred \n\n\n";
        return EXIT_FAILURE;
    }


    return EXIT_SUCCESS;
}
