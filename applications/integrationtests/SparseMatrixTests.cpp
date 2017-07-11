#include <iostream>


#include "math/SparseMatrixCSRVector2General.h"
#include "math/SparseMatrixCSRVector2Symmetric.h"

#include "base/Timer.h"

void ThrowIfNotEqual(const NuTo::SparseMatrix<double>& rSparse, const Eigen::MatrixXd& rFull, double rTolerance = 1e-6)
{
    const Eigen::MatrixXd tmp = rSparse.ConvertToFullMatrix();
    if ((tmp - rFull).cwiseAbs().maxCoeff() > rTolerance)
    {
        std::cout << "############ SPARSE ############ \n" << tmp << std::endl;
        std::cout << "############# FULL ############# \n" << rFull << std::endl;
        throw NuTo::MathException("Matrices not equal");
    }

    if (rSparse.GetNumEntries() == 0)
    {
        throw NuTo::MathException("Sparse matrix is zero...");
    }
}

void SparseMatrixVector2GeneralSymmetricTests()
{
    NuTo::Timer timer("SparseMatrixVector2Tests::Init");
    const int dim = 10;
    const int seedG1 = 42;
    const int seedG2 = 15;
    const int seedS1 = 1337;
    const int seedS2 = 6174;


    const auto G1 = NuTo::SparseMatrixCSRVector2General<double>::Random(dim, dim, .5, seedG1);
    const auto G2 = NuTo::SparseMatrixCSRVector2General<double>::Random(dim, dim, .5, seedG2);
    const auto S1 = NuTo::SparseMatrixCSRVector2Symmetric<double>::Random(dim, .5, seedS1);
    const auto S2 = NuTo::SparseMatrixCSRVector2Symmetric<double>::Random(dim, .5, seedS2);

    const Eigen::MatrixXd fG1 = G1.ConvertToFullMatrix();
    const Eigen::MatrixXd fG2 = G2.ConvertToFullMatrix();
    const Eigen::MatrixXd fS1 = S1.ConvertToFullMatrix();
    const Eigen::MatrixXd fS2 = S2.ConvertToFullMatrix();

    /*
     * General
     */
    timer.Reset("SparseMatrixVector2Tests:: G1 + G2");
    ThrowIfNotEqual(G1 + G2, fG1 + fG2);

    timer.Reset("SparseMatrixVector2Tests:: G1 - G2");
    ThrowIfNotEqual(G1 - G2, fG1 - fG2);

    timer.Reset("SparseMatrixVector2Tests:: G1 * G2");
    ThrowIfNotEqual(G1 * G2, fG1 * fG2);

    timer.Reset("SparseMatrixVector2Tests:: G1 * 42");
    ThrowIfNotEqual(G1 * 42, fG1 * 42);

    timer.Reset("SparseMatrixVector2Tests:: G1.AddScal(G2,2)");
    auto tmpG = G1;
    tmpG.AddScal(G2, 2.);
    ThrowIfNotEqual(tmpG, fG1 + fG2 * 2);

    /*
     * Symmetric
     */

    timer.Reset("SparseMatrixVector2Tests:: S1 + S2");
    ThrowIfNotEqual(S1 + S2, fS1 + fS2);

    timer.Reset("SparseMatrixVector2Tests:: S1 - S2");
    ThrowIfNotEqual(S1 - S2, fS1 - fS2);

    //    timer.Reset("SparseMatrixVector2Tests:: S1 * S2");
    //    ThrowIfNotEqual(S1 * S2, fS1*fS2);

    timer.Reset("SparseMatrixVector2Tests:: S1 * 42");
    ThrowIfNotEqual(S1 * 42, fS1 * 42);

    timer.Reset("SparseMatrixVector2Tests:: S1.AddScal(S2,2)");
    auto tmpS = S1;
    tmpS.AddScal(S2, 2.);
    ThrowIfNotEqual(tmpS, fS1 + fS2 * 2);

    /*
     * mixed
     */

    timer.Reset("SparseMatrixVector2Tests:: G1 + S");
    ThrowIfNotEqual(G1 + S1, fG1 + fS1);

    timer.Reset("SparseMatrixVector2Tests:: G1 - S");
    ThrowIfNotEqual(G1 - S1, fG1 - fS1);

    timer.Reset("SparseMatrixVector2Tests:: S + G1");
    ThrowIfNotEqual(S1 + G1, fS1 + fG1);

    timer.Reset("SparseMatrixVector2Tests:: S - G1");
    ThrowIfNotEqual(S1 - G1, fS1 - fG1);

    timer.Reset("SparseMatrixVector2Tests:: G1.AddScal(S,2)");
    tmpG = G1;
    tmpG.AddScal(S1, 2.);
    ThrowIfNotEqual(tmpG, fG1 + fS1 * 2);
}

void SparseMatrixVector2Tests(int rNumActDofs, int rNumDepDofs, double rDensity, double rScalar)
{
    // general part
    {
        NuTo::Timer timer("SparseMatrixVector2Tests - General::Init");

        auto JJ = NuTo::SparseMatrixCSRVector2General<double>::Random(rNumActDofs, rNumActDofs, rDensity);
        auto JK = NuTo::SparseMatrixCSRVector2General<double>::Random(rNumActDofs, rNumDepDofs, rDensity);
        auto KJ = NuTo::SparseMatrixCSRVector2General<double>::Random(rNumDepDofs, rNumActDofs, rDensity);
        auto KK = NuTo::SparseMatrixCSRVector2General<double>::Random(rNumDepDofs, rNumDepDofs, rDensity);

        auto Cm = NuTo::SparseMatrixCSRVector2General<double>::Random(rNumDepDofs, rNumActDofs, rDensity);

        auto H1 = NuTo::SparseMatrixCSRVector2General<double>(rNumActDofs, rNumActDofs);
        auto H2 = NuTo::SparseMatrixCSRVector2General<double>(rNumActDofs, rNumActDofs);

        // get base references

        const NuTo::SparseMatrixCSRVector2<double>& bKK = KK;
        const NuTo::SparseMatrixCSRVector2<double>& bKJ = KJ;
        const NuTo::SparseMatrixCSRVector2<double>& bJK = JK;

        const NuTo::SparseMatrixCSRVector2<double>& bCm = Cm;

        NuTo::SparseMatrixCSRVector2<double>& bH1 = H1;

        timer.Reset("SparseMatrixVector2Tests - General::methods");

        bH1.AddScal(JJ, rScalar);
        bH1.Sub_TransA_B_Plus_C_D_Scal(bCm, bKJ, bJK, bCm, rScalar);
        bH1.Add_TransA_B_C_Scal(bCm, bKK, bCm, rScalar);

        timer.Reset("SparseMatrixVector2Tests - General::operator");

        H2 += JJ * rScalar;
        H2 -= (Cm.Transpose() * KJ + JK * Cm) * rScalar;
        H2 += (Cm.Transpose() * KK * Cm) * rScalar;


        timer.Reset("SparseMatrixVector2Tests - General::comparison");
        // compare
        auto diff = H2 - bH1.AsSparseMatrixCSRVector2General();
        double tolerance = 1.e-8;
        double diffMaxMin = diff.Max() - diff.Min();
        if (diffMaxMin > tolerance)
        {
            std::cout << diffMaxMin << std::endl;
            std::cout << diff.ConvertToFullMatrix() << std::endl;
            throw NuTo::MathException("[SparseMatrixVector2Tests - General] ApplyCMatrix incorrect.");
        }

        timer.Reset("SparseMatrixVector2Tests - General:: RJ + Cmat.Trans * RK");

        Eigen::VectorXd RJ = Eigen::VectorXd::Random(rNumActDofs);
        Eigen::VectorXd RK = Eigen::VectorXd::Random(rNumDepDofs);

        auto R1 = RJ;
        R1 += bCm.TransMult(RK);

        auto R2 = RJ;
        R2 += Cm.Transpose() * RK;

        Eigen::VectorXd diffVec = R1 - R2;
        double diffVecMaxMin = diffVec.maxCoeff() - diffVec.minCoeff();
        if (diffVecMaxMin > tolerance)
        {
            std::cout << diffVecMaxMin << std::endl;
            std::cout << diffVec << std::endl;
            throw NuTo::MathException("[SparseMatrixVector2Tests - General] RJ + Cmat.Trans * RK incorrect.");
        }


        timer.Reset("SparseMatrixVector2Tests - General::Residual");
    }


    // symmetric part
    {
        NuTo::Timer timer("SparseMatrixVector2Tests - Symmetric::Init");

        auto JJ = NuTo::SparseMatrixCSRVector2Symmetric<double>::Random(rNumActDofs, rDensity);
        auto JK = NuTo::SparseMatrixCSRVector2General<double>::Random(rNumActDofs, rNumDepDofs, rDensity);
        auto KJ = JK.Transpose();
        auto KK = NuTo::SparseMatrixCSRVector2Symmetric<double>::Random(rNumDepDofs, rDensity);
        NuTo::SparseMatrixCSRVector2General<double> KKg(KK);


        auto Cm = NuTo::SparseMatrixCSRVector2General<double>::Random(rNumDepDofs, rNumActDofs, rDensity);

        auto H1 = NuTo::SparseMatrixCSRVector2Symmetric<double>(rNumActDofs, rNumActDofs);
        auto H2 = NuTo::SparseMatrixCSRVector2General<double>(rNumActDofs, rNumActDofs);

        // get base references

        const NuTo::SparseMatrixCSRVector2<double>& bKK = KK;
        const NuTo::SparseMatrixCSRVector2<double>& bKJ = KJ;
        const NuTo::SparseMatrixCSRVector2<double>& bJK = JK;

        const NuTo::SparseMatrixCSRVector2<double>& bCm = Cm;

        NuTo::SparseMatrixCSRVector2<double>& bH1 = H1;

        timer.Reset("SparseMatrixVector2Tests - Symmetric::methods");

        bH1.AddScal(JJ, rScalar);
        bH1.Sub_TransA_B_Plus_C_D_Scal(bCm, bKJ, bJK, bCm, rScalar);
        bH1.Add_TransA_B_C_Scal(bCm, bKK, bCm, rScalar);

        timer.Reset("SparseMatrixVector2Tests - Symmetric::operator");

        H2 += JJ * rScalar;
        H2 -= (Cm.Transpose() * KJ + JK * Cm) * rScalar;
        H2 += (Cm.Transpose() * KKg * Cm) * rScalar;


        timer.Reset("SparseMatrixVector2Tests - Symmetric::cleanup");
        // compare
        auto diff = H2 - bH1.AsSparseMatrixCSRVector2Symmetric();
        double tolerance = 1.e-8;
        double diffMaxMin = diff.Max() - diff.Min();
        if (diffMaxMin > tolerance)
        {
            std::cout << diffMaxMin << std::endl;
            std::cout << diff.ConvertToFullMatrix() << std::endl;
            throw NuTo::MathException("[SparseMatrixVector2Tests - Symmetric] ApplyCMatrix incorrect.");
        }
    }
}


void GaussEliminationTests(int rNumActDofs, int rNumDepDofs, double rDensity)
{
    NuTo::Timer timer("GaussEliminationTests::Init");
    NuTo::SparseMatrixCSRVector2General<double> cMatVector2 =
            NuTo::SparseMatrixCSRVector2General<double>::Random(rNumDepDofs, rNumActDofs + rNumDepDofs, rDensity);
    timer.Reset("GaussEliminationTests::Init");
    NuTo::SparseMatrixCSRVector2General<double> RHSVector2(rNumDepDofs, rNumDepDofs);
    for (int i = 0; i < rNumDepDofs; ++i)
        RHSVector2.AddValue(i, i, 1.);

    NuTo::SparseMatrixCSRVector2<double>& cMatVector2ref = cMatVector2;
    NuTo::SparseMatrixCSRVector2<double>& RHSVector2ref = RHSVector2;

    NuTo::SparseMatrixCSRGeneral<double> cMat(cMatVector2);
    NuTo::SparseMatrixCSRGeneral<double> RHS(RHSVector2);


    std::vector<int> mappingInitialToNewOrderingVector2;
    std::vector<int> mappingNewToInitialOrderingVector2;
    std::vector<int> mappingInitialToNewOrdering;
    std::vector<int> mappingNewToInitialOrdering;

    timer.Reset("GaussEliminationTests::Gauss-CSRGeneral");
    cMat.Gauss(RHS, mappingNewToInitialOrdering, mappingInitialToNewOrdering);

    timer.Reset("GaussEliminationTests::Gauss-CSRVector2");
    cMatVector2ref.Gauss(RHSVector2ref, mappingNewToInitialOrderingVector2, mappingInitialToNewOrderingVector2);

    timer.Reset("GaussEliminationTests::checks");


    // check results
    if (mappingInitialToNewOrdering != mappingInitialToNewOrderingVector2)
        throw NuTo::MathException("[GaussEliminationTests] wrong mappingInitialToNewOrdering");

    if (mappingNewToInitialOrdering != mappingNewToInitialOrderingVector2)
        throw NuTo::MathException("[GaussEliminationTests] wrong mappingInitialToNewOrdering");

    auto diffCMat = cMatVector2 - cMat;
    if (diffCMat.Max() - diffCMat.Min() > 1e-8)
    {
        std::cout << diffCMat.ConvertToFullMatrix() << std::endl;
        throw NuTo::MathException("[GaussEliminationTests] wrong cMat");
    }

    auto diffRHS = RHSVector2 - RHS;
    if (diffRHS.Max() - diffRHS.Min() > 1e-8)
    {
        std::cout << diffRHS.ConvertToFullMatrix() << std::endl;
        throw NuTo::MathException("[GaussEliminationTests] wrong diffRHS");
    }

    std::vector<int> tmpMapping(rNumDepDofs + rNumActDofs);
    std::iota(tmpMapping.begin(), tmpMapping.begin() + rNumDepDofs, rNumActDofs);
    std::iota(tmpMapping.begin() + rNumDepDofs, tmpMapping.end(), 0);

    timer.Reset("GaussEliminationTests::ReorderColumns-CSRGeneral");

    cMat.ReorderColumns(tmpMapping);

    timer.Reset("GaussEliminationTests::ReorderColumns-CSRVector2");

    cMatVector2ref.ReorderColumns(tmpMapping);

    timer.Reset("GaussEliminationTests::cleanup");


    diffCMat = cMatVector2 - cMat;
    if (diffCMat.Max() - diffCMat.Min() > 1e-8)
    {
        std::cout << diffCMat.ConvertToFullMatrix() << std::endl;
        throw NuTo::MathException("[GaussEliminationTests] wrong cMat after renumbering.");
    }


    timer.Reset("GaussEliminationTests::RemoveLastColumns-CSRGeneral");

    cMat.RemoveLastColumns(rNumDepDofs);

    timer.Reset("GaussEliminationTests::RemoveLastColumns-CSRVector2");

    cMatVector2ref.RemoveLastColumns(rNumDepDofs);


    diffCMat = cMatVector2 - cMat;
    if (diffCMat.Max() - diffCMat.Min() > 1e-8)
    {
        std::cout << diffCMat.ConvertToFullMatrix() << std::endl;
        throw NuTo::MathException("[GaussEliminationTests] wrong cMat after removal of the last columns.");
    }
}

void SerializationTest()
{
}

int main()
{
    try
    {
        SparseMatrixVector2GeneralSymmetricTests();
        SparseMatrixVector2Tests(10, 5, 1, 1);

        //        SparseMatrixVector2Tests(1e6, 1e3, 0.0001, 2);

        GaussEliminationTests(12, 3, 1);
        SerializationTest();
    }
    catch (NuTo::MathException& e)
    {
        std::cout << "\n\n\n errors occurred \n\n\n";
        std::cout << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
