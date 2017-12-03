#pragma once

#include <mpi.h>
#include <boost/mpi.hpp>
#include "mechanics/timeIntegration/TimeIntegrationBase.h"
#include "mechanics/timeIntegration/postProcessing/PostProcessor.h"
#include "mechanics/timeIntegration/NewmarkDirect.h"
#include "mechanics/structures/StructureBase.h"
#include "base/Exception.h"
#include "base/Timer.h"

#include "mechanics/structures/StructureBase.h"
#include "mechanics/feti/StructureFeti.h"
#include "mechanics/structures/StructureBaseEnum.h"
#include "mechanics/structures/StructureOutputBlockMatrix.h"
#include "mechanics/structures/StructureOutputDummy.h"

#include "mechanics/feti/FetiSolver.h"

#include "base/CallbackInterface.h"
#include "math/SparseMatrixCSRGeneral.h"
#include <Eigen/Dense>
#include "mechanics/constitutive/ConstitutiveEnum.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveCalculateStaticData.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveIOMap.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveTimeStep.h"
#include "FetiSolver.h"
#include <Eigen/Dense>

#include <Eigen/Sparse>

#include "mechanics/feti/FetiLumpedPreconditioner.h"
#include <cmath>

namespace NuTo
{
template <class EigenSolver = Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>>
class NewmarkFeti : public NewmarkDirect
{
public:
    enum class eIterativeSolver
    {
        ConjugateGradient,
        BiconjugateGradientStabilized,
        ProjectedGmres,
        Direct
    };

    enum class eFetiScaling
    {
        None,
        Multiplicity,
        Superlumped
    };


    using VectorXd = Eigen::VectorXd;
    using MatrixXd = Eigen::MatrixXd;
    using SparseMatrix = Eigen::SparseMatrix<double>;
    using PreconditionerType = std::unique_ptr<FetiPreconditioner>;


    NewmarkFeti(StructureFeti* rStructure)
        : NewmarkDirect(rStructure)
        , mStructureFeti(rStructure)
        , mNumTotalActiveDofs(rStructure->GetNumTotalActiveDofs())
        , mNumTotalDofs(rStructure->GetNumTotalDofs())
        , mPreconditioner(std::make_unique<FetiLumpedPreconditioner>())
    {
    }


    ///
    /// \brief MpiGatherRecvCountAndDispls
    void MpiGatherRecvCountAndDispls(std::vector<int>& recvCount, std::vector<int>& displs, const int numValues)
    {
        boost::mpi::communicator world;
        const int numProcesses = world.size();
        // recvCount:
        // Contais the number of elements that are received from each process.
        recvCount.clear();
        recvCount.resize(numProcesses, 0);

        boost::mpi::all_gather<int>(world, numValues, recvCount);

        // displs:
        // Entry i specifies the displacement (relative to recvbuf) at which to place the incoming data from process i.
        displs.clear();
        displs.resize(numProcesses, 0);
        for (int i = 1; i < numProcesses; ++i)
            displs[i] = displs[i - 1] + recvCount[i - 1];
    }

    //! @brief Projected stabilized Bi-conjugate gradient method (BiCGStab)
    //!
    //! Iterative method to solve the interface problem.
    //! @param projection: Projection matrix that guarantees that the solution satisfies the constraint at every
    //! iteration
    int ProjBiCgStab(const MatrixXd& projection, VectorXd& x, const VectorXd& rhs)
    {


        const SparseMatrix Btrans = mB.transpose();


        VectorXd Ax = CalculateFx(x);


        VectorXd Ay = Ax;
        VectorXd Az = Ax;

        const int n = rhs.rows();
        VectorXd r = rhs - Ax;

        VectorXd rProj = projection * r;
        VectorXd rProj0 = rProj;

        double rProj0_sqnorm = rProj0.squaredNorm();
        double rhs_sqnorm = rhs.squaredNorm();
        double rho = 1.;
        double alpha = 1.;
        double w = 1.;

        VectorXd v = VectorXd::Zero(n);
        VectorXd p = VectorXd::Zero(n);
        VectorXd y(n);
        VectorXd z(n);
        VectorXd s(n);
        VectorXd t(n);

        const double tol2 = mFetiTolerance * mFetiTolerance * rhs_sqnorm;

        int i = 0;
        int restarts = 0;


        while (rProj.squaredNorm() > tol2 and i < mMaxNumFetiIterations)
        {
            const double rho_old = rho;

            rho = rProj0.dot(rProj);

            if (std::fabs(rho) < 1.0e-20)
            {

                mStructure->GetLogger()
                        << "The new residual vector became too orthogonal to the arbitrarily chosen direction r0.\n"
                           "Let's restart with a new r0!"
                        << "\n\n";

                Ax = CalculateFx(x);

                rProj = projection * (rhs - Ax);
                rProj0 = rProj;
                rho = rProj0_sqnorm = rProj.squaredNorm();
                if (restarts++ == 0)
                    i = 0;
            }

            const double beta = (rho / rho_old) * (alpha / w);

            p = rProj + beta * (p - w * v);

            // precondition
            y = mPreconditioner->ApplyOnTheLeft(p);

            // reprojection
            y = projection * y;

            Ay = CalculateFx(y);

            v = projection * Ay;

            alpha = rho / rProj0.dot(v);

            s = rProj - alpha * v;

            // precondition
            z = mPreconditioner->ApplyOnTheLeft(s);

            // reprojection
            z = projection * z;

            Az = CalculateFx(z);

            t = projection * Az;

            const double tSquaredNorm = t.squaredNorm();

            if (tSquaredNorm > 0.0)
                w = t.dot(s) / tSquaredNorm;
            else
                w = 0.0;

            x += alpha * y + w * z;
            rProj = s - w * t;

            mStructure->GetLogger() << "BiCGStab rel. error = " << rProj.squaredNorm() / rhs_sqnorm
                                    << "\t at iteration = " << i << "/" << mMaxNumFetiIterations << "\n";
            ++i;
        }

        return i;
    }


    //! \brief Applies the interface flexibility operator to some vector x
    //!
    //! \param x
    //! \return sum of B_s K_s^+ B_s^T x
    VectorXd CalculateFx(const VectorXd& x)
    {

        VectorXd tmp = VectorXd::Zero(mNumTotalActiveDofs + mNumRigidBodyModes);
        tmp.head(mNumTotalActiveDofs) = mSolver.solve((mB.transpose() * x).head(mNumTotalActiveDofs));
        VectorXd Ax = mB * tmp;

        MPI_Allreduce(MPI_IN_PLACE, Ax.data(), Ax.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        return Ax;
    }


    //! @brief Conjugate projected gradient method (CG)
    //!
    //! Iterative method to solve the interface problem.
    //! @param projection: Projection matrix that guarantees that the solution satisfies the constraint at every
    //! iteration
    int ProjConjugateGradient(const Eigen::MatrixXd& projection, Eigen::VectorXd& lambda, const Eigen::VectorXd& rhs)
    {

        VectorXd tmp = CalculateFx(lambda);

        // initial projected residual
        VectorXd projResidual = projection * (rhs - tmp);

        // apply precondtioner and reproject
        VectorXd precondProjResidual = projection * mPreconditioner->ApplyOnTheLeft(projResidual);

        VectorXd z;

        const double rhs_sqnorm = rhs.squaredNorm();
        const double threshold = mFetiTolerance * mFetiTolerance * rhs_sqnorm;

        double absNew = projResidual.dot(precondProjResidual);
        int iteration = 0;
        while (iteration < mMaxNumFetiIterations)
        {
            // at every iteration i the r has to be recomputed which is quite expensive
            tmp = CalculateFx(precondProjResidual);

            // step size
            const double alpha = absNew / precondProjResidual.dot(tmp);

            // update solution
            lambda += alpha * precondProjResidual;

            // update projected residual
            projResidual -= alpha * projection * tmp;

            if (projResidual.squaredNorm() < threshold)
                break;

            // precondition and reproject
            z = projection * mPreconditioner->ApplyOnTheLeft(projResidual);

            const double absOld = absNew;
            absNew = projResidual.dot(z);
            const double beta = absNew / absOld;
            // update search direction
            precondProjResidual = z + beta * precondProjResidual;

            mStructure->GetLogger() << "ProjConjugateGradient rel. error = " << projResidual.squaredNorm() / rhs_sqnorm
                                    << "\t at iteration = " << iteration << "/" << mMaxNumFetiIterations << "\n";

            ++iteration;
        }

        return iteration;
    }


    //! @brief Projected generalized minimal residual method
    //!
    //! Iterative method to solve the interface problem.
    //! @param projection: Projection matrix that guarantees that the solution satisfies the constraint at every
    //! iteration
    int ProjGmres(const Eigen::MatrixXd& projection, Eigen::VectorXd& x, const Eigen::VectorXd& rhs)
    {

        using MatrixType = Eigen::MatrixXd;
        using VectorType = Eigen::VectorXd;

        VectorType Fx = CalculateFx(x);
        VectorType r = rhs - Fx;

        VectorType projResidual = projection * r;

        // precondition
        VectorType precondProjResidual = mPreconditioner->ApplyOnTheLeft(projResidual);

        // reproject
        precondProjResidual = projection * precondProjResidual;

        //        double projResidualNorm = projResidual.norm();
        double precondProjResidualNorm = precondProjResidual.norm();

        double rhsNorm = rhs.norm();
        if (rhsNorm < 1.e-5)
            rhsNorm = 1.0;

        int n = r.rows();

        // krylovDimension = max number of vectors that form the Krylov subspace
        int krylovDimension = 40;

        // initialize upper Hessenberg matrix
        MatrixType H = MatrixType::Zero(krylovDimension + 1, krylovDimension);

        // initialize Krylov subspace
        MatrixType V = MatrixType::Zero(n, krylovDimension + 1);

        // container for Givens rotation matrices, i.e. a vector of Matrix2d with cosines and sines
        std::vector<Eigen::JacobiRotation<double>> Givens(krylovDimension + 1);

        VectorType e1 = VectorType::Unit(n, 0);

        double tolerance = mFetiTolerance;
        int maxNumRestarts = mMaxNumFetiIterations;

        int numRestarts = 0;

        while (precondProjResidualNorm > tolerance * rhsNorm and numRestarts < maxNumRestarts)
        {
            // the first vector in the Krylov subspace is the normalized residual
            V.col(0) = precondProjResidual / precondProjResidualNorm;

            // initialize the s vector used to estimate the residual
            VectorType s = precondProjResidualNorm * e1;

            for (int i = 0; i < krylovDimension; ++i)
            {
                // calculate w = A * V.col(i)
                VectorType w = CalculateFx(V.col(i));

                // project w
                w = projection * w;

                // precondition
                w = mPreconditioner->ApplyOnTheLeft(w);

                // reproject
                w = projection * w;

                for (int iRow = 0; iRow < i + 1; ++iRow)
                {
                    H(iRow, i) = w.transpose() * V.col(iRow);
                    w = w - H(iRow, i) * V.col(iRow);
                }

                H(i + 1, i) = w.norm();
                V.col(i + 1) = w / H(i + 1, i);

                // Apply the Givens Rotations to ensure that H is an upper triangular matrix.
                // First apply previous rotations to the current matrix
                for (int iRow = 0; iRow < i; ++iRow)
                {
                    H.col(i).applyOnTheLeft(iRow, iRow + 1, Givens[iRow].adjoint());
                }

                // form the i-th rotation matrix
                Givens[i].makeGivens(H(i, i), H(i + 1, i));

                // Apply the new Givens rotation on the
                // new entry in the uppper Hessenberg matrix
                H.col(i).applyOnTheLeft(i, i + 1, Givens[i].adjoint());

                // Finally apply the new Givens rotation on the s vector
                s.applyOnTheLeft(i, i + 1, Givens[i].adjoint());

                const double error = std::abs(s[i + 1]);

                mStructure->GetLogger() << "projected GMRES relative error: \t" << error / (rhsNorm)
                                        << "\t #restarts: \t" << numRestarts << "\n";

                if (error < tolerance * rhsNorm)
                {
                    VectorType y = s.head(i + 1);
                    H.topLeftCorner(i + 1, i + 1).triangularView<Eigen::Upper>().solveInPlace(y);
                    x = x + V.block(0, 0, n, i + 1) * y;

                    // compute preconditioned residual
                    Fx = CalculateFx(x);
                    r = rhs - Fx;

                    mStructure->GetLogger() << "Residual norm: \t" << r.norm() << "\n";

                    projResidual = projection * r;

                    mStructure->GetLogger() << "Projected residual norm: \t" << projResidual.norm() << "\n";

                    // precondition
                    precondProjResidual = mPreconditioner->ApplyOnTheLeft(projResidual);

                    // reproject
                    precondProjResidual = projection * precondProjResidual;

                    precondProjResidualNorm = precondProjResidual.norm();

                    mStructure->GetLogger() << "precond. proj. residual norm: \t" << precondProjResidual.norm() << "\n";

                    WriteGmresInfo(i, numRestarts, krylovDimension);
                    return numRestarts;
                }
            }

            // we have exceeded the number of iterations. Update the approximation and start over
            VectorType y = s.head(krylovDimension);
            H.topLeftCorner(krylovDimension, krylovDimension).triangularView<Eigen::Upper>().solveInPlace(y);
            x = x + V.block(0, 0, n, krylovDimension) * y;

            // compute preconditioned residual
            Fx = CalculateFx(x);
            r = rhs - Fx;
            projResidual = projection * r;

            // precondition
            precondProjResidual = mPreconditioner->ApplyOnTheLeft(projResidual);

            // reproject
            precondProjResidual = projection * precondProjResidual;

            precondProjResidualNorm = precondProjResidual.norm();

            ++numRestarts;
        }


        std::cout << "Not converged! #restarts in projected GMRES: \t" << numRestarts << std::endl;

        return numRestarts;
    }


    //! \brief Writes information of current GMRES iteration to a file
    //! \param numIterations
    //! \param numRestarts
    //! \param krylovDimension
    void WriteGmresInfo(const int numIterations, const int numRestarts, const int krylovDimension)
    {
        std::ofstream file(mPostProcessor->GetResultDirectory() + "/GmresInfo.dat", std::ios::app);
        file << mTime << "\t" << numIterations + (numRestarts * krylovDimension) << "\t" << krylovDimension << "\n";
        file.close();
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void CheckIfGtransDeltaLambdaEqualsE(const MatrixXd& G, const VectorXd& lambda, const VectorXd& e,
                                         const double tolerance = 1.e-6) const
    {
        constexpr double precision = 1.;
        const VectorXd zero = G.transpose() * lambda - e;

        if (not(zero.isMuchSmallerThan(tolerance, precision)))
            throw Exception(__PRETTY_FUNCTION__,
                            "Gtrans * lambda - e = 0 not satisfied. Norm is: " + std::to_string(zero.norm()));
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void CheckIfRigidBodyModesAreOrthogonalToRightHandSide(const MatrixXd& rigidBodyModes, const VectorXd& residual,
                                                           const VectorXd& lambda, const double tolerance = 1.e-6) const
    {
        constexpr double precision = 1.;
        const VectorXd zero = rigidBodyModes.transpose() * (residual - mB.transpose() * lambda);

        if (not(zero.isMuchSmallerThan(tolerance, precision)))
            throw Exception(__PRETTY_FUNCTION__, "Rtrans ( f - Btrans * lambda ) = 0 not satisfied. Norm is: " +
                                                         std::to_string(zero.norm()));
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void CheckIfFtimesLambdaPLusGtimesAlphaEqualsD(const VectorXd& FtimesLambda, const MatrixXd& G,
                                                   const VectorXd& alpha, const VectorXd& d,
                                                   const double tolerance = 1.e-2) const
    {
        constexpr double precision = 1.;
        const VectorXd zero = G * alpha - (d - FtimesLambda);

        if (not(zero.isMuchSmallerThan(tolerance, precision)))
            throw Exception(__PRETTY_FUNCTION__,
                            "G*alpha - (d- F*lambda) = 0 not satisfied. Norm is: " + std::to_string(zero.norm()));
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void CheckIfProjectionOfDminusFtimesLambdaEqualsZero(const VectorXd& d, const VectorXd& FtimesLambda,
                                                         const double tolerance = 1.e-2)
    {
        const MatrixXd& P = mStructureFeti->GetProjectionMatrix();

        constexpr double precision = 1.;
        VectorXd zero = P.transpose() * (d - FtimesLambda);
        if (not(zero.isMuchSmallerThan(tolerance, precision)))
            throw Exception(__PRETTY_FUNCTION__,
                            "Ptrans * (d- F*lambda) = 0 not satisfied. Norm is: " + std::to_string(zero.norm()));
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    VectorXd ExtractLocalAlpha(const VectorXd& alphaGlobal)
    {
        const int rank = static_cast<StructureFeti*>(mStructure)->mRank;

        std::vector<int> recvCount;
        std::vector<int> displs;
        MpiGatherRecvCountAndDispls(recvCount, displs, mNumRigidBodyModes);
        return alphaGlobal.segment(displs[rank], mNumRigidBodyModes);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    VectorXd SolveForActiveDegreesOfFreedom(const VectorXd& residual, const VectorXd& lambda,
                                            const MatrixXd& rigidBodyModes, const VectorXd& alpha)
    {
        VectorXd delta_dof_active = VectorXd::Zero(mNumTotalDofs);
        delta_dof_active.head(mNumTotalActiveDofs) =
                mSolver.solve((residual - mB.transpose() * lambda).head(mNumTotalActiveDofs));

        delta_dof_active = delta_dof_active - rigidBodyModes * alpha;

        return delta_dof_active;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    VectorXd CalculateDisplacementGap(const VectorXd& residual, const double timeStep)
    {

        VectorXd pseudoInvVec = VectorXd::Zero(mNumTotalDofs);
        pseudoInvVec.head(mNumTotalActiveDofs) = mSolver.solve(residual.head(mNumTotalActiveDofs));

        VectorXd displacementGap = mB * pseudoInvVec;


        MPI_Allreduce(MPI_IN_PLACE, displacementGap.data(), displacementGap.size(), MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);


        displacementGap = displacementGap - mStructureFeti->GetPrescribedDofVector() * CalculateLoadFactor(timeStep);

        return displacementGap;
    }

    void WriteFetiSolverInfoToFile(int numIterations, double timeStep)
    {
        std::ofstream filestream(mPostProcessor->GetResultDirectory() + std::string("/FetiSolverInfo.txt"),
                                 std::ofstream::app);
        filestream << "iterations"
                   << "\t"
                   << "time step"
                   << "\n";
        filestream << numIterations << "\t" << timeStep << "\n";
        filestream.close();
    }

    //! @brief Solves the global and local problem
    //!
    //! Calculates the rigid body modes of the structure.
    //! Calculates the initial guess for the projected CG/BiCGStab method
    //! Solves for the Lagrange multipliers at the subdomain interfaces.
    //! Calculates the increment of the free degrees of freedom
    //!
    StructureOutputBlockVector FetiSolve(const VectorXd& residual_mod, const std::set<Node::eDof>& activeDofSet,
                                         VectorXd& deltaLambda, const double timeStep, bool& converged)
    {

        NuTo::Timer timer("Time for FetiSolve", mStructure->GetShowTime());


        const MatrixXd& rigidBodyModes = mStructureFeti->GetRigidBodyModes();

        // G.size() = (number of Lagrange multipliers) x (total number of rigid body modes)
        const MatrixXd& G = mStructureFeti->GetG();

        // GtransGinv.size() = (total number of rigid body modes) x (total number of rigid body modes)
        const MatrixXd GtransGinv = (G.transpose() * G).inverse();

        // R_s^T f
        VectorXd rigidBodyForceVectorLocal = rigidBodyModes.transpose() * residual_mod;


        //
        //     | R_1^T f        |
        // e = |     :          |
        //     | R_{N_s}^T f    |
        //
        VectorXd rigidBodyForceVectorGlobal =
                mFetiSolver.GatherRigidBodyForceVector(rigidBodyForceVectorLocal, mNumRigidBodyModesTotal);


        // initial guess for lambda
        deltaLambda = G * GtransGinv * rigidBodyForceVectorGlobal;


        //*****************************************
        //
        //       | K_11^-1 * r |
        // d = B |  0          |
        //
        //*****************************************
        VectorXd displacementGap = CalculateDisplacementGap(residual_mod, timeStep);

        CheckIfGtransDeltaLambdaEqualsE(G, deltaLambda, rigidBodyForceVectorGlobal);


        mStructure->GetLogger() << "residual_mod.norm() \n" << residual_mod.norm() << "\n \n";


        mStructure->GetLogger() << "*********************************** \n"
                                << "**      Start Interface Problem  ** \n"
                                << "*********************************** \n\n";
        int iterations = 0;

        switch (mIterativeSolver)
        {
        case eIterativeSolver::ConjugateGradient:
        {
            iterations = ProjConjugateGradient(mStructureFeti->GetProjectionMatrix(), deltaLambda, displacementGap);
            break;
        }
        case eIterativeSolver::BiconjugateGradientStabilized:
        {
            iterations = ProjBiCgStab(mStructureFeti->GetProjectionMatrix(), deltaLambda, displacementGap);
            break;
        }
        case eIterativeSolver::ProjectedGmres:
        {
            iterations = ProjGmres(mStructureFeti->GetProjectionMatrix(), deltaLambda, displacementGap);
            break;
        }
        default:
            throw Exception(__PRETTY_FUNCTION__, "Requested solver not yet implemented.");
        }


        if (iterations >= mMaxNumFetiIterations - 1)
        {
            converged = false;
            return StructureOutputBlockVector(mStructure->GetDofStatus());
            //            throw Exception(__PRETTY_FUNCTION__, "Maximum number of iterations exceeded.");
        }

        WriteFetiSolverInfoToFile(iterations, timeStep);

        if (mStructureFeti->mRank == 0)
            std::cout << "Iterative solver converged after iterations = \t" << iterations << "\n\n";


        mStructure->GetLogger() << "*********************************** \n"
                                << "**      End   Interface Problem  ** \n"
                                << "*********************************** \n\n";

        const VectorXd FtimesDeltaLambda = CalculateFx(deltaLambda);
        const VectorXd alphaGlobal = GtransGinv * G.transpose() * (displacementGap - FtimesDeltaLambda);


        // Gtrans * lambda = e?
        CheckIfGtransDeltaLambdaEqualsE(G, deltaLambda, rigidBodyForceVectorGlobal);

        // Rtrans * (residual - Btrans * lambda) = 0?
        CheckIfRigidBodyModesAreOrthogonalToRightHandSide(rigidBodyModes, residual_mod, deltaLambda);

        // F * lambda + G * alpha = d?
        CheckIfFtimesLambdaPLusGtimesAlphaEqualsD(FtimesDeltaLambda, G, alphaGlobal, displacementGap);

        // Ptrans * (d - F * lambda) = 0?
        CheckIfProjectionOfDminusFtimesLambdaEqualsZero(displacementGap, FtimesDeltaLambda);

        VectorXd alphaLocal = ExtractLocalAlpha(alphaGlobal);

        VectorXd delta_dof_active =
                SolveForActiveDegreesOfFreedom(residual_mod, deltaLambda, rigidBodyModes, alphaLocal);

        CheckContinuity(delta_dof_active, mStructureFeti->mNumInterfaceNodesTotal);

        StructureOutputBlockVector delta_dof_dt0(mStructure->GetDofStatus());
        int offset = 0;
        for (const auto& dof : activeDofSet)
        {
            const int numActiveDofs = mStructure->GetNumActiveDofs(dof);

            delta_dof_dt0.J[dof] = delta_dof_active.segment(offset, numActiveDofs);
            offset += numActiveDofs;
        }

        delta_dof_dt0.K[Node::eDof::DISPLACEMENTS] =
                delta_dof_active.segment(offset, mStructure->GetNumDependentDofs(Node::eDof::DISPLACEMENTS));

        return delta_dof_dt0;
    }


    //! \brief Checks if the field variables are continuous across the interfaces
    //! \param delta_dof_active
    //! \param numInterfaceNodesTotal
    void CheckContinuity(const VectorXd& delta_dof_active, const int numInterfaceNodesTotal)
    {
        VectorXd zeroVec = mB * delta_dof_active;

        MPI_Allreduce(MPI_IN_PLACE, zeroVec.data(), zeroVec.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        const int numLagrangeMultipliersDisplacement = numInterfaceNodesTotal * 2;
        if (not(zeroVec.head(numLagrangeMultipliersDisplacement).isMuchSmallerThan(1.e-4, 1.e-1)))
        {
            throw Exception(
                    __PRETTY_FUNCTION__,
                    "Continuity of the displacement field is not satisfied, i.e. B*u = 0 not satisfied. Norm is: " +
                            std::to_string(zeroVec.head(numLagrangeMultipliersDisplacement).norm()));
        }

        //        const int numLagrangeMultipliersCrack = numInterfaceNodesTotal;
        //
        //        if (not(zeroVec.segment(numLagrangeMultipliersDisplacement, numLagrangeMultipliersCrack)
        //                        .isMuchSmallerThan(1.e-2, 1.e-1)))
        //        {
        //            throw Exception(
        //                    __PRETTY_FUNCTION__,
        //                    "Continuity of the crack phase-field is not satisfied, i.e. B*d = 0 not satisfied. Norm
        //                    is: " +
        //                            std::to_string(
        //                                    zeroVec.segment(numLagrangeMultipliersDisplacement,
        //                                    numLagrangeMultipliersCrack)
        //                                            .norm()));
        //        }
    }

    //! \brief Check the rank of a BlockSparseMatrix using a QR factorization
    //! \param blockSparseMatrix
    void CheckRankOfBlockSparseMatrix(const BlockSparseMatrix& blockSparseMatrix)
    {
        if (mStructure->GetVerboseLevel() > 10)
        {
            mStructure->GetLogger() << "Start QR factorization to check the rank of the matrix"
                                    << "\n\n";

            Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> qr(
                    blockSparseMatrix.ExportToEigenSparseMatrix());

            mStructure->GetLogger() << "Number of rows:  \t" << qr.rows() << "\n\n";
            mStructure->GetLogger() << "Rank:            \t" << qr.rank() << "\n\n";
            mStructure->GetLogger() << "Rank deficiency: \t" << qr.rows() - qr.rank() << "\n\n";
        }
        else
        {
            mStructure->GetLogger() << "If you want to check the rank of a matrix, set verbose level >10"
                                    << "\n\n";
        }
    }

    void BuildAndFactorizeMatrix(StructureOutputBlockMatrix& hessian0, StructureOutputBlockMatrix& hessian1)
    {

        if (mStructure->GetNumTimeDerivatives() >= 1)
            hessian0.AddScal(hessian1, mGamma / (mBeta * mTimeStep));

        CheckRankOfBlockSparseMatrix(hessian0.JJ);

        Timer timerBenchmark("Time for factorization", mStructure->GetShowTime());
        // K_{JJ}^{-1}
        mSolver.compute(hessian0.JJ.ExportToEigenSparseMatrix());
    }

    double KeepOrIncreaseTimeStep(const int iteration)
    {
        // maybe increase next time step
        if (mAutomaticTimeStepping && iteration < 0.25 * mMaxNumIterations)
        {
            mTimeStep *= 1.5;
            if (mTimeStep > mMaxTimeStep)
                mTimeStep = mMaxTimeStep;

            return mTimeStep;
        }
        else
        {
            return mTimeStep;
        }
    }

    double ReduceTimeStep(double& curTime)
    {
        // no convergence
        if (mAutomaticTimeStepping)
        {
            // no convergence, reduce the time step and start from scratch
            curTime -= mTimeStep;
            mTimeStep *= 0.5;
            if (mTimeStep < mMinTimeStep)
            {
                mStructure->GetLogger() << "The minimal time step achieved, the actual time step is " << mTimeStep
                                        << "\n";
                throw Exception(__PRETTY_FUNCTION__, "No convergence, the current time step is too short.");
            }
        }
        else
        {
            throw Exception(__PRETTY_FUNCTION__, "No convergence with the current maximum "
                                                 "number of iterations, either use automatic "
                                                 "time stepping, reduce the time step or the "
                                                 "minimal line search cut back factor.");
        };

        return mTimeStep;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void CalculateScalingMatrix()
    {
        mScalingMatrix.resize(mNumLagrangeMultipliers, mNumLagrangeMultipliers);
        switch (mFetiScaling)
        {
        case eFetiScaling::None:
        {
            for (int i = 0; i < mNumLagrangeMultipliers; ++i)
                mScalingMatrix.insert(i, i) = 1;

            break;
        }
        case eFetiScaling::Multiplicity:
        {
            mScalingMatrix = mStructureFeti->MultiplicityScaling();
            break;
        }
        case eFetiScaling::Superlumped:
        {
            mScalingMatrix = mStructureFeti->MultiplicityScaling();
            break;
        }
        default:
            throw Exception(__PRETTY_FUNCTION__, "Scaling method not implemented.");
        }
    }

    //! @brief perform the time integration
    //! @param timeFinal ... length of the simulation
    void Solve(double timeFinal) override
    {

        mTimeControl.SetTimeFinal(timeFinal);

        mStructure->NodeBuildGlobalDofs(__PRETTY_FUNCTION__);

        mToleranceResidual.DefineDefaultValueToIninitializedDofTypes(mToleranceForce);

        mStructure->GetLogger() << "Total number of DOFs: \t" << mStructure->GetNumTotalDofs() << "\n\n";

        mStructure->GetLogger() << "******************************************** \n"
                                << "**        calculate RBMs                  ** \n"
                                << "******************************************** \n\n";

        mStructureFeti->CalculateRigidBodyModesTotalFETI();
        mNumRigidBodyModes = mStructureFeti->GetNumRigidBodyModes();
        mNumRigidBodyModesTotal = mStructureFeti->GetNumRigidBodyModesTotal();

        mStructure->GetLogger() << "******************************************** \n"
                                << "**        calculate connectivity matrix   ** \n"
                                << "******************************************** \n\n";

        mStructureFeti->AssembleConnectivityMatrix();
        mB = mStructureFeti->GetConnectivityMatrix();
        mNumLagrangeMultipliers = mStructureFeti->mNumLagrangeMultipliers;

        mStructure->GetLogger() << "******************************************** \n"
                                << "**        calculate interface RBMs        ** \n"
                                << "******************************************** \n\n";

        mStructureFeti->CalculateInterfaceRigidBodyModes();

        mStructure->GetLogger() << "******************************************** \n"
                                << "**        calculate projection            ** \n"
                                << "******************************************** \n\n";

        mStructureFeti->CalculateG();
        mStructureFeti->CalculateProjectionMatrix();

        mStructureFeti->CheckProjectionMatrix();
        mStructureFeti->CheckProjectionOfCoarseGrid();

        mStructure->GetLogger() << "******************************************** \n"
                                << "**        calculate scaling               ** \n"
                                << "******************************************** \n\n";

        CalculateScalingMatrix();


        // initialize interface forces
        mDeltaLambda.setZero(mB.rows());
        mLambda.setZero(mB.rows());
        mLambdaOld.setZero(mB.rows());

        double curTime = mTime;

        CalculateStaticAndTimeDependentExternalLoad();

        // only necessary for CheckRigidBodyModes
        const DofStatus& dofStatus = mStructure->GetDofStatus();

        if (mStepActiveDofs.empty())
        {
            mStepActiveDofs.push_back(mStructure->DofTypesGetActive());
        }
        else
        {
            for (unsigned int i = 0; i < mStepActiveDofs.size(); ++i)
            {
                if (mStepActiveDofs[i].empty())
                {
                    throw Exception(__PRETTY_FUNCTION__,
                                    "Calculation step " + std::to_string(i) + " has no active DOFs.");
                }
            }
        }


        //********************************************
        //        Allocate Variables
        //********************************************

        StructureOutputBlockMatrix hessian0(dofStatus, true);
        StructureOutputBlockMatrix hessian1(dofStatus);
        StructureOutputBlockVector delta_dof_dt0(dofStatus, true);

        StructureOutputBlockVector dof_dt0(dofStatus, true); // e.g. disp
        StructureOutputBlockVector dof_dt1(dofStatus, true); // e.g. velocity
        StructureOutputBlockVector dof_dt2(dofStatus, true); // e.g. accelerations

        StructureOutputBlockVector extForce(dofStatus, true);
        StructureOutputBlockVector intForce(dofStatus, true);

        BlockScalar normResidual(dofStatus);

        //********************************************
        //        Declare and fill Output maps
        //********************************************

        // Declare output maps

        std::map<NuTo::eStructureOutput, NuTo::StructureOutputBase*> evaluateInternalGradient;
        std::map<NuTo::eStructureOutput, NuTo::StructureOutputBase*> evaluateInternalGradientHessian0Hessian1;
        std::map<NuTo::eStructureOutput, NuTo::StructureOutputBase*> evaluateHessian0Hessian1;

        evaluateInternalGradient[eStructureOutput::INTERNAL_GRADIENT] = &intForce;
        evaluateInternalGradientHessian0Hessian1[eStructureOutput::INTERNAL_GRADIENT] = &intForce;
        evaluateInternalGradientHessian0Hessian1[eStructureOutput::HESSIAN0] = &hessian0;
        evaluateHessian0Hessian1[eStructureOutput::HESSIAN0] = &hessian0;

        if (mStructure->GetNumTimeDerivatives() >= 1 and mMuDampingMass == 0.)
        {
            hessian1.Resize(dofStatus.GetNumActiveDofsMap(), dofStatus.GetNumDependentDofsMap());
            evaluateInternalGradientHessian0Hessian1[eStructureOutput::HESSIAN1] = &hessian1;
            evaluateHessian0Hessian1[eStructureOutput::HESSIAN1] = &hessian1;
        }

        //********************************************
        //        Declare and fill Input map
        //********************************************

        ConstitutiveInputMap inputMap;
        inputMap[Constitutive::eInput::CALCULATE_STATIC_DATA] =
                std::make_unique<ConstitutiveCalculateStaticData>(eCalculateStaticData::EULER_BACKWARD);

        auto lastConvergedDofValues = ExtractDofValues();

        UpdateAndGetAndMergeConstraintRHS(curTime, lastConvergedDofValues[0]);

        PostProcess(mLambda, dofStatus);


        mStructure->GetLogger() << "******************************************** \n"
                                << "**        Start time stepping             ** \n"
                                << "******************************************** \n\n";

        while (curTime < timeFinal)
        {

            mStructure->DofTypeActivateAll();

            curTime += mTimeStep;

            SetTimeAndTimeStep(curTime, mTimeStep, timeFinal);

            extForce = CalculateCurrentExternalLoad(curTime);

            for (const auto& activeDofSet : mStepActiveDofs)
            {
                mStructure->DofTypeSetIsActive(activeDofSet);

                // ******************************************************
                mStructure->Evaluate(inputMap, evaluateInternalGradientHessian0Hessian1);
                // ******************************************************

                if (mFetiScaling == eFetiScaling::Superlumped)
                    mScalingMatrix = mStructureFeti->SuperlumpedScaling(hessian0);

                BuildAndFactorizeMatrix(hessian0, hessian1);

                VectorXd residual = CalculateResidual(extForce, intForce);

                mPreconditioner->Compute(hessian0, mB, mStructureFeti->GetLagrangeMultipliersGlobalIdToLocalId());
                mPreconditioner->AddScaling(mScalingMatrix);

                bool fetiSolverConverged = true;
                delta_dof_dt0 = FetiSolve(residual, activeDofSet, mDeltaLambda, mTimeStep, fetiSolverConverged);
                int numNewtonIterations = 0;
                if (fetiSolverConverged)
                {

                    // calculate trial state
                    dof_dt0 = lastConvergedDofValues[0] + delta_dof_dt0;
                    if (mStructure->GetNumTimeDerivatives() >= 1)
                        dof_dt1 = CalculateDof1(delta_dof_dt0, lastConvergedDofValues[1], lastConvergedDofValues[2],
                                                mTimeStep);
                    if (mStructure->GetNumTimeDerivatives() >= 2)
                        dof_dt2 = CalculateDof2(delta_dof_dt0, lastConvergedDofValues[1], lastConvergedDofValues[2],
                                                mTimeStep);


                    MergeDofValues(dof_dt0, dof_dt1, dof_dt2, false);

                    mLambda = mLambdaOld + mDeltaLambda;


                    // ******************************************************
                    mStructure->Evaluate(inputMap, evaluateInternalGradient);
                    // ******************************************************


                    residual = CalculateResidual(extForce, intForce);

                    normResidual = CalculateResidualNorm(residual, dofStatus, activeDofSet);


                    mStructure->GetLogger() << "Residual norm: \t" << normResidual << "\n\n";

                    while (not(normResidual < mToleranceResidual) and numNewtonIterations < mMaxNumIterations)
                    {


                        // ******************************************************
                        mStructure->Evaluate(inputMap, evaluateHessian0Hessian1);
                        // ******************************************************

                        BuildAndFactorizeMatrix(hessian0, hessian1);

                        delta_dof_dt0 = FetiSolve(residual, activeDofSet, mDeltaLambda, 0., fetiSolverConverged);

                        if (not fetiSolverConverged)
                            break;

                        IncrementDofs(dof_dt0, dof_dt1, dof_dt2, delta_dof_dt0);
                        mLambda += mDeltaLambda;

                        MergeDofValues(dof_dt0, dof_dt1, dof_dt2, false);


                        // ******************************************************
                        mStructure->Evaluate(inputMap, evaluateInternalGradient);
                        // ******************************************************


                        residual = CalculateResidual(extForce, intForce);

                        normResidual = CalculateResidualNorm(residual, dofStatus, activeDofSet);

                        PrintInfoIteration(normResidual, numNewtonIterations);

                        numNewtonIterations++;
                    }
                }


                if (normResidual < mToleranceResidual and fetiSolverConverged)
                {
                    // converged solution
                    mStructure->ElementTotalUpdateStaticData();

                    // store converged step
                    lastConvergedDofValues[0] = dof_dt0;
                    lastConvergedDofValues[1] = dof_dt1;
                    lastConvergedDofValues[2] = dof_dt2;
                    mLambdaOld = mLambda;

                    MergeDofValues(dof_dt0, dof_dt1, dof_dt2, true);

                    mTime += mTimeStep;


                    mStructure->GetLogger() << "Convergence after " << numNewtonIterations
                                            << " newton iterations at time " << mTime << " (time step " << mTimeStep
                                            << ").\n";

                    mTimeStep = KeepOrIncreaseTimeStep(numNewtonIterations);


                    // perform Postprocessing
                    std::ofstream file(mPostProcessor->GetResultDirectory() + "/statistics.dat", std::ios::app);
                    file << "Time: \t" << mTime << "\t # newton iterations: \t" << numNewtonIterations << "\n";
                    file.close();

                    mTimeControl.SetCurrentTime(mTime);
                    PostProcess(mLambda, dofStatus);

                    mStructure->GetLogger() << "normResidual: \n" << normResidual << "\n\n";

                    if (mCallback && mCallback->Exit(*mStructure))
                        return;
                }
                else
                {
                    mStructure->GetLogger() << "No convergence with timestep " << mTimeStep << "\n\n";

                    MergeDofValues(lastConvergedDofValues[0], lastConvergedDofValues[1], lastConvergedDofValues[2],
                                   true);

                    mTimeStep = ReduceTimeStep(curTime);
                }

            } // end for mStepActiveDofs


        } // end while
    }

    StructureOutputBlockVector CalculateDof1(const StructureOutputBlockVector& rDeltaDof_dt0,
                                             const StructureOutputBlockVector& rDof_dt1,
                                             const StructureOutputBlockVector& rDof_dt2, const double timeStep) const
    {
        return rDeltaDof_dt0 * (mGamma / (timeStep * mBeta)) + rDof_dt1 * (1. - mGamma / mBeta) +
               rDof_dt2 * (timeStep * (1. - mGamma / (2. * mBeta)));
    }


    StructureOutputBlockVector CalculateDof2(const StructureOutputBlockVector& rDeltaDof_dt0,
                                             const StructureOutputBlockVector& rDof_dt1,
                                             const StructureOutputBlockVector& rDof_dt2, const double timeStep) const
    {
        return rDeltaDof_dt0 * (1. / (timeStep * timeStep * mBeta)) - rDof_dt1 * (1. / (timeStep * mBeta)) -
               rDof_dt2 * ((0.5 - mBeta) / mBeta);
    }

    void IncrementDofs(StructureOutputBlockVector& dof_dt0, StructureOutputBlockVector& dof_dt1,
                       StructureOutputBlockVector& dof_dt2, const StructureOutputBlockVector& delta_dof_dt0) const
    {
        dof_dt0 += delta_dof_dt0;
        if (mStructure->GetNumTimeDerivatives() >= 1)
            dof_dt1 += delta_dof_dt0 * (mGamma / (mTimeStep * mBeta));
        if (mStructure->GetNumTimeDerivatives() >= 2)
            dof_dt2 += delta_dof_dt0 * (1. / (mTimeStep * mTimeStep * mBeta));
    }

    void PostProcess(const VectorXd& lambda, const DofStatus& dofStatus)
    {
        StructureOutputBlockVector outOfBalance(dofStatus, true);
        outOfBalance.J = BlockFullVector<double>((mB.transpose() * lambda).head(mNumTotalActiveDofs), dofStatus);
        mPostProcessor->PostProcess(outOfBalance);
    }


    BlockScalar CalculateResidualNorm(VectorXd& residual, const DofStatus& dofStatus,
                                      const std::set<Node::eDof>& activeDofSet) const
    {

        BlockScalar normResidual =
                BlockFullVector<double>(residual.head(mNumTotalActiveDofs), dofStatus).CalculateInfNorm();
        for (const auto& dof : activeDofSet)
            MPI_Allreduce(MPI_IN_PLACE, &normResidual[dof], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        return normResidual;
    }


    //! @brief Calculates the right hand side of the system
    //! rhs = extForce - intForce - Btrans*lambda
    VectorXd CalculateResidual(const StructureOutputBlockVector& extForce, const StructureOutputBlockVector& intForce)
    {
        return (extForce.ExportToEigenVector() - intForce.ExportToEigenVector() - mB.transpose() * mLambda);
    }

    double CalculateLoadFactor(double rTime)
    {

        // extract the two data points
        const double loadFactor1 = mTimeDependentLoadFactor(0, 1);
        const double loadFactor2 = mTimeDependentLoadFactor(1, 1);
        const double time1 = mTimeDependentLoadFactor(0, 0);
        const double time2 = mTimeDependentLoadFactor(1, 0);

        return loadFactor1 + (loadFactor2 - loadFactor1) / (time2 - time1) * (rTime - time1);
    }


    void SetIterativeSolver(eIterativeSolver iterativeSolver)
    {
        mIterativeSolver = iterativeSolver;
    }

    void SetToleranceIterativeSolver(const double tolerance)
    {
        mFetiTolerance = tolerance;
    }

    void SetFetiPreconditioner(PreconditionerType fetiPreconditioner)
    {
        mPreconditioner = std::move(fetiPreconditioner);
    }

    void SetFetiScaling(eFetiScaling fetiScaling)
    {
        mFetiScaling = fetiScaling;
    }

    void SetMaxNumberOfFetiIterations(int maxNumFetiIterations)
    {
        mMaxNumFetiIterations = maxNumFetiIterations;
    }

    // temporary fix
    //! @brief sets the  time step for the time integration procedure (initial value)
    virtual void SetTimeStep(double rTimeStep) override
    {
        mTimeStep = rTimeStep;
    }

    //! @brief returns the  time step for the time integration procedure (current value)
    virtual double GetTimeStep() const override
    {
        return mTimeStep;
    }


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //          MEMBER VARIABLES
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

protected:
    double mTime = 0.;
    double mTimeStep = 0.;
    double mMinTimeStep = 0.;
    double mMaxTimeStep = 1.;

private:
    StructureFeti* mStructureFeti;

    // auxilliary member variables
    const int mNumTotalActiveDofs;
    const int mNumTotalDofs;
    int mNumRigidBodyModes = -1337;
    SparseMatrix mB;
    SparseMatrix mScalingMatrix;

    FetiSolver mFetiSolver;
    EigenSolver mSolver;
    SparseMatrix mLocalPreconditioner;
    double mFetiTolerance = 1.0e-12;
    int mMaxNumFetiIterations = 100;
    eIterativeSolver mIterativeSolver = eIterativeSolver::ConjugateGradient;
    eFetiScaling mFetiScaling = eFetiScaling::None;

    // Lagrange multiplier
    VectorXd mDeltaLambda;
    VectorXd mLambda;
    VectorXd mLambdaOld;

    int mNumRigidBodyModesTotal = -1337;
    int mNumLagrangeMultipliers = 0;

    PreconditionerType mPreconditioner;
};
} // namespace NuTo
