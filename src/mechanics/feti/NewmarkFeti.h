#pragma once

#include <mpi.h>
#include <boost/mpi.hpp>
#include "mechanics/timeIntegration/TimeIntegrationBase.h"
#include "mechanics/timeIntegration/NewmarkDirect.h"
#include "mechanics/structures/StructureBase.h"
#include "mechanics/MechanicsException.h"
#include "base/Timer.h"

#include "mechanics/structures/StructureBase.h"
#include "mechanics/feti/StructureFeti.h"
#include "mechanics/structures/StructureBaseEnum.h"
#include "mechanics/structures/StructureOutputBlockMatrix.h"
#include "mechanics/structures/StructureOutputDummy.h"

#include "mechanics/feti/FetiSolver.h"

#include "base/CallbackInterface.h"
#include "math/SparseMatrixCSRGeneral.h"
#include <eigen3/Eigen/Dense>
#include "mechanics/constitutive/ConstitutiveEnum.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveCalculateStaticData.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveIOMap.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveTimeStep.h"
#include "FetiSolver.h"
#include <eigen3/Eigen/Dense>

#include <eigen3/Eigen/Sparse>
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

    enum class eFetiPreconditioner
    {
        None,
        Lumped,
        Dirichlet
    };

    using VectorXd = Eigen::VectorXd;
    using MatrixXd = Eigen::MatrixXd;
    using SparseMatrix = Eigen::SparseMatrix<double>;
    ///
    /// \brief NewmarkFeti
    /// \param rStructure
    ///
    NewmarkFeti(StructureFeti* rStructure)
        : NewmarkDirect(rStructure)
        , mNumTotalActiveDofs(rStructure->GetNumTotalActiveDofs())
        , mNumTotalDofs(rStructure->GetNumTotalDofs())
    {
    }

    ///
    /// \brief HasCriticalTimeStep
    /// \return
    ///
    bool HasCriticalTimeStep() const override
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "Not implemented!");
    }

    ///
    /// \brief CalculateCriticalTimeStep
    /// \return
    ///
    double CalculateCriticalTimeStep() const override
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "Not implemented!");
    }


    ///
    /// \brief GatherInterfaceRigidBodyModes
    /// \param interfaceRigidBodyModes
    /// \param numRigidBodyModesGlobal
    /// \return
    ///
    Eigen::MatrixXd GatherInterfaceRigidBodyModes(Eigen::MatrixXd& interfaceRigidBodyModes,
                                                  const int numRigidBodyModesGlobal)
    {

        std::vector<int> recvCount;
        std::vector<int> displ;
        MpiGatherRecvCountAndDispls(recvCount, displ, interfaceRigidBodyModes.size());

        const int numInterfaceEqs = interfaceRigidBodyModes.rows();
        MatrixXd interfaceRigidBodyModesGlobal = MatrixXd::Zero(numInterfaceEqs, numRigidBodyModesGlobal);

        MPI_Allgatherv(interfaceRigidBodyModes.data(), interfaceRigidBodyModes.size(), MPI_DOUBLE,
                       interfaceRigidBodyModesGlobal.data(), recvCount.data(), displ.data(), MPI_DOUBLE,
                       MPI_COMM_WORLD);

        return interfaceRigidBodyModesGlobal;
    }

    ///
    /// \brief GatherRigidBodyForceVector
    /// \param rigidBodyForceVectorLocal
    /// \param numRigidBodyModesGlobal
    /// \return
    ///
    Eigen::VectorXd GatherRigidBodyForceVector(Eigen::VectorXd& rigidBodyForceVectorLocal,
                                               const int numRigidBodyModesGlobal)
    {

        const int numRigidBodyModesLocal = rigidBodyForceVectorLocal.rows();
        std::vector<int> recvCount;
        std::vector<int> displ;
        MpiGatherRecvCountAndDispls(recvCount, displ, numRigidBodyModesLocal);

        VectorXd rigidBodyForceVectorGlobal = VectorXd::Zero(numRigidBodyModesGlobal);
        MPI_Allgatherv(rigidBodyForceVectorLocal.data(), rigidBodyForceVectorLocal.size(), MPI_DOUBLE,
                       rigidBodyForceVectorGlobal.data(), recvCount.data(), displ.data(), MPI_DOUBLE, MPI_COMM_WORLD);

        return rigidBodyForceVectorGlobal;
    }
    ///
    /// \brief MpiGatherRecvCountAndDispls
    /// \param recvCount
    /// \param displs
    /// \param numValues
    ///
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
            y = ApplyPreconditionerOnTheLeft(p);

            // reprojection
            y = projection * y;

            Ay = CalculateFx(y);

            v = projection * Ay;

            alpha = rho / rProj0.dot(v);

            s = rProj - alpha * v;

            // precondition
            z = ApplyPreconditionerOnTheLeft(s);

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


    VectorXd CalculateFx(const VectorXd& x)
    {
        StructureFeti* structure = static_cast<StructureFeti*>(mStructure);

        const SparseMatrix& B = structure->GetConnectivityMatrix();
        const SparseMatrix Btrans = B.transpose();

        const int numActiveDofs = mStructure->GetNumTotalActiveDofs();
        const int numRigidBodyModes = structure->GetNumRigidBodyModes();

        VectorXd tmp = VectorXd::Zero(numActiveDofs + numRigidBodyModes);
        tmp.head(numActiveDofs) = mSolver.solve((Btrans * x).head(numActiveDofs));
        VectorXd Ax = B * tmp;

        MPI_Allreduce(MPI_IN_PLACE, Ax.data(), Ax.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        return Ax;
    }

    VectorXd ApplyPreconditionerOnTheLeft(const VectorXd& x)
    {
        VectorXd vec = mLocalPreconditioner * x;
        MPI_Allreduce(MPI_IN_PLACE, vec.data(), vec.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        return vec;
    }

    //! @brief Conjugate projected gradient method (CG)
    //!
    //! Iterative method to solve the interface problem.
    //! @param projection: Projection matrix that guarantees that the solution satisfies the constraint at every
    //! iteration
    int ProjConjugateGradient(const Eigen::MatrixXd& projection, Eigen::VectorXd& x, const Eigen::VectorXd& rhs)
    {

        VectorXd tmp = CalculateFx(x);

        // initial projected residual
        VectorXd projResidual = projection * (rhs - tmp);

        // apply precondtioner and reproject
        VectorXd precondProjResidual = projection * ApplyPreconditionerOnTheLeft(projResidual);

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
            x += alpha * precondProjResidual;

            // update projected residual
            projResidual -= alpha * projection * tmp;

            if (projResidual.squaredNorm() < threshold)
                break;

            // precondition and reproject
            z = projection * ApplyPreconditionerOnTheLeft(projResidual);

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

        // calculates Fx = sum(Btrans * K^+ * B * x)
        VectorType Fx = CalculateFx(x);
        VectorType r = rhs - Fx;

        // initial projected search direction
        VectorType projResidual = projection * r;

        // precondition
        VectorType precondProjResidual = ApplyPreconditionerOnTheLeft(projResidual);

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

        double tolerance = 1.e-12;
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
                w = ApplyPreconditionerOnTheLeft(w);

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

                //! TODO: Think about why this is std::abs and not norm
                const double error = std::abs(s[i + 1]);

                mStructure->GetLogger() << "projected GMRES relative error: \t" << error / (rhsNorm)
                                        << "\t #restarts: \t" << numRestarts << "\n";

                if (error < tolerance * rhsNorm)
                {
                    VectorType y = s.head(i + 1);
                    H.topLeftCorner(i + 1, i + 1).triangularView<Eigen::Upper>().solveInPlace(y);
                    x = x + V.block(0, 0, n, i + 1) * y;
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
            precondProjResidual = ApplyPreconditionerOnTheLeft(projResidual);

            // reproject
            precondProjResidual = projection * precondProjResidual;

            precondProjResidualNorm = precondProjResidual.norm();

            ++numRestarts;
        }


        std::cout << "Not converged! #restarts in projected GMRES: \t" << numRestarts << std::endl;


        return numRestarts;
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

        StructureFeti* structure = static_cast<StructureFeti*>(mStructure);


        const SparseMatrix& B = structure->GetConnectivityMatrix();
        const SparseMatrix Btrans = B.transpose();

        const int numRigidBodyModesLocal = structure->GetNumRigidBodyModes();
        const int numRigidBodyModesGlobal = structure->GetNumRigidBodyModesTotal();

        const int numTotalActiveDofs = mStructure->GetNumTotalActiveDofs();
        const int numTotalDofs = mStructure->GetNumTotalDofs();

        const MatrixXd& rigidBodyModes = structure->GetRigidBodyModes();


        // G.size() = (number of Lagrange multipliers) x (total number of rigid body modes)
        const MatrixXd& G = structure->GetG();


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
                mFetiSolver.GatherRigidBodyForceVector(rigidBodyForceVectorLocal, numRigidBodyModesGlobal);


        // initial guess for lambda
        deltaLambda = G * GtransGinv * rigidBodyForceVectorGlobal;


        //*****************************************
        //
        //       | K_11^-1 * r |
        // d = B |  0          |
        //
        //*****************************************
        VectorXd pseudoInvVec;
        pseudoInvVec.setZero(numTotalDofs);


        pseudoInvVec.head(numTotalActiveDofs) = mSolver.solve(residual_mod.head(numTotalActiveDofs));


        VectorXd displacementGap = B * pseudoInvVec;


        MPI_Allreduce(MPI_IN_PLACE, displacementGap.data(), displacementGap.size(), MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);


        //        std::cout << "time step:                       \t" << timeStep << std::endl;
        //        std::cout << "CalculateLoadFactor load factor: \t" << CalculateLoadFactor(timeStep) << std::endl;
        displacementGap = displacementGap - structure->GetPrescribedDofVector() * CalculateLoadFactor(timeStep);


        VectorXd zeroVec = G.transpose() * deltaLambda - rigidBodyForceVectorGlobal;


        if (not(zeroVec.isMuchSmallerThan(1.e-6, 1.e-1)))
            throw MechanicsException(__PRETTY_FUNCTION__, "Gtrans * lambda - e = 0 not satisfied. Norm is: " +
                                                                  std::to_string(zeroVec.norm()));

        mStructure->GetLogger() << "residual_mod.norm() \n" << residual_mod.norm() << "\n \n";


        structure->GetLogger() << "*********************************** \n"
                               << "**      Start Interface Problem  ** \n"
                               << "*********************************** \n\n";

        int iterations = 0;

        switch (mIterativeSolver)
        {
        case eIterativeSolver::ConjugateGradient:
        {
            iterations = ProjConjugateGradient(structure->GetProjectionMatrix(), deltaLambda, displacementGap);
            break;
        }

        case eIterativeSolver::BiconjugateGradientStabilized:
        {
            iterations = ProjBiCgStab(structure->GetProjectionMatrix(), deltaLambda, displacementGap);
            break;
        }
        case eIterativeSolver::ProjectedGmres:
        {
            iterations = ProjGmres(structure->GetProjectionMatrix(), deltaLambda, displacementGap);
            break;
        }
        default:
            throw MechanicsException(__PRETTY_FUNCTION__, "Requested solver not yet implemented.");
        }


        if (iterations >= mMaxNumFetiIterations - 1)
        {
            converged = false;
            return StructureOutputBlockVector(structure->GetDofStatus());
            //            throw MechanicsException(__PRETTY_FUNCTION__, "Maximum number of iterations exceeded.");
        }


        if (structure->mRank == 0)
            std::cout << "Iterative solver converged after iterations = \t" << iterations << "\n\n";

        structure->GetLogger() << "*********************************** \n"
                               << "**      End   Interface Problem  ** \n"
                               << "*********************************** \n\n";


        // MPI_Barrier(MPI_COMM_WORLD);

        zeroVec = G.transpose() * deltaLambda - rigidBodyForceVectorGlobal;

        if (not(zeroVec.isMuchSmallerThan(1.e-6, 1.e-1)))
        {
            throw MechanicsException(__PRETTY_FUNCTION__, "Gtrans * lambda - e = 0 not satisfied. Norm is: " +
                                                                  std::to_string(zeroVec.norm()));
        }


        zeroVec = rigidBodyModes.transpose() * (residual_mod - B.transpose() * deltaLambda);
        if (not(zeroVec.isMuchSmallerThan(1.e-6, 1.e-1)))
            throw MechanicsException(__PRETTY_FUNCTION__, "Rtrans ( f - Btrans*lambda ) = 0 not satisfied. Norm is: " +
                                                                  std::to_string(zeroVec.norm()));


        //*****************************************
        //
        //         | K_11^-1 * Btrans*lambda    |
        // tmp = B |  0                         |
        //
        //*****************************************
        pseudoInvVec.setZero(mStructure->GetNumTotalDofs());
        pseudoInvVec.head(mStructure->GetNumTotalActiveDofs()) =
                mSolver.solve((Btrans * deltaLambda).head(mStructure->GetNumTotalActiveDofs()));

        VectorXd tmp = B * pseudoInvVec;
        MPI_Allreduce(MPI_IN_PLACE, tmp.data(), tmp.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        VectorXd alphaGlobal = VectorXd::Zero(numRigidBodyModesGlobal);
        alphaGlobal = GtransGinv * G.transpose() * (displacementGap - tmp);


        // MPI_Barrier(MPI_COMM_WORLD);
        //        VectorXd zeroVecTmp = G * alphaGlobal - (displacementGap - tmp);
        //        if (not(zeroVecTmp.isMuchSmallerThan(1.e-4, 1.e-1)))
        //        {
        //            structure->GetLogger() << "G*alpha - (d- F*lambda) = 0 \n" << zeroVecTmp << "\n\n";
        //            throw MechanicsException(__PRETTY_FUNCTION__, "G*alpha - (d- F*lambda) = 0 not satisfied");
        //        }
        //
        //        VectorXd zeroVecTmp2 = structure->GetProjectionMatrix().transpose() * (displacementGap - tmp);
        //        if (not(zeroVecTmp2.isMuchSmallerThan(1.e-4, 1.e-1)))
        //            throw MechanicsException(__PRETTY_FUNCTION__, "Ptrans * (d- F*lambda) = 0 not satisfied");


        // MPI_Barrier(MPI_COMM_WORLD);


        std::vector<int> recvCount;
        std::vector<int> displs;
        MpiGatherRecvCountAndDispls(recvCount, displs, numRigidBodyModesLocal);
        VectorXd alphaLocal = alphaGlobal.segment(displs[structure->mRank], numRigidBodyModesLocal);

        //    structure->GetLogger() << "alphaGlobal: \n" << alphaGlobal                << "\n\n";
        //    structure->GetLogger() << "alphaLocal: \n" << alphaLocal                << "\n\n";
        //    structure->GetLogger() << "rigidBodyModes: \n" << rigidBodyModes                << "\n\n";

        VectorXd delta_dof_active;
        delta_dof_active.setZero(mStructure->GetNumTotalDofs());
        delta_dof_active.head(mStructure->GetNumTotalActiveDofs()) =
                mSolver.solve((residual_mod - Btrans * deltaLambda).head(structure->GetNumTotalActiveDofs()));


        delta_dof_active = delta_dof_active - rigidBodyModes * alphaLocal;

        //        zeroVec = B * delta_dof_active;
        //
        //        structure->GetLogger() << "B_s * u_s: \n" << zeroVec << "\n\n";
        //
        //        MPI_Allreduce(MPI_IN_PLACE, zeroVec.data(), zeroVec.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        //
        //        structure->GetLogger() << "Sum B * u: \n" << zeroVec << "\n\n";
        //        structure->GetLogger() << "delta lambda: \n" << deltaLambda << "\n\n";


        StructureOutputBlockVector delta_dof_dt0(structure->GetDofStatus());
        int offset = 0;
        for (const auto& dof : activeDofSet)
        {
            const int numActiveDofs = structure->GetNumActiveDofs(dof);

            delta_dof_dt0.J[dof] = delta_dof_active.segment(offset, numActiveDofs);
            offset += numActiveDofs;
        }

        delta_dof_dt0.K[Node::eDof::DISPLACEMENTS] =
                delta_dof_active.segment(offset, structure->GetNumDependentDofs(Node::eDof::DISPLACEMENTS));
        //    structure->GetLogger() << "delta_dof_active: \n" << delta_dof_active                << "\n\n";
        //    structure->GetLogger() << "delta_dof_active.size(): \n" << delta_dof_active.size()  << "\n\n";
        //    structure->GetLogger() << "delta_dof_dt0.J: \n" << delta_dof_dt0.J                  << "\n\n";
        //    structure->GetLogger() << "delta_dof_dt0.K: \n" << delta_dof_dt0.K                  << "\n\n";


        return delta_dof_dt0;
    }


    void CalculateLocalPreconditioner(const StructureOutputBlockMatrix& hessian)
    {
        mStructure->GetLogger() << "******************************************** \n"
                                << "**        calculate local preconditioner  ** \n"
                                << "******************************************** \n\n";

        Timer timer("Time to calculate the local preconditioner", mStructure->GetShowTime(), mStructure->GetLogger());

        const int numTotalDofs = mStructure->GetNumTotalDofs();
        const auto& B = static_cast<StructureFeti*>(mStructure)->GetConnectivityMatrix();

        mLocalPreconditioner.resize(B.rows(), B.rows());

        switch (mFetiPreconditioner)
        {
        case eFetiPreconditioner::None:
        {
            mLocalPreconditioner.setIdentity();
            break;
        }
        case eFetiPreconditioner::Lumped:
        {
            auto K = hessian.JJ.ExportToEigenSparseMatrix();
            K.conservativeResize(numTotalDofs, numTotalDofs);

            mLocalPreconditioner = B * K * B.transpose();

            break;
        }
        case eFetiPreconditioner::Dirichlet:
        {

            //            SparseMatrix Hkk = hessian.KK.ExportToEigenSparseMatrix();
            //            mStructure->GetLogger() << "Hkk \n" << Hkk << "\n\n";
            SparseMatrix H = hessian.ExportToEigenSparseMatrix();
            //            mStructure->GetLogger() << "H \n" << H << "\n\n";

            //            mStructure->GetLogger() << "Searching for seg fault \n\n";


            std::vector<int> lagrangeMultiplierDofIds =
                    static_cast<StructureFeti*>(mStructure)->CalculateLagrangeMultiplierIds();


            std::set<int> internalDofIds;

            for (int j = 0; j < mStructure->GetNumTotalDofs(); ++j)
            {
                internalDofIds.insert(j);
            }


            for (const auto& id : lagrangeMultiplierDofIds)
            {
                internalDofIds.erase(id);
            }


            //            mStructure->GetLogger() << "LagrangeMultiplierDOfIds    "
            //                                    << "\n";
            //            for (const auto& i : lagrangeMultiplierDofIds)
            //                mStructure->GetLogger() << i << "\n";
            //
            //            mStructure->GetLogger() << "InternalDofIds    "
            //                                    << "\n";
            //
            //            for (const auto& i : internalDofIds)
            //                mStructure->GetLogger() << i << "\n";

            std::vector<int> internalDofIdsVec(internalDofIds.begin(), internalDofIds.end());

            mStructure->GetLogger() << "mLocalPreconditioner.rows() \n" << mLocalPreconditioner.rows() << "\n\n";
            mStructure->GetLogger() << "mLocalPreconditioner.cols() \n" << mLocalPreconditioner.cols() << "\n\n";

            //            mStructure->GetLogger() << "\n\n";
            //
            //            mStructure->GetLogger() << "InternalDofIdsVec    "
            //                                    << "\n";
            //
            //            for (const auto& i : internalDofIdsVec)
            //                mStructure->GetLogger() << i << "\n";


            Timer extrac("Time to extract");
            mStructure->GetLogger() << "Extract Kbb"
                                    << "\n\n";
            SparseMatrix Kbb = ExtractSubMatrix(H, lagrangeMultiplierDofIds, lagrangeMultiplierDofIds);
            mStructure->GetLogger() << "Extract Ki"
                                    << "\n\n";
            SparseMatrix Kii = ExtractSubMatrix(H, internalDofIdsVec, internalDofIdsVec);
            mStructure->GetLogger() << "Extract Kbi"
                                    << "\n\n";
            SparseMatrix Kbi = ExtractSubMatrix(H, lagrangeMultiplierDofIds, internalDofIdsVec);
            mStructure->GetLogger() << "Extract Kib"
                                    << "\n\n";
            SparseMatrix Kib = ExtractSubMatrix(H, internalDofIdsVec, lagrangeMultiplierDofIds);

            extrac.Reset();

            //            mStructure->GetLogger() << "Searching for seg fault \n\n";
            //            mStructure->GetLogger() << "Kbb \n" << Kbb << "\n\n";
            //            mStructure->GetLogger() << "Kii \n" << Kii << "\n\n";
            //            mStructure->GetLogger() << "Kbi \n" << Kbi << "\n\n";
            //            mStructure->GetLogger() << "Kib \n" << Kib << "\n\n";

            Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> luKii(Kii);
            SparseMatrix Sbb = Kbb - Kbi * luKii.solve(Kib);
            //            mStructure->GetLogger() << "Sbb \n" << Sbb << "\n\n";

            //
            //     | 0  0   |
            // S = | 0  Sbb |
            //
            SparseMatrix S(mStructure->GetNumTotalDofs(), mStructure->GetNumTotalDofs());


            for (size_t rowId = 0; rowId < lagrangeMultiplierDofIds.size(); ++rowId)
                for (size_t colId = 0; colId < lagrangeMultiplierDofIds.size(); ++colId)
                    S.insert(lagrangeMultiplierDofIds[rowId], lagrangeMultiplierDofIds[colId]) =
                            Sbb.coeff(rowId, colId);


            mLocalPreconditioner = B * S * B.transpose();

            //            mStructure->GetLogger() << "mLocalPrecon \n" << mLocalPreconditioner << "\n\n";


            MPI_Barrier(MPI_COMM_WORLD);

            //            throw MechanicsException(__PRETTY_FUNCTION__, "Dirichlet preconditioner is not implemented
            //            yet.");

            break;
        }
        default:
            throw MechanicsException(__PRETTY_FUNCTION__, "Preconditioner is not implemented yet.");
        }
    }

    template <class A, class B>
    Eigen::SparseMatrix<double> ExtractSubMatrix(const Eigen::SparseMatrix<double>& mat, const A& rowIds,
                                                 const B& colIds)
    {
        Eigen::SparseMatrix<double> subMatrix(rowIds.size(), colIds.size());

        mStructure->GetLogger() << "rowIds.size() \n"
                                << rowIds.size() << "\ncolIds.size() \n"
                                << colIds.size() << "\n\n";

        //        for (size_t rowId = 0; rowId < rowIds.size(); ++rowId)
        //            for (size_t colId = 0; colId < colIds.size(); ++colId)
        //                subMatrix.insert(rowId, colId) = mat.coeff(rowIds[rowId], colIds[colId]);


        //        mStructure->GetLogger() << "Lgrange multiplier ids" << "\n";
        //        for (const auto& i : rowIds)
        //            mStructure->GetLogger() << i << "\n";


        for (int k = 0; k < mat.outerSize(); ++k)
        {
            auto colIt = std::find(colIds.begin(), colIds.end(), k);
            if (colIt != colIds.end())
            {
                for (typename Eigen::SparseMatrix<double>::InnerIterator it(mat, k); it; ++it)
                {
                    auto rowIt = std::find(rowIds.begin(), rowIds.end(), it.row());
                    if (rowIt != rowIds.end())
                    {
                        int rowId = std::distance(rowIds.begin(), rowIt);
                        int colId = std::distance(colIds.begin(), colIt);
                        //                        mStructure->GetLogger() << "rowId \t" << rowId
                        //                                                << "\t it.row() \t" << it.row() <<"\n"
                        //                                                << "colId \t" << colId
                        //                                                << "\t it.col() \t" << it.col() <<"\n";
                        subMatrix.insert(rowId, colId) = it.value();

                        ++rowId;
                    }
                }
            }
        }


        //        mStructure->GetLogger() << "submatrix2 \n" << subMatrix2 << "\n\n";

        return subMatrix;
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
                throw MechanicsException(__PRETTY_FUNCTION__, "No convergence, the current time step is too short.");
            }
        }
        else
        {
            throw MechanicsException(__PRETTY_FUNCTION__, "No convergence with the current maximum "
                                                          "number of iterations, either use automatic "
                                                          "time stepping, reduce the time step or the "
                                                          "minimal line search cut back factor.");
        };

        return mTimeStep;
    }

    //! @brief perform the time integration
    //! @param rTimeDelta ... length of the simulation
    void Solve(double rTimeDelta) override
    {

        try
        {
            StructureFeti* structure = static_cast<StructureFeti*>(mStructure);
            mStructure->NodeBuildGlobalDofs(__PRETTY_FUNCTION__);

            mStructure->GetLogger() << "Total number of DOFs: \t" << mStructure->GetNumTotalDofs() << "\n\n";

            mStructure->GetLogger() << "******************************************** \n"
                                    << "**        calculate RBMs                  ** \n"
                                    << "******************************************** \n\n";

            structure->CalculateRigidBodyModesTotalFETI();
            mNumRigidBodyModes = structure->GetNumRigidBodyModes();

            mStructure->GetLogger() << "******************************************** \n"
                                    << "**        calculate connectivity matrix   ** \n"
                                    << "******************************************** \n\n";

            structure->AssembleConnectivityMatrix();
            mB = structure->GetConnectivityMatrix();

            mStructure->GetLogger() << "******************************************** \n"
                                    << "**        calculate interface RBMs        ** \n"
                                    << "******************************************** \n\n";

            structure->CalculateInterfaceRigidBodyModes();

            mStructure->GetLogger() << "******************************************** \n"
                                    << "**        calculate projection            ** \n"
                                    << "******************************************** \n\n";

            structure->CalculateG();
            structure->CalculateProjectionMatrix();

            structure->CheckProjectionMatrix();
            structure->CheckProjectionOfCoarseGrid();

            // initialize interface forces
            mDeltaLambda.setZero(mB.rows());
            mLambda.setZero(mB.rows());
            mLambdaOld.setZero(mB.rows());

            double curTime = mTime;

            CalculateStaticAndTimeDependentExternalLoad();

            // sets hasInteractingConstraints to true to also assemble K_KJ and K_KK
            // only necessary for CheckRigidBodyModes
            mStructure->DofStatusSetHasInteractingConstraints(true);
            const DofStatus& dofStatus = mStructure->GetDofStatus();

            assert(dofStatus.HasInteractingConstraints() == true);

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
                        throw MechanicsException(__PRETTY_FUNCTION__,
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

            StructureOutputBlockVector lastConverged_dof_dt0(dofStatus, true); // e.g. disp
            StructureOutputBlockVector lastConverged_dof_dt1(dofStatus, true); // e.g. velocity
            StructureOutputBlockVector lastConverged_dof_dt2(dofStatus, true); // e.g. accelerations

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

            ExtractDofValues(lastConverged_dof_dt0, lastConverged_dof_dt1, lastConverged_dof_dt2);

            UpdateAndGetAndMergeConstraintRHS(curTime, lastConverged_dof_dt0);

            PostProcess(mLambda, dofStatus);


            mStructure->GetLogger() << "******************************************** \n"
                                    << "**        Start time stepping             ** \n"
                                    << "******************************************** \n\n";

            while (curTime < rTimeDelta)
            {

                mStructure->DofTypeActivateAll();

                curTime += mTimeStep;

                SetTimeAndTimeStep(curTime, mTimeStep, rTimeDelta);

                extForce = CalculateCurrentExternalLoad(curTime);

                for (const auto& activeDofSet : mStepActiveDofs)
                {
                    mStructure->DofTypeSetIsActive(activeDofSet);

                    mStructure->DofStatusSetHasInteractingConstraints(true);

                    // ******************************************************
                    mStructure->Evaluate(inputMap, evaluateInternalGradientHessian0Hessian1);
                    // ******************************************************

                    BuildAndFactorizeMatrix(hessian0, hessian1);

                    VectorXd residual = CalculateResidual(extForce, intForce);

                    CalculateLocalPreconditioner(hessian0);

                    bool fetiSolverConverged = true;
                    delta_dof_dt0 = FetiSolve(residual, activeDofSet, mDeltaLambda, mTimeStep, fetiSolverConverged);
                    int iteration = 0;
                    if (fetiSolverConverged)
                    {

                        // calculate trial state
                        dof_dt0 = lastConverged_dof_dt0 + delta_dof_dt0;
                        if (mStructure->GetNumTimeDerivatives() >= 1)
                            dof_dt1 = CalculateDof1(delta_dof_dt0, lastConverged_dof_dt1, lastConverged_dof_dt2,
                                                    mTimeStep);
                        if (mStructure->GetNumTimeDerivatives() >= 2)
                            dof_dt2 = CalculateDof2(delta_dof_dt0, lastConverged_dof_dt1, lastConverged_dof_dt2,
                                                    mTimeStep);


                        MergeDofValues(dof_dt0, dof_dt1, dof_dt2, false);

                        mLambda = mLambdaOld + mDeltaLambda;


                        // ******************************************************
                        mStructure->Evaluate(inputMap, evaluateInternalGradient);
                        // ******************************************************


                        residual = CalculateResidual(extForce, intForce);

                        normResidual = CalculateResidualNorm(residual, dofStatus, activeDofSet);

                        mStructure->GetLogger() << "Residual norm: \t" << normResidual << "\n\n";

                        while (not(normResidual < mToleranceResidual) and iteration < mMaxNumIterations)
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

                            PrintInfoIteration(normResidual, iteration);

                            iteration++;

                        } // end of while(normResidual<mToleranceForce && iteration<mMaxNumIterations)
                    }


                    if (normResidual < mToleranceResidual and fetiSolverConverged)
                    {
                        // converged solution
                        mStructure->ElementTotalUpdateStaticData();

                        // store converged step
                        lastConverged_dof_dt0 = dof_dt0;
                        lastConverged_dof_dt1 = dof_dt1;
                        lastConverged_dof_dt2 = dof_dt2;
                        mLambdaOld = mLambda;

                        MergeDofValues(dof_dt0, dof_dt1, dof_dt2, true);

                        mTime += mTimeStep;


                        mStructure->GetLogger() << "Convergence after " << iteration << " iterations at time " << mTime
                                                << " (time step " << mTimeStep << ").\n";

                        mTimeStep = KeepOrIncreaseTimeStep(iteration);


                        // perform Postprocessing
                        PostProcess(mLambda, dofStatus);

                        mStructure->GetLogger() << "normResidual: \n" << normResidual << "\n\n";

                        if (mCallback && mCallback->Exit(*mStructure))
                            return;
                    }
                    else
                    {
                        mStructure->GetLogger() << "No convergence with timestep " << mTimeStep << "\n\n";

                        MergeDofValues(lastConverged_dof_dt0, lastConverged_dof_dt1, lastConverged_dof_dt2, true);

                        mTimeStep = ReduceTimeStep(curTime);
                    }

                } // end for mStepActiveDofs


            } // end while
        }
        catch (MechanicsException& e)
        {
            e.AddMessage(__PRETTY_FUNCTION__, " ERROR performing FETI.");
            throw;
        }
    }

    void IncrementDofs(StructureOutputBlockVector &dof_dt0, StructureOutputBlockVector &dof_dt1,
                       StructureOutputBlockVector &dof_dt2, const StructureOutputBlockVector &delta_dof_dt0) const
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
        TimeIntegrationBase::PostProcess(outOfBalance);
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
    VectorXd CalculateResidual(const StructureOutputBlockVector &extForce, const StructureOutputBlockVector &intForce)
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

    void SetFetiPreconditioner(eFetiPreconditioner fetiPreconditioner)
    {
        mFetiPreconditioner = fetiPreconditioner;
    }

private:
    // auxilliary member variables
    const int mNumTotalActiveDofs;
    const int mNumTotalDofs;
    int mNumRigidBodyModes = -1337;
    SparseMatrix mB;
    FetiSolver mFetiSolver;
    EigenSolver mSolver;
    SparseMatrix mLocalPreconditioner;
    double mFetiTolerance = 1.0e-4;
    const int mMaxNumFetiIterations = 100;
    eIterativeSolver mIterativeSolver = eIterativeSolver::ConjugateGradient;
    eFetiPreconditioner mFetiPreconditioner = eFetiPreconditioner::Lumped;

    // Lagrange multiplier
    VectorXd mDeltaLambda;
    VectorXd mLambda;
    VectorXd mLambdaOld;
};
} // namespace NuTo
