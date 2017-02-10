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
template <class EigenSolver = Eigen::SparseLU<Eigen::SparseMatrix<double>,Eigen::COLAMDOrdering<int>>>
class NewmarkFeti : public NewmarkDirect
{
public:

    using VectorXd      = Eigen::VectorXd;
    using MatrixXd      = Eigen::MatrixXd;
    using SparseMatrix  = Eigen::SparseMatrix<double>;
    ///
    /// \brief NewmarkFeti
    /// \param rStructure
    ///
    NewmarkFeti(StructureFeti* rStructure) : NewmarkDirect (rStructure)
    {

    }

    ///
    /// \brief GetTypeId
    /// \return
    ///
    std::string GetTypeId() const override
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "Not implemented!");
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
    Eigen::MatrixXd GatherInterfaceRigidBodyModes(Eigen::MatrixXd& interfaceRigidBodyModes, const int numRigidBodyModesGlobal)
    {

        std::vector<int> recvCount;
        std::vector<int> displ;
        MpiGatherRecvCountAndDispls(recvCount, displ, interfaceRigidBodyModes.size());

        const int numInterfaceEqs               = interfaceRigidBodyModes.rows();
        MatrixXd interfaceRigidBodyModesGlobal   = MatrixXd::Zero(numInterfaceEqs,numRigidBodyModesGlobal);

        MPI_Allgatherv(interfaceRigidBodyModes.data(),
                       interfaceRigidBodyModes.size(),
                       MPI_DOUBLE,
                       interfaceRigidBodyModesGlobal.data(),
                       recvCount.data(),
                       displ.data(),
                       MPI_DOUBLE,
                       MPI_COMM_WORLD);

        return interfaceRigidBodyModesGlobal;
    }

    ///
    /// \brief GatherRigidBodyForceVector
    /// \param rigidBodyForceVectorLocal
    /// \param numRigidBodyModesGlobal
    /// \return
    ///
    Eigen::VectorXd GatherRigidBodyForceVector(Eigen::VectorXd &rigidBodyForceVectorLocal, const int numRigidBodyModesGlobal)
    {

        const int numRigidBodyModesLocal        = rigidBodyForceVectorLocal.rows();
        std::vector<int> recvCount;
        std::vector<int> displ;
        MpiGatherRecvCountAndDispls(recvCount, displ, numRigidBodyModesLocal);

        VectorXd rigidBodyForceVectorGlobal      =VectorXd::Zero(numRigidBodyModesGlobal);
        MPI_Allgatherv(rigidBodyForceVectorLocal.data(),
                       rigidBodyForceVectorLocal.size(),
                       MPI_DOUBLE,
                       rigidBodyForceVectorGlobal.data(),
                       recvCount.data(),
                       displ.data(),
                       MPI_DOUBLE,
                       MPI_COMM_WORLD);

        return rigidBodyForceVectorGlobal;
    }
    ///
    /// \brief MpiGatherRecvCountAndDispls
    /// \param recvCount
    /// \param displs
    /// \param numValues
    ///
    void MpiGatherRecvCountAndDispls(std::vector<int> &recvCount, std::vector<int> &displs, const int numValues)
    {
        boost::mpi::communicator world;
        const int numProcesses = world.size();
        // recvCount:
        // Contais the number of elements that are received from each process.
        recvCount.clear();
        recvCount.resize(numProcesses, 0);

        boost::mpi::all_gather<int>(world,numValues,recvCount);

        // displs:
        // Entry i specifies the displacement (relative to recvbuf) at which to place the incoming data from process i.
        displs.clear();
        displs.resize(numProcesses, 0);
        for (int i = 1; i < numProcesses; ++i)
            displs[i] = displs[i-1] + recvCount[i-1];

    }

    //! @brief Projected stabilized Bi-conjugate gradient method (BiCGStab)
    //!
    //! Iterative method to solve the interface problem.
    //! @param projection: Projection matrix that guarantees that the solution satisfies the constraint at every iteration
    int BiCgStab(const MatrixXd &projection, VectorXd &x, const VectorXd &rhs)
    {

        boost::mpi::communicator world;

        StructureFeti* structure = static_cast<StructureFeti*>(mStructure);
        const SparseMatrix& B      = structure->GetConnectivityMatrix();
        const SparseMatrix Btrans = B.transpose();

        world.barrier();


        const int numActiveDofs     = mStructure->GetNumTotalActiveDofs();
        const int numRigidBodyModes = structure->GetNumRigidBodyModes();

        // solve the linear system Ax = b
        VectorXd tmp;
        tmp.setZero(numActiveDofs + numRigidBodyModes);
        tmp.head(numActiveDofs) = mSolver.solve( (Btrans*x).head(numActiveDofs));
        VectorXd Ax = B * tmp;
        MPI_Allreduce(MPI_IN_PLACE, Ax.data(), Ax.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        VectorXd& Ay = Ax;
        VectorXd& Az = Ax;

        const int n = rhs.rows();
        VectorXd r  = rhs - Ax;
//    VectorXd r0 = r;

        VectorXd rProj   = projection * r;
        VectorXd rProj0  = rProj;

        double rProj0_sqnorm = rProj0.squaredNorm();
        double rhs_sqnorm = rhs.squaredNorm();

        double rho    = 1.;
        double alpha  = 1.;
        double w      = 1.;

        VectorXd v = VectorXd::Zero(n);
        VectorXd p = VectorXd::Zero(n);
        VectorXd y(n);
        VectorXd z(n);
        VectorXd s(n);
        VectorXd t(n);

//    const double tol    = mCpgTolerance;
        //        const double tol2   = tol*tol*rhs_sqnorm;
        const double tol2   = 1.0e-7*rhs_sqnorm; // the phase-field model is worse conditioned it seems. Lower tolerance required?

        //const double eps    = std::numeric_limits<double>::epsilon();
        //const double eps2   = eps*eps;

        int i = 0;
        int restarts = 0;



        const int maxIters = mCpgMaxIterations;


        while (rProj.squaredNorm() > tol2 and i<maxIters )
        {
            double rho_old = rho;

            rho = rProj0.dot(rProj);



            if (std::fabs(rho) < 1.0e-20)
            {

                structure->GetLogger() << "The new residual vector became too orthogonal to the arbitrarily chosen direction r0.\n"
                        "Let's restart with a new r0!" << "\n\n";

                // The new residual vector became too orthogonal to the arbitrarily chosen direction r0
                // Let's restart with a new r0:

                tmp.head(numActiveDofs) = mSolver.solve( (Btrans*x).head(numActiveDofs));
                Ax.noalias() = B * tmp;

                MPI_Allreduce(MPI_IN_PLACE, Ax.data(), Ax.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);


                rProj  = projection*(rhs - Ax);
                rProj0 = rProj;
                rho = rProj0_sqnorm = rProj.squaredNorm();
                if(restarts++ == 0)
                    i = 0;
            }


            double beta = (rho/rho_old) * (alpha / w);

            p = rProj + beta * (p - w * v);


            //      y = precond.solve(p);
            y = p;

            tmp.head(numActiveDofs) = mSolver.solve( (Btrans*y).head(numActiveDofs));
            Ay.noalias() = B * tmp;

            MPI_Allreduce(MPI_IN_PLACE, Ay.data(), Ay.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

            v.noalias() = projection * Ay;

            alpha = rho / rProj0.dot(v);

            s = rProj - alpha * v;

            //      z = precond.solve(s);
            z = s;

            tmp.head(numActiveDofs) = mSolver.solve( (Btrans*z).head(numActiveDofs));
            Az.noalias() = B * tmp;

            MPI_Allreduce(MPI_IN_PLACE, Az.data(), Az.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

            t.noalias() = projection * Az;

            double tmp = t.dot(t);

            if(tmp>0.0)
                w = t.dot(s) / tmp;
            else
                w = 0.0;

            x += alpha * y + w * z;
            rProj = z - w * t;

            structure->GetLogger() << "BiCGStab rel. error = "   << rProj.squaredNorm()/rhs_sqnorm
                                   << "\t at iteration = "  << i << "/" << maxIters <<"\n";
            ++i;

        }
        //    tol_error = sqrt(r.squaredNorm()/rhs_sqnorm);
        return i;

    }



    //! @brief Conjugate projected gradient method (CG)
    //!
    //! Iterative method to solve the interface problem.
    //! @param projection: Projection matrix that guarantees that the solution satisfies the constraint at every iteration
    int CPG(const Eigen::MatrixXd &projection, Eigen::VectorXd &x, const Eigen::VectorXd &rhs)
    {

        boost::mpi::communicator world;
        VectorXd        Ap;

        StructureFeti* structure = static_cast<StructureFeti*>(mStructure);

        const SparseMatrix& B       = structure->GetConnectivityMatrix();
        const SparseMatrix Btrans   = B.transpose();

        world.barrier();

        const int numActiveDofs     = mStructure->GetNumTotalActiveDofs();
        const int numRigidBodyModes = structure->GetNumRigidBodyModes();


        VectorXd tmp;
        tmp.setZero(numActiveDofs + numRigidBodyModes);
        tmp.head(numActiveDofs) = mSolver.solve( (Btrans*x).head(numActiveDofs));
        VectorXd Ax = B * tmp;

        world.barrier();


        MPI_Allreduce(MPI_IN_PLACE, Ax.data(), Ax.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        world.barrier();
        // initial residual
        VectorXd r = rhs - Ax;
        VectorXd z;

        // initial projected search direction
        VectorXd w = projection * r;
        VectorXd p = w;

        // precondition
        structure->GetLogger() << "projection.cols() = \t"   << projection.cols()
                               << "\nmLocalPreconditioner.rows() = \t"  << mLocalPreconditioner.rows()
                               << "\nmLocalPreconditioner.cols() = \t"  << mLocalPreconditioner.cols()
                               << "\nw.rows() = \t"  << w.rows();

        p =  projection * mLocalPreconditioner * w;
        boost::mpi::all_reduce(world,boost::mpi::inplace(p.data()), p.size(),std::plus<double>());
//    p = w;

        const double rhs_sqnorm = rhs.squaredNorm();
        const double threshold = mCpgTolerance * mCpgTolerance * rhs_sqnorm;

        double absNew = w.dot(p);
        int iteration = 0;
        while(iteration < mCpgMaxIterations)
        {
            // at every iteration i the r has to be recomputd which is quite expensive
            world.barrier();

            tmp.head(numActiveDofs) = mSolver.solve( (Btrans*p).head(numActiveDofs));
            Ap.noalias() = B * tmp;

            world.barrier();
            MPI_Allreduce(MPI_IN_PLACE, Ap.data(), Ap.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

            // step size
            double alpha = absNew / p.dot(Ap);

            // update solution
            x += alpha * p;

            // update projected residual
            w -= alpha * projection * Ap;

            // precondition
            z =  projection * mLocalPreconditioner * w;
            boost::mpi::all_reduce(world,boost::mpi::inplace(z.data()), z.size(),std::plus<double>());
//        z = w;

            if (w.squaredNorm() < threshold)
                break;

            double absOld = absNew;
            absNew = w.dot(z);                        // update the absolute value of r
            double beta = absNew / absOld;            // calculate the Gram-Schmidt value used to create the new search direction
            p = z + beta * p;                         // update search direction

            structure->GetLogger() << "CPG rel. error = "   << w.squaredNorm()/rhs_sqnorm
                                   << "\t at iteration = "  << iteration << "/" << mCpgMaxIterations <<"\n";


            ++iteration;


        }

        structure->GetLogger() << "\nConverged! \n"
                               << "CPG rel. error = "   << w.squaredNorm()/rhs_sqnorm
                               << "\t at iteration = "  << iteration << "/" << mCpgMaxIterations <<"\n";

        return iteration;
    }


    //! @brief Solves the global and local problem
    //!
    //! Calculates the rigid body modes of the structure.
    //! Calculates the initial guess for the projected CG/BiCGStab method
    //! Solves for the Lagrange multipliers at the subdomain interfaces.
    //! Calculates the increment of the free degrees of freedom
    //!
    StructureOutputBlockVector FetiSolve(const VectorXd& residual_mod, const std::set<Node::eDof> &activeDofSet, VectorXd &deltaLambda, const double timeStep)
    {


        boost::mpi::communicator world;

        StructureFeti* structure = static_cast<StructureFeti*>(mStructure);



        const SparseMatrix& B                    = structure->GetConnectivityMatrix();
        const SparseMatrix  Btrans               = B.transpose();

        const int numRigidBodyModesLocal        = structure->GetNumRigidBodyModes();
        const int numRigidBodyModesGlobal       = structure->GetNumRigidBodyModesTotal();

        const int numTotalActiveDofs            = mStructure->GetNumTotalActiveDofs();
        const int numTotalDofs                  = mStructure->GetNumTotalDofs();

        const MatrixXd& rigidBodyModes              = structure->GetRigidBodyModes();




        world.barrier();

        // G.size() = (number of Lagrange multipliers) x (total number of rigid body modes)
        const MatrixXd& G          =   structure->GetG();

        // GtransGinv.size() = (total number of rigid body modes) x (total number of rigid body modes)
        const MatrixXd GtransGinv =   (G.transpose() * G).inverse();



        // R_s^T f
        VectorXd rigidBodyForceVectorLocal = rigidBodyModes.transpose() * residual_mod;

        //
        //     | R_1^T f        |
        // e = |     :          |
        //     | R_{N_s}^T f    |
        //
        VectorXd rigidBodyForceVectorGlobal = mFetiSolver.GatherRigidBodyForceVector(rigidBodyForceVectorLocal, numRigidBodyModesGlobal);

        // initial guess for lambda
        deltaLambda = G * GtransGinv * rigidBodyForceVectorGlobal;

        world.barrier();

        //*****************************************
        //
        //       | K_11^-1 * r |
        // d = B |  0          |
        //
        //*****************************************
        VectorXd pseudoInvVec;
        pseudoInvVec.setZero(numTotalDofs);
        pseudoInvVec.head(numTotalActiveDofs)   = mSolver.solve(residual_mod.head(numTotalActiveDofs));
        VectorXd displacementGap            = B * pseudoInvVec;


        boost::mpi::all_reduce(world,boost::mpi::inplace(displacementGap.data()), displacementGap.size(),std::plus<double>());

        std::cout << "time step:                       \t" << timeStep << std::endl;
        std::cout << "CalculateLoadFactor load factor: \t" << CalculateLoadFactor(timeStep) << std::endl;
        displacementGap = displacementGap - structure->GetPrescribedDofVector() * CalculateLoadFactor(timeStep);



        MPI_Barrier(MPI_COMM_WORLD);

        VectorXd zeroVec = G.transpose() * deltaLambda - rigidBodyForceVectorGlobal;


        if( not (zeroVec.isMuchSmallerThan(1.e-4, 1.e-1)) )
            throw MechanicsException(__PRETTY_FUNCTION__, "Gtrans * lambda - e = 0 not satisfied. Norm is: " + std::to_string(zeroVec.norm()));

        mStructure->GetLogger() << "residual_mod.norm() \n" << residual_mod.norm() << "\n \n";

        MPI_Barrier(MPI_COMM_WORLD);

        structure->GetLogger() << "\n*************************************\n";
        structure->GetLogger() << "       Start Interface Problem           ";
        structure->GetLogger() << "\n*************************************\n\n";

        MPI_Barrier(MPI_COMM_WORLD);

        int iteration = CPG(structure->GetProjectionMatrix(),deltaLambda,displacementGap);
//    int iteration = BiCgStab(structure->GetProjectionMatrix(),deltaLambda,displacementGap);

        if (iteration >= mCpgMaxIterations)
            throw MechanicsException(__PRETTY_FUNCTION__,"Maximum number of iterations exceeded.");

        structure->GetLogger() << "\n*************************************\n";
        structure->GetLogger() << "       End Interface Problem             ";
        structure->GetLogger() << "\n*************************************\n\n";



        MPI_Barrier(MPI_COMM_WORLD);

        zeroVec = G.transpose() * deltaLambda - rigidBodyForceVectorGlobal;

        if( not (zeroVec.isMuchSmallerThan(1.e-4, 1.e-1)) )
            throw MechanicsException(__PRETTY_FUNCTION__, "Gtrans * lambda - e = 0 not satisfied. Norm is: " + std::to_string(zeroVec.norm()));


        zeroVec = rigidBodyModes.transpose() * ( residual_mod - B.transpose() *deltaLambda );
        if( not (zeroVec.isMuchSmallerThan(1.e-4, 1.e-1)) )
            throw MechanicsException(__PRETTY_FUNCTION__, "Rtrans ( f - Btrans*lambda ) = 0 not satisfied. Norm is: " + std::to_string(zeroVec.norm()));

        MPI_Barrier(MPI_COMM_WORLD);

        //*****************************************
        //
        //         | K_11^-1 * Btrans*lambda    |
        // tmp = B |  0                         |
        //
        //*****************************************
        pseudoInvVec.setZero(mStructure->GetNumTotalDofs());
        pseudoInvVec.head(mStructure->GetNumTotalActiveDofs()) = mSolver.solve((Btrans*deltaLambda).head(mStructure->GetNumTotalActiveDofs()));

        VectorXd tmp = B * pseudoInvVec;
        boost::mpi::all_reduce(world,boost::mpi::inplace(tmp.data()), tmp.size(),std::plus<double>());

        VectorXd alphaGlobal =VectorXd::Zero(numRigidBodyModesGlobal);
        alphaGlobal = GtransGinv * G.transpose() * (displacementGap - tmp);


        MPI_Barrier(MPI_COMM_WORLD);
        VectorXd zeroVecTmp = G * alphaGlobal - (displacementGap - tmp);
        if( not(zeroVecTmp.isMuchSmallerThan(1.e-4, 1.e-1)) )
            throw MechanicsException(__PRETTY_FUNCTION__, "G*alpha - (d- F*lambda) = 0 not satisfied");

        VectorXd zeroVecTmp2 = structure->GetProjectionMatrix().transpose() * (displacementGap - tmp);
        if( not(zeroVecTmp2.isMuchSmallerThan(1.e-4, 1.e-1)) )
            throw MechanicsException(__PRETTY_FUNCTION__, "Ptrans * (d- F*lambda) = 0 not satisfied");


        MPI_Barrier(MPI_COMM_WORLD);



        std::vector<int> recvCount;
        std::vector<int> displs;
        MpiGatherRecvCountAndDispls(recvCount, displs, numRigidBodyModesLocal);
        VectorXd alphaLocal = alphaGlobal.segment(displs[structure->mRank], numRigidBodyModesLocal);

//    structure->GetLogger() << "alphaGlobal: \n" << alphaGlobal                << "\n\n";
//    structure->GetLogger() << "alphaLocal: \n" << alphaLocal                << "\n\n";
//    structure->GetLogger() << "rigidBodyModes: \n" << rigidBodyModes                << "\n\n";

        VectorXd delta_dof_active;
        delta_dof_active.setZero(mStructure->GetNumTotalDofs());
        delta_dof_active.head(mStructure->GetNumTotalActiveDofs())    = mSolver.solve( (residual_mod - Btrans * deltaLambda).head(structure->GetNumTotalActiveDofs()));



        delta_dof_active   = delta_dof_active - rigidBodyModes * alphaLocal;


        StructureOutputBlockVector  delta_dof_dt0(structure->GetDofStatus());
        int offset = 0;
        for (const auto& dof : activeDofSet)
        {
            const int numActiveDofs       = structure->GetNumActiveDofs(dof);

            delta_dof_dt0.J[dof]    = delta_dof_active.segment(offset, numActiveDofs);
            offset += numActiveDofs;


        }

        delta_dof_dt0.K[Node::eDof::DISPLACEMENTS]    = delta_dof_active.segment(offset, structure->GetNumDependentDofs(Node::eDof::DISPLACEMENTS));
        //    structure->GetLogger() << "delta_dof_active: \n" << delta_dof_active                << "\n\n";
        //    structure->GetLogger() << "delta_dof_active.size(): \n" << delta_dof_active.size()  << "\n\n";
        //    structure->GetLogger() << "delta_dof_dt0.J: \n" << delta_dof_dt0.J                  << "\n\n";
        //    structure->GetLogger() << "delta_dof_dt0.K: \n" << delta_dof_dt0.K                  << "\n\n";


        return delta_dof_dt0;


    }


    //! @brief perform the time integration
    //! @param rTimeDelta ... length of the simulation
    void Solve(double rTimeDelta)
    {

        try
        {
            boost::mpi::communicator world;
            StructureFeti* structure = static_cast<StructureFeti*>(mStructure);
            mStructure->NodeBuildGlobalDofs(__PRETTY_FUNCTION__);

            mStructure->GetLogger() << "********************************************" << "\n";
            mStructure->GetLogger() << "**        Calculate rigid body modes      **" << "\n";
            mStructure->GetLogger() << "********************************************" << "\n\n";

            structure->CalculateRigidBodyModesTotalFETI();
//        structure->CalculateRigidBodyModes();

            const int numRigidBodyModes = structure->GetNumRigidBodyModes();



            mStructure->GetLogger() << "************************************************" << "\n";
            mStructure->GetLogger() << "**        Assemble connectivity matrix B      **" << "\n";
            mStructure->GetLogger() << "************************************************" << "\n\n";

            structure->AssembleConnectivityMatrix();
            const SparseMatrix& B       = structure->GetConnectivityMatrix();
            const SparseMatrix Btrans   = B.transpose();

            mStructure->GetLogger() << "******************************************************" << "\n";
            mStructure->GetLogger() << "**        Calculate interface rigid body modes      **" << "\n";
            mStructure->GetLogger() << "******************************************************" << "\n\n";



            structure->CalculateInterfaceRigidBodyModes();

            mStructure->GetLogger() << "********************************************" << "\n";
            mStructure->GetLogger() << "**        Calculate projection            **" << "\n";
            mStructure->GetLogger() << "********************************************" << "\n\n";

            structure->CalculateG();
            structure->CalculateProjectionMatrix();

            structure->CheckProjectionMatrix();
            structure->CheckProjectionOfCoarseGrid();

            const int numActiveDofs     = mStructure->GetNumTotalActiveDofs();
            const int numTotalDofs      = mStructure->GetNumTotalDofs();


            VectorXd deltaLambda            = VectorXd::Zero(B.rows());
            VectorXd lambda                 = VectorXd::Zero(B.rows());
            VectorXd lastConverged_lambda   = VectorXd::Zero(B.rows());

            double curTime  = mTime;
            double timeStep = mTimeStep;
            mStructure->SetPrevTime(curTime);
            mStructure->SetTime(curTime);


            CalculateStaticAndTimeDependentExternalLoad();

            // sets hasInteractingConstraints to true to also assemble K_KJ and K_KK
            // only necessary for CheckRigidBodyModes
            mStructure->DofStatusSetHasInteractingConstraints(true);
            const DofStatus& dofStatus = mStructure->GetDofStatus();

            assert(dofStatus.HasInteractingConstraints()==true);

            if(mStepActiveDofs.empty())
            {
                mStepActiveDofs.push_back(mStructure->DofTypesGetActive());
            }
            else
            {
                for(unsigned int i=0; i<mStepActiveDofs.size();++i)
                {
                    if(mStepActiveDofs[i].empty())
                    {
                        throw MechanicsException(__PRETTY_FUNCTION__, "Calculation step " +std::to_string(i)+ " has no active DOFs.");
                    }
                }
            }


            //********************************************
            //        Allocate Variables
            //********************************************

            StructureOutputBlockMatrix  hessian0(dofStatus, true);

            StructureOutputBlockVector  delta_dof_dt0(dofStatus, true);

            StructureOutputBlockVector  dof_dt0(dofStatus, true); // e.g. disp

            StructureOutputBlockVector  lastConverged_dof_dt0(dofStatus, true); // e.g. disp
            StructureOutputBlockVector  lastConverged_dof_dt1(dofStatus, true); // e.g. velocity
            StructureOutputBlockVector  lastConverged_dof_dt2(dofStatus, true); // e.g. accelerations

            StructureOutputBlockVector  extForce(dofStatus, true);
            StructureOutputBlockVector  intForce(dofStatus, true);
            StructureOutputBlockVector  residual(dofStatus, true);

            StructureOutputBlockVector  prevResidual(dofStatus, true);


            // for constraints
            // ---------------
//        const auto& cmat = mStructure->GetConstraintMatrix();

            BlockScalar normResidual(dofStatus);
            BlockScalar normResidualK(dofStatus);

            //********************************************
            //        Declare and fill Output maps
            //********************************************

            // Declare output maps

            std::map<NuTo::eStructureOutput, NuTo::StructureOutputBase*> evalInternalGradient;
            std::map<NuTo::eStructureOutput, NuTo::StructureOutputBase*> evalInternalGradientAndHessian0;
            std::map<NuTo::eStructureOutput, NuTo::StructureOutputBase*> evalHessian0;

            evalInternalGradient                [eStructureOutput::INTERNAL_GRADIENT]   = &intForce;

            evalInternalGradientAndHessian0     [eStructureOutput::INTERNAL_GRADIENT]   = &intForce;
            evalInternalGradientAndHessian0     [eStructureOutput::HESSIAN0]            = &hessian0;

            evalHessian0                        [eStructureOutput::HESSIAN0]            = &hessian0;

            //********************************************
            //        Declare and fill Input map
            //********************************************

            ConstitutiveInputMap inputMap;
            inputMap[Constitutive::eInput::CALCULATE_STATIC_DATA] = std::make_unique<ConstitutiveCalculateStaticData>(
                    eCalculateStaticData::EULER_BACKWARD);

            ExtractDofValues(lastConverged_dof_dt0, lastConverged_dof_dt1, lastConverged_dof_dt2);

            UpdateAndGetAndMergeConstraintRHS(curTime, lastConverged_dof_dt0);

            PostProcess(residual);





            mStructure->GetLogger() << "********************************************" << "\n";
            mStructure->GetLogger() << "**        Start time stepping             **" << "\n";
            mStructure->GetLogger() << "********************************************" << "\n\n";

            while (curTime < rTimeDelta)
            {


                mStructure->DofTypeActivateAll();

                // calculate Delta_BRhs and Delta_ExtForce
//            auto bRHS           = UpdateAndGetAndMergeConstraintRHS(mTime, lastConverged_dof_dt0);

                curTime += timeStep;
                SetTimeAndTimeStep(curTime, timeStep, rTimeDelta);     //check whether harmonic excitation, check whether curTime is too close to the time data
                mStructure->SetTime(curTime);

//            auto deltaBRHS = UpdateAndGetConstraintRHS(curTime) - bRHS;
                extForce = CalculateCurrentExternalLoad(curTime);

                for (const auto& activeDofSet : mStepActiveDofs)
                {
                    mStructure->DofTypeSetIsActive(activeDofSet);

                    // ******************************************************
                    mStructure->Evaluate(inputMap, evalInternalGradientAndHessian0);
                    // ******************************************************

//                structure->CheckRigidBodyModes(hessian0);
                    structure->CheckStiffnessPartitioning(hessian0);

                    mStructure->GetLogger() << "extForce norm: \t" << extForce.J.Export().norm() + extForce.K.Export().norm() << "\n\n";
                    mStructure->GetLogger() << "intForce norm: \t" << intForce.J.Export().norm() + intForce.K.Export().norm() << "\n\n";
                    VectorXd rhs;
                    residual =  extForce - intForce;
                    residual.J = BlockFullVector<double>(residual.J.Export() - (Btrans * lambda).head(numActiveDofs), dofStatus);


                    rhs.setZero(numTotalDofs);
                    rhs.head(numActiveDofs)     = residual.J.Export();
                    rhs.tail(numRigidBodyModes) = residual.K.Export() - (Btrans * lambda).tail(numRigidBodyModes);


                    mStructure->GetLogger() << "Rhs norm: \t" << rhs.norm() << "\n\n";

                    // K_{JJ}^{-1}
                    mSolver.compute(hessian0.JJ.ExportToEigenSparseMatrix());

                    mLocalPreconditioner.resize(B.rows(), B.rows());
//                mLocalPreconditioner.setIdentity();

                    VectorXd hessianJJdiag = hessian0.JJ.ExportToEigenSparseMatrix().diagonal();
                    VectorXd hessianKKdiag = hessian0.KK.ExportToEigenSparseMatrix().diagonal();
                    SparseMatrix Kdiag(B.cols(),B.cols());

                    mStructure->GetLogger() << "mLocalPreconditioner.rows(): \t" << mLocalPreconditioner.rows() << "\n\n";
                    mStructure->GetLogger() << "hessianJJdiag.rows(): \t" << hessianJJdiag.rows() << "\n\n";
                    mStructure->GetLogger() << "hessianKKdiag.rows(): \t" << hessianKKdiag.rows() << "\n\n";
                    for (int i = 0; i < hessianJJdiag.rows(); ++i)
                        Kdiag.insert(i,i) = hessianJJdiag[i];

                    for (int i = 0; i < hessianKKdiag.rows(); ++i)
                        Kdiag.insert(i + hessianJJdiag.rows(),i + hessianJJdiag.rows()) = hessianKKdiag[i];

                    mStructure->GetLogger() << "mLOca: \t" << Kdiag << "\n\n";

                    MPI_Barrier(MPI_COMM_WORLD);

                    mLocalPreconditioner = B * Kdiag * Btrans;


                    delta_dof_dt0   = FetiSolve(rhs, activeDofSet, deltaLambda, timeStep);


                    dof_dt0 = lastConverged_dof_dt0 + delta_dof_dt0;

                    lambda  = lastConverged_lambda + deltaLambda;

//                mStructure->GetLogger() << "active dofs disp: \n" << mStructure->GetNumActiveDofs(Node::eDof::DISPLACEMENTS) << "\n\n";
//                mStructure->GetLogger() << "depe dofs disp: \n" << mStructure->GetNumDependentDofs(Node::eDof::DISPLACEMENTS) << "\n\n";
//                mStructure->GetLogger() << "active dofs nl eq strain: \n" << mStructure->GetNumActiveDofs(Node::eDof::NONLOCALEQSTRAIN) << "\n\n";
//                mStructure->GetLogger() << "depe dofs nl eq strain: \n" << mStructure->GetNumDependentDofs(Node::eDof::NONLOCALEQSTRAIN) << "\n\n";




                    mStructure->NodeMergeDofValues(dof_dt0);



                    // ******************************************************
                    mStructure->Evaluate(inputMap, evalInternalGradient);
                    // ******************************************************

                    residual =  extForce - intForce;
                    residual.J = BlockFullVector<double>(residual.J.Export() - (Btrans * lambda).head(numActiveDofs), dofStatus);

                    rhs.setZero(numTotalDofs);
                    rhs.head(numActiveDofs)     = residual.J.Export();
                    rhs.tail(numRigidBodyModes) = residual.K.Export() - (Btrans * lambda).tail(numRigidBodyModes);

                    normResidual = residual.J.CalculateInfNorm();
                    for (const auto& dof : activeDofSet)
                        boost::mpi::all_reduce(world,boost::mpi::inplace(normResidual[dof]),std::plus<double>());

                    mStructure->GetLogger() << "Residual J: \t" << normResidual << "\n\n";
                    mStructure->GetLogger() << "Rhs: \t" << rhs.norm() << "\n\n";


                    int iteration = 0;
                    while( not(normResidual < mToleranceResidual) and iteration < mMaxNumIterations)
                    {

                        // ******************************************************
                        mStructure->Evaluate(inputMap, evalHessian0);
                        // ******************************************************

                        mSolver.compute(hessian0.JJ.ExportToEigenSparseMatrix());

                        delta_dof_dt0   = FetiSolve(rhs, activeDofSet, deltaLambda, 0.);

                        dof_dt0 += delta_dof_dt0;
                        lambda += deltaLambda;
                        mStructure->NodeMergeDofValues(dof_dt0);

                        // ******************************************************
                        mStructure->Evaluate(inputMap, evalInternalGradient);
                        // ******************************************************


                        residual =  extForce - intForce;
                        residual.J = BlockFullVector<double>(residual.J.Export() - (Btrans * lambda).head(numActiveDofs), dofStatus);

                        rhs.setZero(numTotalDofs);
                        rhs.head(numActiveDofs)     = residual.J.Export();
                        rhs.tail(numRigidBodyModes) = residual.K.Export() - (Btrans * lambda).tail(numRigidBodyModes);

                        normResidual = residual.J.CalculateInfNorm();
                        for (const auto& dof : activeDofSet)
                            boost::mpi::all_reduce(world,boost::mpi::inplace(normResidual[dof]),std::plus<double>());

                        PrintInfoIteration(normResidual,iteration);
                        mStructure->GetLogger() << "delta_dof_dt0 J: \n" << delta_dof_dt0.J.Export().norm() << "\n\n";
                        mStructure->GetLogger() << "delta_dof_dt0 K: \n" << delta_dof_dt0.K.Export().norm() << "\n\n";

                        iteration++;

                    } // end of while(normResidual<mToleranceForce && iteration<mMaxNumIterations)


                    world.barrier();

                    if (normResidual < mToleranceResidual)
                    {
                        //converged solution
                        mStructure->ElementTotalUpdateStaticData();

                        //store converged step
                        lastConverged_dof_dt0 = dof_dt0;
                        lastConverged_lambda  = lambda;

                        prevResidual = residual;

                        mStructure->NodeMergeDofValues(dof_dt0);

                        //update structure time
                        mStructure->SetPrevTime(curTime);

                        mTime+=timeStep;


                        mStructure->GetLogger() << "Convergence after " << iteration << " iterations at time " << mTime << " (timestep " << timeStep << ").\n";
                        mStructure->GetLogger() << "Residual: \t" << normResidual << "\n\n";


                        //eventually increase next time step
                        if (mAutomaticTimeStepping && iteration<0.25*mMaxNumIterations)
                        {
                            timeStep*=1.5;
                            if (timeStep>mMaxTimeStep)
                                timeStep = mMaxTimeStep;
                        }

                        //perform Postprocessing
                        residual.J = BlockFullVector<double>(residual.J.Export() + (Btrans * lambda).head(numActiveDofs), dofStatus);
                        PostProcess(residual);

                        mStructure->GetLogger() << "normResidual: \n" << normResidual << "\n\n";

                        if (mCallback && mCallback->Exit(*mStructure))
                            return;

                    }
                    else
                    {
                        if (structure->mRank == 0)
                            mStructure->GetLogger() << "No convergence with timestep " << timeStep << "\n\n";

                        //no convergence
                        if (mAutomaticTimeStepping)
                        {
                            //no convergence, reduce the time step and start from scratch
                            curTime -= timeStep;
                            timeStep *= 0.5;
                            if (timeStep < mMinTimeStep) {
                                mStructure->GetLogger() << "The minimal time step achieved, the actual time step is " << timeStep << "\n";
                                throw MechanicsException(__PRETTY_FUNCTION__, "No convergence, the current time step is too short.");
                            }
                        }
                        else
                        {
                            throw MechanicsException(__PRETTY_FUNCTION__, "No convergence with the current maximum number of iterations, either use automatic time stepping, reduce the time step or the minimal line search cut back factor.");
                        }

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


    double CalculateLoadFactor(double rTime)
    {

        //extract the two data points
        const double loadFactor1  = mTimeDependentLoadFactor(0,1);
        const double loadFactor2  = mTimeDependentLoadFactor(1,1);
        const double time1        = mTimeDependentLoadFactor(0,0);
        const double time2        = mTimeDependentLoadFactor(1,0);

        return loadFactor1 + (loadFactor2-loadFactor1)/(time2-time1) * (rTime-time1);

    }
private:

    FetiSolver mFetiSolver;
    EigenSolver mSolver;
    SparseMatrix mLocalPreconditioner;
    SparseMatrix mTangentStiffnessMatrix;
    const double    mCpgTolerance     = 1.0e-6;
    const int       mCpgMaxIterations = 1000;

};
}// namespace NuTo
