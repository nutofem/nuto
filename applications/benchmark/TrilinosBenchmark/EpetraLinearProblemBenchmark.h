#include "../LinearElasticBenchmarkStructure.h"

#include <mpi.h>
#include <Epetra_Map.h>
#include <Epetra_LinearProblem.h>
//#include <Epetra_RowMatrix.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_MultiVector.h>
#include <Epetra_Vector.h>
#include <Epetra_MpiComm.h>
#include <Amesos.h>
#include <Amesos_BaseSolver.h>
#include <Amesos_Mumps.h>
#include <Amesos_ConfigDefs.h>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <AztecOO.h>
#include <BelosLinearProblem.hpp>
#include <BelosBlockGmresSolMgr.hpp>
#include <BelosPseudoBlockGmresSolMgr.hpp>
#include <BelosEpetraAdapter.hpp>
#include <BelosPseudoBlockCGSolMgr.hpp>
#include <Ifpack.h>
#include <Ifpack_AdditiveSchwarz.h>
#include <eigen3/Eigen/Core>
#include <EpetraExt_RowMatrixOut.h>

#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/structures/StructureOutputBlockMatrix.h"
#include "mechanics/structures/StructureOutputBlockVector.h"

#include "../../applications/custom/TestClasses/ConversionTools.h"

namespace NuTo
{

namespace Benchmark
{


class EpetraLinearProblemBenchmark
{
public:
//    EpetraLinearProblemBenchmark(Epetra_RowMatrix& rA, Epetra_MultiVector& rLhs, Epetra_MultiVector& rRhs) : mProb()
    EpetraLinearProblemBenchmark(LinearElasticBenchmarkStructure& rS, Epetra_MpiComm rComm, bool rGlobalStructure = true) : mProb()
    {
        ConversionTools converter(rComm);
        mA = new Epetra_CrsMatrix(converter.convertEigen2EpetraCrsMatrix(rS.GetStructure().BuildGlobalHessian0().JJ.ExportToEigenSparseMatrix(), rGlobalStructure, true));
        Eigen::VectorXd grad = rS.GetStructure().BuildGlobalInternalGradient().J.Export();

//        int n = grad.rows();
//        int j = 0;
//        for (int i = 0; i < n; ++i)
//        {
//            j = std::rand() % n;
//            grad[j] = pow(-1, j) * 2 * i;
//        }

        mRhs = new Epetra_Vector(converter.convertEigen2EpetraVector(grad, rGlobalStructure));
        mLhs = new Epetra_Vector(*mRhs);
//        mLhs = new Epetra_Vector(mA->ColMap());
        mLhs->Random();

        mProb.SetOperator(mA);
        mProb.SetRHS(mRhs);
        mProb.SetLHS(mLhs);
    }


    void solveDirect(std::string rSolverType, bool printStatus = false, bool printTiming = false)
    {
        Amesos Factory;
        bool solverAvail = Factory.Query(rSolverType);
        if (solverAvail)
        {
            Teuchos::ParameterList params;
            params.set("PrintStatus", printStatus);
            params.set("PrintTiming", printTiming);
            params.set("MaxProcs", -3); //all processes in communicator will be used
            if (rSolverType == "Mumps")
            {
                Amesos_Mumps solverMumps(mProb);
                Teuchos::ParameterList mumpsList = params.sublist("mumps");
                int* icntl = new int[40];
                icntl[14] = 50;
                mumpsList.set("ICNTL", icntl);
                EpetraExt::RowMatrixToMatlabFile("mA.mat", *mA);
                solverMumps.SetParameters(params);
                solverMumps.Solve();
            }
            else
            {
                Amesos_BaseSolver* solver;
                solver = Factory.Create(rSolverType, mProb);
                solver->SetParameters(params);
                solver->Solve();
                delete solver;
            }

        }
        else
        {
            throw NuTo::Exception("Solver '" + rSolverType + "' not available");
        }
    }


    void solveIterative_AztecOO(std::string rSolverType, bool printStatus = false)
    {
        int method = str2AZSolverType(rSolverType);
        bool solverAvail = (method > -1);
        if (solverAvail)
        {
            /*METHODS
             * AZ_cg               0 --> preconditioned conjugate gradient method
             * AZ_gmres            1 --> preconditioned gmres method
             * AZ_cgs              2 --> preconditioned cg squared method
             * AZ_tfqmr            3 --> preconditioned transpose-free qmr method
             * AZ_bicgstab         4 --> preconditioned stabilized bi-cg method
             * AZ_fixed_pt         8 --> fixed point iteration
             * AZ_cg_condnum      11
             * AZ_gmres_condnum   12
             */

            int output = (printStatus ? AZ_all : AZ_none);

            AztecOO Solver(mProb);
            Solver.SetAztecOption(AZ_solver, method);
    //        Solver.SetAztecOption(AZ_output, output);
            Solver.SetAztecOption(AZ_diagnostics, output);
            Solver.SetAztecOption(AZ_precond, AZ_dom_decomp);
            Solver.SetAztecOption(AZ_subdomain_solve, AZ_ilut);
            Solver.SetAztecOption(AZ_overlap, 1);
            Solver.SetAztecOption(AZ_orthog, AZ_classic);
            Solver.SetAztecOption(AZ_kspace, 50);
            Solver.Iterate(1000, 1e-8);
        }
        else
        {
            throw NuTo::Exception("Solver '" + rSolverType + "' not available");
        }
    }


    void solveIterative_Belos(std::string rSolverType, bool usePreconditioner = false)
    {
        int method = str2BelosSolverType(rSolverType);
        bool solverAvail = (method > -1);
        if (solverAvail)
        {
            using Teuchos::RCP;
            using Teuchos::rcp;


            RCP<Epetra_CrsMatrix> A = rcp(new Epetra_CrsMatrix(*mA));
            RCP<Epetra_MultiVector> LHS = rcp (new Epetra_MultiVector (*mLhs));
            RCP<Epetra_MultiVector> RHS = rcp (new Epetra_MultiVector (*mRhs));

            RCP<Belos::LinearProblem<double, Epetra_MultiVector, Epetra_Operator> > belosProblem = rcp (new Belos::LinearProblem<double, Epetra_MultiVector, Epetra_Operator>(A, LHS, RHS));

            if (usePreconditioner)
            {
                Teuchos::ParameterList paramList;
                Ifpack Factory;

                std::string PrecType = "ILU"; // incomplete LU
                int OverlapLevel = 1;

                RCP<Ifpack_Preconditioner> prec = rcp(Factory.Create(PrecType, &*A, OverlapLevel));
                TEUCHOS_TEST_FOR_EXCEPTION(prec == Teuchos::null, std::runtime_error,
                         "IFPACK failed to create a preconditioner of type \""
                         << PrecType << "\" with overlap level "
                         << OverlapLevel << ".");

                paramList.set("fact: drop tolerance", 1e-9);
                paramList.set("fact: level-of-fill", 1);

                paramList.set("schwarz: combine mode", "Add");
                prec->SetParameters(paramList);

    //            IFPACK_CHK_ERR(prec->Initialize());
                prec->Initialize();

                prec->Compute();

                RCP<Belos::EpetraPrecOp> belosPrec = rcp (new Belos::EpetraPrecOp(prec));
                belosProblem->setRightPrec (belosPrec);
            }


            bool set = belosProblem->setProblem();
            TEUCHOS_TEST_FOR_EXCEPTION( ! set,
                      std::runtime_error,
                      "*** Belos::LinearProblem failed to set up correctly! ***");

            // Create a parameter list to define the Belos solver.
            RCP<Teuchos::ParameterList> belosList = rcp (new Teuchos::ParameterList ());
            belosList->set ("Block Size", 1);              // Blocksize to be used by iterative solver
            belosList->set ("Num Blocks", 30);              //Krylov dimension
            belosList->set ("Maximum Restarts", 20);
            belosList->set ("Maximum Iterations", 1000);   // Maximum number of iterations allowed
            belosList->set ("Convergence Tolerance", 1e-8);// Relative convergence tolerance requested
            belosList->set ("Verbosity", Belos::Errors+Belos::Warnings+Belos::TimingDetails+Belos::FinalSummary );

            // Create an iterative solver manager.
            if (method == 0)
            {
                Belos::PseudoBlockCGSolMgr<double, Epetra_MultiVector, Epetra_Operator> belosCGSolver(belosProblem, belosList);
                Belos::ReturnType retCG = belosCGSolver.solve();
            }
            else if (method == 1)
            {
                Belos::PseudoBlockGmresSolMgr<double, Epetra_MultiVector, Epetra_Operator> belosSolver(belosProblem, belosList);
                Belos::ReturnType ret = belosSolver.solve();
            }

//            RCP<const Epetra_MultiVector> solu = belosSolver.getProblem().getLHS();
//            solu->Print(std::cout);

        }
        else
        {
            throw NuTo::Exception("Solver '" + rSolverType + "' not available");
        }
    }


    Epetra_LinearProblem getProblemEquation()
    {
        return mProb;
    }

    Epetra_CrsMatrix getSystemMatrix()
    {
        return *mA;
    }

    Epetra_Vector getRhs()
    {
        return *mRhs;
    }

    Epetra_Vector getLhs()
    {
        return *mLhs;
    }



private:
    Epetra_LinearProblem mProb;
    Epetra_CrsMatrix* mA;
    Epetra_Vector* mRhs;
    Epetra_Vector* mLhs;


    int str2AZSolverType(std::string rSolverName)
    {

        if (rSolverName == "CG")
            return AZ_cg;

        if (rSolverName == "CGS")
            return AZ_cgs;

        if (rSolverName == "CG_CondNum")
            return AZ_cg_condnum;

        if (rSolverName == "GMRES")
            return AZ_gmres;

        if (rSolverName == "GMRES_CondNum")
            return AZ_gmres_condnum;

        if (rSolverName == "BiCGStab")
            return AZ_bicgstab;

        return -1;
    }


    int str2BelosSolverType(std::string rSolverName)
    {
        if (rSolverName == "CG")
            return 0;

        if (rSolverName == "GMRES")
            return 1;

        return -1;
    }

};



} // Benchmark

} // NuTo
