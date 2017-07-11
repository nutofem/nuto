#include <iostream>
#include <mpi.h>

#include <eigen3/Eigen/Core>
#include <boost/filesystem.hpp>
#include <Epetra_DataAccess.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_MultiVector.h>
#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include <Epetra_MpiComm.h>
#include <Epetra_LinearProblem.h>
#include <Epetra_Export.h>
#include <AztecOO.h>
#include <Amesos.h>
#include <Amesos_BaseSolver.h>
#include <Amesos_Mumps.h>
#include <Amesos_ConfigDefs.h>
#include <Teuchos_GlobalMPISession.hpp>
#include <EpetraExt_VectorOut.h>
#include <EpetraExt_MultiVectorOut.h>
#include <EpetraExt_RowMatrixOut.h>
#include <BelosLinearProblem.hpp>
#include <BelosBlockGmresSolMgr.hpp>
#include <BelosPseudoBlockGmresSolMgr.hpp>
#include <BelosEpetraAdapter.hpp>
#include <Ifpack.h>
#include <Ifpack_AdditiveSchwarz.h>
//#include <Ifpack2_AdditiveSchwarz.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/structures/StructureBase.h"
#include "mechanics/sections/SectionTruss.h"
#include "mechanics/sections/SectionPlane.h"
#include "mechanics/nodes/NodeDof.h"
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/MechanicsEnums.h"
#include "mechanics/structures/StructureOutputBlockMatrix.h"
#include "mechanics/structures/StructureOutputBlockVector.h"
#include "visualize/VisualizeEnum.h"
#include "mechanics/timeIntegration/NewmarkDirect.h"
#include "mechanics/dofSubMatrixStorage/BlockFullMatrix.h"


#include "ConversionTools.h"
#include "PrintTools.h"
#include "TrilinosUtils.h"


//#define SHOW_INTERMEDIATE_RESULTS
//#define SHOW_SOLUTION


void checkDirectories(const std::string rLogDirectory, const std::string rResultDirectory)
{
    if (!boost::filesystem::exists(rLogDirectory))
        boost::filesystem::create_directory(rLogDirectory);

    if (!boost::filesystem::exists(rResultDirectory))
        boost::filesystem::create_directory(rResultDirectory);
}


Epetra_MultiVector solveSystem(Epetra_CrsMatrix rA, Epetra_MultiVector rLhs, Epetra_MultiVector rRhs, bool iterative = true, bool useAztecOO = true)
{
    Epetra_LinearProblem problem(&rA, &rLhs, &rRhs);

    if (iterative)
    {
        if (useAztecOO)
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
            int method = AZ_gmres;

            AztecOO Solver(problem);
            Solver.SetAztecOption(AZ_solver, method);
    //        Solver.SetAztecOption(AZ_output, AZ_all);
            Solver.SetAztecOption(AZ_diagnostics, AZ_all);
            Solver.SetAztecOption(AZ_precond, AZ_dom_decomp);
    //        Solver.SetAztecOption(AZ_precond, AZ_Jacobi);
            Solver.SetAztecOption(AZ_subdomain_solve, AZ_ilut);
            Solver.SetAztecOption(AZ_overlap, 1);
            Solver.SetAztecOption(AZ_orthog, AZ_classic);
            Solver.SetAztecOption(AZ_kspace, 50);
            Solver.Iterate(1000,1e-8);
        }
        else
        {
            //++++++++++++Ifpack preconditioner+++++++++++++
            Teuchos::ParameterList paramList;
            Ifpack Factory;

            // Create the preconditioner. For the list of PrecType values check the IFPACK documentation.
            string PrecType = "ILU"; // incomplete LU
            int OverlapLevel = 1; // must be >= 0. If Comm.NumProc() == 1,
                                // it is ignored.

            RCP<Epetra_CrsMatrix> A = rcp (new Epetra_CrsMatrix(rA));
            RCP<Ifpack_Preconditioner> prec = rcp(Factory.Create(PrecType, &*A, OverlapLevel));
            TEUCHOS_TEST_FOR_EXCEPTION(prec == null, std::runtime_error,
                     "IFPACK failed to create a preconditioner of type \""
                     << PrecType << "\" with overlap level "
                     << OverlapLevel << ".");

            // Specify parameters for ILU.  ILU is local to each MPI process.
            paramList.set("fact: drop tolerance", 1e-9);
            paramList.set("fact: level-of-fill", 1);

            // how to combine overlapping results:
            // "Add", "Zero", "Insert", "InsertAdd", "Average", "AbsMax"
            paramList.set("schwarz: combine mode", "Add");
//            IFPACK_CHK_ERR(prec->SetParameters(paramList));
            prec->SetParameters(paramList);

            // Initialize the preconditioner. At this point the matrix must have
            // been FillComplete()'d, but actual values are ignored.
//            IFPACK_CHK_ERR(prec->Initialize());
            prec->Initialize();

            // Build the preconditioner, by looking at the values of the matrix.
//            IFPACK_CHK_ERR(prec->Compute());
            prec->Compute();

            // Create the Belos preconditioned operator from the Ifpack preconditioner.
            // NOTE:  This is necessary because Belos expects an operator to apply the
            //        preconditioner with Apply() NOT ApplyInverse().
            RCP<Belos::EpetraPrecOp> belosPrec = rcp (new Belos::EpetraPrecOp(prec));

            //++++++++++++end preconditioner definition+++++++++++++

            //++++++++++++perform solve+++++++++++++++++++++

            RCP<Epetra_MultiVector> LHS = rcp (new Epetra_MultiVector (rLhs));
            RCP<Epetra_MultiVector> RHS = rcp (new Epetra_MultiVector (rRhs));

            // Need a Belos::LinearProblem to define a Belos solver
            typedef Epetra_MultiVector                MV;
            typedef Epetra_Operator                   OP;
            RCP<Belos::LinearProblem<double,MV,OP> > belosProblem
            = rcp (new Belos::LinearProblem<double,MV,OP>(A, LHS, RHS));

            belosProblem->setRightPrec (belosPrec);

            bool set = belosProblem->setProblem();
            TEUCHOS_TEST_FOR_EXCEPTION( ! set,
                      std::runtime_error,
                      "*** Belos::LinearProblem failed to set up correctly! ***");

            // Create a parameter list to define the Belos solver.
            RCP<ParameterList> belosList = rcp (new ParameterList ());
            belosList->set ("Block Size", 1);              // Blocksize to be used by iterative solver
            belosList->set ("Num Blocks", 30);              //Krylov dimension
            belosList->set ("Maximum Restarts", 20);
            belosList->set ("Maximum Iterations", 1000);   // Maximum number of iterations allowed
            belosList->set ("Convergence Tolerance", 1e-8);// Relative convergence tolerance requested
            belosList->set ("Verbosity", Belos::Errors+Belos::Warnings+Belos::TimingDetails+Belos::FinalSummary );

            // Create an iterative solver manager.
            Belos::PseudoBlockGmresSolMgr<double, MV, OP> belosSolver(belosProblem, belosList);

            // Perform solve.
            Belos::ReturnType ret = belosSolver.solve();
            return *belosSolver.getProblem().getLHS().get();

            //++++++++++++end of solve++++++++++++++++++++

        }
    }
    else
    {
        /*METHODS
         * Klu
         * Mumps
         * Lapack
         * Scalapack
         * Umfpack
         * Superlu
         * Superludist
         * Dscpack
         * Taucs
         */
        Amesos Factory;
        std::string solverType = "Mumps";
        bool solverAvail = Factory.Query(solverType);
        if (rA.Comm().MyPID() == 0)
            std::cout << "Direct Solver '" << solverType << "' available: " << (solverAvail ? "oh yeah" : "oh no") << std::endl;
        if (solverAvail)
        {
            Teuchos::ParameterList params;
            params.set("PrintStatus", true);
            params.set("PrintTiming", true);
            params.set("MaxProcs", -3); //all processes in communicator will be used
//            Amesos_BaseSolver* solver;
//            solver = Factory.Create(solverType, problem);
            Amesos_Mumps* solver;
            solver = new Amesos_Mumps(problem);
            Teuchos::ParameterList mumpsList = params.sublist("mumps");
            int* icntl = new int[40];
            double* cntl = new double[5];
            icntl[0] = 0;
            icntl[1] = 0;
            icntl[2] = 0;
            icntl[3] = 0;
            icntl[5] = 7;
            icntl[6] = 7;
            icntl[7] = 0;
            icntl[28] = 2;
            icntl[29] = 0;
            mumpsList.set("ICNTL", icntl);
            solver->SetParameters(params);
            solver->Solve();
//            for (int i = 0; i < 40; ++i)
//            {
//                std::cout << solver->GetINFO()[i] << std::endl;
//            }
            delete[] icntl;
            delete[] cntl;
            delete solver;
        }
        else
        {
            Epetra_MultiVector zeroVec(*(problem.GetLHS()));
            zeroVec.PutScalar(0.);
            problem.SetLHS(&zeroVec);
        }
    }

    Epetra_MultiVector* lhs = problem.GetLHS();
    return *lhs;
}


void visualizeResults(NuTo::Structure* rS, const std::string& rResultDir, double rTime, int rTimeStep)
{
    //plot the solution vtk file
    std::stringstream ssTimeStepVTK;
    ssTimeStepVTK << rTimeStep;
    boost::filesystem::path resultFile(rResultDir);

//    if (mExportDataFileNodes==true)
//    {
//        resultFile /= std::string("Nodes") + ssTimeStepVTK.str() + std::string(".vtu");
//        rS->ExportVtkDataFileNodes(resultFile.string(), true);
//    }

    std::stringstream timeFormatted;
    timeFormatted.width(15);
    timeFormatted.precision(12);
    timeFormatted << rTime;

    //plot all groups separately
    for (auto const & iVisualizePair : rS->GetGroupVisualizeComponentsMap())
    {
        //plot all elements
        resultFile = rResultDir;
        resultFile /= std::string("Group") + std::to_string(iVisualizePair.first) + std::string("_Elements") + ssTimeStepVTK.str() + std::string(".vtu");
        rS->ElementGroupExportVtkDataFile(iVisualizePair.first, resultFile.string(), true);

        //write an additional pvd file
        resultFile = rResultDir;
        resultFile /= std::string("Group") + std::to_string(iVisualizePair.first) + std::string("_ElementsAll") + std::string(".pvd");

        std::fstream file;
        if (rTimeStep == 0)
        {
            file.open(resultFile.string(), std::fstream::out);
        } else
        {
            file.open(resultFile.string(), std::fstream::out | std::fstream::in | std::ios_base::ate);
        }
        if (!file.is_open())
        {
            throw NuTo::MechanicsException(std::string("[NuTo::TimeIntegrationBase::ExportVisualizationFiles] Error opening file ") + resultFile.string());
        }
        std::stringstream endOfXML;
        endOfXML << "</Collection>" << std::endl;
        endOfXML << "</VTKFile>" << std::endl;
        if (rTimeStep == 0)
        {
            // header /////////////////////////////////////////////////////////////////
            file << "<?xml version=\"1.0\"?>" << std::endl;
            file << "<VTKFile type=\"Collection\">" << std::endl;
            file << "<Collection>" << std::endl;
        } else
        {
            //delete the last part of the xml file
            file.seekp(-endOfXML.str().length(), std::ios_base::end);
        }
        file << "<DataSet timestep=\"" << timeFormatted.str() << "\" file=\"Group" << iVisualizePair.first << "_Elements" << rTimeStep << ".vtu\"/>" << std::endl;
        file << endOfXML.str();
        file.close();
    }
}

void run_Test_2D_GivenMesh()
{
    Epetra_MpiComm Comm(MPI_COMM_WORLD);
    int rank = Comm.MyPID();
    int numProc = Comm.NumProc();


    //CREATE DIRECTORIES FOR LOGS AND RESULTS
    const std::string logDirectory = "2D_Test_GivenMesh_Logs";
    const std::string resultDirectory = "2D_Test_GivenMesh_Results";
    checkDirectories(logDirectory, resultDirectory);

    // one-dimensional problem
    const int dim = 2;

    // --------------------
    // | DOMAIN VARIABLES |
    // --------------------
    int numSubDomains_X = (numProc == 1 ? 1 : 2);
    int numSubDomains_Y = 1;
    int numSubDomains = numSubDomains_X * numSubDomains_Y;
    int* numElementsPerSubDomain_X = new int[numSubDomains_X];
    int* numElementsPerSubDomain_Y = new int[numSubDomains_Y];
    int* numActiveDofsPerSubDomain = new int[numSubDomains];
    int* numActiveDofsTillSubDomain_inclusive = new int[numSubDomains];
    int numElementsTotal_X = 0;
    int numElementsTotal_Y = 0;
    bool setTotalLength = true;
    double* lengthPerSubDomain_X = new double[numSubDomains_X];
    double lengthTotal_X = (setTotalLength ? 8. : 0.);
    double* lengthPerSubDomain_Y = new double[numSubDomains_Y];
    double lengthTotal_Y = (setTotalLength ? 4. : 0.);

    for (int i = 0; i < numSubDomains_X; ++i)
    {
        numElementsPerSubDomain_X[i] = (numProc == 1 ? 8 : 4);
        numElementsTotal_X += numElementsPerSubDomain_X[i];

        if (setTotalLength)
        {
            lengthPerSubDomain_X[i] = lengthTotal_X/numSubDomains_X;
        }
        else
        {
            lengthPerSubDomain_X[i] = 4;
            lengthTotal_X += lengthPerSubDomain_X[i];
        }
    }

    for (int i = 0; i < numSubDomains_Y; ++i)
    {
        numElementsPerSubDomain_Y[i] = 4;
        numElementsTotal_Y += numElementsPerSubDomain_Y[i];

        if (setTotalLength)
        {
            lengthPerSubDomain_Y[i] = lengthTotal_Y/numSubDomains_Y;
        }
        else
        {
            lengthPerSubDomain_Y[i] = 4;
            lengthTotal_Y += lengthPerSubDomain_Y[i];
        }
    }

    // ------------------------
    // | MECHANICAL VARIABLES |
    // ------------------------
    double thickness = 10.;
    double YoungsModulus = 2.e4;
    double PoissonRatio = 0.3;
    double force = 1.e5;

    double fixedDisplacement = 0.;
    Eigen::VectorXd directionX(dim);
    directionX(0) = 1;
    directionX(1) = 0;
    Eigen::VectorXd directionY(dim);
    directionY(0) = 0;
    directionY(1) = 1;
    bool enableDisplacementControl = false;


    //CONSTRUCT STRUCTURE BY IMPORT
    std::map<int, int> newNodes;
    std::map<int, int> gmshNodes;
    NuTo::Structure structure(dim);
    std::vector<std::pair<int, int>> importContainer;
    if (numProc == 1)
    {
        importContainer = structure.ImportFromGmsh("trilinos.msh", true, newNodes, gmshNodes);
    }
    else
    {
        importContainer = structure.ImportFromGmsh("trilinos.msh_00000" + std::to_string(rank+1), true, newNodes, gmshNodes);
    }



    structure.SetVerboseLevel(10);
    structure.GetLogger().OpenFile(logDirectory + "/output_" + std::to_string(rank));
    structure.GetLogger().SetQuiet(false);

    // assign interpolation type
    const int InterpolationType = importContainer[0].second;
    structure.InterpolationTypeAdd(InterpolationType, NuTo::Node::eDof::DISPLACEMENTS,
            NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    structure.ElementTotalConvertToInterpolationType();


    auto section = NuTo::SectionPlane::Create(thickness, true);
    structure.ElementTotalSetSection(section);

    auto material = structure.ConstitutiveLawCreate("Linear_Elastic_Engineering_Stress");
    structure.ConstitutiveLawSetParameterDouble(material, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, YoungsModulus);
    structure.ConstitutiveLawSetParameterDouble(material, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, PoissonRatio);
    structure.ElementTotalSetConstitutiveLaw(material);



    // ----------------------------------------
    // | DEFINE BOUNDARY AND OVERLAPPING AREAS|
    // ----------------------------------------
    int nodesBCLeft = structure.GroupCreate(NuTo::eGroupId::Nodes);
    int nodesBCRight = structure.GroupCreate(NuTo::eGroupId::Nodes);
    structure.GroupAddNodeCoordinateRange(nodesBCLeft, 0, -1e-6, 1e-6);
//    structure.GroupAddNodeCoordinateRange(nodesBCRight, 0, lengthTotal_X-1e-6, lengthTotal_X+1e-6);
    Eigen::Vector2d nodeCoordinates;
    nodeCoordinates(0) = lengthTotal_X;
    nodeCoordinates(1) = 0.;
    structure.GroupAddNodeRadiusRange(nodesBCRight, nodeCoordinates, 0, 1.e-6);

    //COMPUTE OVERLAPPING AREAS
    int numOverlapping = numSubDomains_X - 1;
    int* nodeGroupsOverlap = new int[numOverlapping];
    for (int i = 0; i < numOverlapping; ++i)
        nodeGroupsOverlap[i] = structure.GroupCreate(NuTo::eGroupId::Nodes);

    double currLength = lengthPerSubDomain_X[0];
    for (int i = 0; i < numOverlapping; ++i)
    {
        structure.GroupAddNodeCoordinateRange(nodeGroupsOverlap[i], 0, currLength - 1e-6, currLength + 1e-6);
//        currLength += lengthPerSubDomain_X[i+1];
    }


    //FIX LEFT BOUNDARY
    structure.ConstraintLinearSetDisplacementNodeGroup(nodesBCLeft, directionX, fixedDisplacement);
    structure.ConstraintLinearSetDisplacementNodeGroup(nodesBCLeft, directionY, fixedDisplacement);

    //SET LOAD ON RIGHT BOUNDARY
    structure.SetNumLoadCases(1);
    if (enableDisplacementControl)
    {
//        structure.ConstraintLinearSetDisplacementNode(numElements, directionX, displacement);
    }
    else
    {
//        structure.LoadCreateNodeForce(0, numElements, directionX, force);
        structure.LoadCreateNodeGroupForce(0, nodesBCRight, directionX, force);
    }

    //COMPUTE (LOCAL) HESSIAN AND RESIDUAL
    NuTo::StructureOutputBlockMatrix hessian0 = structure.BuildGlobalHessian0();
    NuTo::StructureOutputBlockVector dofs(structure.GetDofStatus(), true);
    dofs.J.SetZero();
    dofs.K.SetZero();

    NuTo::StructureOutputBlockVector residual = hessian0*dofs - structure.BuildGlobalExternalLoadVector(0) + structure.BuildGlobalInternalGradient();

    Eigen::SparseMatrix<double> hessian0_eigen = hessian0.JJ.ExportToEigenSparseMatrix();
    Eigen::MatrixXd residual_eigen = residual.J.Export();


    int ownedNumActiveDOFs = -1;
    ownedNumActiveDOFs = numElementsPerSubDomain_X[rank]*(numElementsPerSubDomain_Y[0] +1)*dim;

    int localNumActiveDOFs = -1;
    if (rank == 0)
    {
        localNumActiveDOFs = numElementsPerSubDomain_X[rank]*(numElementsPerSubDomain_Y[0] + 1)*dim;
    }
    else
    {
        localNumActiveDOFs = (numElementsPerSubDomain_X[rank] + 1)*(numElementsPerSubDomain_Y[0] + 1)*dim;
    }

    //GET IDS OF OVERLAPPING NODES AND LEFT (FIXED) BOUNDARY
    bool leftNodeIDFound = false;
    bool overlappingNodeIDFound = false;
    std::vector<int> overlapNodeIDs;
    if (numProc > 1)
    {
        if (rank == 0)
            overlapNodeIDs = structure.GroupGetMemberIds(nodeGroupsOverlap[rank]);
        else
            overlapNodeIDs = structure.GroupGetMemberIds(nodeGroupsOverlap[rank-1]);
    }

    int subDomainCounter = 0;

    for (int k = 0; k < numSubDomains_Y; ++k)
    {
        numActiveDofsPerSubDomain[subDomainCounter] = (numElementsPerSubDomain_X[0]) * (numElementsPerSubDomain_Y[k] + 1)*dim;
        numActiveDofsTillSubDomain_inclusive[subDomainCounter] = numActiveDofsPerSubDomain[subDomainCounter];
        ++subDomainCounter;
    }


    for (int j = 1; j < numSubDomains_X; ++j)
    {
        for (int k = 0; k < numSubDomains_Y; ++k)
        {
            numActiveDofsPerSubDomain[subDomainCounter] = (numElementsPerSubDomain_X[j] + 1) * (numElementsPerSubDomain_Y[k] + 1)*dim;
            numActiveDofsTillSubDomain_inclusive[subDomainCounter] = numActiveDofsTillSubDomain_inclusive[subDomainCounter-1] + numActiveDofsPerSubDomain[subDomainCounter];
            ++subDomainCounter;
        }
    }

    std::vector<int> leftNodeIDs = structure.GroupGetMemberIds(nodesBCLeft);
    int leftNodeCount = leftNodeIDs.size();
    int overlappingNodeCount = overlapNodeIDs.size();

    int internDofCounter = -1;
    if (rank == 0)
        internDofCounter = 0;
    else
        internDofCounter = numActiveDofsTillSubDomain_inclusive[rank-1];

    //DEFINE LOCAL-TO-GLOBAL DOF MAPPING
    int* local2GlobalMapping = new int[localNumActiveDOFs];
    for (int i = 0; i < structure.GetNumNodes(); ++i)
    {
        leftNodeIDFound = false;
        overlappingNodeIDFound = false;

        for (int j = 0; j < leftNodeCount; ++j)
        {
            if (leftNodeIDs[j] == i)
            {
                leftNodeIDFound = true;
            }
        }

        for (int j = 0; j < overlappingNodeCount; ++j)
        {
            if (overlapNodeIDs[j] == i)
            {
                overlappingNodeIDFound = true;
            }
        }

        if (!leftNodeIDFound && !overlappingNodeIDFound)
        {
            std::vector<int> dofIDs = structure.NodeGetDofIds(i, NuTo::Node::eDof::DISPLACEMENTS);
            for (int dofID : dofIDs)
            {
                local2GlobalMapping[dofID] = internDofCounter;
                ++internDofCounter;
            }
        }
    }

    internDofCounter = 0;
    if (numProc > 1)
    {
        for (int nodeID : overlapNodeIDs)
        {
            std::vector<int> overlapDofIDs = structure.NodeGetDofIds(nodeID, NuTo::Node::eDof::DISPLACEMENTS);
            for (int dofID : overlapDofIDs)
            {
                if (rank == 0)
                {
                    local2GlobalMapping[dofID] = numActiveDofsTillSubDomain_inclusive[0] - overlappingNodeCount*dim + internDofCounter;
                }
                else
                {
                    local2GlobalMapping[dofID] = numActiveDofsTillSubDomain_inclusive[rank-1] - overlappingNodeCount*dim + internDofCounter;
                }

                ++internDofCounter;
            }
        }
    }

    //DEFINE MAPS FOR LOCAL OWNED DOFS AND OVERLAPPING DOFS
    Epetra_Map owningMap(-1, ownedNumActiveDOFs, 0, Comm);
    Epetra_Map overlappingMap(-1, localNumActiveDOFs, local2GlobalMapping, 0, Comm);
    Epetra_Export exporter(overlappingMap, owningMap);

    //CREATE CORRESPONDING GRAPHS (FOR MATRIX STRUCTURE)
    Epetra_CrsGraph overlappingGraph(Epetra_DataAccess::Copy, overlappingMap, 0);
    Epetra_CrsGraph owningGraph(Epetra_DataAccess::Copy, owningMap, 0);

    for (int i = 0; i < localNumActiveDOFs; ++i)
    {
        for (int j = 0; j < localNumActiveDOFs; ++j)
        {
            overlappingGraph.InsertGlobalIndices(local2GlobalMapping[i], 1, &local2GlobalMapping[j]);
        }
    }

    overlappingGraph.FillComplete();
    owningGraph.Export(overlappingGraph, exporter, Insert); //inter-process communication of matrix structure
    owningGraph.FillComplete();

    //CREATE (GLOBAL) MATRIX AND VECTOR WITH KNOWN GRAPH STRUCTURE
    Epetra_CrsMatrix globalMatrix(Epetra_DataAccess::Copy, owningGraph);
    globalMatrix.FillComplete();
    Epetra_Vector globalRhsVector(owningMap);
    globalRhsVector.PutScalar(0.0);

    //CONVERT EIGEN::MATRIX AND EIGEN::VECTOR TO EPETRA FORMAT
    ConversionTools converter;
    Epetra_CrsMatrix localMatrix = converter.convertEigen2EpetraCrsMatrix(hessian0_eigen, overlappingGraph);
    Epetra_Vector localRhsVector = converter.convertEigen2EpetraVector(residual_eigen, overlappingMap);

    globalMatrix.PutScalar(0.0);
    globalMatrix.Export(localMatrix, exporter, Add);    //inter-process communication of matrix entries
    globalRhsVector.Export(localRhsVector, exporter, Add);    //inter-process communication of vector entries

#ifdef SHOW_INTERMEDIATE_RESULTS
    //PRINT SOME INTERMEDIATE RESULTS
    std::ostream& os = std::cout;

    os << "Hessian0 by NuTo:\n------------------\n" << hessian0.JJ.ExportToFullMatrix() << "\n\n";
    os << "residual by NuTo:\n------------------\n" << residual_eigen << "\n\n";
    os << "Owning Map:\n-------------\n";
    owningMap.Print(os);
    os << "Overlap Map:\n-------------\n";
    overlappingMap.Print(os);
    os << "Local matrix on proc " << rank << ":\n-----------------\n";
    localMatrix.Print(std::cout);
    os << "Global matrix:\n-------------\n";
    globalMatrix.Print(std::cout);
    os << "Local RHS on proc " << rank << ":\n-----------------\n";
    localRhsVector.Print(std::cout);
    os << "Global RHS:\n-------------\n";
    globalRhsVector.Print(std::cout);
#endif

    // --------------------------------
    // | SOLVE GLOBAL SYSTEM PARALLEL |
    // --------------------------------
    Epetra_Vector lhs(owningMap);
    lhs.PutScalar(0.0);
    Epetra_MultiVector sol = solveSystem(globalMatrix, lhs, globalRhsVector, false);
    sol.Scale(-1);
    //save solution
    EpetraExt::MultiVectorToMatlabFile((resultDirectory + "/solution.mat").c_str(), sol);
    EpetraExt::RowMatrixToMatlabFile((resultDirectory + "/globalMatrix.mat").c_str(), globalMatrix);
    EpetraExt::VectorToMatlabFile((resultDirectory + "/globalRHS.mat").c_str(), globalRhsVector);

#ifdef SHOW_SOLUTION
    std::cout << "Solution:\n---------------\n";
    sol.Print(std::cout);
#endif

    //solve system serial (for comparison only)
    if (numProc == 1)
    {
        Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
        solver.compute(hessian0_eigen);
        Eigen::VectorXd displ = solver.solve(residual_eigen);
        displ *= -1;

        //save direct solution
        std::ofstream file(resultDirectory + "/solution_direct.mat");
        file << displ;
        file.close();
    }
}


void separateDomain(int numProc, int& subDomains_X, int& subDomains_Y)
{
    double rootedValue = std::sqrt(numProc);
    if (abs(std::round(rootedValue) - rootedValue) <= 1e-12)
    {
        subDomains_X = int(rootedValue);
        subDomains_Y = subDomains_X;
    }
    else if (numProc % 2 == 0)
    {
        subDomains_X = numProc / 2;
        subDomains_Y = 2;
    }
    else
    {
        subDomains_X = numProc;
        subDomains_Y = 1;
    }
}


std::vector<int> vectorDifference(std::vector<int> rVec_1, std::vector<int> rVec_2)
{
    int n_1 = rVec_1.size();
    int n_2 = rVec_2.size();
    std::vector<int> diffVector;
    bool identValue = false;

    if ((n_1 == 0) || (n_2 == 0))
    {
        diffVector = rVec_1;
    }
    else
    {
        for (int i_1 = 0; i_1 < n_1; ++i_1)
        {
            identValue = false;
            for (int i_2 = 0; i_2 < n_2; ++i_2)
            {
                if (rVec_1[i_1] == rVec_2[i_2])
                {
                    identValue = true;
                    break;
                }
            }
            if (!identValue)
            {
                diffVector.push_back(rVec_1[i_1]);
            }
        }
    }

    return diffVector;
}


std::vector<int> vectorIntersection(std::vector<int> rVec_1, std::vector<int> rVec_2)
{
    int n_1 = rVec_1.size();
    int n_2 = rVec_2.size();
    std::vector<int> intersectVector;

    if ((n_1 == 0) || (n_2 == 0))
    {
//        intersectVector = rVec_1;
    }
    else
    {
        for (int i_1 = 0; i_1 < n_1; ++i_1)
        {
            for (int i_2 = 0; i_2 < n_2; ++i_2)
            {
                if (rVec_1[i_1] == rVec_2[i_2])
                {
                    intersectVector.push_back(rVec_1[i_1]);
                    break;
                }
            }
        }
    }

    return intersectVector;
}


void checkErrorCode(int rErrorCode, std::string rDescription)
{
    if (rErrorCode != 0)
    {
        std::cout << "ErrorCode " << rErrorCode << ":\t" << rDescription << std::endl;
    }
}


void generateSubDomainIdentifiers(int rNumSubDomains_X, int rNumSubDomains_Y, std::vector<int>& subDomainIdentifier_X, std::vector<int>& subDomainIdentifier_Y)
{
//    int numSubDomains = rNumSubDomains_X*rNumSubDomains_Y;
//    std::vector<int> ident_X(numSubDomains);
//    std::vector<int> ident_Y(numSubDomains);
    std::vector<int> ident_X;
    std::vector<int> ident_Y;
    for (int ix = 0; ix < rNumSubDomains_X; ++ix)
    {
        for (int iy = 0; iy < rNumSubDomains_Y; ++iy)
        {
            ident_X.push_back(ix);
            ident_Y.push_back(iy);
        }
    }

    subDomainIdentifier_X = ident_X;
    subDomainIdentifier_Y = ident_Y;
}





void run_Test_2D(Epetra_MpiComm rComm, double rTotalLength_X, double rTotalLength_Y, double rThickness, int rNumSubdomains_X, int rNumSubdomains_Y, int rNumElementsPerSubdomain_X, int rNumElementsPerSubdomain_Y, double rTotalForce, bool rSaveSolution = false)
{
//    Epetra_MpiComm Comm(MPI_COMM_WORLD);
    Epetra_MpiComm Comm(rComm);
    int rank = Comm.MyPID();
    int numProc = Comm.NumProc();

    //CREATE DIRECTORIES FOR LOGS AND RESULTS
    const std::string logDirectory = "2D_Test_Logs";
    const std::string resultDirectory = "2D_Test_Results";
    checkDirectories(logDirectory, resultDirectory);


    // two-dimensional problem
    const int dim = 2;

    // --------------------
    // | DOMAIN VARIABLES |
    // --------------------
    int numSubDomains_X = 0;
    int numSubDomains_Y = 0;
    if ((rNumSubdomains_X == 0) || (rNumSubdomains_Y == 0) || (rNumSubdomains_X*rNumSubdomains_Y != numProc))
        separateDomain(numProc, numSubDomains_X, numSubDomains_Y);
    else
    {
        numSubDomains_X = rNumSubdomains_X;
        numSubDomains_Y = rNumSubdomains_Y;
    }

    if (rank == 0)
    std::cout << "DOMAIN SEPARATION\n---------------\nNum Subdomains X: " << numSubDomains_X << "\nNum Subdomains Y: " << numSubDomains_Y
              << "\nNum Elements per Subdomain X: " << rNumElementsPerSubdomain_X << "\nNum Elements per Subdomain Y: " << rNumElementsPerSubdomain_Y << std::endl;

    std::vector<int> subDomains_X;
    std::vector<int> subDomains_Y;
    generateSubDomainIdentifiers(numSubDomains_X, numSubDomains_Y, subDomains_X, subDomains_Y);

    int numSubDomains = numSubDomains_X * numSubDomains_Y;
//    int* numElementsPerSubDomain_X = new int[numSubDomains_X];
//    int* numElementsPerSubDomain_Y = new int[numSubDomains_Y];
    int numElementsPerSubDomain_X = rNumElementsPerSubdomain_X;
    int numElementsPerSubDomain_Y = rNumElementsPerSubdomain_Y;
    int* numActiveDofsPerSubDomain = new int[numSubDomains];
    int* numActiveDofsTillSubDomain_inclusive = new int[numSubDomains];
    int numElementsTotal_X = 0;
    int numElementsTotal_Y = 0;
    bool setTotalLength = true;
    double* lengthPerSubDomain_X = new double[numSubDomains_X];
    double lengthTotal_X = (setTotalLength ? rTotalLength_X : 0.);
    double* lengthPerSubDomain_Y = new double[numSubDomains_Y];
    double lengthTotal_Y = (setTotalLength ? rTotalLength_Y : 0.);

    for (int i = 0; i < numSubDomains_X; ++i)
    {
//        numElementsPerSubDomain_X[i] = 1;
//        numElementsTotal_X += numElementsPerSubDomain_X[i];
        numElementsTotal_X += numElementsPerSubDomain_X;

        if (setTotalLength)
        {
            lengthPerSubDomain_X[i] = lengthTotal_X/numSubDomains_X;
        }
        else
        {
            lengthPerSubDomain_X[i] = 5;
            lengthTotal_X += lengthPerSubDomain_X[i];
        }
    }

    for (int i = 0; i < numSubDomains_Y; ++i)
    {
//        numElementsPerSubDomain_Y[i] = 1;
//        numElementsTotal_Y += numElementsPerSubDomain_Y[i];
        numElementsTotal_Y += numElementsPerSubDomain_Y;

        if (setTotalLength)
        {
            lengthPerSubDomain_Y[i] = lengthTotal_Y/numSubDomains_Y;
        }
        else
        {
            lengthPerSubDomain_Y[i] = 2.5;
            lengthTotal_Y += lengthPerSubDomain_Y[i];
        }
    }

    // ------------------------
    // | MECHANICAL VARIABLES |
    // ------------------------
    double thickness = rThickness;
    double YoungsModulus = 2.e4;
    double PoissonRatio = 0.3;
    double force = rTotalForce;

//    double fixedDisplacement = 0.;
    Eigen::VectorXd directionX(dim);
    directionX(0) = 1;
    directionX(1) = 0;
    Eigen::VectorXd directionY(dim);
    directionY(0) = 0;
    directionY(1) = 1;
    bool enableDisplacementControl = false;
    bool loadOnNode = false;


    //DEFINE STRUCTURE
    NuTo::Structure structure(dim);

    structure.SetVerboseLevel(10);
    structure.GetLogger().OpenFile(logDirectory + "/output_" + std::to_string(rank));
    structure.GetLogger().SetQuiet(false);

    auto section = NuTo::SectionPlane::Create(thickness, true);
    structure.ElementTotalSetSection(section);

    auto material = structure.ConstitutiveLawCreate("Linear_Elastic_Engineering_Stress");
    structure.ConstitutiveLawSetParameterDouble(material, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, YoungsModulus);
    structure.ConstitutiveLawSetParameterDouble(material, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, PoissonRatio);
    structure.ElementTotalSetConstitutiveLaw(material);

    // DEFINE INTERPOLATION
    const int interpolationType = structure.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::QUAD2D);
    structure.InterpolationTypeAdd(interpolationType, NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    structure.InterpolationTypeAdd(interpolationType, NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);


    //CREATE NODES
    Eigen::Vector2d nodeCoords;
    double currLength_X = 0.;
    double currLength_Y = 0.;
    int currSubDomain = 0;
    int currSubDomain_X = subDomains_X[rank];
    int currSubDomain_Y = subDomains_Y[rank];
    int currNode = 0;
    bool leaveLoops = false;
    for (int ix = 0; ix < numSubDomains_X; ++ix)
    {
        for (int iy = 0; iy < numSubDomains_Y; ++iy)
        {
            if (rank == currSubDomain)
            {
//                for (int jx = 0; jx < numElementsPerSubDomain_X[currSubDomain_X]+1; ++jx)
                for (int jx = 0; jx < numElementsPerSubDomain_X+1; ++jx)
                {
//                    nodeCoords(0) = currLength_X + jx*lengthPerSubDomain_X[currSubDomain_X]/numElementsPerSubDomain_X[currSubDomain_X];
                    nodeCoords(0) = currLength_X + jx*lengthPerSubDomain_X[currSubDomain_X]/numElementsPerSubDomain_X;
//                    for (int jy = 0; jy < numElementsPerSubDomain_Y[currSubDomain_Y]+1; ++jy)
                    for (int jy = 0; jy < numElementsPerSubDomain_Y+1; ++jy)
                    {
//                        nodeCoords(1) = currLength_Y + jy*lengthPerSubDomain_Y[currSubDomain_Y]/numElementsPerSubDomain_Y[currSubDomain_Y];
                        nodeCoords(1) = currLength_Y + jy*lengthPerSubDomain_Y[currSubDomain_Y]/numElementsPerSubDomain_Y;
                        structure.NodeCreate(currNode, nodeCoords);
                        ++currNode;
                    }
                }
                leaveLoops = true;
            }
            else
            {
                ++currSubDomain;
            }

            currLength_Y += lengthPerSubDomain_Y[iy];
            if (leaveLoops)
                break;
        }
        currLength_Y = 0;
        currLength_X += lengthPerSubDomain_X[ix];
        if (leaveLoops)
            break;
    }


    //CREATE ELEMENTS
    std::vector<int> elementNodes(4);
    int currElement = 0;
    currSubDomain = 0;
    currNode = 0;
    for (int ix = 0; ix < numSubDomains_X; ++ix)
    {
        for (int iy = 0; iy < numSubDomains_Y; ++iy)
        {
            if (rank == currSubDomain)
            {
//                for (int jx = 0; jx < numElementsPerSubDomain_X[currSubDomain_X]; ++jx)
                for (int jx = 0; jx < numElementsPerSubDomain_X; ++jx)
                {
//                    for (int jy = 0; jy < numElementsPerSubDomain_Y[currSubDomain_Y]; ++jy)
                    for (int jy = 0; jy < numElementsPerSubDomain_Y; ++jy)
                    {
                        elementNodes[0] = currNode;
                        elementNodes[1] = currNode + 1;
//                        elementNodes[2] = currNode + (numElementsPerSubDomain_Y[currSubDomain_Y]+1) + 1;
//                        elementNodes[3] = currNode + (numElementsPerSubDomain_Y[currSubDomain_Y]+1);
                        elementNodes[2] = currNode + (numElementsPerSubDomain_Y+1) + 1;
                        elementNodes[3] = currNode + (numElementsPerSubDomain_Y+1);
                        structure.ElementCreate(interpolationType, elementNodes);
                        ++currNode;
                        ++currElement;
                    }
                    ++currNode;
                }
            }
            ++currSubDomain;

        }
    }

    structure.ElementTotalSetSection(section);
    structure.ElementTotalSetConstitutiveLaw(material);
    structure.ElementTotalConvertToInterpolationType();


    // ----------------------------------------
    // | DEFINE BOUNDARY AND OVERLAPPING AREAS|
    // ----------------------------------------
    int nodesBCLeft = structure.GroupCreate(NuTo::eGroupId::Nodes);
    int nodesBCRight = structure.GroupCreate(NuTo::eGroupId::Nodes);
    structure.GroupAddNodeCoordinateRange(nodesBCLeft, 0, -1e-6, 1e-6);
    if (loadOnNode)
    {
        Eigen::Vector2d nodeCoordinates;
        nodeCoordinates(0) = lengthTotal_X;
        nodeCoordinates(1) = lengthTotal_Y;
        structure.GroupAddNodeRadiusRange(nodesBCRight, nodeCoordinates, 0, 1.e-6);
    }
    else
        structure.GroupAddNodeCoordinateRange(nodesBCRight, 0, lengthTotal_X-1e-6, lengthTotal_X+1e-6);


    //COMPUTE OVERLAPPING AREAS
    int numOverlapping_X = numSubDomains_X - 1;
    int numOverlapping_Y = numSubDomains_Y - 1;
    int* nodeGroupsOverlap_X = new int[numOverlapping_X];
    int* nodeGroupsOverlap_Y = new int[numOverlapping_Y];

    for (int ix = 0; ix < numOverlapping_X; ++ix)
    {
        nodeGroupsOverlap_X[ix] = structure.GroupCreate(NuTo::eGroupId::Nodes);
    }

    for (int iy = 0; iy < numOverlapping_Y; ++iy)
    {
        nodeGroupsOverlap_Y[iy] = structure.GroupCreate(NuTo::eGroupId::Nodes);
    }

    currLength_X = lengthPerSubDomain_X[0];
    for (int ix = 0; ix < numOverlapping_X; ++ix)
    {
        structure.GroupAddNodeCoordinateRange(nodeGroupsOverlap_X[ix], 0, currLength_X - 1e-6, currLength_X + 1e-6);
        currLength_X += lengthPerSubDomain_X[ix+1];
    }

    currLength_Y = lengthPerSubDomain_Y[0];
    for (int iy = 0; iy < numOverlapping_Y; ++iy)
    {
        structure.GroupAddNodeCoordinateRange(nodeGroupsOverlap_Y[iy], 1, currLength_Y - 1e-6, currLength_Y + 1e-6);
        currLength_Y += lengthPerSubDomain_Y[iy+1];
    }

    delete[] lengthPerSubDomain_X;
    delete[] lengthPerSubDomain_Y;

    //FIX LEFT BOUNDARY
//    structure.ConstraintLinearSetDisplacementNodeGroup(nodesBCLeft, directionX, fixedDisplacement);
//    structure.ConstraintLinearSetDisplacementNodeGroup(nodesBCLeft, directionY, fixedDisplacement);

    //SET LOAD ON RIGHT BOUNDARY
    structure.SetNumLoadCases(1);
    int loadID = -1;
    if (enableDisplacementControl)
    {
//        loadID = structure.ConstraintLinearSetDisplacementNode(numElements, directionX, displacement);
    }
    else
    {
        loadID = structure.LoadCreateNodeGroupForce(0, nodesBCRight, directionX, force);
    }

    //COMPUTE (LOCAL) HESSIAN AND RESIDUAL
    NuTo::StructureOutputBlockMatrix hessian0 = structure.BuildGlobalHessian0();
    NuTo::StructureOutputBlockVector dofs(structure.GetDofStatus(), true);
    dofs.J.SetZero();
    dofs.K.SetZero();

    NuTo::StructureOutputBlockVector residual = hessian0*dofs - structure.BuildGlobalExternalLoadVector(0) + structure.BuildGlobalInternalGradient();
//    structure.Info();
//    std::cout << "NuTo - Solve starts..." << std::endl;
//    structure.SolveGlobalSystemStaticElastic();
//    std::cout << "NuTo - Solve finished" << std::endl;

    Eigen::SparseMatrix<double> hessian0_eigen = hessian0.JJ.ExportToEigenSparseMatrix();
    Eigen::MatrixXd residual_eigen = residual.J.Export();


    std::vector<int> structNodes = structure.GroupGetMemberIds(nodesBCLeft);
    for (int currNode : structNodes)
    {
        std::vector<int> structDofIDs = structure.NodeGetDofIds(currNode, NuTo::Node::eDof::DISPLACEMENTS);
        for (int dofID : structDofIDs)
        {
            hessian0_eigen.coeffRef(dofID, dofID) = 1.e15;
        }
    }

    int ownedNumActiveDOFs = -1;
//    if (currSubDomain_Y == 0)
//    {
////        ownedNumActiveDOFs = numElementsPerSubDomain_X[currSubDomain_X]*(numElementsPerSubDomain_Y[currSubDomain_Y]+1)*dim;
//        ownedNumActiveDOFs = numElementsPerSubDomain_X*(numElementsPerSubDomain_Y+1)*dim;
//    }
//    else
//    {
////        ownedNumActiveDOFs = numElementsPerSubDomain_X[currSubDomain_X]*numElementsPerSubDomain_Y[currSubDomain_Y]*dim;
//        ownedNumActiveDOFs = numElementsPerSubDomain_X*numElementsPerSubDomain_Y*dim;
//    }
    if (subDomains_X[rank] == 0)
    {
        if (subDomains_Y[rank] == 0)
        {
            ownedNumActiveDOFs = (numElementsPerSubDomain_X+1)*(numElementsPerSubDomain_Y+1)*dim;
        }
        else
        {
            ownedNumActiveDOFs = (numElementsPerSubDomain_X+1)*numElementsPerSubDomain_Y*dim;
        }
    }
    else
    {
        if (subDomains_Y[rank] == 0)
        {
            ownedNumActiveDOFs = numElementsPerSubDomain_X*(numElementsPerSubDomain_Y+1)*dim;
        }
        else
        {
            ownedNumActiveDOFs = numElementsPerSubDomain_X*numElementsPerSubDomain_Y*dim;
        }
    }


    //GET IDS OF OVERLAPPING NODES AND LEFT (FIXED) BOUNDARY
    std::vector<int> overlapNodeIDs_X;
    std::vector<int> overlapNodeIDs_Y_total;
    std::vector<int> overlapNodeIDs_Corner;

    if (numSubDomains_X > 1)
    {
        if (currSubDomain_X == 0)
        {
            overlapNodeIDs_X = structure.GroupGetMemberIds(nodeGroupsOverlap_X[currSubDomain_X]);
        }
        else
        {
            overlapNodeIDs_X = structure.GroupGetMemberIds(nodeGroupsOverlap_X[currSubDomain_X-1]);
        }
    }

    if (numSubDomains_Y > 1)
    {
        if (currSubDomain_Y == 0)
        {
            overlapNodeIDs_Y_total = structure.GroupGetMemberIds(nodeGroupsOverlap_Y[currSubDomain_Y]);
        }
        else
        {
            overlapNodeIDs_Y_total = structure.GroupGetMemberIds(nodeGroupsOverlap_Y[currSubDomain_Y-1]);
        }
    }

    delete[] nodeGroupsOverlap_X;
    delete[] nodeGroupsOverlap_Y;

    std::vector<int> leftNodeIDs = structure.GroupGetMemberIds(nodesBCLeft);
//    int leftNodeCount = leftNodeIDs.size();
    int leftNodeCount = 0;
    std::vector<int> overlapNodeIDs_Y;
    overlapNodeIDs_Corner = vectorIntersection(overlapNodeIDs_X, overlapNodeIDs_Y_total);
//    overlapNodeIDs_Y = vectorDifference(overlapNodeIDs_Y_total, leftNodeIDs);
    overlapNodeIDs_Y = overlapNodeIDs_Y_total;
    overlapNodeIDs_Y = vectorDifference(overlapNodeIDs_Y, overlapNodeIDs_Corner);
    overlapNodeIDs_X = vectorDifference(overlapNodeIDs_X, overlapNodeIDs_Corner);
    int overlappingNodeCount_X = overlapNodeIDs_X.size();
    int overlappingNodeCount_Y = overlapNodeIDs_Y.size();
    int overlappingNodeCount_Corner = overlapNodeIDs_Corner.size();

    int subDomainCounter = 0;
    for (int k = 0; k < numSubDomains_Y; ++k)
    {
//        numActiveDofsPerSubDomain[subDomainCounter] = (numElementsPerSubDomain_X[0]) * (numElementsPerSubDomain_Y[k] + 1)*dim;
//        numActiveDofsPerSubDomain[subDomainCounter] = (numElementsPerSubDomain_X) * (numElementsPerSubDomain_Y + 1)*dim;
        numActiveDofsPerSubDomain[subDomainCounter] = (numElementsPerSubDomain_X + 1) * (numElementsPerSubDomain_Y + 1)*dim;
        if (subDomainCounter == 0)
            numActiveDofsTillSubDomain_inclusive[subDomainCounter] = numActiveDofsPerSubDomain[subDomainCounter];
        else
        {
            if ((subDomains_X[rank] == 0) && (numSubDomains_X > 1))
                numActiveDofsTillSubDomain_inclusive[subDomainCounter] = numActiveDofsTillSubDomain_inclusive[subDomainCounter-1]
                                                                        + numActiveDofsPerSubDomain[subDomainCounter] - dim*(overlappingNodeCount_Y+overlappingNodeCount_Corner);
            else
                numActiveDofsTillSubDomain_inclusive[subDomainCounter] = numActiveDofsTillSubDomain_inclusive[subDomainCounter-1]
                                                                        + numActiveDofsPerSubDomain[subDomainCounter] - dim*(overlappingNodeCount_Y + overlappingNodeCount_Corner);

        }
        ++subDomainCounter;
    }


    for (int j = 1; j < numSubDomains_X; ++j)
    {
        for (int k = 0; k < numSubDomains_Y; ++k)
        {
//            numActiveDofsPerSubDomain[subDomainCounter] = (numElementsPerSubDomain_X[j] + 1) * (numElementsPerSubDomain_Y[k] + 1)*dim;
            numActiveDofsPerSubDomain[subDomainCounter] = (numElementsPerSubDomain_X + 1) * (numElementsPerSubDomain_Y + 1)*dim;
//            if (subDomains_X[rank] == 0)
            if ((subDomains_X[rank] == 0) && (numSubDomains_X > 1))
            {
                if (subDomains_Y[subDomainCounter] == 0)
                    numActiveDofsTillSubDomain_inclusive[subDomainCounter] = numActiveDofsTillSubDomain_inclusive[subDomainCounter-1]
                                                                        + numActiveDofsPerSubDomain[subDomainCounter] - dim*(overlappingNodeCount_X+overlappingNodeCount_Corner);
                else
//                    numActiveDofsTillSubDomain_inclusive[subDomainCounter] = numActiveDofsTillSubDomain_inclusive[subDomainCounter-1]
//                                                                        + numActiveDofsPerSubDomain[subDomainCounter] - dim*(overlappingNodeCount_X+overlappingNodeCount_Corner+overlappingNodeCount_Y+1);
                    numActiveDofsTillSubDomain_inclusive[subDomainCounter] = numActiveDofsTillSubDomain_inclusive[subDomainCounter-1]
                                                                        + numActiveDofsPerSubDomain[subDomainCounter] - dim*(overlappingNodeCount_X+overlappingNodeCount_Corner+overlappingNodeCount_Y);
            }
            else
            {
                if (subDomains_Y[subDomainCounter] == 0)
                    numActiveDofsTillSubDomain_inclusive[subDomainCounter] = numActiveDofsTillSubDomain_inclusive[subDomainCounter-1]
                                                                        + numActiveDofsPerSubDomain[subDomainCounter] - dim*(overlappingNodeCount_X+overlappingNodeCount_Corner);
                else
                    numActiveDofsTillSubDomain_inclusive[subDomainCounter] = numActiveDofsTillSubDomain_inclusive[subDomainCounter-1]
                                                                        + numActiveDofsPerSubDomain[subDomainCounter] - dim*(overlappingNodeCount_X+overlappingNodeCount_Corner+overlappingNodeCount_Y);
            }

            ++subDomainCounter;
        }
    }
    int localNumActiveDOFs = numActiveDofsPerSubDomain[rank];
    delete[] numActiveDofsPerSubDomain;

    int internDofCounter = -1;
    if (rank == 0)
        internDofCounter = 0;
    else
        internDofCounter = numActiveDofsTillSubDomain_inclusive[rank-1];

    //DEFINE LOCAL-TO-GLOBAL DOF MAPPING
    bool leftNodeIDFound = false;
    bool overlappingNodeIDFound_X = false;
    bool overlappingNodeIDFound_Y = false;
    bool overlappingNodeIDFound_Corner = false;
    int* local2GlobalMapping = new int[localNumActiveDOFs];
    for (int i = 0; i < structure.GetNumNodes(); ++i)
    {
        leftNodeIDFound = false;
        overlappingNodeIDFound_X = false;
        overlappingNodeIDFound_Y = false;
        overlappingNodeIDFound_Corner = false;

        for (int j = 0; j < leftNodeCount; ++j)
        {
            if (leftNodeIDs[j] == i)
            {
                leftNodeIDFound = true;
            }
        }

        for (int j = 0; j < overlappingNodeCount_X; ++j)
        {
            if (overlapNodeIDs_X[j] == i)
            {
                overlappingNodeIDFound_X = true;
            }
        }

        for (int j = 0; j < overlappingNodeCount_Y; ++j)
        {
            if (overlapNodeIDs_Y[j] == i)
            {
                overlappingNodeIDFound_Y = true;
            }
        }

        if (overlappingNodeCount_Corner > 0)
        {
            if (overlapNodeIDs_Corner[0] == i)
                overlappingNodeIDFound_Corner = true;
        }

        if (!leftNodeIDFound && !overlappingNodeIDFound_X && !overlappingNodeIDFound_Y && !overlappingNodeIDFound_Corner)
        {
            std::vector<int> dofIDs = structure.NodeGetDofIds(i, NuTo::Node::eDof::DISPLACEMENTS);
            for (int dofID : dofIDs)
            {
                local2GlobalMapping[dofID] = internDofCounter;
                ++internDofCounter;
            }
        }
    }

    internDofCounter = 0;
    int newID = -1;
    if (numProc > 1)
    {
        if ((subDomains_X[rank] == 0) && (numSubDomains_X > 1))
        {
            if (subDomains_Y[rank] == 0)
            {
                newID = numActiveDofsTillSubDomain_inclusive[rank] - (overlappingNodeCount_X+overlappingNodeCount_Corner+overlappingNodeCount_Y)*dim;
                for (int nodeID_X : overlapNodeIDs_X)
                {
                    std::vector<int> overlapDofIDs_X = structure.NodeGetDofIds(nodeID_X, NuTo::Node::eDof::DISPLACEMENTS);
                    for (int dofID : overlapDofIDs_X)
                    {
                        local2GlobalMapping[dofID] = newID + internDofCounter;
                        ++internDofCounter;
                    }
                }

                for (int nodeID_Y : overlapNodeIDs_Y)
                {
                    std::vector<int> overlapDofIDs_Y = structure.NodeGetDofIds(nodeID_Y, NuTo::Node::eDof::DISPLACEMENTS);
                    for (int dofID : overlapDofIDs_Y)
                    {
                        local2GlobalMapping[dofID] = newID + internDofCounter;
                        ++internDofCounter;
                    }
                }

                if (overlappingNodeCount_Corner > 0)
                {
                    int nodeID_Corner = overlapNodeIDs_Corner[0];
                    std::vector<int> overlapDofIDs_Corner = structure.NodeGetDofIds(nodeID_Corner, NuTo::Node::eDof::DISPLACEMENTS);
                    for (int dofID : overlapDofIDs_Corner)
                    {
                        local2GlobalMapping[dofID] = newID + internDofCounter;
                        ++internDofCounter;
                    }
                }

            }
            else if (subDomains_Y[rank] == 1)
            {
                newID = numActiveDofsTillSubDomain_inclusive[rank] - (overlappingNodeCount_X)*dim;
                for (int nodeID_X : overlapNodeIDs_X)
                {
                    std::vector<int> overlapDofIDs_X = structure.NodeGetDofIds(nodeID_X, NuTo::Node::eDof::DISPLACEMENTS);
                    for (int dofID : overlapDofIDs_X)
                    {
                        local2GlobalMapping[dofID] = newID + internDofCounter;
                        ++internDofCounter;
                    }
                }

                newID = numActiveDofsTillSubDomain_inclusive[rank-1] - (overlappingNodeCount_X+overlappingNodeCount_Corner+overlappingNodeCount_Y)*dim;
//                newID = numActiveDofsTillSubDomain_inclusive[rank-1] - (overlappingNodeCount_X+overlappingNodeCount_Corner+overlappingNodeCount_X)*dim;
                for (int nodeID_Y : overlapNodeIDs_Y)
                {
                    std::vector<int> overlapDofIDs_Y = structure.NodeGetDofIds(nodeID_Y, NuTo::Node::eDof::DISPLACEMENTS);
                    for (int dofID : overlapDofIDs_Y)
                    {
                        local2GlobalMapping[dofID] = newID + internDofCounter;
                        ++internDofCounter;
                    }
                }

//                newID = numActiveDofsTillSubDomain_inclusive[rank-1] - (overlappingNodeCount_X+overlappingNodeCount_Corner+overlappingNodeCount_Y)*dim;
                if (overlappingNodeCount_Corner > 0)
                {
                    int nodeID_Corner = overlapNodeIDs_Corner[0];
                    std::vector<int> overlapDofIDs_Corner = structure.NodeGetDofIds(nodeID_Corner, NuTo::Node::eDof::DISPLACEMENTS);
                    for (int dofID : overlapDofIDs_Corner)
                    {
                        local2GlobalMapping[dofID] = newID + internDofCounter;
                        ++internDofCounter;
                    }
                }
            }
            else
            {
                newID = numActiveDofsTillSubDomain_inclusive[rank] - (overlappingNodeCount_X)*dim;
                for (int nodeID_X : overlapNodeIDs_X)
                {
                    std::vector<int> overlapDofIDs_X = structure.NodeGetDofIds(nodeID_X, NuTo::Node::eDof::DISPLACEMENTS);
                    for (int dofID : overlapDofIDs_X)
                    {
                        local2GlobalMapping[dofID] = newID + internDofCounter;
                        ++internDofCounter;
                    }
                }

//                newID = numActiveDofsTillSubDomain_inclusive[rank-1] - (overlappingNodeCount_X+overlappingNodeCount_Corner+overlappingNodeCount_Y)*dim;
//                newID = numActiveDofsTillSubDomain_inclusive[rank-1] - (overlappingNodeCount_X+overlappingNodeCount_Corner+overlappingNodeCount_X)*dim;
                newID = numActiveDofsTillSubDomain_inclusive[rank-2] + (overlappingNodeCount_X-1)*dim - overlappingNodeCount_X*dim;
                int offset = 0;
                for (int nodeID_Y : overlapNodeIDs_Y)
                {
                    std::vector<int> overlapDofIDs_Y = structure.NodeGetDofIds(nodeID_Y, NuTo::Node::eDof::DISPLACEMENTS);
                    for (int dofID : overlapDofIDs_Y)
                    {
//                        local2GlobalMapping[dofID] = newID + internDofCounter;
                        local2GlobalMapping[dofID] = newID + internDofCounter + offset*dim*(numElementsPerSubDomain_Y-1);
                        ++internDofCounter;
                    }
                    ++offset;
                }

                newID = numActiveDofsTillSubDomain_inclusive[rank-1] - (overlappingNodeCount_X+overlappingNodeCount_Corner+overlappingNodeCount_Y)*dim;
                if (overlappingNodeCount_Corner > 0)
                {
                    int nodeID_Corner = overlapNodeIDs_Corner[0];
                    std::vector<int> overlapDofIDs_Corner = structure.NodeGetDofIds(nodeID_Corner, NuTo::Node::eDof::DISPLACEMENTS);
                    for (int dofID : overlapDofIDs_Corner)
                    {
                        local2GlobalMapping[dofID] = newID + internDofCounter;
                        ++internDofCounter;
                    }
                }


            }
        }
        else if ((subDomains_X[rank] == 0) && (numSubDomains_X == 1))
        {
            if (subDomains_Y[rank] == 0)
            {
                newID = numActiveDofsTillSubDomain_inclusive[rank] - (overlappingNodeCount_Y)*dim;
                for (int nodeID_Y : overlapNodeIDs_Y)
                {
                    std::vector<int> overlapDofIDs_Y = structure.NodeGetDofIds(nodeID_Y, NuTo::Node::eDof::DISPLACEMENTS);
                    for (int dofID : overlapDofIDs_Y)
                    {
                        local2GlobalMapping[dofID] = newID + internDofCounter;
                        ++internDofCounter;
                    }
                }
            }
            else if (subDomains_Y[rank] == 1)
            {
                newID = numActiveDofsTillSubDomain_inclusive[rank-1] - (overlappingNodeCount_Y)*dim;
                for (int nodeID_Y : overlapNodeIDs_Y)
                {
                    std::vector<int> overlapDofIDs_Y = structure.NodeGetDofIds(nodeID_Y, NuTo::Node::eDof::DISPLACEMENTS);
                    for (int dofID : overlapDofIDs_Y)
                    {
                        local2GlobalMapping[dofID] = newID + internDofCounter;
                        ++internDofCounter;
                    }
                }
            }
            else
            {
//                newID = numActiveDofsTillSubDomain_inclusive[rank-1] - (overlappingNodeCount_Y)*dim;
                newID = numActiveDofsTillSubDomain_inclusive[rank-1] - (overlappingNodeCount_Y+(numElementsPerSubDomain_X)*(numElementsPerSubDomain_Y-1))*dim;
                int offset = 0;
                for (int nodeID_Y : overlapNodeIDs_Y)
                {
                    std::vector<int> overlapDofIDs_Y = structure.NodeGetDofIds(nodeID_Y, NuTo::Node::eDof::DISPLACEMENTS);
                    for (int dofID : overlapDofIDs_Y)
                    {
                        local2GlobalMapping[dofID] = newID + internDofCounter + offset*dim*(numElementsPerSubDomain_Y-1);;
                        ++internDofCounter;
                    }
                    ++offset;
                }
            }
        }
        else if (subDomains_X[rank] == 1)
        {
            if (subDomains_Y[rank] == 0)
            {
//                newID = numActiveDofsTillSubDomain_inclusive[rank-numSubDomains_Y] - (overlappingNodeCount_X+overlappingNodeCount_Corner)*dim;
                newID = numActiveDofsTillSubDomain_inclusive[rank-numSubDomains_Y] - (overlappingNodeCount_X+overlappingNodeCount_Corner+overlappingNodeCount_Y)*dim;
                for (int nodeID_X : overlapNodeIDs_X)
                {
                    std::vector<int> overlapDofIDs_X = structure.NodeGetDofIds(nodeID_X, NuTo::Node::eDof::DISPLACEMENTS);
                    for (int dofID : overlapDofIDs_X)
                    {
//                        if (overlappingNodeCount_Y > 0)
////                            local2GlobalMapping[dofID] = newID - (overlappingNodeCount_Y-1)*dim + internDofCounter;
//                            local2GlobalMapping[dofID] = newID - (overlappingNodeCount_Y)*dim + internDofCounter;
//                        else
                            local2GlobalMapping[dofID] = newID + internDofCounter;

                        ++internDofCounter;
                    }
                }

//                newID = numActiveDofsTillSubDomain_inclusive[rank] - (overlappingNodeCount_X+overlappingNodeCount_Corner+overlappingNodeCount_Y)*dim;
                newID = numActiveDofsTillSubDomain_inclusive[rank] - (overlappingNodeCount_X+overlappingNodeCount_Y)*dim;
                for (int nodeID_Y : overlapNodeIDs_Y)
                {
                    std::vector<int> overlapDofIDs_Y = structure.NodeGetDofIds(nodeID_Y, NuTo::Node::eDof::DISPLACEMENTS);
                    for (int dofID : overlapDofIDs_Y)
                    {
                        local2GlobalMapping[dofID] = newID + internDofCounter;
                        ++internDofCounter;
                    }
                }

                newID = numActiveDofsTillSubDomain_inclusive[rank-numSubDomains_Y] - (overlappingNodeCount_X+overlappingNodeCount_Corner+overlappingNodeCount_Y)*dim;
                if (overlappingNodeCount_Corner > 0)
                {
                    int nodeID_Corner = overlapNodeIDs_Corner[0];
                    std::vector<int> overlapDofIDs_Corner = structure.NodeGetDofIds(nodeID_Corner, NuTo::Node::eDof::DISPLACEMENTS);
                    for (int dofID : overlapDofIDs_Corner)
                    {
                        local2GlobalMapping[dofID] = newID + internDofCounter;

                        ++internDofCounter;
                    }
                }

            }
            else if (subDomains_Y[rank] == 1)
            {
                newID = numActiveDofsTillSubDomain_inclusive[rank-numSubDomains_Y] - (overlappingNodeCount_X)*dim;
                for (int nodeID_X : overlapNodeIDs_X)
                {
                    std::vector<int> overlapDofIDs_X = structure.NodeGetDofIds(nodeID_X, NuTo::Node::eDof::DISPLACEMENTS);
                    for (int dofID : overlapDofIDs_X)
                    {
                        local2GlobalMapping[dofID] = newID + internDofCounter;
                        ++internDofCounter;
                    }
                }

                newID = numActiveDofsTillSubDomain_inclusive[rank-1] - (overlappingNodeCount_X+overlappingNodeCount_Y)*dim;
                for (int nodeID_Y : overlapNodeIDs_Y)
                {
                    std::vector<int> overlapDofIDs_Y = structure.NodeGetDofIds(nodeID_Y, NuTo::Node::eDof::DISPLACEMENTS);
                    for (int dofID : overlapDofIDs_Y)
                    {
                        local2GlobalMapping[dofID] = newID + internDofCounter;
                        ++internDofCounter;
                    }
                }

                if (overlappingNodeCount_Corner > 0)
                {
//                    newID = numActiveDofsTillSubDomain_inclusive[rank-numSubDomains_Y-1] - (overlappingNodeCount_X+overlappingNodeCount_Corner)*dim;
                    newID = numActiveDofsTillSubDomain_inclusive[rank-numSubDomains_Y-1] - (overlappingNodeCount_X+overlappingNodeCount_Corner+overlappingNodeCount_Y)*dim;
                    int nodeID_Corner = overlapNodeIDs_Corner[0];
                    std::vector<int> overlapDofIDs_Corner = structure.NodeGetDofIds(nodeID_Corner, NuTo::Node::eDof::DISPLACEMENTS);
                    for (int dofID : overlapDofIDs_Corner)
                    {
//                        if (overlappingNodeCount_Y > 0)
////                            local2GlobalMapping[dofID] = newID - (overlappingNodeCount_Y-1)*dim + internDofCounter;
//                            local2GlobalMapping[dofID] = newID - (overlappingNodeCount_Y)*dim + internDofCounter;
//                        else
                        local2GlobalMapping[dofID] = newID + internDofCounter;
                        ++internDofCounter;
                    }
                }

            }
            else
            {
                newID = numActiveDofsTillSubDomain_inclusive[rank-numSubDomains_Y] - (overlappingNodeCount_X)*dim;
                for (int nodeID_X : overlapNodeIDs_X)
                {
                    std::vector<int> overlapDofIDs_X = structure.NodeGetDofIds(nodeID_X, NuTo::Node::eDof::DISPLACEMENTS);
                    for (int dofID : overlapDofIDs_X)
                    {
                        local2GlobalMapping[dofID] = newID + internDofCounter;
                        ++internDofCounter;
                    }
                }

//                newID = numActiveDofsTillSubDomain_inclusive[rank-1] - (overlappingNodeCount_X+overlappingNodeCount_Y)*dim;
                newID = numActiveDofsTillSubDomain_inclusive[rank-1] - (overlappingNodeCount_X+overlappingNodeCount_Y+(numElementsPerSubDomain_X-1)*(numElementsPerSubDomain_Y-1))*dim;
                int offset = 0;
                for (int nodeID_Y : overlapNodeIDs_Y)
                {
                    std::vector<int> overlapDofIDs_Y = structure.NodeGetDofIds(nodeID_Y, NuTo::Node::eDof::DISPLACEMENTS);
                    for (int dofID : overlapDofIDs_Y)
                    {
                        local2GlobalMapping[dofID] = newID + internDofCounter + offset*dim*(numElementsPerSubDomain_Y-1);
                        ++internDofCounter;
                    }
                    ++offset;
                }

                if (overlappingNodeCount_Corner > 0)
                {
//                    newID = numActiveDofsTillSubDomain_inclusive[rank-numSubDomains_Y-1] - (overlappingNodeCount_X+overlappingNodeCount_Corner)*dim;
                    newID = numActiveDofsTillSubDomain_inclusive[rank-numSubDomains_Y-1] - (overlappingNodeCount_X+overlappingNodeCount_Corner+overlappingNodeCount_Y)*dim;
                    int nodeID_Corner = overlapNodeIDs_Corner[0];
                    std::vector<int> overlapDofIDs_Corner = structure.NodeGetDofIds(nodeID_Corner, NuTo::Node::eDof::DISPLACEMENTS);
                    for (int dofID : overlapDofIDs_Corner)
                    {
//                        if (overlappingNodeCount_Y > 0)
////                            local2GlobalMapping[dofID] = newID - (overlappingNodeCount_Y-1)*dim + internDofCounter;
//                            local2GlobalMapping[dofID] = newID - (overlappingNodeCount_Y)*dim + internDofCounter;
//                        else
                            local2GlobalMapping[dofID] = newID + internDofCounter;
                        ++internDofCounter;
                    }
                }

            }
        }
        else
        {
            if (subDomains_Y[rank] == 0)
            {
//                newID = numActiveDofsTillSubDomain_inclusive[rank-numSubDomains_Y] - (overlappingNodeCount_X+overlappingNodeCount_Corner)*dim;
//                newID = numActiveDofsTillSubDomain_inclusive[rank-numSubDomains_Y] - (overlappingNodeCount_X+overlappingNodeCount_Corner)*dim;
                newID = numActiveDofsTillSubDomain_inclusive[rank-numSubDomains_Y] - (overlappingNodeCount_X+overlappingNodeCount_Y)*dim;
                for (int nodeID_X : overlapNodeIDs_X)
                {
                    std::vector<int> overlapDofIDs_X = structure.NodeGetDofIds(nodeID_X, NuTo::Node::eDof::DISPLACEMENTS);
                    for (int dofID : overlapDofIDs_X)
                    {
//                        if (overlappingNodeCount_Y > 0)
//                            local2GlobalMapping[dofID] = newID - (overlappingNodeCount_Y-1)*dim + internDofCounter;
////                            local2GlobalMapping[dofID] = newID - (overlappingNodeCount_Y)*dim + internDofCounter;
//                        else
                            local2GlobalMapping[dofID] = newID + internDofCounter;

                        ++internDofCounter;
                    }
                }

                newID = numActiveDofsTillSubDomain_inclusive[rank] - (overlappingNodeCount_X+overlappingNodeCount_Y)*dim;
                for (int nodeID_Y : overlapNodeIDs_Y)
                {
                    std::vector<int> overlapDofIDs_Y = structure.NodeGetDofIds(nodeID_Y, NuTo::Node::eDof::DISPLACEMENTS);
                    for (int dofID : overlapDofIDs_Y)
                    {
                        local2GlobalMapping[dofID] = newID + internDofCounter;
                        ++internDofCounter;
                    }
                }

                newID = numActiveDofsTillSubDomain_inclusive[rank-numSubDomains_Y] - (overlappingNodeCount_X+overlappingNodeCount_Corner+overlappingNodeCount_Y)*dim;
                if (overlappingNodeCount_Corner > 0)
                {
                    int nodeID_Corner = overlapNodeIDs_Corner[0];
                    std::vector<int> overlapDofIDs_Corner = structure.NodeGetDofIds(nodeID_Corner, NuTo::Node::eDof::DISPLACEMENTS);
                    for (int dofID : overlapDofIDs_Corner)
                    {
//                        if (overlappingNodeCount_Y > 0)
//                            local2GlobalMapping[dofID] = newID - (overlappingNodeCount_Y-1)*dim + internDofCounter;
////                            local2GlobalMapping[dofID] = newID - (overlappingNodeCount_Y)*dim + internDofCounter;
//                        else
                            local2GlobalMapping[dofID] = newID + internDofCounter;

                        ++internDofCounter;
                    }
                }

            }
            else if (subDomains_Y[rank] == 1)
            {
                newID = numActiveDofsTillSubDomain_inclusive[rank-numSubDomains_Y] - (overlappingNodeCount_X)*dim;
                for (int nodeID_X : overlapNodeIDs_X)
                {
                    std::vector<int> overlapDofIDs_X = structure.NodeGetDofIds(nodeID_X, NuTo::Node::eDof::DISPLACEMENTS);
                    for (int dofID : overlapDofIDs_X)
                    {
                        local2GlobalMapping[dofID] = newID + internDofCounter;
                        ++internDofCounter;
                    }
                }

                newID = numActiveDofsTillSubDomain_inclusive[rank-1] - (overlappingNodeCount_X+overlappingNodeCount_Y)*dim;
                for (int nodeID_Y : overlapNodeIDs_Y)
                {
                    std::vector<int> overlapDofIDs_Y = structure.NodeGetDofIds(nodeID_Y, NuTo::Node::eDof::DISPLACEMENTS);
                    for (int dofID : overlapDofIDs_Y)
                    {
                        local2GlobalMapping[dofID] = newID + internDofCounter;
                        ++internDofCounter;
                    }
                }

                if (overlappingNodeCount_Corner > 0)
                {
                    newID = numActiveDofsTillSubDomain_inclusive[rank-numSubDomains_Y-1] - (overlappingNodeCount_X+overlappingNodeCount_Corner+overlappingNodeCount_Y)*dim;
                    int nodeID_Corner = overlapNodeIDs_Corner[0];
                    std::vector<int> overlapDofIDs_Corner = structure.NodeGetDofIds(nodeID_Corner, NuTo::Node::eDof::DISPLACEMENTS);
                    for (int dofID : overlapDofIDs_Corner)
                    {
//                        if (overlappingNodeCount_Y > 0)
//                            local2GlobalMapping[dofID] = newID - (overlappingNodeCount_Y-1)*dim + internDofCounter;
////                            local2GlobalMapping[dofID] = newID - (overlappingNodeCount_Y)*dim + internDofCounter;
//                        else
                            local2GlobalMapping[dofID] = newID + internDofCounter;
                        ++internDofCounter;
                    }
                }

            }
            else
            {
                newID = numActiveDofsTillSubDomain_inclusive[rank-numSubDomains_Y] - (overlappingNodeCount_X)*dim;
                for (int nodeID_X : overlapNodeIDs_X)
                {
                    std::vector<int> overlapDofIDs_X = structure.NodeGetDofIds(nodeID_X, NuTo::Node::eDof::DISPLACEMENTS);
                    for (int dofID : overlapDofIDs_X)
                    {
                        local2GlobalMapping[dofID] = newID + internDofCounter;
                        ++internDofCounter;
                    }
                }

//                newID = numActiveDofsTillSubDomain_inclusive[rank-1] - (overlappingNodeCount_X+overlappingNodeCount_Y)*dim;
                newID = numActiveDofsTillSubDomain_inclusive[rank-1] - (overlappingNodeCount_X+overlappingNodeCount_Y+(numElementsPerSubDomain_X-1)*(numElementsPerSubDomain_Y-1))*dim;
                int offset = 0;
                for (int nodeID_Y : overlapNodeIDs_Y)
                {
                    std::vector<int> overlapDofIDs_Y = structure.NodeGetDofIds(nodeID_Y, NuTo::Node::eDof::DISPLACEMENTS);
                    for (int dofID : overlapDofIDs_Y)
                    {
                        local2GlobalMapping[dofID] = newID + internDofCounter + offset*dim*(numElementsPerSubDomain_Y-1);
                        ++internDofCounter;
                    }
                    ++offset;
                }

                if (overlappingNodeCount_Corner > 0)
                {
                    newID = numActiveDofsTillSubDomain_inclusive[rank-numSubDomains_Y-1] - (overlappingNodeCount_X+overlappingNodeCount_Corner+overlappingNodeCount_Y)*dim;
                    int nodeID_Corner = overlapNodeIDs_Corner[0];
                    std::vector<int> overlapDofIDs_Corner = structure.NodeGetDofIds(nodeID_Corner, NuTo::Node::eDof::DISPLACEMENTS);
                    for (int dofID : overlapDofIDs_Corner)
                    {
//                        if (overlappingNodeCount_Y > 0)
//                            local2GlobalMapping[dofID] = newID - (overlappingNodeCount_Y-1)*dim + internDofCounter;
////                            local2GlobalMapping[dofID] = newID - (overlappingNodeCount_Y)*dim + internDofCounter;
//                        else
                            local2GlobalMapping[dofID] = newID + internDofCounter;
                        ++internDofCounter;
                    }
                }

            }
        }
    }

    delete[] numActiveDofsTillSubDomain_inclusive;

    //DEFINE MAPS FOR LOCAL OWNED DOFS AND OVERLAPPING DOFS
    Epetra_Map owningMap(-1, ownedNumActiveDOFs, 0, Comm);
    Epetra_Map overlappingMap(-1, localNumActiveDOFs, local2GlobalMapping, 0, Comm);
//    Epetra_Map owningMap = TrilinosUtils::createLinearMap(Comm, ownedNumActiveDOFs, 0);
//    Epetra_Map overlappingMap = TrilinosUtils::createSpecificMap(Comm, localNumActiveDOFs, 0, local2GlobalMapping);
    Epetra_Export exporter(overlappingMap, owningMap);

    //CREATE CORRESPONDING GRAPHS (FOR MATRIX STRUCTURE)
    int maxNonZeros = 18;   //got by interpolation type (order), e.g. QUAD2 + EQUIDISTANT1 => 18
//    Epetra_CrsGraph owningGraph(Epetra_DataAccess::Copy, owningMap, 0, false);
    Epetra_CrsGraph owningGraph(Epetra_DataAccess::Copy, owningMap, maxNonZeros, false);
//    Epetra_CrsGraph overlappingGraph(Epetra_DataAccess::Copy, overlappingMap, 0, false);
    Epetra_CrsGraph overlappingGraph(Epetra_DataAccess::Copy, overlappingMap, maxNonZeros, false);

    int errCode = -1;
//    for (int i = 0; i < localNumActiveDOFs; ++i)
//    {
//        for (int j = 0; j < localNumActiveDOFs; ++j)
//        {

//            errCode = overlappingGraph.InsertGlobalIndices(local2GlobalMapping[i], 1, &local2GlobalMapping[j]);
//            checkErrorCode(errCode, "i = " + std::to_string(i) + ", j = " + std::to_string(j));

//        }
//    }

    for (int k=0; k<hessian0_eigen.outerSize(); ++k)
    {
      for (Eigen::SparseMatrix<double>::InnerIterator it(hessian0_eigen,k); it; ++it)
      {
//        std::cout << it.value() << std::endl;
//        std::cout << it.row() << std::endl;   // row index
//        std::cout << it.col() << std::endl;   // col index (here it is equal to k)
//        std::cout << it.index() << std::endl; // inner index, here it is equal to it.row()
        errCode = overlappingGraph.InsertGlobalIndices(local2GlobalMapping[it.row()], 1, &local2GlobalMapping[it.col()]);
        checkErrorCode(errCode, "k = " + std::to_string(k) + ", it = " + std::to_string(it.row()));
      }
    }

    delete[] local2GlobalMapping;

    errCode = overlappingGraph.FillComplete();
    checkErrorCode(errCode, "overlappingGraph.FillComplete()");
    errCode = owningGraph.Export(overlappingGraph, exporter, Insert); //inter-process communication of matrix structure
    checkErrorCode(errCode, "owningGraph.Export(...)");
    errCode = owningGraph.FillComplete();
    checkErrorCode(errCode, "owningGraph.FillComplete()");

    //CREATE (GLOBAL) MATRIX AND VECTOR WITH KNOWN GRAPH STRUCTURE
    Epetra_CrsMatrix globalMatrix(Epetra_DataAccess::Copy, owningGraph);
//    globalMatrix.FillComplete();
//    globalMatrix.FillComplete(false);
    Epetra_Vector globalRhsVector(owningMap);
    globalRhsVector.PutScalar(0.0);

    //CONVERT EIGEN::MATRIX AND EIGEN::VECTOR TO EPETRA FORMAT
    ConversionTools converter;
    Eigen::SparseMatrix<double, Eigen::RowMajor> hess(hessian0_eigen);
//    Epetra_CrsMatrix localMatrix = converter.convertEigen2EpetraCrsMatrix(hessian0_eigen, overlappingGraph);
    Epetra_CrsMatrix localMatrix = converter.convertEigen2EpetraCrsMatrix(hess, overlappingGraph);
    Epetra_Vector localRhsVector = converter.convertEigen2EpetraVector(residual_eigen, overlappingMap);
    if (!loadOnNode)
        localRhsVector.Scale(1/(double(numSubDomains_Y*numElementsPerSubDomain_Y+1)));

    //globalMatrix.PutScalar(0.0);
    errCode = globalMatrix.Export(localMatrix, exporter, Add);    //inter-process communication of matrix entries
    checkErrorCode(errCode, "globalMatrix.Export(...)");
//    globalRhsVector.Export(localRhsVector, exporter, Add);    //inter-process communication of vector entries
    errCode = globalRhsVector.Export(localRhsVector, exporter, Insert);
    checkErrorCode(errCode, "globalRhsVector.Export(...)");

#ifdef SHOW_INTERMEDIATE_RESULTS
    //PRINT SOME INTERMEDIATE RESULTS
    std::ostream& os = std::cout;

    os << "Hessian0 by NuTo on proc " << rank << ":\n----------------------\n" << hessian0_eigen.toDense() << "\n\n";
    os << "residual by NuTo on proc " << rank << ":\n----------------------\n" << residual_eigen << "\n\n";
    os << "Owning Map:\n-------------\n";
    owningMap.Print(os);
    os << "Overlap Map:\n-------------\n";
    overlappingMap.Print(os);
    os << "Local matrix on proc " << rank << ":\n-----------------\n";
    localMatrix.Print(os);
    os << "Global matrix:\n-------------\n";
    globalMatrix.Print(os);
    os << "Local RHS on proc " << rank << ":\n-----------------\n";
    localRhsVector.Print(os);
    os << "Global RHS:\n-------------\n";
    globalRhsVector.Print(os);
#endif

    // --------------------------------
    // | SOLVE GLOBAL SYSTEM PARALLEL |
    // --------------------------------
//    Epetra_Vector lhs(owningMap);
//    lhs.PutScalar(0.0);
    Epetra_Vector lhs(globalRhsVector);
    Epetra_MultiVector sol = solveSystem(globalMatrix, lhs, globalRhsVector, true, true);
    sol.Scale(-1);
    //save solution
    if (rSaveSolution)
    {
        EpetraExt::MultiVectorToMatlabFile((resultDirectory + "/solution.mat").c_str(), sol);
        EpetraExt::RowMatrixToMatlabFile((resultDirectory + "/globalMatrix.mat").c_str(), globalMatrix);
        EpetraExt::VectorToMatlabFile((resultDirectory + "/globalRHS.mat").c_str(), globalRhsVector);
    }

#ifdef SHOW_SOLUTION
    std::cout << "Solution:\n------------\n";
    sol.Print(std::cout);
#else
    if (rank == 0)
        std::cout << "Number of DOFs: " << overlappingMap.MaxAllGID() + 1 << std::endl;
#endif

    //solve system serial (for comparison only)
//    if (numProc == 1)
//    {
//        Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
//        solver.compute(hessian0_eigen);
//        Eigen::VectorXd displ = solver.solve(residual_eigen);
//        displ *= -1;
//        std::cout << rank << ":\n" << displ << std::endl;
//        //save direct solution
//        std::ofstream file(resultDirectory + "/solution_direct.mat");
//        file << displ;
//        file.close();
//    }
}

void run_Test_2D(int argc, char** argv)
{
    Epetra_MpiComm Comm(MPI_COMM_WORLD);
    double totalLength_X = 1000.;
    double totalLength_Y = 150.;
    double thickness = 10.;
    double totalForce = 1.e5;
    int numSubDomains_X = 0;
    int numSubDomains_Y = 0;
    int numElementsPerSubDomain_X = 1;  //200
    int numElementsPerSubDomain_Y = 1;  //100
    bool saveSolution = true;

    if (argc >= 2 && (strcmp(argv[1], "help")== 0))
    {
        if (Comm.MyPID() == 0)
        {
            std::cout << "Abbreviations:\n--------------"
                      << "\n- numElX -> number of elements per subdomain in direction X\n- numElY -> ... in direction Y"
                      << "\n- numSubX -> number of subdomains in direction X\n- numSubY -> ... in direction Y"
                      << "\n- lenX -> length in direction X\n- lenY -> ... Y"
                      << "\n- thick -> thickness of plane"
                      << "\n\nRemark: numSubX * numSubY = <number of processes>"
                      << "\n\nDefault values:\n---------------"
                      << "\n- numElX = numElY = 1"
                      << "\n- lenX = 1000"
                      << "\n- lenY = 150"
                      << "\n- thick = 10"
                      << "\n\n\nUSAGE: 2D_Test help --> show this help"
                      << "\nUSAGE: 2D_Test [numElX] [numElY] [numSubX] [numSubY] [lenX] [lenY] [thick]\n\n";
        }
    }
    else
    {

        if (argc >= 3)
        {
            numElementsPerSubDomain_X = atoi(argv[1]);
            numElementsPerSubDomain_Y = atoi(argv[2]);
        }
        if (argc >= 5)
        {
            numSubDomains_X = atoi(argv[3]);
            numSubDomains_Y = atoi(argv[4]);
        }
        if (argc >= 7)
        {
            totalLength_X = atoi(argv[5]);
            totalLength_Y = atoi(argv[6]);
        }
        if (argc >= 8)
            thickness = atoi(argv[7]);

        run_Test_2D(Comm, totalLength_X, totalLength_Y, thickness, numSubDomains_X, numSubDomains_Y, numElementsPerSubDomain_X, numElementsPerSubDomain_Y, totalForce, saveSolution);
    }
}


int main(int argc, char** argv)
{
    Teuchos::GlobalMPISession mpiSession(&argc, &argv);

//    run_Test_2D_GivenMesh();

    run_Test_2D(argc, argv);

    return 0;
}
