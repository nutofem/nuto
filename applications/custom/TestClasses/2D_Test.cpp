#include <iostream>
#include <mpi.h>

#include <eigen3/Eigen/Core>
#include <boost/filesystem.hpp>

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


//#define SHOW_INTERMEDIATE_RESULTS
#define SHOW_SOLUTION


void checkDirectories(const std::string rLogDirectory, const std::string rResultDirectory)
{
    if (!boost::filesystem::exists(rLogDirectory))
        boost::filesystem::create_directory(rLogDirectory);

    if (!boost::filesystem::exists(rResultDirectory))
        boost::filesystem::create_directory(rResultDirectory);
}


Epetra_MultiVector solveSystem(Epetra_CrsMatrix rA, Epetra_MultiVector rLhs, Epetra_MultiVector rRhs, bool iterative = true)
{
    Epetra_LinearProblem problem(&rA, &rLhs, &rRhs);

    if (iterative)
    {
        AztecOO Solver(problem);
        Solver.Iterate(1000,1e-8);
    }
    else
    {
        Amesos Factory;
        Amesos_BaseSolver* solver;
//        std::string solverType = "Mumps";
        std::string solverType = "Klu";
        solver = Factory.Create(solverType, problem);
        solver->Solve();
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
    Epetra_CrsGraph overlappingGraph(Copy, overlappingMap, 0);
    Epetra_CrsGraph owningGraph(Copy, owningMap, 0);

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
    Epetra_CrsMatrix globalMatrix(Copy, owningGraph);
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


void run_Test_2D()
{
    Epetra_MpiComm Comm(MPI_COMM_WORLD);
    int rank = Comm.MyPID();
    int numProc = Comm.NumProc();

    //CREATE DIRECTORIES FOR LOGS AND RESULTS
    const std::string logDirectory = "2D_Test_Logs";
    const std::string resultDirectory = "2D_Test_Results";
    checkDirectories(logDirectory, resultDirectory);


    // one-dimensional problem
    const int dim = 2;

    // --------------------
    // | DOMAIN VARIABLES |
    // --------------------
    int numSubDomains_X = 0;
    int numSubDomains_Y = 0;
    separateDomain(numProc, numSubDomains_X, numSubDomains_Y);
    std::cout << "SEPARATE DOMAIN: " << numProc << ", " << numSubDomains_X << ", " << numSubDomains_Y << std::endl;

    std::vector<int> subDomains_X;
    std::vector<int> subDomains_Y;
    generateSubDomainIdentifiers(numSubDomains_X, numSubDomains_Y, subDomains_X, subDomains_Y);

    int numSubDomains = numSubDomains_X * numSubDomains_Y;
    int* numElementsPerSubDomain_X = new int[numSubDomains_X];
    int* numElementsPerSubDomain_Y = new int[numSubDomains_Y];
    int* numActiveDofsPerSubDomain = new int[numSubDomains];
    int* numActiveDofsTillSubDomain_inclusive = new int[numSubDomains];
    int numElementsTotal_X = 0;
    int numElementsTotal_Y = 0;
    bool setTotalLength = true;
    double* lengthPerSubDomain_X = new double[numSubDomains_X];
    double lengthTotal_X = (setTotalLength ? 10. : 0.);
    double* lengthPerSubDomain_Y = new double[numSubDomains_Y];
    double lengthTotal_Y = (setTotalLength ? 5. : 0.);

    for (int i = 0; i < numSubDomains_X; ++i)
    {
        numElementsPerSubDomain_X[i] = 1;
        numElementsTotal_X += numElementsPerSubDomain_X[i];

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
        numElementsPerSubDomain_Y[i] = 1;
        numElementsTotal_Y += numElementsPerSubDomain_Y[i];

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
    bool loadOnNode = true;


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
                for (int jx = 0; jx < numElementsPerSubDomain_X[currSubDomain_X]+1; ++jx)
                {
                    nodeCoords(0) = currLength_X + jx*lengthPerSubDomain_X[currSubDomain_X]/numElementsPerSubDomain_X[currSubDomain_X];
                    for (int jy = 0; jy < numElementsPerSubDomain_Y[currSubDomain_Y]+1; ++jy)
                    {
                        nodeCoords(1) = currLength_Y + jy*lengthPerSubDomain_Y[currSubDomain_Y]/numElementsPerSubDomain_Y[currSubDomain_Y];
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
                for (int jx = 0; jx < numElementsPerSubDomain_X[currSubDomain_X]; ++jx)
                {
                    for (int jy = 0; jy < numElementsPerSubDomain_Y[currSubDomain_Y]; ++jy)
                    {
                        elementNodes[0] = currNode;
                        elementNodes[1] = currNode + 1;
                        elementNodes[2] = currNode + (numElementsPerSubDomain_Y[currSubDomain_Y]+1) + 1;
                        elementNodes[3] = currNode + (numElementsPerSubDomain_Y[currSubDomain_Y]+1);
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

    //FIX LEFT BOUNDARY
    structure.ConstraintLinearSetDisplacementNodeGroup(nodesBCLeft, directionX, fixedDisplacement);
    structure.ConstraintLinearSetDisplacementNodeGroup(nodesBCLeft, directionY, fixedDisplacement);

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

    Eigen::SparseMatrix<double> hessian0_eigen = hessian0.JJ.ExportToEigenSparseMatrix();
    Eigen::MatrixXd residual_eigen = residual.J.Export();

    int ownedNumActiveDOFs = -1;
    if (currSubDomain_Y == 0)
    {
        ownedNumActiveDOFs = numElementsPerSubDomain_X[currSubDomain_X]*(numElementsPerSubDomain_Y[currSubDomain_Y]+1)*dim;
    }
    else
    {
        ownedNumActiveDOFs = numElementsPerSubDomain_X[currSubDomain_X]*numElementsPerSubDomain_Y[currSubDomain_Y]*dim;
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

    std::vector<int> leftNodeIDs = structure.GroupGetMemberIds(nodesBCLeft);
    int leftNodeCount = leftNodeIDs.size();
    std::vector<int> overlapNodeIDs_Y;
    overlapNodeIDs_Corner = vectorIntersection(overlapNodeIDs_X, overlapNodeIDs_Y_total);
    overlapNodeIDs_Y = vectorDifference(overlapNodeIDs_Y_total, leftNodeIDs);
    overlapNodeIDs_Y = vectorDifference(overlapNodeIDs_Y, overlapNodeIDs_Corner);
    overlapNodeIDs_X = vectorDifference(overlapNodeIDs_X, overlapNodeIDs_Corner);
    int overlappingNodeCount_X = overlapNodeIDs_X.size();
    int overlappingNodeCount_Y = overlapNodeIDs_Y.size();
    int overlappingNodeCount_Corner = overlapNodeIDs_Corner.size();

    int subDomainCounter = 0;
    for (int k = 0; k < numSubDomains_Y; ++k)
    {
        numActiveDofsPerSubDomain[subDomainCounter] = (numElementsPerSubDomain_X[0]) * (numElementsPerSubDomain_Y[k] + 1)*dim;
        if (subDomainCounter == 0)
            numActiveDofsTillSubDomain_inclusive[subDomainCounter] = numActiveDofsPerSubDomain[subDomainCounter];
        else
        {
            if ((subDomains_X[rank] == 0) && (numSubDomains_X > 1))
                numActiveDofsTillSubDomain_inclusive[subDomainCounter] = numActiveDofsTillSubDomain_inclusive[subDomainCounter-1]
                                                                        + numActiveDofsPerSubDomain[subDomainCounter] - dim*(overlappingNodeCount_Y+overlappingNodeCount_Corner);
            else
                numActiveDofsTillSubDomain_inclusive[subDomainCounter] = numActiveDofsTillSubDomain_inclusive[subDomainCounter-1]
                                                                        + numActiveDofsPerSubDomain[subDomainCounter] - dim*(overlappingNodeCount_Y);

        }
        ++subDomainCounter;
    }


    for (int j = 1; j < numSubDomains_X; ++j)
    {
        for (int k = 0; k < numSubDomains_Y; ++k)
        {
            numActiveDofsPerSubDomain[subDomainCounter] = (numElementsPerSubDomain_X[j] + 1) * (numElementsPerSubDomain_Y[k] + 1)*dim;
//            if (subDomains_X[rank] == 0)
            if ((subDomains_X[rank] == 0) && (numSubDomains_X > 1))
            {
                if (subDomains_Y[subDomainCounter] == 0)
                    numActiveDofsTillSubDomain_inclusive[subDomainCounter] = numActiveDofsTillSubDomain_inclusive[subDomainCounter-1]
                                                                        + numActiveDofsPerSubDomain[subDomainCounter] - dim*(overlappingNodeCount_X+overlappingNodeCount_Corner);
                else
                    numActiveDofsTillSubDomain_inclusive[subDomainCounter] = numActiveDofsTillSubDomain_inclusive[subDomainCounter-1]
                                                                        + numActiveDofsPerSubDomain[subDomainCounter] - dim*(overlappingNodeCount_X+overlappingNodeCount_Corner+overlappingNodeCount_Y+1);
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
    if (numProc > 1)
    {
        if ((subDomains_X[rank] == 0) && (numSubDomains_X > 1))
        {
            if (subDomains_Y[rank] == 0)
            {
                for (int nodeID_X : overlapNodeIDs_X)
                {
                    std::vector<int> overlapDofIDs_X = structure.NodeGetDofIds(nodeID_X, NuTo::Node::eDof::DISPLACEMENTS);
                    for (int dofID : overlapDofIDs_X)
                    {
                        local2GlobalMapping[dofID] = numActiveDofsTillSubDomain_inclusive[rank] - (overlappingNodeCount_X+overlappingNodeCount_Corner+overlappingNodeCount_Y)*dim + internDofCounter;
                        ++internDofCounter;
                    }
                }

                if (overlappingNodeCount_Corner > 0)
                {
                    int nodeID_Corner = overlapNodeIDs_Corner[0];
                    std::vector<int> overlapDofIDs_Corner = structure.NodeGetDofIds(nodeID_Corner, NuTo::Node::eDof::DISPLACEMENTS);
                    for (int dofID : overlapDofIDs_Corner)
                    {
                        local2GlobalMapping[dofID] = numActiveDofsTillSubDomain_inclusive[rank] - (overlappingNodeCount_X+overlappingNodeCount_Corner+overlappingNodeCount_Y)*dim + internDofCounter;
                        ++internDofCounter;
                    }
                }

                for (int nodeID_Y : overlapNodeIDs_Y)
                {
                    std::vector<int> overlapDofIDs_Y = structure.NodeGetDofIds(nodeID_Y, NuTo::Node::eDof::DISPLACEMENTS);
                    for (int dofID : overlapDofIDs_Y)
                    {
                        local2GlobalMapping[dofID] = numActiveDofsTillSubDomain_inclusive[rank] - (overlappingNodeCount_X+overlappingNodeCount_Corner+overlappingNodeCount_Y)*dim + internDofCounter;
                        ++internDofCounter;
                    }
                }
            }
            else if (subDomains_Y[rank] == 1)
            {
                for (int nodeID_X : overlapNodeIDs_X)
                {
                    std::vector<int> overlapDofIDs_X = structure.NodeGetDofIds(nodeID_X, NuTo::Node::eDof::DISPLACEMENTS);
                    for (int dofID : overlapDofIDs_X)
                    {
                        local2GlobalMapping[dofID] = numActiveDofsTillSubDomain_inclusive[rank] - (overlappingNodeCount_X)*dim + internDofCounter;
                        ++internDofCounter;
                    }
                }

                if (overlappingNodeCount_Corner > 0)
                {
                    int nodeID_Corner = overlapNodeIDs_Corner[0];
                    std::vector<int> overlapDofIDs_Corner = structure.NodeGetDofIds(nodeID_Corner, NuTo::Node::eDof::DISPLACEMENTS);
                    for (int dofID : overlapDofIDs_Corner)
                    {
                        local2GlobalMapping[dofID] = numActiveDofsTillSubDomain_inclusive[rank-1] - (overlappingNodeCount_X+overlappingNodeCount_Corner+overlappingNodeCount_Y)*dim + internDofCounter;
                        ++internDofCounter;
                    }
                }

                for (int nodeID_Y : overlapNodeIDs_Y)
                {
                    std::vector<int> overlapDofIDs_Y = structure.NodeGetDofIds(nodeID_Y, NuTo::Node::eDof::DISPLACEMENTS);
                    for (int dofID : overlapDofIDs_Y)
                    {
                        local2GlobalMapping[dofID] = numActiveDofsTillSubDomain_inclusive[rank-1] - (overlappingNodeCount_X+overlappingNodeCount_Corner+overlappingNodeCount_Y)*dim + internDofCounter;
                        ++internDofCounter;
                    }
                }
            }
            else
            {
                for (int nodeID_X : overlapNodeIDs_X)
                {
                    std::vector<int> overlapDofIDs_X = structure.NodeGetDofIds(nodeID_X, NuTo::Node::eDof::DISPLACEMENTS);
                    for (int dofID : overlapDofIDs_X)
                    {
                        local2GlobalMapping[dofID] = numActiveDofsTillSubDomain_inclusive[rank] - (overlappingNodeCount_X)*dim + internDofCounter;
                        ++internDofCounter;
                    }
                }

                if (overlappingNodeCount_Corner > 0)
                {
                    int nodeID_Corner = overlapNodeIDs_Corner[0];
                    std::vector<int> overlapDofIDs_Corner = structure.NodeGetDofIds(nodeID_Corner, NuTo::Node::eDof::DISPLACEMENTS);
                    for (int dofID : overlapDofIDs_Corner)
                    {
                        local2GlobalMapping[dofID] = numActiveDofsTillSubDomain_inclusive[rank-1] - (overlappingNodeCount_X+overlappingNodeCount_Corner+overlappingNodeCount_Y)*dim + internDofCounter;
                        ++internDofCounter;
                    }
                }

                for (int nodeID_Y : overlapNodeIDs_Y)
                {
                    std::vector<int> overlapDofIDs_Y = structure.NodeGetDofIds(nodeID_Y, NuTo::Node::eDof::DISPLACEMENTS);
                    for (int dofID : overlapDofIDs_Y)
                    {
                        local2GlobalMapping[dofID] = numActiveDofsTillSubDomain_inclusive[rank-1] - (overlappingNodeCount_X+overlappingNodeCount_Corner+overlappingNodeCount_Y)*dim + internDofCounter;
                        ++internDofCounter;
                    }
                }
            }
        }
        else if ((subDomains_X[rank] == 0) && (numSubDomains_X == 1))
        {
            if (subDomains_Y[rank] == 0)
            {
                for (int nodeID_Y : overlapNodeIDs_Y)
                {
                    std::vector<int> overlapDofIDs_Y = structure.NodeGetDofIds(nodeID_Y, NuTo::Node::eDof::DISPLACEMENTS);
                    for (int dofID : overlapDofIDs_Y)
                    {
                        local2GlobalMapping[dofID] = numActiveDofsTillSubDomain_inclusive[rank] - (overlappingNodeCount_Y)*dim + internDofCounter;
                        ++internDofCounter;
                    }
                }
            }
            else if (subDomains_Y[rank] == 1)
            {
                for (int nodeID_Y : overlapNodeIDs_Y)
                {
                    std::vector<int> overlapDofIDs_Y = structure.NodeGetDofIds(nodeID_Y, NuTo::Node::eDof::DISPLACEMENTS);
                    for (int dofID : overlapDofIDs_Y)
                    {
                        local2GlobalMapping[dofID] = numActiveDofsTillSubDomain_inclusive[rank-1] - (overlappingNodeCount_Y)*dim + internDofCounter;
                        ++internDofCounter;
                    }
                }
            }
            else
            {
                for (int nodeID_Y : overlapNodeIDs_Y)
                {
                    std::vector<int> overlapDofIDs_Y = structure.NodeGetDofIds(nodeID_Y, NuTo::Node::eDof::DISPLACEMENTS);
                    for (int dofID : overlapDofIDs_Y)
                    {
                        local2GlobalMapping[dofID] = numActiveDofsTillSubDomain_inclusive[rank-1] - (overlappingNodeCount_Y)*dim + internDofCounter;
                        ++internDofCounter;
                    }
                }
            }
        }
        else if (subDomains_X[rank] == 1)
        {
            if (subDomains_Y[rank] == 0)
            {
                for (int nodeID_X : overlapNodeIDs_X)
                {
                    std::vector<int> overlapDofIDs_X = structure.NodeGetDofIds(nodeID_X, NuTo::Node::eDof::DISPLACEMENTS);
                    for (int dofID : overlapDofIDs_X)
                    {
                        if (overlappingNodeCount_Y > 0)
                            local2GlobalMapping[dofID] = numActiveDofsTillSubDomain_inclusive[rank-numSubDomains_Y] - (overlappingNodeCount_X+overlappingNodeCount_Corner+overlappingNodeCount_Y-1)*dim + internDofCounter;
                        else
                            local2GlobalMapping[dofID] = numActiveDofsTillSubDomain_inclusive[rank-numSubDomains_Y] - (overlappingNodeCount_X+overlappingNodeCount_Corner)*dim + internDofCounter;

                        ++internDofCounter;
                    }
                }

                if (overlappingNodeCount_Corner > 0)
                {
                    int nodeID_Corner = overlapNodeIDs_Corner[0];
                    std::vector<int> overlapDofIDs_Corner = structure.NodeGetDofIds(nodeID_Corner, NuTo::Node::eDof::DISPLACEMENTS);
                    for (int dofID : overlapDofIDs_Corner)
                    {
                        if (overlappingNodeCount_Y > 0)
                            local2GlobalMapping[dofID] = numActiveDofsTillSubDomain_inclusive[rank-numSubDomains_Y] - (overlappingNodeCount_X+overlappingNodeCount_Corner+overlappingNodeCount_Y-1)*dim + internDofCounter;
                        else
                            local2GlobalMapping[dofID] = numActiveDofsTillSubDomain_inclusive[rank-numSubDomains_Y] - (overlappingNodeCount_X+overlappingNodeCount_Corner)*dim + internDofCounter;

                        ++internDofCounter;
                    }
                }

                for (int nodeID_Y : overlapNodeIDs_Y)
                {
                    std::vector<int> overlapDofIDs_Y = structure.NodeGetDofIds(nodeID_Y, NuTo::Node::eDof::DISPLACEMENTS);
                    for (int dofID : overlapDofIDs_Y)
                    {
                        local2GlobalMapping[dofID] = numActiveDofsTillSubDomain_inclusive[rank] - (overlappingNodeCount_X+overlappingNodeCount_Corner+overlappingNodeCount_Y)*dim + internDofCounter;
                        ++internDofCounter;
                    }
                }
            }
            else if (subDomains_Y[rank] == 1)
            {
                for (int nodeID_X : overlapNodeIDs_X)
                {
                    std::vector<int> overlapDofIDs_X = structure.NodeGetDofIds(nodeID_X, NuTo::Node::eDof::DISPLACEMENTS);
                    for (int dofID : overlapDofIDs_X)
                    {
                        local2GlobalMapping[dofID] = numActiveDofsTillSubDomain_inclusive[rank-numSubDomains_Y] - (overlappingNodeCount_X)*dim + internDofCounter;
                        ++internDofCounter;
                    }
                }

                if (overlappingNodeCount_Corner > 0)
                {
                    int nodeID_Corner = overlapNodeIDs_Corner[0];
                    std::vector<int> overlapDofIDs_Corner = structure.NodeGetDofIds(nodeID_Corner, NuTo::Node::eDof::DISPLACEMENTS);
                    for (int dofID : overlapDofIDs_Corner)
                    {
                        if (overlappingNodeCount_Y > 0)
                            local2GlobalMapping[dofID] = numActiveDofsTillSubDomain_inclusive[rank-numSubDomains_Y-1] - (overlappingNodeCount_X+overlappingNodeCount_Corner+overlappingNodeCount_Y-1)*dim + internDofCounter;
                        else
                            local2GlobalMapping[dofID] = numActiveDofsTillSubDomain_inclusive[rank-numSubDomains_Y-1] - (overlappingNodeCount_X+overlappingNodeCount_Corner)*dim + internDofCounter;

                        ++internDofCounter;
                    }
                }

                for (int nodeID_Y : overlapNodeIDs_Y)
                {
                    std::vector<int> overlapDofIDs_Y = structure.NodeGetDofIds(nodeID_Y, NuTo::Node::eDof::DISPLACEMENTS);
                    for (int dofID : overlapDofIDs_Y)
                    {
                        local2GlobalMapping[dofID] = numActiveDofsTillSubDomain_inclusive[rank-1] - (overlappingNodeCount_X+overlappingNodeCount_Corner+overlappingNodeCount_Y)*dim + internDofCounter;
                        ++internDofCounter;
                    }
                }
            }
            else
            {
                for (int nodeID_X : overlapNodeIDs_X)
                {
                    std::vector<int> overlapDofIDs_X = structure.NodeGetDofIds(nodeID_X, NuTo::Node::eDof::DISPLACEMENTS);
                    for (int dofID : overlapDofIDs_X)
                    {
                        local2GlobalMapping[dofID] = numActiveDofsTillSubDomain_inclusive[rank-numSubDomains_Y] - (overlappingNodeCount_X)*dim + internDofCounter;
                        ++internDofCounter;
                    }
                }

                if (overlappingNodeCount_Corner > 0)
                {
                    int nodeID_Corner = overlapNodeIDs_Corner[0];
                    std::vector<int> overlapDofIDs_Corner = structure.NodeGetDofIds(nodeID_Corner, NuTo::Node::eDof::DISPLACEMENTS);
                    for (int dofID : overlapDofIDs_Corner)
                    {
                        if (overlappingNodeCount_Y > 0)
                            local2GlobalMapping[dofID] = numActiveDofsTillSubDomain_inclusive[rank-numSubDomains_Y-1] - (overlappingNodeCount_X+overlappingNodeCount_Corner+overlappingNodeCount_Y-1)*dim + internDofCounter;
                        else
                            local2GlobalMapping[dofID] = numActiveDofsTillSubDomain_inclusive[rank-numSubDomains_Y-1] - (overlappingNodeCount_X+overlappingNodeCount_Corner)*dim + internDofCounter;
                        ++internDofCounter;
                    }
                }

                for (int nodeID_Y : overlapNodeIDs_Y)
                {
                    std::vector<int> overlapDofIDs_Y = structure.NodeGetDofIds(nodeID_Y, NuTo::Node::eDof::DISPLACEMENTS);
                    for (int dofID : overlapDofIDs_Y)
                    {
                        local2GlobalMapping[dofID] = numActiveDofsTillSubDomain_inclusive[rank-1] - (overlappingNodeCount_X+overlappingNodeCount_Corner+overlappingNodeCount_Y)*dim + internDofCounter;
                        ++internDofCounter;
                    }
                }
            }
        }
        else
        {
            if (subDomains_Y[rank] == 0)
            {
                for (int nodeID_X : overlapNodeIDs_X)
                {
                    std::vector<int> overlapDofIDs_X = structure.NodeGetDofIds(nodeID_X, NuTo::Node::eDof::DISPLACEMENTS);
                    for (int dofID : overlapDofIDs_X)
                    {
                        local2GlobalMapping[dofID] = numActiveDofsTillSubDomain_inclusive[rank-numSubDomains_Y] - (overlappingNodeCount_X+overlappingNodeCount_Corner+overlappingNodeCount_Y)*dim + internDofCounter;
                        ++internDofCounter;
                    }
                }

                if (overlappingNodeCount_Corner > 0)
                {
                    int nodeID_Corner = overlapNodeIDs_Corner[0];
                    std::vector<int> overlapDofIDs_Corner = structure.NodeGetDofIds(nodeID_Corner, NuTo::Node::eDof::DISPLACEMENTS);
                    for (int dofID : overlapDofIDs_Corner)
                    {
                        local2GlobalMapping[dofID] = numActiveDofsTillSubDomain_inclusive[rank-numSubDomains_Y] - (overlappingNodeCount_X+overlappingNodeCount_Corner+overlappingNodeCount_Y)*dim + internDofCounter;
                        ++internDofCounter;
                    }
                }

                for (int nodeID_Y : overlapNodeIDs_Y)
                {
                    std::vector<int> overlapDofIDs_Y = structure.NodeGetDofIds(nodeID_Y, NuTo::Node::eDof::DISPLACEMENTS);
                    for (int dofID : overlapDofIDs_Y)
                    {
                        local2GlobalMapping[dofID] = numActiveDofsTillSubDomain_inclusive[rank] - (overlappingNodeCount_X+overlappingNodeCount_Corner+overlappingNodeCount_Y)*dim + internDofCounter;
                        ++internDofCounter;
                    }
                }
            }
            else if (subDomains_Y[rank] == 1)
            {
                for (int nodeID_X : overlapNodeIDs_X)
                {
                    std::vector<int> overlapDofIDs_X = structure.NodeGetDofIds(nodeID_X, NuTo::Node::eDof::DISPLACEMENTS);
                    for (int dofID : overlapDofIDs_X)
                    {
                        local2GlobalMapping[dofID] = numActiveDofsTillSubDomain_inclusive[rank-numSubDomains_Y] - (overlappingNodeCount_X)*dim + internDofCounter;
                        ++internDofCounter;
                    }
                }

                if (overlappingNodeCount_Corner > 0)
                {
                    int nodeID_Corner = overlapNodeIDs_Corner[0];
                    std::vector<int> overlapDofIDs_Corner = structure.NodeGetDofIds(nodeID_Corner, NuTo::Node::eDof::DISPLACEMENTS);
                    for (int dofID : overlapDofIDs_Corner)
                    {
                        local2GlobalMapping[dofID] = numActiveDofsTillSubDomain_inclusive[rank-numSubDomains_Y-1] - (overlappingNodeCount_X+overlappingNodeCount_Corner+overlappingNodeCount_Y)*dim + internDofCounter;
                        ++internDofCounter;
                    }
                }

                for (int nodeID_Y : overlapNodeIDs_Y)
                {
                    std::vector<int> overlapDofIDs_Y = structure.NodeGetDofIds(nodeID_Y, NuTo::Node::eDof::DISPLACEMENTS);
                    for (int dofID : overlapDofIDs_Y)
                    {
                        local2GlobalMapping[dofID] = numActiveDofsTillSubDomain_inclusive[rank-1] - (overlappingNodeCount_X+overlappingNodeCount_Corner+overlappingNodeCount_Y)*dim + internDofCounter;
                        ++internDofCounter;
                    }
                }
            }
            else
            {
                for (int nodeID_X : overlapNodeIDs_X)
                {
                    std::vector<int> overlapDofIDs_X = structure.NodeGetDofIds(nodeID_X, NuTo::Node::eDof::DISPLACEMENTS);
                    for (int dofID : overlapDofIDs_X)
                    {
                        local2GlobalMapping[dofID] = numActiveDofsTillSubDomain_inclusive[rank-numSubDomains_Y] - (overlappingNodeCount_X)*dim + internDofCounter;
                        ++internDofCounter;
                    }
                }

                if (overlappingNodeCount_Corner > 0)
                {
                    int nodeID_Corner = overlapNodeIDs_Corner[0];
                    std::vector<int> overlapDofIDs_Corner = structure.NodeGetDofIds(nodeID_Corner, NuTo::Node::eDof::DISPLACEMENTS);
                    for (int dofID : overlapDofIDs_Corner)
                    {
                        local2GlobalMapping[dofID] = numActiveDofsTillSubDomain_inclusive[rank-numSubDomains_Y-1] - (overlappingNodeCount_X+overlappingNodeCount_Corner+overlappingNodeCount_Y)*dim + internDofCounter;
                        ++internDofCounter;
                    }
                }

                for (int nodeID_Y : overlapNodeIDs_Y)
                {
                    std::vector<int> overlapDofIDs_Y = structure.NodeGetDofIds(nodeID_Y, NuTo::Node::eDof::DISPLACEMENTS);
                    for (int dofID : overlapDofIDs_Y)
                    {
                        local2GlobalMapping[dofID] = numActiveDofsTillSubDomain_inclusive[rank-1] - (overlappingNodeCount_X+overlappingNodeCount_Corner+overlappingNodeCount_Y)*dim + internDofCounter;
                        ++internDofCounter;
                    }
                }
            }
        }
    }

    //DEFINE MAPS FOR LOCAL OWNED DOFS AND OVERLAPPING DOFS
    Epetra_Map owningMap(-1, ownedNumActiveDOFs, 0, Comm);
    Epetra_Map overlappingMap(-1, localNumActiveDOFs, local2GlobalMapping, 0, Comm);
    Epetra_Export exporter(overlappingMap, owningMap);

    //CREATE CORRESPONDING GRAPHS (FOR MATRIX STRUCTURE)
    Epetra_CrsGraph overlappingGraph(Copy, overlappingMap, 0);
    Epetra_CrsGraph owningGraph(Copy, owningMap, 0);

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
    Epetra_CrsMatrix globalMatrix(Copy, owningGraph);
    globalMatrix.FillComplete();
    Epetra_Vector globalRhsVector(owningMap);
    globalRhsVector.PutScalar(0.0);

    //CONVERT EIGEN::MATRIX AND EIGEN::VECTOR TO EPETRA FORMAT
    ConversionTools converter;
    Epetra_CrsMatrix localMatrix = converter.convertEigen2EpetraCrsMatrix(hessian0_eigen, overlappingGraph);
//    if (rank == 2)
//        std::cout << hessian0_eigen.toDense() << std::endl;
//    if (rank == numProc-2)
//        std::cout << hessian0_eigen.toDense() << std::endl;
    Epetra_Vector localRhsVector = converter.convertEigen2EpetraVector(residual_eigen, overlappingMap);

    globalMatrix.PutScalar(0.0);
    globalMatrix.Export(localMatrix, exporter, Add);    //inter-process communication of matrix entries
    globalRhsVector.Export(localRhsVector, exporter, Add);    //inter-process communication of vector entries

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
    Epetra_MultiVector sol = solveSystem(globalMatrix, lhs, globalRhsVector, true);
    sol.Scale(-1);
    //save solution
    EpetraExt::MultiVectorToMatlabFile((resultDirectory + "/solution.mat").c_str(), sol);
    EpetraExt::RowMatrixToMatlabFile((resultDirectory + "/globalMatrix.mat").c_str(), globalMatrix);
    EpetraExt::VectorToMatlabFile((resultDirectory + "/globalRHS.mat").c_str(), globalRhsVector);

#ifdef SHOW_SOLUTION
    std::cout << "Solution:\n------------\n";
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


int main(int argc, char** argv)
{
    Teuchos::GlobalMPISession mpiSession(&argc, &argv);

//    run_Test_2D_GivenMesh();
    run_Test_2D();

    return 0;
}
