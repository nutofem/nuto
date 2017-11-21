#include <iostream>
#include <mpi.h>

#include <eigen3/Eigen/Core>
#include <boost/filesystem.hpp>

//#include <Epetra_CrsMatrix.h>
//#include <Epetra_MultiVector.h>
//#include <Epetra_Vector.h>
//#include <Epetra_Map.h>
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
#include "mechanics/nodes/NodeDof.h"
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/MechanicsEnums.h"
#include "mechanics/structures/StructureOutputBlockMatrix.h"
#include "mechanics/structures/StructureOutputBlockVector.h"
#include "visualize/VisualizeEnum.h"
#include "mechanics/groups/Group.h"
#include "mechanics/constraints/ConstraintCompanion.h"

#include "ConversionTools.h"
#include "PrintTools.h"


//#define SHOW_INTERMEDIATE_RESULTS
//#define SHOW_SOLUTION


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

//        Solver.SetAztecOption(AZ_precond, AZ_dom_decomp);
//        Solver.SetAztecOption(AZ_subdomain_solve, AZ_ilu);
//        Solver.SetAztecOption(AZ_graph_fill, 1);
//        Solver.SetAztecOption(AZ_overlap, 1);
        Solver.Iterate(1000,1e-8);

//        double condest = 1e5;
//        Solver.ConstructPreconditioner(condest);
//        Solver.AdaptiveIterate(1000, 30, 1e-8);
    }
    else
    {
        Amesos Factory;
        Amesos_BaseSolver* solver;
        std::string solverType = "Mumps";
//        std::string solverType = "Klu";
        solver = Factory.Create(solverType, problem);
        solver->Solve();
    }

    Epetra_MultiVector* lhs = problem.GetLHS();
    return *lhs;
}


void saveMultiVector(Epetra_MultiVector rVector, const std::string rFileName, const std::string description)
{
    std::ofstream file(rFileName.c_str());

    if (!file.is_open())
    {
        std::cout << "ERROR: Unable to open file \"" + rFileName + "\"";
    }
    else
    {
        file << description;
        file << "\n\n";
        file << "\tLID\t\tGID\t\tValue";
        rVector.Print(std::cout);
        rVector.Print(file);
    }
    file.close();
}


//!
//! \brief Calculates strain of a rod with one fixed (left) and  one loaded (right) boundary
//!
int main(int argc, char** argv)
{
    Teuchos::GlobalMPISession mpiSession(&argc, &argv);

    int rank = -1;
    int size = -1;

    Epetra_MpiComm Comm(MPI_COMM_WORLD);
    rank = Comm.MyPID();
    size = Comm.NumProc();

    //CREATE DIRECTORIES FOR LOGS AND RESULTS
    const std::string logDirectory = "1D_Test_Logs";
    const std::string resultDirectory = "1D_Test_Results";
    checkDirectories(logDirectory, resultDirectory);


    // one-dimensional problem
    const int dim = 1;

    // --------------------
    // | DOMAIN VARIABLES |
    // --------------------
    int numSubDomains = size;   //number of sub domains equals number of processes
    int* numElementsPerSubDomain = new int[numSubDomains];
    int* numElementsTillSubDomain_inclusive = new int[numSubDomains];
    int numElementsTotal = 0;
    bool setTotalLength = true;
    double* lengthPerSubDomain = new double[numSubDomains];
    double lengthTotal = (setTotalLength ? 10. : 0.);

    for (int i = 0; i < numSubDomains; ++i)
    {
        numElementsPerSubDomain[i] = i+1;
        numElementsTotal += numElementsPerSubDomain[i];
        numElementsTillSubDomain_inclusive[i] = numElementsTotal;
        if (setTotalLength)
        {
            lengthPerSubDomain[i] = lengthTotal/numSubDomains;
        }
        else
        {
            lengthPerSubDomain[i] = 5;
            lengthTotal += lengthPerSubDomain[i];
        }
    }

    // ------------------------
    // | MECHANICAL VARIABLES |
    // ------------------------
    double area = 10.;
    double YoungsModulus = 2.e4;
    double force = 1000.;
    std::cout << YoungsModulus << std::endl;
    double fixedDisplacement = 0.;
    Eigen::VectorXd directionX(1);
    directionX(0) = 1;
    bool enableDisplacementControl = false;


    NuTo::Structure structure(dim);

    structure.SetVerboseLevel(10);
    structure.GetLogger().OpenFile(logDirectory + "/output_" + std::to_string(rank));
    structure.GetLogger().SetQuiet(false);

    auto section = NuTo::SectionTruss::Create(area);

    auto material = structure.ConstitutiveLawCreate("Linear_Elastic_Engineering_Stress");
    structure.ConstitutiveLawSetParameterDouble(material, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, YoungsModulus);

    //CREATE NODES
    Eigen::VectorXd nodeCoords(1);
    double currLength = 0.;
    for (int i = 0; i < numSubDomains; ++i)
    {
        if (rank == i)
        {
            for (int currNode = 0; currNode < numElementsPerSubDomain[i] + 1; ++currNode)
            {
                nodeCoords(0) = currNode * lengthPerSubDomain[i]/numElementsPerSubDomain[i] + currLength;
                structure.NodeCreate(currNode, nodeCoords);
            }
        }
        else
        {
            currLength += lengthPerSubDomain[i];
        }
    }

    //DEFINE INTERPOLATION
    int interpolationType = structure.InterpolationTypeCreate("Truss1D");
    structure.InterpolationTypeAdd(interpolationType, NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    structure.InterpolationTypeAdd(interpolationType, NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);

    //CREATE ELEMENTS
    std::vector<int> elementNodes(2);
    for (int i = 0; i < numSubDomains; ++i)
    {
        if (rank == i)
        {
            for (int currElement = 0; currElement < numElementsPerSubDomain[i]; ++currElement)
            {
                elementNodes[0] = currElement;
                elementNodes[1] = currElement + 1;
                structure.ElementCreate(interpolationType, elementNodes);
            }
        }
    }

    structure.ElementTotalSetSection(section);
    structure.ElementTotalSetConstitutiveLaw(material);
    structure.ElementTotalConvertToInterpolationType();

    // -----------------------------------------
    // | DEFINE BOUNDARY AND OVERLAPPING AREAS |
    // -----------------------------------------
    int nodeBCLeft = structure.GroupCreate(NuTo::eGroupId::Nodes);
    int nodeBCRight = structure.GroupCreate(NuTo::eGroupId::Nodes);
    structure.GroupAddNodeCoordinateRange(nodeBCLeft, 0, -1e-6, 1e-6);
    structure.GroupAddNodeCoordinateRange(nodeBCRight, 0, lengthTotal-1e-6, lengthTotal+1e-6);

    //COMPUTE OVERLAPPING AREAS
    int numOverlapping = numSubDomains - 1;
    int* nodeGroupsOverlap = new int[numOverlapping];
    for (int i = 0; i < numOverlapping; ++i)
        nodeGroupsOverlap[i] = structure.GroupCreate(NuTo::eGroupId::Nodes);

    currLength = lengthPerSubDomain[0];
    for (int i = 0; i < numOverlapping; ++i)
    {
        structure.GroupAddNodeCoordinateRange(nodeGroupsOverlap[i], 0, currLength - 1e-6, currLength + 1e-6);
        currLength += lengthPerSubDomain[i+1];
    }


    //FIX LEFT BOUNDARY
//    structure.ConstraintLinearSetDisplacementNodeGroup(nodeBCLeft, directionX, fixedDisplacement);
    auto leftBoundary = structure.GroupGetNodeCoordinateRange(NuTo::eDirection::X, -1e-6, 1e-6);
    std::vector<NuTo::eDirection> directionsXY(2);
    directionsXY[0] = NuTo::eDirection::X;
    directionsXY[1] = NuTo::eDirection::Y;
    structure.Constraints().Add(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Constraint::Component(leftBoundary, directionsXY, fixedDisplacement));

    //SET LOAD ON RIGHT BOUNDARY
    if (enableDisplacementControl)
    {
//        structure.ConstraintLinearSetDisplacementNode(numElements, directionX, displacement);
    }
    else
    {
//        structure.LoadCreateNodeForce(0, numElements, directionX, force);
        structure.LoadCreateNodeGroupForce(nodeBCRight, directionX, force);
    }

    //COMPUTE (LOCAL) HESSIAN AND RESIDUAL
    NuTo::StructureOutputBlockMatrix hessian0 = structure.BuildGlobalHessian0();
    NuTo::StructureOutputBlockVector dofs(structure.GetDofStatus(), true);
    dofs.J.SetZero();
    dofs.K.SetZero();

    NuTo::StructureOutputBlockVector residual = hessian0*dofs - structure.BuildGlobalExternalLoadVector() + structure.BuildGlobalInternalGradient();

    Eigen::SparseMatrix<double> hessian0_eigen = hessian0.JJ.ExportToEigenSparseMatrix();
    Eigen::MatrixXd residual_eigen = residual.J.Export();


    int ownedNumActiveDOFs = -1;
    ownedNumActiveDOFs = numElementsPerSubDomain[rank];

    int localNumActiveDOFs = -1;
    if (rank == 0)
    {
        localNumActiveDOFs = numElementsPerSubDomain[rank];
    }
    else
    {
        localNumActiveDOFs = numElementsPerSubDomain[rank] + 1;
    }


    //GET IDS OF OVERLAPPING NODES AND LEFT (FIXED) BOUNDARY
    bool leftNodeIDFound = false;
    bool overlappingNodeIDFound = false;
    std::vector<int> overlapNodeIDs;
    if (size > 1)
    {
        if (rank == 0)
            overlapNodeIDs = structure.GroupGetMemberIds(nodeGroupsOverlap[rank]);
        else
            overlapNodeIDs = structure.GroupGetMemberIds(nodeGroupsOverlap[rank-1]);
    }

    std::vector<int> leftNodeIDs = structure.GroupGetMemberIds(nodeBCLeft);
    int leftNodeCount = leftNodeIDs.size();
    int overlappingNodeCount = overlapNodeIDs.size();

    int internDofCounter = -1;
    if (rank == 0)
        internDofCounter = 0;
    else
        internDofCounter = numElementsTillSubDomain_inclusive[rank-1];


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
                if (rank == 0)
                {
                    local2GlobalMapping[dofID] = internDofCounter;  //nothing happend here
                }
                else
                {
                    local2GlobalMapping[dofID] = internDofCounter;
                }
                ++internDofCounter;
            }
        }
    }

    if (size > 1)
    {
        //Just one DOF ID
        std::vector<int> overlapDofIDs = structure.NodeGetDofIds(overlapNodeIDs[0], NuTo::Node::eDof::DISPLACEMENTS);
        if (rank == 0)
            local2GlobalMapping[overlapDofIDs[0]] = numElementsTillSubDomain_inclusive[0]-1;
        else
            local2GlobalMapping[overlapDofIDs[0]] = numElementsTillSubDomain_inclusive[rank-1]-1;
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
    owningGraph.Export(overlappingGraph, exporter, Epetra_CombineMode::Insert); //inter-process communication of matrix structure
    owningGraph.FillComplete();

    //CREATE (GLOBAL) MATRIX AND VECTOR WITH KNOWN GRAPH STRUCTURE
    Epetra_CrsMatrix globalMatrix(Epetra_DataAccess::Copy, owningGraph);
    globalMatrix.FillComplete();
    Epetra_Vector globalRhsVector(owningMap);
    globalRhsVector.PutScalar(0.0);

    //CONVERT EIGEN::MATRIX AND EIGEN::VECTOR TO EPETRA FORMAT
    ConversionTools converter(Comm);
    Epetra_CrsMatrix localMatrix = converter.convertEigen2EpetraCrsMatrix(hessian0_eigen, overlappingGraph);
    Epetra_Vector localRhsVector = converter.convertEigen2EpetraVector(residual_eigen, overlappingMap);
    globalMatrix.PutScalar(0.0);
    globalMatrix.Export(localMatrix, exporter, Epetra_CombineMode::Add);    //inter-process communication of matrix entries
    globalRhsVector.Export(localRhsVector, exporter, Epetra_CombineMode::Add);    //inter-process communication of vector entries
#ifdef SHOW_INTERMEDIATE_RESULTS
    //PRINT SOME INTERMEDIATE RESULTS
    std::ostream& os = std::cout;

    os << "Hessian0 by NuTo:\n------------------\n" << hessian0.JJ.ExportToFullMatrix() << "\n\n";
    os << "residual by NuTo:\n------------------\n" << residual << "\n\n";
    os << "Overlap Map:\n-------------\n";
    overlappingMap.Print(os);
    os << "Owning Map:\n-------------\n";
    owningMap.Print(os);
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
    Epetra_MultiVector sol = solveSystem(globalMatrix, lhs, globalRhsVector);
    sol.Scale(-1);
    //save solution
//    saveMultiVector(sol, resultDirectory + "/solution_vector", "Solution, describes displacements");
    EpetraExt::MultiVectorToMatlabFile((resultDirectory + "/solution.mat").c_str(), sol);
    EpetraExt::RowMatrixToMatlabFile((resultDirectory + "/globalMatrix.mat").c_str(), globalMatrix);
    EpetraExt::VectorToMatlabFile((resultDirectory + "/globalRHS.mat").c_str(), globalRhsVector);

#ifdef SHOW_SOLUTION
    std::cout << "Solution:\n---------------\n";
    sol.Print(std::cout);
#endif

    //solve system serial (for comparison only)
    if (size == 1)
    {
        Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
        solver.compute(hessian0_eigen);
        Eigen::VectorXd displ = solver.solve(residual.J.Export());
        displ *= -1;

        //save direct solution
        std::ofstream file(resultDirectory + "/solution_direct.mat");
        file << displ;
        file.close();
    }


    return 0;
}
