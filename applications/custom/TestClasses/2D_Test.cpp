#include <iostream>
#include <mpi.h>

#include <eigen3/Eigen/Core>

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

#include "ConversionTools.h"
#include "PrintTools.h"


//#define SHOW_INTERMEDIATE_RESULTS
#define SHOW_SOLUTION



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

void run_Test_2D_GivenMesh()
{
    Epetra_MpiComm Comm(MPI_COMM_WORLD);
    int rank = Comm.MyPID();
    int numProc = Comm.NumProc();

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
    structure.GetLogger().OpenFile("2D_Test_Logs/output_" + std::to_string(rank));
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
    EpetraExt::MultiVectorToMatlabFile("2D_Test_Results/solution.mat", sol);
    EpetraExt::RowMatrixToMatlabFile("2D_Test_Results/globalMatrix.mat", globalMatrix);
    EpetraExt::VectorToMatlabFile("2D_Test_Results/globalRHS.mat", globalRhsVector);

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
        std::ofstream file("2D_Test_Results/solution_direct.mat");
        file << displ;
        file.close();
    }
}

int main(int argc, char** argv)
{
    Teuchos::GlobalMPISession mpiSession(&argc, &argv);

    run_Test_2D_GivenMesh();

    return 0;
}
