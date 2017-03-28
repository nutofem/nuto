#include <iostream>
#include <mpi.h>

#include <eigen3/Eigen/Core>

#include <Epetra_CrsMatrix.h>
#include <Epetra_MultiVector.h>
#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include <Epetra_MpiComm.h>
#include <Epetra_LinearProblem.h>
#include <AztecOO.h>
#include <Amesos.h>
#include <Amesos_BaseSolver.h>
#include <Amesos_Mumps.h>
#include <Amesos_ConfigDefs.h>
#include <Teuchos_GlobalMPISession.hpp>

#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/structures/StructureBase.h"
#include "mechanics/sections/SectionTruss.h"
#include "mechanics/nodes/NodeDof.h"
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/MechanicsEnums.h"
#include "mechanics/structures/StructureOutputBlockMatrix.h"
#include "mechanics/structures/StructureOutputBlockVector.h"
#include "visualize/VisualizeEnum.h"

#include "ConversionTools.h"
#include "PrintTools.h"




Epetra_MultiVector solveSystem(Epetra_CrsMatrix rA, Epetra_Vector rLhs, Epetra_Vector rRhs)
{
    Epetra_LinearProblem problem(&rA, &rLhs, &rRhs);
    AztecOO Solver(problem);
    Solver.Iterate(1000,1e-8);
//    double condest = 1e5;
//    Solver.ConstructPreconditioner(condest);
//    Solver.AdaptiveIterate(1000, 30, 1e-8);
//    Amesos Factory;
//    Amesos_BaseSolver* solver;
//    std::string solverType = "Mumps";
//    std::string solverType = "Klu";
//    solver = Factory.Create(solverType, problem);
//    solver->Solve();


    Epetra_MultiVector* lhs = problem.GetLHS();
    return *lhs;
}



int main(int argc, char** argv)
{
    Teuchos::GlobalMPISession mpiSession(&argc, &argv);

    int rank = -1;
    int size = -1;

    Epetra_MpiComm Comm(MPI_COMM_WORLD);
    rank = Comm.MyPID();
    size = Comm.NumProc();

    int dim = 1;
    double area = 10.;
    int numLeftElements = 1;
    int numRightElements = 1;
    int numElements = numLeftElements + numRightElements;
    double length = 10.;
    double leftLength = numLeftElements*length/numElements;
    double rightLength = numRightElements*length/numElements;
    double YoungsModulus = 20000.;
    double force = 1000.;
    double displacement = 0.1;
    double fixedDisplacement = 0.;
    Eigen::VectorXd directionX(1);
    directionX(0) = 1;
    bool enableDisplacementControl = false;

    NuTo::Structure structure(dim);

    structure.SetVerboseLevel(10);
    structure.GetLogger().OpenFile("1D_Test_Results/output_" + std::to_string(rank));
    structure.GetLogger().SetQuiet(false);

    auto section = NuTo::SectionTruss::Create(area);

    auto material = structure.ConstitutiveLawCreate("Linear_Elastic_Engineering_Stress");
    structure.ConstitutiveLawSetParameterDouble(material, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, YoungsModulus);


    Eigen::VectorXd nodeCoords(1);
    if (Comm.MyPID() == 0)
    {
        for (int currNode = 0; currNode < numLeftElements + 1; ++currNode)
        {
            nodeCoords(0) = currNode * leftLength/numLeftElements;
            structure.NodeCreate(currNode, nodeCoords);
        }
    }
    else
    {
        for (int currNode = 0; currNode < numRightElements + 1; ++currNode)
        {
            nodeCoords(0) = currNode * rightLength/numRightElements + leftLength;
            structure.NodeCreate(currNode, nodeCoords);
        }
    }


    int interpolationType = structure.InterpolationTypeCreate("Truss1D");
    structure.InterpolationTypeAdd(interpolationType, NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    structure.InterpolationTypeAdd(interpolationType, NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);

    std::vector<int> elementNodes(2);
    if (Comm.MyPID() == 0)
    {
        for (int currElement = 0; currElement < numLeftElements; ++currElement)
        {
            elementNodes[0]= currElement;
            elementNodes[1] = currElement + 1;
            structure.ElementCreate(interpolationType, elementNodes);
            structure.ElementSetSection(currElement, section);
//            structure.ElementSetConstitutiveLaw(currElement, material);
        }
    }
    else
    {
        for (int currElement = 0; currElement < numRightElements; ++currElement)
        {
            elementNodes[0]= currElement;
            elementNodes[1] = currElement + 1;
            structure.ElementCreate(interpolationType, elementNodes);
            structure.ElementSetSection(currElement, section);
//            structure.ElementSetConstitutiveLaw(currElement, material);
        }
    }

    structure.ElementTotalSetConstitutiveLaw(material);
    structure.ElementTotalConvertToInterpolationType();


    int nodeBCLeft = structure.GroupCreate(NuTo::eGroupId::Nodes);
    int nodeCenter = structure.GroupCreate(NuTo::eGroupId::Nodes);
    int nodeBCRight = structure.GroupCreate(NuTo::eGroupId::Nodes);
    structure.GroupAddNodeCoordinateRange(nodeBCLeft, 0, -1e-6, 1e-6);
    structure.GroupAddNodeCoordinateRange(nodeCenter, 0, leftLength-1e-6, leftLength+1e-6);
    structure.GroupAddNodeCoordinateRange(nodeBCRight, 0, length-1e-6, length+1e-6);

//    structure.ConstraintLinearSetDisplacementNode(0, directionX, 0.0);
    structure.ConstraintLinearSetDisplacementNodeGroup(nodeBCLeft, directionX, fixedDisplacement);

    structure.SetNumLoadCases(1);
    if (enableDisplacementControl)
    {
//        structure.ConstraintLinearSetDisplacementNode(numElements, directionX, displacement);
    }
    else
    {
//        structure.LoadCreateNodeForce(0, numElements, directionX, force);
        structure.LoadCreateNodeGroupForce(0, nodeBCRight, directionX, force);
    }

    NuTo::StructureOutputBlockMatrix hessian0 = structure.BuildGlobalHessian0();
    NuTo::StructureOutputBlockVector dofs(structure.GetDofStatus(), true);
    dofs.J.SetZero();
    dofs.K.SetZero();

    NuTo::StructureOutputBlockVector residual = hessian0*dofs - structure.BuildGlobalExternalLoadVector(0) + structure.BuildGlobalInternalGradient();

//    structure.Info();

    Eigen::SparseMatrix<double> hessian0_eigen = hessian0.JJ.ExportToEigenSparseMatrix();
    std::cout << "Hessian0: \n" << hessian0.JJ.ExportToFullMatrix() << std::endl;
    std::cout << "residual: \n" << residual << std::endl;



    int currNumDOFs;
    if (Comm.MyPID() == 0)
    {
        currNumDOFs = numLeftElements;
    }
    else
    {
        currNumDOFs = numRightElements+1;
    }

    int* local2GlobalMapping = new int[currNumDOFs];
    int internDofCounter = 1;   //0 is reserved for overlapping node

    bool leftNodeIDFound = false;
    bool centerNodeIDFound = false;
    std::vector<int> leftNodeIDs = structure.GroupGetMemberIds(nodeBCLeft);
    std::vector<int> centerNodeIDs = structure.GroupGetMemberIds(nodeCenter);
    int leftNodeCount = leftNodeIDs.size();
    int centerNodeCount = centerNodeIDs.size();

    for (int i = 0; i < structure.GetNumNodes(); ++i)
    {
        leftNodeIDFound = false;
        centerNodeIDFound = false;

        for (int j = 0; j < leftNodeCount; ++j)
        {
            if (leftNodeIDs[j] == i)
            {
                leftNodeIDFound = true;
            }
        }

        for (int j = 0; j < centerNodeCount; ++j)
        {
            if (centerNodeIDs[j] == i)
            {
                centerNodeIDFound = true;
            }
        }


        if (!leftNodeIDFound && !centerNodeIDFound)
        {
            std::vector<int> dofIDs = structure.NodeGetDofIds(i, NuTo::Node::eDof::DISPLACEMENTS);
            for (int dofID : dofIDs)
            {
                if (Comm.MyPID() == 0)
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

    std::vector<int> centerDofIDs = structure.NodeGetDofIds(centerNodeIDs[0], NuTo::Node::eDof::DISPLACEMENTS);
    local2GlobalMapping[centerDofIDs[0]] = 0;


    PrintTools printer;
    printer.printArray_int(local2GlobalMapping, currNumDOFs , "local2Global", Comm.MyPID());

    Epetra_Map rowMap(-1, currNumDOFs, 0, Comm);
    Epetra_Map colMap(-1, currNumDOFs, local2GlobalMapping, 0, Comm);
//    colMap.Print(std::cout);

    ConversionTools converter;
    Epetra_CrsMatrix mat = converter.convertEigen2EpetraCrsMatrix(hessian0_eigen, rowMap, colMap);
    mat.Print(std::cout);

    Epetra_Vector vec = converter.convertEigen2EpetraVector(residual.J.Export(), colMap);
    vec.Print(std::cout);
    Epetra_Vector lhs(colMap);
    Epetra_MultiVector sol = solveSystem(mat, lhs, vec);
    sol.Scale(-1);
    sol.Print(std::cout);


    Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
    solver.compute(hessian0_eigen);
    Eigen::VectorXd displ = solver.solve(residual.J.Export());
    displ *= -1;

    std::cout << "displ:\n" << displ << std::endl;

    Eigen::VectorXd newCoordinates(numElements+1);
    std::vector<NuTo::NodeBase*> nodes;
    structure.GetNodesTotal(nodes);

    //not numElements + 1, because one node (the 0th) is constrained to the left side
//    newCoordinates(0) = nodes[0]->Get(NuTo::Node::eDof::COORDINATES)(0);
//    for (int i = 0; i < numElements; ++i)
//    {
//        newCoordinates(i+1) = displ(i) + nodes[i+1]->Get(NuTo::Node::eDof::COORDINATES)(0);
//    }
//    std::cout << "newCoordinates:\n" << newCoordinates << std::endl;


    int visualGroup = structure.GroupCreate(NuTo::eGroupId::Elements);
    structure.GroupAddElementsTotal(visualGroup);

    structure.AddVisualizationComponent(visualGroup, NuTo::eVisualizeWhat::DISPLACEMENTS);
    structure.AddVisualizationComponent(visualGroup, NuTo::eVisualizeWhat::ENGINEERING_STRAIN);
    structure.AddVisualizationComponent(visualGroup, NuTo::eVisualizeWhat::ENGINEERING_STRESS);
    structure.ExportVtkDataFileElements("1D_Test_Results/result.vtk");

    return 0;
}
