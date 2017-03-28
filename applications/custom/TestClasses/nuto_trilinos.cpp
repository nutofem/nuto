#include <iostream>
#include <mpi.h>
#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/structures/StructureBase.h"
#include "mechanics/MechanicsEnums.h"
#include "math/SparseMatrixCSRGeneral.h"
#include "mechanics/structures/StructureOutputBlockMatrix.h"
#include "mechanics/structures/StructureOutputBlockVector.h"
#include "visualize/VisualizeEnum.h"
#include "mechanics/timeIntegration/NewmarkDirect.h"
#include "mechanics/nodes/NodeDof.h"
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/sections/SectionPlane.h"
//#include "mechanics/feti/StructureFeti.h"

#include <Epetra_CrsMatrix.h>
#include <Epetra_MultiVector.h>
#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include <Epetra_MpiComm.h>
#include <Epetra_LinearProblem.h>
#include <Epetra_Import.h>
#include <Epetra_CrsGraph.h>
#include <AztecOO.h>
#include <Amesos.h>
#include <Amesos_BaseSolver.h>
#include <Amesos_Mumps.h>
#include <Amesos_ConfigDefs.h>
#include <Teuchos_GlobalMPISession.hpp>

#include <eigen3/Eigen/Core>



constexpr int dim = 2;


int* map2Array_Int(std::map<int, int> rMap)
{
    int* newArray = new int[rMap.size()];

    for (int i = 0; i < rMap.size(); ++i)
    {
        newArray[i] = rMap[i];
    }

    return newArray;
}


std::map<int, int> invertMap_int(std::map<int, int> rMap)
{
    std::map<int, int> newMap;

    for (int i = 0; i < rMap.size(); ++i)
    {
        newMap[rMap[i]] = i;
    }

    return newMap;
}



Epetra_CrsMatrix convertEigen2EpetraCrsMatrix(Eigen::SparseMatrix<double> rEigenMatrix, Epetra_Map rRowMap, Epetra_Map rColMap)
{
    int columnCount = rEigenMatrix.cols();

    int numMyElements_Row = rRowMap.NumMyElements();
//    int* myGlobalNodeIndices_Row = rRowMap.MyGlobalElements();


    Epetra_CrsMatrix convMatrix(Copy, rRowMap, rColMap, 0);


    std::vector<int> indicesVector;
    int* indices;
    int numEntries = 0;
    std::vector<double> entriesVector;
    double* entries;
    double eigenValue = 0.;


    for (int i = 0; i < numMyElements_Row; ++i)
    {
        numEntries = 0;
        for (int j = 0; j < columnCount; ++j)
        {
            eigenValue = rEigenMatrix.coeff(i, j);
            if (eigenValue != 0)
            {
                ++numEntries;
                entriesVector.push_back(eigenValue);
                indicesVector.push_back(j);
            }
        }

        entries = &entriesVector[0];
        indices = &indicesVector[0];

//        convMatrix.InsertGlobalValues(myGlobalNodeIndices_Row[i], numEntries, entries, indices);
        convMatrix.InsertMyValues(i, numEntries, entries, indices);

        entriesVector.clear();
        indicesVector.clear();
        numEntries = 0;
    }

//    convMatrix.FillComplete();
    convMatrix.FillComplete(rColMap, rRowMap);



    return convMatrix;
}


Epetra_Vector convertEigen2EpetraVector(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> rVector, Epetra_Map rMap)
{
    int oldVectorSize = int(rVector.rows());
    Epetra_Vector newVector(rMap);
//    int* myGlobalNodeIndices = rMap.MyGlobalElements();
    int numMyElements_Row = rMap.NumMyElements();


    double* vals = new double[1];
    int* indices = new int[1];
    for (int i = 0; i < oldVectorSize; ++i)
    {
        vals[0] = rVector(i, 0);
//        indices[0] = myGlobalNodeIndices[i];
        indices[0] = i;
//        newVector.ReplaceGlobalValues(1, vals, indices);
        newVector.ReplaceMyValues(1, vals, indices);
    }

    delete[] vals;
    delete[] indices;

    return newVector;

}


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


double computeNorm1Between(Eigen::VectorXd rVec1, Epetra_MultiVector rVec2)
{
    double normValue = 0.;

    for (int i = 0; i < int(rVec1.rows()); ++i)
    {
        normValue += fabs(rVec1(i,0) - rVec2[0][i]);
    }

    return normValue;
}


double computeNorm2Between(Eigen::VectorXd rVec1, Epetra_MultiVector rVec2)
{
    double tailValue = 0.;

    for (int i = 0; i < int(rVec1.rows()); ++i)
    {
        tailValue += rVec1(i,0)*rVec2[0][i];
    }

    double normVec2;
    rVec2.Norm2(&normVec2);
    double normValue = pow(rVec1.norm(),2) + pow(normVec2,2) - 2*tailValue;

    normValue = sqrt(normValue);


    return normValue;
}


void printArray_int(int* rArray, int rSize, char* rTitle, int rPID)
{
    std::cout << rTitle << "\n" << std::endl;
    for (int j = 0; j < rSize; ++j)
    {
        std::cout << "PID: " << rPID << ":  " << j << " --> " << rArray[j]  << std::endl;
    }
}


void printMap_int_int(std::map<int, int> rMap, char* rTitle, int rPID)
{
    std::cout << rTitle << "\n" << std::endl;
    for (int j = 0; j < rMap.size(); ++j)
    {
        std::cout << "PID: " << rPID << ":  " << j << " --> [" << rMap[j] << "]" << std::endl;
    }
}



int main(int argc, char** argv)
{
//    MPI_Init(&argc, &argv);
    Teuchos::GlobalMPISession mpiSession(&argc, &argv);

    int rank = -1; 
//    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int size = -1;
//    MPI_Comm_size(MPI_COMM_WORLD, &size);

    Epetra_MpiComm Comm(MPI_COMM_WORLD);
    rank = Comm.MyPID();
    size = Comm.NumProc();

    std::map<int, int> newNodes;
    std::map<int, int> gmshNodes;

    NuTo::Structure s(dim);
//    auto importContainer = s.ImportFromGmsh("trilinos.msh_00000" + std::to_string(rank+1));
    auto importContainer = s.ImportFromGmsh("trilinos.msh_00000" + std::to_string(rank+1), true, newNodes, gmshNodes);
//    auto importContainer = s.ImportFromGmsh("trilinos.msh_0", true, newNodes, gmshNodes);

//    printMap_int_int(gmshNodes, "Gmsh Nodes", Comm.MyPID());


    s.SetVerboseLevel(10);
//    s.GetLogger().OpenFile("output_" + std::to_string(0))
    s.GetLogger().OpenFile("output_" + std::to_string(rank));
    s.GetLogger().SetQuiet(false);

    // assign interpolation type
    const int InterpolationType = importContainer[0].second;
    s.InterpolationTypeAdd(InterpolationType, NuTo::Node::eDof::DISPLACEMENTS,
            NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);

    s.ElementTotalConvertToInterpolationType();

    const int nodeDOFCount = 2;


    Eigen::VectorXd directionX(dim);
    directionX(0) = 1;
    directionX(1) = 0;

    Eigen::VectorXd directionY(dim);
    directionY(0) = 0;
    directionY(1) = 1;

    //Create Load on Node
    double load = 1.e5;
    int nodeID = 2; //global node ID

    bool nodeFound = false;
    int localNodeID = -1;
    for (int i = 0; i < gmshNodes.size(); ++i)
    {
        if (gmshNodes[i] == nodeID)
        {
            nodeFound = true;
            localNodeID = i;
        }
    }

    int loadId = -1;
    if (nodeFound)
    {
        int groupLoad = s.GroupCreate(NuTo::eGroupId::Nodes);
//        s.GroupAddNode(groupLoad, nodeID);
        s.GroupAddNode(groupLoad, localNodeID);
        loadId = s.LoadCreateNodeGroupForce(0,groupLoad, directionX, load);
    }


    double thickness = 10.;
//    int section = s.SectionCreate("Plane_Strain");
//    s.SectionSetThickness(section, thickness);
    auto section = NuTo::SectionPlane::Create(thickness, true);


    s.ElementTotalSetSection(section);

    double YoungsModulus = 20000.;
    double PoissonRatio = 0.3;
    auto material = s.ConstitutiveLawCreate("Linear_Elastic_Engineering_Stress");
    s.ConstitutiveLawSetParameterDouble(material, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, YoungsModulus);
    s.ConstitutiveLawSetParameterDouble(material, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, PoissonRatio);


    s.ElementTotalSetConstitutiveLaw(material);

    //visualize
    int visualizationGroup = s.GroupCreate(NuTo::eGroupId::Elements);
    s.GroupAddElementsTotal(visualizationGroup);

    s.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::DISPLACEMENTS);
    s.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ENGINEERING_STRAIN);
    s.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ENGINEERING_STRESS);

    s.ExportVtkDataFileElements("nuto_trilinos_test.vtk");


    std::cout << "***********************************" << std::endl;
    std::cout << "**      Boundary Conditions      **" << std::endl;
    std::cout << "***********************************" << std::endl;

    int groupNodeBCLeft = s.GroupCreate(NuTo::eGroupId::Nodes);
    int groupNodeCentre = s.GroupCreate(NuTo::eGroupId::Nodes);
//    s.GroupAddNodeCoordinateRange(groupNodeBCLeft, 0, 4-1e-6, 4+1e-6);
    s.GroupAddNodeCoordinateRange(groupNodeBCLeft, 0, -1e-6, 1e-6);
    s.GroupAddNodeCoordinateRange(groupNodeCentre, 0, 4-1e-6, 4+1e-6);


    double displacement = 0;
//    int loadId = s.ConstraintLinearSetDisplacementNode(0,directionX, displacement);
    s.ConstraintLinearSetDisplacementNodeGroup(groupNodeBCLeft, directionX, displacement);
    s.ConstraintLinearSetDisplacementNodeGroup(groupNodeBCLeft, directionY, displacement);

    // solve system
//    NuTo::NewmarkDirect newmark(&s);
//    double simulationTime = 1.0;
//    double loadFactor     = 1.0;
//    newmark.SetTimeStep(0.5);
//    bool deleteDirectory = true;
//    newmark.SetResultDirectory("trilinosResults", deleteDirectory);

//    Eigen::MatrixXd timeDependentLoad(2,2);
//    timeDependentLoad(0, 0) = 0;
//    timeDependentLoad(1, 0) = simulationTime;
//    timeDependentLoad(0, 1) = 0;
//    timeDependentLoad(1, 1) = loadFactor;

//    newmark.SetTimeDependentLoadCase(loadId, timeDependentLoad);
////    newmark.AddTimeDependentConstraint(loadId, timeDependentLoad);
//    newmark.Solve(simulationTime);


    NuTo::StructureOutputBlockMatrix hessian0 = s.BuildGlobalHessian0();
    NuTo::StructureOutputBlockVector dofs(s.GetDofStatus(), true);
    dofs.J.SetZero();
    dofs.K.SetZero();

    NuTo::StructureOutputBlockVector residual = hessian0 * dofs - s.BuildGlobalExternalLoadVector(0) + s.BuildGlobalInternalGradient();


    s.Info();


    //++++++++++++++++++++++ CREATE [NODE <--> DOF]-MAPPINGS ++++++++++++++++++++++

    std::vector<std::pair<int, NuTo::NodeBase*>> nodes;
    s.GetNodesTotal(nodes);


    std::map<int, std::pair<int, int>> nodeToDOFMapping;
    std::map<int, int> dofToNodeMapping;

    for (int i = 0; i < nodes.size(); ++i)
    {
        nodeToDOFMapping[i].first = nodes[i].second->GetDof(NuTo::Node::eDof::DISPLACEMENTS, 0);
        nodeToDOFMapping[i].second = nodes[i].second->GetDof(NuTo::Node::eDof::DISPLACEMENTS, 1);

        dofToNodeMapping[int(nodes[i].second->GetDof(NuTo::Node::eDof::DISPLACEMENTS, 0))] = i;
        dofToNodeMapping[int(nodes[i].second->GetDof(NuTo::Node::eDof::DISPLACEMENTS, 1))] = i;
    }


    //++++++++++++++++++++++ GET HESSIAN MATRIX FROM NUTO ++++++++++++++++++++++


    Eigen::SparseMatrix<double> hessian0_eigen = hessian0.JJ.ExportToEigenSparseMatrix();
    std::map<NuTo::Node::eDof, int> activeDOFMap =  s.GetDofStatus().GetNumActiveDofsMap();
    std::set<NuTo::Node::eDof> dofTypes =  s.GetDofStatus().GetActiveDofTypes();


    int numLocalDOFs = activeDOFMap[NuTo::Node::eDof::DISPLACEMENTS];




    int* gmshNodeArray = map2Array_Int(gmshNodes);
//    printArray_int(gmshNodeArray, gmshNodes.size(), "GMSH_NODES\n", Comm.MyPID());
    Epetra_Map targetMapNodes(-1, gmshNodes.size(), gmshNodeArray, 0, Comm);
//    targetMapNodes.Print(std::cout);


    int numLocalNodes = targetMapNodes.NumMyElements();

    int* lIDs = new int[numLocalNodes];
    int* pIDs = new int[numLocalNodes];

    targetMapNodes.RemoteIDList(numLocalNodes, gmshNodeArray, pIDs, lIDs);



    //++++++++++++++++++++++ GLOBAL NUMBERING OF ACTIVE DOFS ++++++++++++++++++++++




    int localDOFCount = nodeDOFCount*numLocalNodes;
    int* myGlobalNodeIndices = targetMapNodes.MyGlobalElements();
    int* local2GlobalDOFMapping = new int[localDOFCount];

    int* resizedLocal2GlobalDOFMapping = new int[numLocalDOFs];
    int internDofCounter = 0;
    int subtractCounter = 0;



    std::map<int, int> old2NewDofsMapping;
    std::vector<int> leftNodeIds = s.GroupGetMemberIds(groupNodeBCLeft);
    std::vector<int> centreNodeIds = s.GroupGetMemberIds(groupNodeCentre);
    bool centreNodeIDFound = false;
    bool leftNodeIDFound = false;
    int centreNodeCount = centreNodeIds.size();
    int leftNodeCount = leftNodeIds.size();
    internDofCounter = centreNodeCount*2;
    int internCentreDofCounter = 0;
    int currDofID = centreNodeCount*2;
//    int currDofID = 50;
    for (int i = 0; i < s.GetNumNodes(); ++i)
    {
        centreNodeIDFound = false;
        for (int j = 0; j < centreNodeCount; ++j)
        {
            if (centreNodeIds[j] == i)
            {
                centreNodeIDFound = true;
            }
        }

        leftNodeIDFound = false;
        for (int j = 0; j < leftNodeCount; ++j)
        {
            if (leftNodeIds[j] == i)
            {
                leftNodeIDFound = true;
            }
        }

        if (!centreNodeIDFound && !leftNodeIDFound)
        {
            std::vector<int> dofIDs = s.NodeGetDofIds(i, NuTo::Node::eDof::DISPLACEMENTS);

            for (int dofID : dofIDs)
            {
                if (Comm.MyPID() == 0)
                {
//                    old2NewDofsMapping[dofID] = currDofID;
                    old2NewDofsMapping[dofID] = internDofCounter;
//                    local2GlobalDOFMapping[currDofID] = internDofCounter;
                    local2GlobalDOFMapping[dofID] = internDofCounter;
//                    local2GlobalDOFMapping[dofID] = dofID;
                }
                else
                {
//                    old2NewDofsMapping[dofID] = currDofID;
                    old2NewDofsMapping[dofID] = internDofCounter + 40;
//                    local2GlobalDOFMapping[currDofID] = internDofCounter + 40;
                    local2GlobalDOFMapping[dofID] = internDofCounter + 40;
//                    local2GlobalDOFMapping[dofID] = dofID + 40;
                }
                ++currDofID;
                ++internDofCounter;
            }


        }


    }

    currDofID = 0;

//    std::vector<int> globalCenterNodeIDs;
//    for (int i = 0; i < centreNodeCount; ++i)
//    {
//        globalCenterNodeIDs.push_back(gmshNodes[centreNodeIds[i]]);


//    }
//    std::sort(globalCenterNodeIDs.begin(), globalCenterNodeIDs.end());

//    std::map<int, int> global2localNodes = invertMap_int(gmshNodes);
//    int localNodeID = -1;

//    for (int i = 0; i < centreNodeCount; ++i)
//    {
//        localNodeID = global2localNodes[globalCenterNodeIDs[i]];

//        for (int k = 0; k < 2; ++k)
//        {
//            local2GlobalDOFMapping[currDofID] = currDofID;
//            ++currDofID;
//        }

//    }


    for (int i = 0; i < centreNodeCount; ++i)
    {
//        std::vector<int> centerDofIDs = s.NodeGetDofIds(i, NuTo::Node::eDof::DISPLACEMENTS);
        std::vector<int> centerDofIDs = s.NodeGetDofIds(centreNodeIds[i], NuTo::Node::eDof::DISPLACEMENTS);

        for (int dofID : centerDofIDs)
        {
//            old2NewDofsMapping[dofID] = currDofID;
            old2NewDofsMapping[dofID] = currDofID;
//            local2GlobalDOFMapping[currDofID] = currDofID;
            local2GlobalDOFMapping[dofID] = currDofID;
            ++currDofID;
        }
    }
//    printMap_int_int(old2NewDofsMapping, "old to new DOFs", Comm.MyPID());
//    printArray_int(local2GlobalDOFMapping, numLocalDOFs, "local 2 global DOFs", Comm.MyPID());



    Epetra_Map targetMapDOFs(-1, numLocalDOFs, local2GlobalDOFMapping, 0, Comm);
//    Epetra_Map targetMapDOFs(-1, numLocalDOFs, old2NewDofsMapping, 0, Comm);
    Epetra_Map sourceMapDOFs(-1, numLocalDOFs, 0, Comm);
    targetMapDOFs.Print(std::cout);
    sourceMapDOFs.Print(std::cout);


    Epetra_Map rangeMap(targetMapDOFs);
//    Epetra_Map domainMap(targetMapDOFs);
    Epetra_Map domainMap(sourceMapDOFs);

    //print first line of hessian0
    Eigen::MatrixXd hess = hessian0.JJ.ExportToFullMatrix();
    int l = 0;
    std::cout << "PID = " << Comm.MyPID() << " : " << hess(0, l) << ", ";
    for (l = 1; l < hess.cols() - 1; ++l)
    {
        std::cout << hess(0, l) << ", ";
    }
    std::cout << hess(0, l) << std::endl;


    Epetra_CrsMatrix hessian0_epetra_old = convertEigen2EpetraCrsMatrix(hessian0_eigen, rangeMap, domainMap);
    Epetra_CrsMatrix hessian0_epetra(hessian0_epetra_old);
    hessian0_epetra.Print(std::cout);
    std::cout << "----------------Conversion Eigen2EpetraCrsMatrix successful-------------------" << std::endl;
    Epetra_Vector residual_epetra = convertEigen2EpetraVector(residual.J.Export(), rangeMap);
//    residual_epetra.Print(std::cout);
    std::cout << "----------------Conversion Eigen2EpetraVector successful-------------------" << std::endl;
    Epetra_Vector lhs(domainMap);
    Epetra_MultiVector sol = solveSystem(hessian0_epetra, lhs, residual_epetra);
    sol.Scale(-1);
    std::cout << "PID: " << Comm.MyPID() << "  " << "solution_epetra:" << std::endl;
    sol.Print(std::cout);



    Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
    solver.compute(hessian0_eigen);
    Eigen::VectorXd disp = solver.solve(residual.J.Export());
    disp *= -1;


    std::cout << "disp \n" << disp << std::endl;
//    std::cout << "PID: " << Comm.MyPID() << "  " << "Norm_2(solution_nuto) = " << disp.norm() << std::endl;
    double norm2;
    sol.Norm2(&norm2);
//    std::cout << "PID: " << Comm.MyPID() << "  " << "Norm_2(solution_trilinos) = " << norm2 << std::endl;
    std::cout << "PID: " << Comm.MyPID() << "  " << "Norm_1(solution_nuto - solution_trilinos) = " << computeNorm1Between(disp, sol) << std::endl;
    std::cout << "PID: " << Comm.MyPID() << "  " << "Norm_2(solution_nuto - solution_trilinos) = " << computeNorm2Between(disp, sol) << std::endl;


//    MPI_Finalize();

}


