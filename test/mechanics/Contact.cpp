#include <iostream>
#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/structures/StructureOutputBlockMatrix.h"
#include "nuto/mechanics/structures/StructureOutputBlockVector.h"

#include <eigen3/Eigen/Core>

#include "nuto/math/FullMatrix.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"
#include "nuto/math/SparseMatrixCSRSymmetric.h"
#include "nuto/math/SparseDirectSolverMUMPS.h"
#include "nuto/math/SparseDirectSolverPardiso.h"

#include "nuto/mechanics/integrationtypes/IntegrationType1D2NLobatto3Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType1D2NLobatto4Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType1D2NLobatto5Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType2D4NLobatto9Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType2D4NLobatto16Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType2D4NLobatto25Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType1D2NGauss12Ip.h"

#include "nuto/mechanics/constitutive/ConstitutiveEnum.h"
#include "nuto/mechanics/interpolationtypes/InterpolationTypeEnum.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeEnum.h"
#include "nuto/mechanics/elements/ElementDataEnum.h"
#include "nuto/mechanics/elements/IpDataEnum.h"

#include "nuto/mechanics/IGA/NURBSCurve.h"
#include "nuto/mechanics/IGA/NURBSSurface.h"

#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/nodes/NodeEnum.h"

#include "nuto/mechanics/groups/GroupEnum.h"

#include <boost/filesystem.hpp>

#include "nuto/mechanics/elements/ElementShapeFunctions.h"

#include "nuto/mechanics/timeIntegration/RungeKutta4.h"

#ifdef ENABLE_VISUALIZE
#include "nuto/visualize/VisualizeEnum.h"
#endif

#define PRINTRESULT true


int buildStructure2D(NuTo::Interpolation::eTypeOrder rElementTypeIdent,
                      int rNumNodesPerElementInOneDir,
                      NuTo::FullVector<double, Eigen::Dynamic>& nodeCoordinatesFirstElement,
                      int NumElementsX,
                      int NumElementsY,
                      double Height,
                      double Length,
                      double startx,
                      double starty,
                      int startNode,
                      NuTo::Structure* myStructure,
                      int rNodeGroup,
                      int rElementGroup,
                      const std::set<NuTo::Node::eDof> &setOfDOFS)
{
    /** parameters **/
    double elementSize = nodeCoordinatesFirstElement(rNumNodesPerElementInOneDir-1) - nodeCoordinatesFirstElement(0);
    // reference element Length is 2.
    double  factorX =  Length/(NumElementsX*elementSize);
    double  factorY =  Height/(NumElementsY*elementSize);

    NuTo::FullVector<double,Eigen::Dynamic> nodeStart(2);
    nodeStart(0) = startx;
    nodeStart(1) = starty;

    /** Nodes **/
    NuTo::FullVector<double,Eigen::Dynamic> nodeCoordinates(2);
    int node = startNode;
    double elementBeginX = 0.;
    double elementBeginY = 0.;

    nodeCoordinates(0) = factorX*nodeCoordinatesFirstElement(0);
    nodeCoordinates(1) = factorY*nodeCoordinatesFirstElement(0);
    nodeCoordinates += nodeStart;

    myStructure->NodeCreateDOFs(node, setOfDOFS, nodeCoordinates);
//    myStructure->NodeCreate(node, nodeCoordinates);
    myStructure->GroupAddNode(rNodeGroup, node);
    node++;

    for(int y = 0; y < NumElementsY; y++)
    {
        for (int i = ((y==0)?0:1); i < nodeCoordinatesFirstElement.size(); i++)
        {
            if (node != startNode + 1)
            {
                nodeCoordinates(0) = factorX*nodeCoordinatesFirstElement(0);
                nodeCoordinates(1) = factorY*(nodeCoordinatesFirstElement(i) + elementBeginY);
                nodeCoordinates += nodeStart;

                myStructure->NodeCreateDOFs(node, setOfDOFS, nodeCoordinates);
                myStructure->GroupAddNode(rNodeGroup, node);
//                myStructure->NodeCreate(node, nodeCoordinates);
                node++;
            }

            elementBeginX = 0.;
            for(int x = 0; x < NumElementsX; x++)
            {
                for (int j = 1; j < nodeCoordinatesFirstElement.size(); j++)
                {
                    nodeCoordinates(0) = factorX*(nodeCoordinatesFirstElement(j) + elementBeginX);
                    nodeCoordinates(1) = factorY*(nodeCoordinatesFirstElement(i) + elementBeginY);
                    nodeCoordinates += nodeStart;

                    myStructure->NodeCreateDOFs(node, setOfDOFS, nodeCoordinates);
                    myStructure->GroupAddNode(rNodeGroup, node);
//                    myStructure->NodeCreate(node, nodeCoordinates);
                    node++;
                }
                elementBeginX += elementSize;
            }
        }
        elementBeginY += elementSize;
    }

    int interpolationType = myStructure->InterpolationTypeCreate("QUAD2D");
    myStructure->InterpolationTypeAdd(interpolationType, NuTo::Node::eDof::DISPLACEMENTS, rElementTypeIdent);
    myStructure->InterpolationTypeAdd(interpolationType, NuTo::Node::eDof::COORDINATES, rElementTypeIdent);

    /** Elements **/
    NuTo::FullVector<int,Eigen::Dynamic> elementIncidence(rNumNodesPerElementInOneDir*rNumNodesPerElementInOneDir);
    int numNodesInRow = NumElementsX*(rNumNodesPerElementInOneDir-1) + 1;
    for(int j = 0; j < NumElementsY; j++)
    {
        for (int i = 0; i < NumElementsX; i++)
        {
            // one element
            for (int k = 0; k < rNumNodesPerElementInOneDir; k++)
                for(int l = 0; l < rNumNodesPerElementInOneDir; l++)
                    elementIncidence(k + l*rNumNodesPerElementInOneDir) = startNode + i*(rNumNodesPerElementInOneDir - 1) + k + l*numNodesInRow + j*(rNumNodesPerElementInOneDir-1)*numNodesInRow;

            int myElement = myStructure->ElementCreate(interpolationType, elementIncidence);
            myStructure->GroupAddElement(rElementGroup, myElement);
        }
    }
    return node;
}

void dynamicExplicitContactTest(const std::string &resultDir,
                                const std::string &path,
                                const std::string &fileNameSlave,
                                const std::string &fileNameMaster,
                                double rPenalty,
                                int rContactAlgo,
                                int numElXSlave, int numElYSlave,
                                int numElXMaster, int numElYMaster)
{
    NuTo::Structure myStructure(2);
    myStructure.SetNumTimeDerivatives(2);

#ifdef _OPENMP
    int numThreads = 4;
    myStructure.SetNumProcessors(numThreads);
#endif

    double E = 1.e5;
    double nue = 0.3;
    double rho = 1.;

    ///////////////////////////////////////////////////////////////////////
    // ====> create SLAVE mesh from gmsh (smth. impacting a rectangle)//
    ///////////////////////////////////////////////////////////////////////

    NuTo::IntegrationType1D2NLobatto3Ip Lobatto1D2N3Ip;

    NuTo::FullVector<double, Eigen::Dynamic> ones(3); ones.fill(1);
    NuTo::FullVector<double, Eigen::Dynamic> nodeCoordinates(3);
    for (int i = 0; i < 3; i++) Lobatto1D2N3Ip.GetLocalIntegrationPointCoordinates1D(i, nodeCoordinates(i));
    nodeCoordinates += ones;

    double startySlave = 0.00001;

    std::set<NuTo::Node::eDof> setOfDOFSSlave;
    setOfDOFSSlave.insert(NuTo::Node::eDof::COORDINATES);
    setOfDOFSSlave.insert(NuTo::Node::eDof::DISPLACEMENTS);

    int groupNodesSlave = myStructure.GroupCreate(NuTo::eGroupId::Nodes);
    int groupElementsSlave = myStructure.GroupCreate(NuTo::eGroupId::Elements);
    int node = buildStructure2D(NuTo::Interpolation::eTypeOrder::LOBATTO2, 3, nodeCoordinates, numElXSlave, numElYSlave, 5., 10., 0., startySlave, 0, &myStructure, groupNodesSlave, groupElementsSlave, setOfDOFSSlave);

    // ===> build contact elements slave elements
    auto LambdaGetSlaveNodesLower = [startySlave](NuTo::NodeBase* rNodePtr) -> bool
    {
        double Tol = 1.e-6;
        if (rNodePtr->GetNum(NuTo::Node::eDof::COORDINATES)>0)
        {
            double y = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[1];
            if ((y >= startySlave - Tol && y <= startySlave + Tol))
            {
                return true;
            }
        }
        return false;
    };  // LambdaGetSlaveNodesLower

    int groupNodesSlaveLower = myStructure.GroupCreate(NuTo::eGroupId::Nodes);
    myStructure.GroupAddNodeFunction(groupNodesSlaveLower, LambdaGetSlaveNodesLower);
    NuTo::FullVector<int, Eigen::Dynamic> members = myStructure.GroupGetMemberIds(groupNodesSlaveLower);

    int groupElementsSlaveLower = myStructure.GroupCreate(NuTo::eGroupId::Elements);
    myStructure.GroupAddElementsFromNodes(groupElementsSlaveLower, groupNodesSlaveLower, false);
    members = myStructure.GroupGetMemberIds(groupElementsSlaveLower);

    //////////////////////////////
    // ====> create MASTER mesh //
    //////////////////////////////

    std::set<NuTo::Node::eDof> setOfDOFSMaster;
    setOfDOFSMaster.insert(NuTo::Node::eDof::COORDINATES);
    setOfDOFSMaster.insert(NuTo::Node::eDof::DISPLACEMENTS);

    int groupNodesMaster = myStructure.GroupCreate(NuTo::eGroupId::Nodes);
    int groupElementsMaster = myStructure.GroupCreate(NuTo::eGroupId::Elements);
    buildStructure2D(NuTo::Interpolation::eTypeOrder::LOBATTO2, 3, nodeCoordinates, numElXMaster, numElYMaster, 5., 20., -5., -5., node, &myStructure, groupNodesMaster, groupElementsMaster, setOfDOFSMaster);

    /////////////////////////////////////////////////////
    // ====> create FE slave mesh from gmsh (rectange) //
    /////////////////////////////////////////////////////

    auto LambdaGetMasterNodesLower = [](NuTo::NodeBase* rNodePtr) -> bool
    {
        double Tol = 1.e-6;
        if (rNodePtr->GetNum(NuTo::Node::eDof::COORDINATES)>0)
        {
            double y = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[1];
            if ((y >= -5. - Tol && y <= -5. + Tol))
            {
                return true;
            }
        }
        return false;
    };  // LambdaGetMasterNodesLower

    // DBC for the bottom of the master
    int groupNodesMasterLower = myStructure.GroupCreate(NuTo::eGroupId::Nodes);
    myStructure.GroupAddNodeFunction(groupNodesMasterLower, LambdaGetMasterNodesLower);

    NuTo::FullVector<int,Eigen::Dynamic> rMembersMasterDBC;
    myStructure.NodeGroupGetMembers(groupNodesMasterLower, rMembersMasterDBC);

    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> coordinatesMasterLower;
    myStructure.NodeGroupGetCoordinates(groupNodesMasterLower, coordinatesMasterLower);

    std::cout << coordinatesMasterLower << std::endl;
    int countDBC = 0;

    NuTo::FullVector<double,Eigen::Dynamic> direction(2);
    direction << 0, 1;
    for(int i = 0; i < rMembersMasterDBC.rows(); i++)
    {
        countDBC = myStructure.ConstraintLinearSetDisplacementNode(rMembersMasterDBC(i), direction, 0.0);
    }
    direction << 1, 0;
    countDBC = myStructure.ConstraintLinearSetDisplacementNode(rMembersMasterDBC(0), direction, 0.0);
    countDBC++;

    //////////////////////////
    // ===> build iga layer //
    //////////////////////////

    auto LambdaGetMasterNodesUpper = [](NuTo::NodeBase* rNodePtr) -> bool
    {
        double Tol = 1.e-6;
        if (rNodePtr->GetNum(NuTo::Node::eDof::COORDINATES)>0)
        {
            double y = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[1];
            if ((y >= 0. - Tol && y <= 0. + Tol))
            {
                return true;
            }
        }
        return false;
    };  // LambdaGetMasterNodesUpper

    int groupNodesMasterUpper = myStructure.GroupCreate(NuTo::eGroupId::Nodes);
    myStructure.GroupAddNodeFunction(groupNodesMasterUpper, LambdaGetMasterNodesUpper);
    NuTo::FullVector<int,Eigen::Dynamic> rMembersMasterContactBoundary;
    myStructure.NodeGroupGetMembers(groupNodesMasterUpper, rMembersMasterContactBoundary);

    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> coordinates;
    myStructure.NodeGroupGetCoordinates(groupNodesMasterUpper, coordinates);

    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> coordinatesAndIDs;
    coordinatesAndIDs.resize(coordinates.rows(), coordinates.cols() + 1);
    coordinatesAndIDs.block(0,0,coordinates.rows(), coordinates.cols()) = coordinates;

    for(int i = 0; i < coordinates.rows(); i++)
        coordinatesAndIDs(i,coordinates.cols()) = rMembersMasterContactBoundary(i);

    coordinatesAndIDs.SortRow(0);

    // create IGA curve
    int degree = 2;

    Eigen::MatrixXd A;

    NuTo::NURBSCurve curve(degree, coordinatesAndIDs.block(0, 0, coordinates.rows(), coordinates.cols()), A);

    for(int i = 0; i < 0; i++)
    {
        curve.DuplicateKnots();
    }

    std::set<NuTo::Node::eDof> setOfDOFS;
    setOfDOFS.insert(NuTo::Node::eDof::COORDINATES);
    setOfDOFS.insert(NuTo::Node::eDof::DISPLACEMENTS);

    int groupNodesIGAlayer    = myStructure.GroupCreate("Nodes");
    int groupElementsIGAlayer = myStructure.GroupCreate("Elements");

    Eigen::Matrix<std::pair<int, int>, Eigen::Dynamic, Eigen::Dynamic> elementsMaster = curve.buildIGAStructure(myStructure, setOfDOFS, groupElementsIGAlayer, groupNodesIGAlayer, "IGA1DLAYER");


    ///////////////////
    // ===> material //
    ///////////////////

    double Thickness = 1.;
    int section = myStructure.SectionCreate("PLANE_STRESS");
    myStructure.SectionSetThickness(section, Thickness);

    int constitutiveLaw = myStructure.ConstitutiveLawCreate("Linear_Elastic_Engineering_Stress");
    myStructure.ConstitutiveLawSetParameterDouble(constitutiveLaw, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, E);
    myStructure.ConstitutiveLawSetParameterDouble(constitutiveLaw, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, nue);
    myStructure.ConstitutiveLawSetParameterDouble(constitutiveLaw, NuTo::Constitutive::eConstitutiveParameter::DENSITY, rho);
    myStructure.ElementTotalSetSection(section);
    myStructure.ElementTotalSetConstitutiveLaw(constitutiveLaw);

    ////////////////////////////
    // ===> CONTACT ELEMENTS  //
    ////////////////////////////

    myStructure.NuTo::Structure::ContactElementsCreate<2,1>(groupElementsSlaveLower, groupNodesSlaveLower, elementsMaster, NuTo::eIntegrationType::IntegrationType1D2NGauss12Ip, rPenalty, rContactAlgo);

    ////////////////////////////
    // ===> Constraints       //
    ////////////////////////////

    // coordinates
    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> coordinatesLayer;
    myStructure.NodeGroupGetCoordinates(groupNodesIGAlayer, coordinatesLayer);

    // IDs
    NuTo::FullVector<int,Eigen::Dynamic> rMembersMasterContactBoundaryLayer;
    myStructure.NodeGroupGetMembers(groupNodesIGAlayer, rMembersMasterContactBoundaryLayer);

    // all together
    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> coordinatesAndIDsLayer;
    coordinatesAndIDsLayer.resize(coordinatesLayer.rows(), coordinatesLayer.cols() + 1);
    coordinatesAndIDsLayer.block(0,0,coordinatesLayer.rows(), coordinatesLayer.cols()) = coordinatesLayer;

    for(int i = 0; i < coordinatesLayer.rows(); i++)
        coordinatesAndIDsLayer(i,coordinatesLayer.cols()) = rMembersMasterContactBoundaryLayer(i);

    coordinatesAndIDsLayer.SortRow(0);

    int dim = coordinatesAndIDs.cols() - 1;

    for(int nodeMaster = 0; nodeMaster < coordinatesAndIDs.rows(); nodeMaster++)
    {
        for(int dof = 0; dof < dim; dof++)
        {
            myStructure.ConstraintLinearEquationCreate (countDBC, coordinatesAndIDs(nodeMaster, dim), NuTo::Node::eDof::DISPLACEMENTS, dof,  1., 0.);
            for(int controlPoint = 0; controlPoint < coordinatesAndIDsLayer.rows(); controlPoint++)
            {
                myStructure.ConstraintLinearEquationAddTerm(countDBC, coordinatesAndIDsLayer(controlPoint, dim), NuTo::Node::eDof::DISPLACEMENTS, dof, -A(nodeMaster, controlPoint));
            }
            countDBC++;
        }
    }

    ///////////////////
    // ===> Solution //
    ///////////////////

    myStructure.NodeInfo(10);

    myStructure.CalculateMaximumIndependentSets();
    myStructure.NodeBuildGlobalDofs();


    NuTo::FullVector<double,Eigen::Dynamic>  velocities(2);
    velocities(0) =   0.;
    velocities(1) = -20.;
    myStructure.NodeGroupSetDisplacements(groupNodesSlave, 1, velocities);

    NuTo::RungeKutta4 myIntegrationScheme(&myStructure);

    double simulationTime(0.00005);

//    myIntegrationScheme.SetMinTimeStep(0.0001 * myIntegrationScheme.GetMaxTimeStep());
    double h = 1.e-8;
    myIntegrationScheme.SetTimeStep(h);
    int numOutputs = 700;
    myIntegrationScheme.SetMinTimeStepPlot(simulationTime/numOutputs);

    //set output during the simulation to false
    myStructure.SetShowTime(false);

    myIntegrationScheme.AddResultTime("Time");

    //set result directory
    bool deleteResultDirectoryFirst(true);
    myIntegrationScheme.SetResultDirectory(resultDir, deleteResultDirectoryFirst);

#ifdef ENABLE_VISUALIZE
    myStructure.AddVisualizationComponent(groupElementsIGAlayer, NuTo::eVisualizeWhat::DISPLACEMENTS);
    // master + slave fe only
    int visualizationGroup = myStructure.GroupUnion(groupElementsSlave, groupElementsMaster);
    myStructure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::DISPLACEMENTS);
    myStructure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ENGINEERING_STRAIN);
    myStructure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ENGINEERING_STRESS);
    // master + layer
    visualizationGroup = myStructure.GroupUnion(groupElementsIGAlayer, groupElementsMaster);
    myStructure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::DISPLACEMENTS);
#endif

    //solve (perform Newton raphson iteration
    myIntegrationScheme.Solve(simulationTime);
}

void patchTest(const std::string &resultDir,
               const std::string &path,
               const std::string &fileNameSlave,
               const std::string &fileNameMaster,
               int rContactAlgo,
               int numElXSlave, int numElYSlave,
               int numElXMaster, int numElYMaster)
{
    NuTo::Structure myStructure(2);
    myStructure.SetNumTimeDerivatives(2);

#ifdef _OPENMP
    int numThreads = 4;
    myStructure.SetNumProcessors(numThreads);
#endif

    double E = 1.e5;
    double nue = 0.3;
    double rho = 1.;

    ///////////////////////////////////////////////////////////////////////
    // ====> create SLAVE mesh from gmsh (smth. impacting a rectangle)//
    ///////////////////////////////////////////////////////////////////////

    NuTo::IntegrationType1D2NLobatto3Ip Lobatto1D2N3Ip;

    NuTo::FullVector<double, Eigen::Dynamic> ones(3); ones.fill(1);
    NuTo::FullVector<double, Eigen::Dynamic> nodeCoordinates(3);
    for (int i = 0; i < 3; i++) Lobatto1D2N3Ip.GetLocalIntegrationPointCoordinates1D(i, nodeCoordinates(i));
    nodeCoordinates += ones;

    double startySlave = 0.0;

    std::set<NuTo::Node::eDof> setOfDOFSSlave;
    setOfDOFSSlave.insert(NuTo::Node::eDof::COORDINATES);
    setOfDOFSSlave.insert(NuTo::Node::eDof::DISPLACEMENTS);
    setOfDOFSSlave.insert(NuTo::Node::eDof::LAGRANGEMULTIPLIER);

    int groupNodesSlave = myStructure.GroupCreate(NuTo::eGroupId::Nodes);
    int groupElementsSlave = myStructure.GroupCreate(NuTo::eGroupId::Elements);
    int node = buildStructure2D(NuTo::Interpolation::eTypeOrder::LOBATTO2, 3, nodeCoordinates, numElXSlave, numElYSlave, 5., 10., 0., startySlave, 0, &myStructure, groupNodesSlave, groupElementsSlave, setOfDOFSSlave);

    // ===> build contact elements slave elements
    auto LambdaGetSlaveNodesLower = [startySlave](NuTo::NodeBase* rNodePtr) -> bool
    {
        double Tol = 1.e-6;
        if (rNodePtr->GetNum(NuTo::Node::eDof::COORDINATES)>0)
        {
            double y = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[1];
            if ((y >= startySlave - Tol && y <= startySlave + Tol))
            {
                return true;
            }
        }
        return false;
    };  // LambdaGetSlaveNodesLower

    int groupNodesSlaveLower = myStructure.GroupCreate(NuTo::eGroupId::Nodes);
    myStructure.GroupAddNodeFunction(groupNodesSlaveLower, LambdaGetSlaveNodesLower);
    NuTo::FullVector<int, Eigen::Dynamic> members = myStructure.GroupGetMemberIds(groupNodesSlaveLower);

    int groupElementsSlaveLower = myStructure.GroupCreate(NuTo::eGroupId::Elements);
    myStructure.GroupAddElementsFromNodes(groupElementsSlaveLower, groupNodesSlaveLower, false);
    members = myStructure.GroupGetMemberIds(groupElementsSlaveLower);

    //////////////////////////////
    // ====> create MASTER mesh //
    //////////////////////////////

    std::set<NuTo::Node::eDof> setOfDOFSMaster;
    setOfDOFSMaster.insert(NuTo::Node::eDof::COORDINATES);
    setOfDOFSMaster.insert(NuTo::Node::eDof::DISPLACEMENTS);

    int groupNodesMaster = myStructure.GroupCreate(NuTo::eGroupId::Nodes);
    int groupElementsMaster = myStructure.GroupCreate(NuTo::eGroupId::Elements);
    buildStructure2D(NuTo::Interpolation::eTypeOrder::LOBATTO2, 3, nodeCoordinates, numElXMaster, numElYMaster, 5., 20., -5., -5., node, &myStructure, groupNodesMaster, groupElementsMaster, setOfDOFSMaster);

    /////////////////////////////////////////////////////
    // ====> create FE slave mesh from gmsh (rectange) //
    /////////////////////////////////////////////////////

    auto LambdaGetMasterNodesLower = [](NuTo::NodeBase* rNodePtr) -> bool
    {
        double Tol = 1.e-6;
        if (rNodePtr->GetNum(NuTo::Node::eDof::COORDINATES)>0)
        {
            double y = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[1];
            if ((y >= -5. - Tol && y <= -5. + Tol))
            {
                return true;
            }
        }
        return false;
    };  // LambdaGetMasterNodesLower

    // DBC for the bottom of the master
    int groupNodesMasterLower = myStructure.GroupCreate(NuTo::eGroupId::Nodes);
    myStructure.GroupAddNodeFunction(groupNodesMasterLower, LambdaGetMasterNodesLower);

    NuTo::FullVector<int,Eigen::Dynamic> rMembersMasterDBC;
    myStructure.NodeGroupGetMembers(groupNodesMasterLower, rMembersMasterDBC);

    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> coordinatesMasterLower;
    myStructure.NodeGroupGetCoordinates(groupNodesMasterLower, coordinatesMasterLower);

    std::cout << coordinatesMasterLower << std::endl;
    int countDBC = 0;

    NuTo::FullVector<double,Eigen::Dynamic> direction(2);
    direction << 0, 1;
    for(int i = 0; i < rMembersMasterDBC.rows(); i++)
    {
        countDBC = myStructure.ConstraintLinearSetDisplacementNode(rMembersMasterDBC(i), direction, 0.0);
    }
    direction << 1, 0;
    countDBC = myStructure.ConstraintLinearSetDisplacementNode(rMembersMasterDBC(0), direction, 0.0);
    countDBC++;

    //////////////////////////
    // ===> build iga layer //
    //////////////////////////

    auto LambdaGetMasterNodesUpper = [](NuTo::NodeBase* rNodePtr) -> bool
    {
        double Tol = 1.e-6;
        if (rNodePtr->GetNum(NuTo::Node::eDof::COORDINATES)>0)
        {
            double y = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[1];
            if ((y >= 0. - Tol && y <= 0. + Tol))
            {
                return true;
            }
        }
        return false;
    };  // LambdaGetMasterNodesUpper

    int groupNodesMasterUpper = myStructure.GroupCreate(NuTo::eGroupId::Nodes);
    myStructure.GroupAddNodeFunction(groupNodesMasterUpper, LambdaGetMasterNodesUpper);
    NuTo::FullVector<int,Eigen::Dynamic> rMembersMasterContactBoundary;
    myStructure.NodeGroupGetMembers(groupNodesMasterUpper, rMembersMasterContactBoundary);

    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> coordinates;
    myStructure.NodeGroupGetCoordinates(groupNodesMasterUpper, coordinates);

    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> coordinatesAndIDs;
    coordinatesAndIDs.resize(coordinates.rows(), coordinates.cols() + 1);
    coordinatesAndIDs.block(0,0,coordinates.rows(), coordinates.cols()) = coordinates;

    for(int i = 0; i < coordinates.rows(); i++)
        coordinatesAndIDs(i,coordinates.cols()) = rMembersMasterContactBoundary(i);

    coordinatesAndIDs.SortRow(0);

    // create IGA curve
    int degree = 2;

    Eigen::MatrixXd A;

    NuTo::NURBSCurve curve(degree, coordinatesAndIDs.block(0, 0, coordinates.rows(), coordinates.cols()), A);

    for(int i = 0; i < 0; i++)
    {
        curve.DuplicateKnots();
    }

    std::set<NuTo::Node::eDof> setOfDOFS;
    setOfDOFS.insert(NuTo::Node::eDof::COORDINATES);
    setOfDOFS.insert(NuTo::Node::eDof::DISPLACEMENTS);

    int groupNodesIGAlayer    = myStructure.GroupCreate("Nodes");
    int groupElementsIGAlayer = myStructure.GroupCreate("Elements");

    Eigen::Matrix<std::pair<int, int>, Eigen::Dynamic, Eigen::Dynamic> elementsMaster = curve.buildIGAStructure(myStructure, setOfDOFS, groupElementsIGAlayer, groupNodesIGAlayer, "IGA1DLAYER");

    ///////////////////
    // ===> material //
    ///////////////////

    double Thickness = 1.;
    int section = myStructure.SectionCreate("PLANE_STRESS");
    myStructure.SectionSetThickness(section, Thickness);

    int constitutiveLaw = myStructure.ConstitutiveLawCreate("Linear_Elastic_Engineering_Stress");
    myStructure.ConstitutiveLawSetParameterDouble(constitutiveLaw, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, E);
    myStructure.ConstitutiveLawSetParameterDouble(constitutiveLaw, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, nue);
    myStructure.ConstitutiveLawSetParameterDouble(constitutiveLaw, NuTo::Constitutive::eConstitutiveParameter::DENSITY, rho);
    myStructure.ElementTotalSetSection(section);
    myStructure.ElementTotalSetConstitutiveLaw(constitutiveLaw);

    ////////////////////////////
    // ===> CONTACT ELEMENTS  //
    ////////////////////////////

    double rPenalty = 1.e9;
    myStructure.NuTo::Structure::ContactElementsCreate<2,1>(groupElementsSlaveLower, groupNodesSlaveLower, elementsMaster, NuTo::eIntegrationType::IntegrationType1D2NGauss12Ip, rPenalty, rContactAlgo);

    ////////////////////////////
    // ===> Constraints       //
    ////////////////////////////

    // coordinates
    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> coordinatesLayer;
    myStructure.NodeGroupGetCoordinates(groupNodesIGAlayer, coordinatesLayer);

    // IDs
    NuTo::FullVector<int,Eigen::Dynamic> rMembersMasterContactBoundaryLayer;
    myStructure.NodeGroupGetMembers(groupNodesIGAlayer, rMembersMasterContactBoundaryLayer);

    // all together
    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> coordinatesAndIDsLayer;
    coordinatesAndIDsLayer.resize(coordinatesLayer.rows(), coordinatesLayer.cols() + 1);
    coordinatesAndIDsLayer.block(0,0,coordinatesLayer.rows(), coordinatesLayer.cols()) = coordinatesLayer;

    for(int i = 0; i < coordinatesLayer.rows(); i++)
        coordinatesAndIDsLayer(i,coordinatesLayer.cols()) = rMembersMasterContactBoundaryLayer(i);

    coordinatesAndIDsLayer.SortRow(0);

    int dim = coordinatesAndIDs.cols() - 1;

    for(int nodeMaster = 0; nodeMaster < coordinatesAndIDs.rows(); nodeMaster++)
    {
        for(int dof = 0; dof < dim; dof++)
        {
            myStructure.ConstraintLinearEquationCreate (countDBC, coordinatesAndIDs(nodeMaster, dim), NuTo::Node::eDof::DISPLACEMENTS, dof,  1., 0.);
            for(int controlPoint = 0; controlPoint < coordinatesAndIDsLayer.rows(); controlPoint++)
            {
                myStructure.ConstraintLinearEquationAddTerm(countDBC, coordinatesAndIDsLayer(controlPoint, dim), NuTo::Node::eDof::DISPLACEMENTS, dof, -A(nodeMaster, controlPoint));
            }
            countDBC++;
        }
    }

    ///////////////////
    // ===> Solution //
    ///////////////////

    myStructure.NodeInfo(10);

    myStructure.CalculateMaximumIndependentSets();
    myStructure.NodeBuildGlobalDofs();

#ifdef ENABLE_VISUALIZE
    myStructure.AddVisualizationComponent(groupElementsIGAlayer, NuTo::eVisualizeWhat::DISPLACEMENTS);
    // master + slave fe only
    int visualizationGroup = myStructure.GroupUnion(groupElementsSlave, groupElementsMaster);
    myStructure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::DISPLACEMENTS);
    myStructure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ENGINEERING_STRAIN);
    myStructure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ENGINEERING_STRESS);
    // master + layer
    visualizationGroup = myStructure.GroupUnion(groupElementsIGAlayer, groupElementsMaster);
    myStructure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::DISPLACEMENTS);
#endif
}


int main(int argc, char* argv[])
{
    Eigen::VectorXd knots(7);
    knots << 0, 0, 0, 0.5, 1, 1, 1;

    Eigen::VectorXd weights(4);
    weights << 1, 1, 1, 1;

    Eigen::VectorXd shapeFunctions = NuTo::ShapeFunctionsIGA::BasisFunctionsAndDerivativesRat(2, 0.24, 2, 2, knots, weights);
    shapeFunctions = NuTo::ShapeFunctionsIGA::BasisFunctionsAndDerivativesRat(1, 0.24, 2, 2, knots, weights);
    shapeFunctions = NuTo::ShapeFunctionsIGA::BasisFunctionsAndDerivativesRat(0, 0.24, 2, 2, knots, weights);


    std::string path = "./";

    std::string fileNameSlave = "Slave.msh";
    std::string fileNameMaster = "Master.msh";

    std::string resultDir = "";
    double penalty;
    int    contactAlgo; //mortar = 0, non-mortar = 1

    /************************************************************************************************************/

    resultDir = "./ResultsContactNM";

    if (boost::filesystem::exists(resultDir))
    {
        if (boost::filesystem::is_directory(resultDir))
        {
            boost::filesystem::remove_all(resultDir);
        }
    }

    // create result directory
    boost::filesystem::create_directory(resultDir);

    penalty = 1.e8;
    contactAlgo = 1; //mortar = 0, non-mortar = 1
    dynamicExplicitContactTest(resultDir, path, fileNameSlave, fileNameMaster, penalty, contactAlgo, 4, 4, 4, 4);

    /************************************************************************************************************/

    resultDir = "./ResultsContactNMFine";

    if (boost::filesystem::exists(resultDir))
    {
        if (boost::filesystem::is_directory(resultDir))
        {
            boost::filesystem::remove_all(resultDir);
        }
    }

    // create result directory
    boost::filesystem::create_directory(resultDir);
    penalty = 1.e8;
    contactAlgo = 1; //mortar = 0, non-mortar = 1
    dynamicExplicitContactTest(resultDir, path, fileNameSlave, fileNameMaster, penalty, contactAlgo, 20, 20, 20, 20);
}
