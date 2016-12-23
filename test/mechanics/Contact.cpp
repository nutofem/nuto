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

#include "nuto/mechanics/dofSubMatrixStorage/BlockScalar.h"
#include "nuto/mechanics/dofSubMatrixStorage/BlockScalar.h"
#include "nuto/mechanics/dofSubMatrixStorage/BlockScalar.h"

#include "nuto/mechanics/timeIntegration/NewmarkDirect.h"

#include "nuto/mechanics/groups/GroupBase.h"

#define _USE_MATH_DEFINES
#include <cmath>

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
                                NuTo::eIntegrationType rIntegrationType,
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

    double E = 1.e5;
    double nue = 0.3;
    double rho = 1.;

    int constitutiveLaw = myStructure.ConstitutiveLawCreate("Linear_Elastic_Engineering_Stress");
    myStructure.ConstitutiveLawSetParameterDouble(constitutiveLaw, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, E);
    myStructure.ConstitutiveLawSetParameterDouble(constitutiveLaw, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, nue);
    myStructure.ConstitutiveLawSetParameterDouble(constitutiveLaw, NuTo::Constitutive::eConstitutiveParameter::DENSITY, rho);
    myStructure.ElementTotalSetSection(section);
    myStructure.ElementTotalSetConstitutiveLaw(constitutiveLaw);

    ////////////////////////////
    // ===> CONTACT ELEMENTS  //
    ////////////////////////////

    myStructure.NuTo::Structure::ContactElementsCreate<2,1>(groupElementsSlaveLower, groupNodesSlaveLower, elementsMaster, rIntegrationType, rPenalty, rContactAlgo);

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

void Elastic(const std::string &path,
             const std::string &fileNameMaster,
             const std::string &resultDir)
{
     NuTo::Structure *myStructure = new NuTo::Structure(2);
     myStructure->SetNumTimeDerivatives(0);
#ifdef _OPENMP
    int numThreads = 4;
    myStructure->SetNumProcessors(numThreads);
#endif

    double E = 1.e5;
    double nue = 0.3;
    double rho = 1.;

    //////////////////////////////
    // ====> create MASTER mesh //
    //////////////////////////////

    NuTo::FullMatrix<int, Eigen::Dynamic, Eigen::Dynamic> groupIndicesMaster = myStructure->ImportFromGmsh(path + fileNameMaster,
                                                                                                          NuTo::ElementData::eElementDataType::CONSTITUTIVELAWIP,
                                                                                                          NuTo::IpData::eIpDataType::NOIPDATA);

    int masterInterpolationType = myStructure->InterpolationTypeCreate(NuTo::Interpolation::eShapeType::QUAD2D);
    myStructure->InterpolationTypeAdd(masterInterpolationType, NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::LOBATTO2);
    myStructure->InterpolationTypeAdd(masterInterpolationType, NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::LOBATTO2);
    int masterElementsGroupId = groupIndicesMaster.GetValue(0, 0);
    myStructure->ElementGroupSetInterpolationType(masterElementsGroupId, masterInterpolationType);
    myStructure->ElementConvertToInterpolationType(masterElementsGroupId);

    int masterNodesGroupId = myStructure->GroupCreate(NuTo::eGroupId::Nodes);
    myStructure->GroupAddNodesFromElements(masterNodesGroupId, masterElementsGroupId);

    auto LambdaGetMasterNodesLower = [](NuTo::NodeBase* rNodePtr) -> bool
    {
        double Tol = 1.e-6;
        if (rNodePtr->GetNum(NuTo::Node::eDof::COORDINATES)>0)
        {
            double y = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[1];
            if ((y >= -1.3 - Tol && y <= -1.3 + Tol))
            {
                return true;
            }
        }
        return false;
    };  // LambdaGetMasterNodesLower

    // DBC for the bottom of the master
    int groupNodesMasterLower = myStructure->GroupCreate(NuTo::eGroupId::Nodes);
    myStructure->GroupAddNodeFunction(groupNodesMasterLower, masterNodesGroupId, LambdaGetMasterNodesLower);

    NuTo::FullVector<int,Eigen::Dynamic> rMembersMasterDBC;
    myStructure->NodeGroupGetMembers(groupNodesMasterLower, rMembersMasterDBC);

    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> coordinatesMasterLower;
    myStructure->NodeGroupGetCoordinates(groupNodesMasterLower, coordinatesMasterLower);

    std::cout << coordinatesMasterLower << std::endl;
    int countDBC(0);

    NuTo::FullVector<double,Eigen::Dynamic> direction(2);
    direction << 0, 1;
    for(int i = 0; i < rMembersMasterDBC.rows(); i++)
    {
        countDBC = myStructure->ConstraintLinearSetDisplacementNode(rMembersMasterDBC(i), direction, 0.0);
    }
    direction << 1, 0;
    countDBC = myStructure->ConstraintLinearSetDisplacementNode(rMembersMasterDBC(0), direction, 0.0);

    auto LambdaGetMasterNodesUpper = [](NuTo::NodeBase* rNodePtr) -> bool
    {
        double Tol = 1.e-6;
        if (rNodePtr->GetNum(NuTo::Node::eDof::COORDINATES)>0)
        {
            double y = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[1];
            if ((y >= -1. - Tol && y <= -1. + Tol))
            {
                return true;
            }
        }
        return false;
    };  // LambdaGetMasterNodesUpper

    int groupNodesMasterUpper = myStructure->GroupCreate(NuTo::eGroupId::Nodes);
    myStructure->GroupAddNodeFunction(groupNodesMasterUpper, masterNodesGroupId, LambdaGetMasterNodesUpper);
    NuTo::FullVector<int,Eigen::Dynamic> rMembersMasterContactBoundary;
    myStructure->NodeGroupGetMembers(groupNodesMasterUpper, rMembersMasterContactBoundary);

    int groupElementsMasterUpper = myStructure->GroupCreate(NuTo::eGroupId::Elements);
    myStructure->GroupAddElementsFromNodes(groupElementsMasterUpper, groupNodesMasterUpper, false);

    myStructure->SetNumLoadCases(1);
    double Stress = 1.e3;
    myStructure->LoadSurfacePressureCreate2D(0, groupElementsMasterUpper, groupNodesMasterUpper, Stress);

    ///////////////////
    // ===> material //
    ///////////////////

    double Thickness = 1.;
    int section = myStructure->SectionCreate("PLANE_STRESS");
    myStructure->SectionSetThickness(section, Thickness);

    int constitutiveLaw = myStructure->ConstitutiveLawCreate("Linear_Elastic_Engineering_Stress");
    myStructure->ConstitutiveLawSetParameterDouble(constitutiveLaw, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, E);
    myStructure->ConstitutiveLawSetParameterDouble(constitutiveLaw, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, nue);
    myStructure->ConstitutiveLawSetParameterDouble(constitutiveLaw, NuTo::Constitutive::eConstitutiveParameter::DENSITY, rho);
    myStructure->ElementTotalSetSection(section);
    myStructure->ElementTotalSetConstitutiveLaw(constitutiveLaw);

    ////////////////
    // ===> solve //
    ////////////////

    myStructure->CalculateMaximumIndependentSets();
    myStructure->NodeBuildGlobalDofs();

//    myStructure->SolveGlobalSystemStaticElastic(0);


    NuTo::BlockScalar tol(myStructure->GetDofStatus());
    NuTo::BlockScalar error(myStructure->GetDofStatus());
    tol.DefineDefaultValueToIninitializedDofTypes(1.e-10);
    error.DefineDefaultValueToIninitializedDofTypes(1.);
    myStructure->SolveGlobalSystemStaticElasticContact(tol, error, 10, 0);



//    NuTo::NewmarkDirect myIntegrationScheme(myStructure);
//    double timeStep = 1.;
//    double simulationTime = 1.;
//    myIntegrationScheme.SetTimeStep(timeStep);
//    myIntegrationScheme.SetPerformLineSearch(false);
//    myIntegrationScheme.SetResultDirectory(resultDir, false);
//    myIntegrationScheme.Solve(simulationTime);

#ifdef ENABLE_VISUALIZE

    //set result directory
    if (boost::filesystem::exists(resultDir))
    {
        if (boost::filesystem::is_directory(resultDir))
        {
            boost::filesystem::remove_all(resultDir);
        }
    }

    // create result directory
    boost::filesystem::create_directory(resultDir);

    myStructure->AddVisualizationComponent(masterElementsGroupId, NuTo::eVisualizeWhat::DISPLACEMENTS);
    myStructure->AddVisualizationComponent(masterElementsGroupId, NuTo::eVisualizeWhat::ENGINEERING_STRAIN);
    myStructure->AddVisualizationComponent(masterElementsGroupId, NuTo::eVisualizeWhat::ENGINEERING_STRESS);
    myStructure->ElementGroupExportVtkDataFile(masterElementsGroupId, resultDir+"/ElementsMasterNewton.vtu", true);
#endif

}


void ContactHertz(const std::string &path,
                  const std::string &fileNameSlave,
                  const std::string &fileNameMaster,
                  int rContactAlgo,
                  NuTo::Structure *myStructure,
                  double rPenalty,
                  NuTo::eIntegrationType rIntegrationType,
                  int &groupElementsIGAlayer, int &slaveElementsGroupId, int &masterElementsGroupId, int &groupElementsSlaveUpper, int &groupNodesSlaveUpper, int &groupElementsSlaveLower, int &groupElementsMasterUpper)
{
#ifdef _OPENMP
    int numThreads = 4;
    myStructure->SetNumProcessors(numThreads);
#endif

    double E = 1.e5;
    double nue = 0.3;
    double rho = 1.;

    ///////////////////////////////////////////////////////////////////////
    // ====> create SLAVE mesh from gmsh (smth. impacting a rectangle)   //
    ///////////////////////////////////////////////////////////////////////

    NuTo::FullMatrix<int, Eigen::Dynamic, Eigen::Dynamic> groupIndicesSlave = myStructure->ImportFromGmsh(path + fileNameSlave,
                                                                                                          NuTo::ElementData::eElementDataType::CONSTITUTIVELAWIP,
                                                                                                          NuTo::IpData::eIpDataType::NOIPDATA);

    int slaveInterpolationType = myStructure->InterpolationTypeCreate(NuTo::Interpolation::eShapeType::QUAD2D);
    myStructure->InterpolationTypeAdd(slaveInterpolationType, NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::LOBATTO2);
    myStructure->InterpolationTypeAdd(slaveInterpolationType, NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::LOBATTO2);
    slaveElementsGroupId = groupIndicesSlave.GetValue(0, 0);
    myStructure->ElementGroupSetInterpolationType(slaveElementsGroupId, slaveInterpolationType);
    myStructure->ElementConvertToInterpolationType(slaveElementsGroupId);

    int slaveNodesGroupId = myStructure->GroupCreate(NuTo::eGroupId::Nodes);
    myStructure->GroupAddNodesFromElements(slaveNodesGroupId, slaveElementsGroupId);

    // ===> build contact elements slave elements
    auto LambdaGetSlaveNodesLower = [](NuTo::NodeBase* rNodePtr) -> bool
    {
        double Tol = 1.e-4;
        //double angleincr = (15. * M_PI)/180.0;
        //double anglemin = 3.*M_PI/2. - angleincr;
        //double anglemax = 3.*M_PI/2. + angleincr;
        double radius = 1.;
        if ( rNodePtr->GetNum(NuTo::Node::eDof::COORDINATES)>0 )
        {
            double x = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[0];
            double y = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[1];
            double r = sqrt(x*x+y*y);
            //double angle = std::atan2(y,x);
            if (r <= radius + Tol && r >= radius - Tol)// && x >= -0.07 && x <= 0.07)//angle <= anglemax && angle >= anglemin)
            {
                return true;
            }
        }
        return false;
    };  // LambdaGetSlaveNodesLower

    myStructure->GroupGetNodesTotal();

    int groupNodesSlaveLower = myStructure->GroupCreate(NuTo::eGroupId::Nodes);
    myStructure->GroupAddNodeFunction(groupNodesSlaveLower, slaveNodesGroupId, LambdaGetSlaveNodesLower);
    NuTo::FullVector<int, Eigen::Dynamic> members = myStructure->GroupGetMemberIds(groupNodesSlaveLower);

    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> coordinatesSlaveC;
    myStructure->NodeGroupGetCoordinates(groupNodesSlaveLower, coordinatesSlaveC);

    groupElementsSlaveLower = myStructure->GroupCreate(NuTo::eGroupId::Elements);
    myStructure->GroupAddElementsFromNodes(groupElementsSlaveLower, groupNodesSlaveLower, false);
    members = myStructure->GroupGetMemberIds(groupElementsSlaveLower);

    auto LambdaGetSlaveNodesUpper = [](NuTo::NodeBase* rNodePtr) -> bool
    {
        double Tol = 1.e-6;
        if (rNodePtr->GetNum(NuTo::Node::eDof::COORDINATES)>0)
        {
            double y = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[1];
            if ((y >= -0.7 - Tol && y <= -0.7 + Tol))
            {
                return true;
            }
        }
        return false;
    };  // LambdaGetSlaveNodesUpper


    groupNodesSlaveUpper = myStructure->GroupCreate(NuTo::eGroupId::Nodes);
    myStructure->GroupAddNodeFunction(groupNodesSlaveUpper, slaveNodesGroupId, LambdaGetSlaveNodesUpper);
    members = myStructure->GroupGetMemberIds(groupNodesSlaveUpper);

    groupElementsSlaveUpper = myStructure->GroupCreate(NuTo::eGroupId::Elements);
    myStructure->GroupAddElementsFromNodes(groupElementsSlaveUpper, groupNodesSlaveUpper, false);
    members = myStructure->GroupGetMemberIds(groupElementsSlaveUpper);

    //////////////////////////////
    // ====> create MASTER mesh //
    //////////////////////////////

    NuTo::FullMatrix<int, Eigen::Dynamic, Eigen::Dynamic> groupIndicesMaster = myStructure->ImportFromGmsh(path + fileNameMaster,
                                                                                                          NuTo::ElementData::eElementDataType::CONSTITUTIVELAWIP,
                                                                                                          NuTo::IpData::eIpDataType::NOIPDATA);

    int masterInterpolationType = myStructure->InterpolationTypeCreate(NuTo::Interpolation::eShapeType::QUAD2D);
    myStructure->InterpolationTypeAdd(masterInterpolationType, NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::LOBATTO2);
    myStructure->InterpolationTypeAdd(masterInterpolationType, NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::LOBATTO2);
    masterElementsGroupId = groupIndicesMaster.GetValue(0, 0);
    myStructure->ElementGroupSetInterpolationType(masterElementsGroupId, masterInterpolationType);
    myStructure->ElementConvertToInterpolationType(masterElementsGroupId);

    int masterNodesGroupId = myStructure->GroupCreate(NuTo::eGroupId::Nodes);
    myStructure->GroupAddNodesFromElements(masterNodesGroupId, masterElementsGroupId);

    auto LambdaGetMasterNodesLower = [](NuTo::NodeBase* rNodePtr) -> bool
    {
        double Tol = 1.e-6;
        if (rNodePtr->GetNum(NuTo::Node::eDof::COORDINATES)>0)
        {
            double y = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[1];
            if ((y >= -1.3 - Tol && y <= -1.3 + Tol))
            {
                return true;
            }
        }
        return false;
    };  // LambdaGetMasterNodesLower

    // DBC for the bottom of the master
    int groupNodesMasterLower = myStructure->GroupCreate(NuTo::eGroupId::Nodes);
    myStructure->GroupAddNodeFunction(groupNodesMasterLower, masterNodesGroupId, LambdaGetMasterNodesLower);

    NuTo::FullVector<int,Eigen::Dynamic> rMembersMasterDBC;
    myStructure->NodeGroupGetMembers(groupNodesMasterLower, rMembersMasterDBC);

    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> coordinatesMasterLower;
    myStructure->NodeGroupGetCoordinates(groupNodesMasterLower, coordinatesMasterLower);

    int countDBC = 0;

    NuTo::FullVector<double,Eigen::Dynamic> direction(2);
    direction << 0, 1;
    for(int i = 0; i < rMembersMasterDBC.rows(); i++)
    {
        countDBC = myStructure->ConstraintLinearSetDisplacementNode(rMembersMasterDBC(i), direction, 0.0);
    }
    direction << 1, 0;
    countDBC = myStructure->ConstraintLinearSetDisplacementNode(rMembersMasterDBC(0), direction, 0.0);

    auto LambdaGetNodesLeft = [](NuTo::NodeBase* rNodePtr) -> bool
    {
        double Tol = 1.e-6;
        if (rNodePtr->GetNum(NuTo::Node::eDof::COORDINATES)>0)
        {
            double x = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[0];
            if ((x >= -Tol && x <= Tol))
            {
                return true;
            }
        }
        return false;
    };  // LambdaGetNodesLeft

//     DBC for the left side
    int groupNodesLeft = myStructure->GroupCreate(NuTo::eGroupId::Nodes);
    myStructure->GroupAddNodeFunction(groupNodesLeft, LambdaGetNodesLeft);

    NuTo::FullVector<int,Eigen::Dynamic> rMembersLeft;
    myStructure->NodeGroupGetMembers(groupNodesLeft, rMembersLeft);

    direction << 1, 0;
    for(int i = 0; i < rMembersLeft.rows(); i++)
    {
        countDBC = myStructure->ConstraintLinearSetDisplacementNode(rMembersLeft(i), direction, 0.0);
    }
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
            if ((y >= -1. - Tol && y <= -1. + Tol))
            {
                return true;
            }
        }
        return false;
    };  // LambdaGetMasterNodesUpper

    int groupNodesMasterUpper = myStructure->GroupCreate(NuTo::eGroupId::Nodes);
    myStructure->GroupAddNodeFunction(groupNodesMasterUpper, masterNodesGroupId, LambdaGetMasterNodesUpper);
    NuTo::FullVector<int,Eigen::Dynamic> rMembersMasterContactBoundary;
    myStructure->NodeGroupGetMembers(groupNodesMasterUpper, rMembersMasterContactBoundary);


    groupElementsMasterUpper = myStructure->GroupCreate(NuTo::eGroupId::Elements);
    myStructure->GroupAddElementsFromNodes(groupElementsMasterUpper, groupNodesMasterUpper, false);

    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> coordinates;
    myStructure->NodeGroupGetCoordinates(groupNodesMasterUpper, coordinates);

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

    int groupNodesIGAlayer = myStructure->GroupCreate("Nodes");
    groupElementsIGAlayer  = myStructure->GroupCreate("Elements");

    Eigen::Matrix<std::pair<int, int>, Eigen::Dynamic, Eigen::Dynamic> elementsMaster = curve.buildIGAStructure(*myStructure, setOfDOFS, groupElementsIGAlayer, groupNodesIGAlayer, "IGA1DLAYER");

    ///////////////////
    // ===> material //
    ///////////////////

    double Thickness = 1.;
    int section = myStructure->SectionCreate("PLANE_STRESS");
    myStructure->SectionSetThickness(section, Thickness);

    int constitutiveLaw = myStructure->ConstitutiveLawCreate("Linear_Elastic_Engineering_Stress");
    myStructure->ConstitutiveLawSetParameterDouble(constitutiveLaw, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, E);
    myStructure->ConstitutiveLawSetParameterDouble(constitutiveLaw, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, nue);
    myStructure->ConstitutiveLawSetParameterDouble(constitutiveLaw, NuTo::Constitutive::eConstitutiveParameter::DENSITY, rho);
    myStructure->ElementTotalSetSection(section);
    myStructure->ElementTotalSetConstitutiveLaw(constitutiveLaw);

    ////////////////////////////
    // ===> CONTACT ELEMENTS  //
    ////////////////////////////

    myStructure->NuTo::Structure::ContactElementsCreate<2,1>(groupElementsSlaveLower, groupNodesSlaveLower, elementsMaster, rIntegrationType, rPenalty, rContactAlgo);

    ////////////////////////////
    // ===> Constraints       //
    ////////////////////////////

    // coordinates
    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> coordinatesLayer;
    myStructure->NodeGroupGetCoordinates(groupNodesIGAlayer, coordinatesLayer);

    // IDs
    NuTo::FullVector<int,Eigen::Dynamic> rMembersMasterContactBoundaryLayer;
    myStructure->NodeGroupGetMembers(groupNodesIGAlayer, rMembersMasterContactBoundaryLayer);

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
            myStructure->ConstraintLinearEquationCreate (countDBC, coordinatesAndIDs(nodeMaster, dim), NuTo::Node::eDof::DISPLACEMENTS, dof,  1., 0.);
            for(int controlPoint = 0; controlPoint < coordinatesAndIDsLayer.rows(); controlPoint++)
            {
                myStructure->ConstraintLinearEquationAddTerm(countDBC, coordinatesAndIDsLayer(controlPoint, dim), NuTo::Node::eDof::DISPLACEMENTS, dof, -A(nodeMaster, controlPoint));
            }
            countDBC++;
        }
    }
}

void ContactHertzQuarterCircle(const std::string &path,
                               const std::string &fileNameSlave,
                               const std::string &fileNameMaster,
                               int rContactAlgo,
                               NuTo::Structure *myStructure,
                               double rPenalty,
                               double offset,
                               NuTo::eIntegrationType rIntegrationType,
                               int &groupElementsIGAlayer, int &slaveElementsGroupId, int &masterElementsGroupId, int &groupElementsSlaveUpper, int &groupNodesSlaveUpper, int &groupElementsSlaveLower, int &groupElementsMasterUpper)
{
#ifdef _OPENMP
    int numThreads = 4;
    myStructure->SetNumProcessors(numThreads);
#endif

    double E = 1.e5;
    double nue = 0.3;
    double rho = 1.;

    ///////////////////////////////////////////////////////////////////////
    // ====> create SLAVE mesh from gmsh (smth. impacting a rectangle)   //
    ///////////////////////////////////////////////////////////////////////

    NuTo::FullMatrix<int, Eigen::Dynamic, Eigen::Dynamic> groupIndicesSlave = myStructure->ImportFromGmsh(path + fileNameSlave,
                                                                                                          NuTo::ElementData::eElementDataType::CONSTITUTIVELAWIP,
                                                                                                          NuTo::IpData::eIpDataType::NOIPDATA);

    int slaveInterpolationType = myStructure->InterpolationTypeCreate(NuTo::Interpolation::eShapeType::QUAD2D);

    myStructure->InterpolationTypeAdd(slaveInterpolationType, NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::LOBATTO2);
    myStructure->InterpolationTypeAdd(slaveInterpolationType, NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::LOBATTO2);
    slaveElementsGroupId = groupIndicesSlave.GetValue(0, 0);
    myStructure->ElementGroupSetInterpolationType(slaveElementsGroupId, slaveInterpolationType);
    myStructure->ElementConvertToInterpolationType(slaveElementsGroupId, 1.e-6, 1.e-3);
    //myStructure->InterpolationTypeSetIntegrationType(slaveInterpolationType,NuTo::eIntegrationType::IntegrationType2D4NLobatto25Ip,NuTo::IpData::eIpDataType::NOIPDATA);

    int slaveNodesGroupId = myStructure->GroupCreate(NuTo::eGroupId::Nodes);
    myStructure->GroupAddNodesFromElements(slaveNodesGroupId, slaveElementsGroupId);

    // ===> build contact elements slave elements
    auto LambdaGetSlaveNodesLower = [offset](NuTo::NodeBase* rNodePtr) -> bool
    {
        double Tol = 1.e-4;
//        double angleincr = (15. * M_PI)/180.0;
//        double anglemin = 3.*M_PI/2. - angleincr;
//        double anglemax = 3.*M_PI/2. + angleincr;
        double radius = 1.;
        if ( rNodePtr->GetNum(NuTo::Node::eDof::COORDINATES)>0 )
        {
            double x = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[0];
            double y = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[1] - offset;
            double r = sqrt(x*x+y*y);
//            double angle = std::atan2(y,x);
            if (r <= radius + Tol && r >= radius - Tol && x >= 0.0 && x <= 0.12)//angle <= anglemax && angle >= anglemin)
            {
                return true;
            }
        }
        return false;
    };  // LambdaGetSlaveNodesLower

    myStructure->GroupGetNodesTotal();

    int groupNodesSlaveLower = myStructure->GroupCreate(NuTo::eGroupId::Nodes);
    myStructure->GroupAddNodeFunction(groupNodesSlaveLower, slaveNodesGroupId, LambdaGetSlaveNodesLower);
    NuTo::FullVector<int, Eigen::Dynamic> members = myStructure->GroupGetMemberIds(groupNodesSlaveLower);

    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> coordinatesSlaveC;
    myStructure->NodeGroupGetCoordinates(groupNodesSlaveLower, coordinatesSlaveC);

    groupElementsSlaveLower = myStructure->GroupCreate(NuTo::eGroupId::Elements);
    myStructure->GroupAddElementsFromNodes(groupElementsSlaveLower, groupNodesSlaveLower, false);
    members = myStructure->GroupGetMemberIds(groupElementsSlaveLower);

    auto LambdaGetSlaveNodesUpper = [offset](NuTo::NodeBase* rNodePtr) -> bool
    {
        double Tol = 1.e-6;
        if (rNodePtr->GetNum(NuTo::Node::eDof::COORDINATES)>0)
        {
            double y = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[1];
            if ((y >= 0. + offset - Tol && y <= 0. + offset + Tol))
            {
                return true;
            }
        }
        return false;
    };  // LambdaGetSlaveNodesUpper


    groupNodesSlaveUpper = myStructure->GroupCreate(NuTo::eGroupId::Nodes);
    myStructure->GroupAddNodeFunction(groupNodesSlaveUpper, slaveNodesGroupId, LambdaGetSlaveNodesUpper);
    members = myStructure->GroupGetMemberIds(groupNodesSlaveUpper);

    groupElementsSlaveUpper = myStructure->GroupCreate(NuTo::eGroupId::Elements);
    myStructure->GroupAddElementsFromNodes(groupElementsSlaveUpper, groupNodesSlaveUpper, false);
    members = myStructure->GroupGetMemberIds(groupElementsSlaveUpper);

    //////////////////////////////
    // ====> create MASTER mesh //
    //////////////////////////////

    NuTo::FullMatrix<int, Eigen::Dynamic, Eigen::Dynamic> groupIndicesMaster = myStructure->ImportFromGmsh(path + fileNameMaster,
                                                                                                          NuTo::ElementData::eElementDataType::CONSTITUTIVELAWIP,
                                                                                                          NuTo::IpData::eIpDataType::NOIPDATA);

    int masterInterpolationType = myStructure->InterpolationTypeCreate(NuTo::Interpolation::eShapeType::QUAD2D);

    myStructure->InterpolationTypeAdd(masterInterpolationType, NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::LOBATTO2);
    myStructure->InterpolationTypeAdd(masterInterpolationType, NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::LOBATTO3);
    masterElementsGroupId = groupIndicesMaster.GetValue(0, 0);
    myStructure->ElementGroupSetInterpolationType(masterElementsGroupId, masterInterpolationType);
    myStructure->ElementConvertToInterpolationType(masterElementsGroupId, 1.e-6, 1.e-3);
    //myStructure->InterpolationTypeSetIntegrationType(masterInterpolationType,NuTo::eIntegrationType::IntegrationType2D4NLobatto25Ip,NuTo::IpData::eIpDataType::STATICDATA);

    int masterNodesGroupId = myStructure->GroupCreate(NuTo::eGroupId::Nodes);
    myStructure->GroupAddNodesFromElements(masterNodesGroupId, masterElementsGroupId);

    auto LambdaGetMasterNodesLower = [](NuTo::NodeBase* rNodePtr) -> bool
    {
        double Tol = 1.e-6;
        if (rNodePtr->GetNum(NuTo::Node::eDof::COORDINATES)>0 && rNodePtr->IsDof(NuTo::Node::eDof::DISPLACEMENTS))
        {
            double y = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[1];
            if ((y >= -2. - Tol && y <= -2. + Tol))
            {
                return true;
            }
        }
        return false;
    };  // LambdaGetMasterNodesLower

    // DBC for the bottom of the master
    int groupNodesMasterLower = myStructure->GroupCreate(NuTo::eGroupId::Nodes);
    myStructure->GroupAddNodeFunction(groupNodesMasterLower, masterNodesGroupId, LambdaGetMasterNodesLower);

    NuTo::FullVector<int,Eigen::Dynamic> rMembersMasterDBC;
    myStructure->NodeGroupGetMembers(groupNodesMasterLower, rMembersMasterDBC);

    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> coordinatesMasterLower;
    myStructure->NodeGroupGetCoordinates(groupNodesMasterLower, coordinatesMasterLower);

    int countDBC = 0;

    NuTo::FullVector<double,Eigen::Dynamic> direction(2);
    direction << 0, 1;
    for(int i = 0; i < rMembersMasterDBC.rows(); i++)
    {
        countDBC = myStructure->ConstraintLinearSetDisplacementNode(rMembersMasterDBC(i), direction, 0.0);
    }
    direction << 1, 0;
    countDBC = myStructure->ConstraintLinearSetDisplacementNode(rMembersMasterDBC(0), direction, 0.0);

    auto LambdaGetNodesLeft = [](NuTo::NodeBase* rNodePtr) -> bool
    {
        double Tol = 1.e-6;
        if (rNodePtr->GetNum(NuTo::Node::eDof::COORDINATES)>0 && rNodePtr->IsDof(NuTo::Node::eDof::DISPLACEMENTS))
        {
            double x = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[0];
            if ((x >= -Tol && x <= Tol))
            {
                return true;
            }
        }
        return false;
    };  // LambdaGetNodesLeft

//     DBC for the left side
    int groupNodesLeft = myStructure->GroupCreate(NuTo::eGroupId::Nodes);
    myStructure->GroupAddNodeFunction(groupNodesLeft, LambdaGetNodesLeft);

    NuTo::FullVector<int,Eigen::Dynamic> rMembersLeft;
    myStructure->NodeGroupGetMembers(groupNodesLeft, rMembersLeft);

    direction << 1, 0;
    for(int i = 0; i < rMembersLeft.rows(); i++)
    {
        countDBC = myStructure->ConstraintLinearSetDisplacementNode(rMembersLeft(i), direction, 0.0);
    }
    countDBC++;

    //////////////////////////
    // ===> build iga layer //
    //////////////////////////

    auto LambdaGetMasterNodesUpper = [](NuTo::NodeBase* rNodePtr) -> bool
    {
        double Tol = 1.e-6;
        if (rNodePtr->GetNum(NuTo::Node::eDof::COORDINATES)>0 && rNodePtr->IsDof(NuTo::Node::eDof::DISPLACEMENTS))
        {
            double x = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[0];
            double y = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[1];
            if (y >= -1. - Tol && y <= -1. + Tol && x >= 0. - Tol && x <= 0.2 + Tol)
            {
                return true;
            }
        }
        return false;
    };  // LambdaGetMasterNodesUpper

    int groupNodesMasterUpper = myStructure->GroupCreate(NuTo::eGroupId::Nodes);
    myStructure->GroupAddNodeFunction(groupNodesMasterUpper, masterNodesGroupId, LambdaGetMasterNodesUpper);
    NuTo::FullVector<int,Eigen::Dynamic> rMembersMasterContactBoundary;
    myStructure->NodeGroupGetMembers(groupNodesMasterUpper, rMembersMasterContactBoundary);

    groupElementsMasterUpper = myStructure->GroupCreate(NuTo::eGroupId::Elements);
    myStructure->GroupAddElementsFromNodes(groupElementsMasterUpper, groupNodesMasterUpper, false);

    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> coordinates;
    myStructure->NodeGroupGetCoordinates(groupNodesMasterUpper, coordinates);

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

    int groupNodesIGAlayer = myStructure->GroupCreate("Nodes");
    groupElementsIGAlayer  = myStructure->GroupCreate("Elements");

    Eigen::Matrix<std::pair<int, int>, Eigen::Dynamic, Eigen::Dynamic> elementsMaster = curve.buildIGAStructure(*myStructure, setOfDOFS, groupElementsIGAlayer, groupNodesIGAlayer, "IGA1DLAYER");

    ///////////////////
    // ===> material //
    ///////////////////

    double Thickness = 1.;
    int section = myStructure->SectionCreate("PLANE_STRESS");
    myStructure->SectionSetThickness(section, Thickness);

    int constitutiveLaw = myStructure->ConstitutiveLawCreate("Linear_Elastic_Engineering_Stress");
    myStructure->ConstitutiveLawSetParameterDouble(constitutiveLaw, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, E);
    myStructure->ConstitutiveLawSetParameterDouble(constitutiveLaw, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, nue);
    myStructure->ConstitutiveLawSetParameterDouble(constitutiveLaw, NuTo::Constitutive::eConstitutiveParameter::DENSITY, rho);
    myStructure->ElementTotalSetSection(section);
    myStructure->ElementTotalSetConstitutiveLaw(constitutiveLaw);

    ////////////////////////////
    // ===> CONTACT ELEMENTS  //
    ////////////////////////////

    myStructure->NuTo::Structure::ContactElementsCreate<2,1>(groupElementsSlaveLower, groupNodesSlaveLower, elementsMaster, rIntegrationType, rPenalty, rContactAlgo);

    ////////////////////////////
    // ===> Constraints       //
    ////////////////////////////

    // coordinates
    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> coordinatesLayer;
    myStructure->NodeGroupGetCoordinates(groupNodesIGAlayer, coordinatesLayer);

    // IDs
    NuTo::FullVector<int,Eigen::Dynamic> rMembersMasterContactBoundaryLayer;
    myStructure->NodeGroupGetMembers(groupNodesIGAlayer, rMembersMasterContactBoundaryLayer);

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
            myStructure->ConstraintLinearEquationCreate (countDBC, coordinatesAndIDs(nodeMaster, dim), NuTo::Node::eDof::DISPLACEMENTS, dof,  1., 0.);
            for(int controlPoint = 0; controlPoint < coordinatesAndIDsLayer.rows(); controlPoint++)
            {
                myStructure->ConstraintLinearEquationAddTerm(countDBC, coordinatesAndIDsLayer(controlPoint, dim), NuTo::Node::eDof::DISPLACEMENTS, dof, -A(nodeMaster, controlPoint));
            }
            countDBC++;
        }
    }
}


void staticSolve(NuTo::Structure *myStructure, const std::string &resultDir, double offset, int groupElementsIGAlayer, int groupElementsSlave, int groupElementsMaster, int groupElementsSlaveUpper, int groupNodesSlaveUpper, int groupElementsSlaveLower, int groupElementsMasterUpper)
{
    ////////////////////////////
    // ===> Neumann BCs       //
    ////////////////////////////

    myStructure->SetNumLoadCases(1);
    double Stress = 10.;
    myStructure->LoadSurfacePressureCreate2D(0, groupElementsSlaveUpper, groupNodesSlaveUpper, Stress);

    /////////////////////
    // ===> DBCs       //
    /////////////////////

    double disp = -offset - 0.0008;

//    NuTo::FullVector<int,Eigen::Dynamic> rMembersSlaveDBC;
//    myStructure->NodeGroupGetMembers(groupNodesSlaveUpper, rMembersSlaveDBC);

//    NuTo::FullVector<double,Eigen::Dynamic> direction(2);
//    direction << 0, 1;
//    for(int i = 0; i < rMembersSlaveDBC.rows(); i++)
//    {
//        myStructure->ConstraintLinearSetDisplacementNode(rMembersSlaveDBC(i), direction, disp);
//    }

    // initial values
    disp = -offset - 0.0005;
    int slaveNodesGroupId = myStructure->GroupCreate(NuTo::eGroupId::Nodes);
    myStructure->GroupAddNodesFromElements(slaveNodesGroupId, groupElementsSlave, NuTo::Node::eDof::DISPLACEMENTS);
    NuTo::FullVector<double,Eigen::Dynamic>  dispVec(2);
    dispVec(0) =   0.;
    dispVec(1) = disp;
    myStructure->NodeGroupSetDisplacements(slaveNodesGroupId, 0, dispVec);

    myStructure->CalculateMaximumIndependentSets();
    myStructure->NodeBuildGlobalDofs();

//    NuTo::BlockScalar tol(myStructure->GetDofStatus());
//    NuTo::BlockScalar error(myStructure->GetDofStatus());
//    tol.DefineDefaultValueToIninitializedDofTypes(1.e-10);
//    error.DefineDefaultValueToIninitializedDofTypes(1.);
//    myStructure->SolveGlobalSystemStaticElasticContact(tol, error, 15, 0);

    NuTo::NewmarkDirect myIntegrationScheme(myStructure);
    double timeStep = 1.;
    double simulationTime = 1.;

    myIntegrationScheme.SetResultDirectory(resultDir, true);

    // speactral elements => #ip = #nodes
    std::vector<int> IPIds;
    IPIds.push_back(3);
    myIntegrationScheme.AddResultElementGroupIpData("ContactStress0",  groupElementsMasterUpper, 0, IPIds, NuTo::IpData::eIpStaticDataType::ENGINEERING_STRESS);
    myIntegrationScheme.AddResultElementGroupIpData("ContactStress1",  groupElementsMasterUpper, 1, IPIds, NuTo::IpData::eIpStaticDataType::ENGINEERING_STRESS);
    myIntegrationScheme.AddResultElementGroupIpData("ContactStress5",  groupElementsMasterUpper, 5, IPIds, NuTo::IpData::eIpStaticDataType::ENGINEERING_STRESS);

//    NuTo::FullMatrix<int, Eigen::Dynamic, Eigen::Dynamic> membersElementsMasterUpper =  myStructure->GroupGetMemberIds(groupElementsMasterUpper);
//    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> coordinates;

//    for(int i = 0; i < membersElementsMasterUpper.rows(); i++)
//    {
//        coordinates = myStructure->ElementGetIntegrationPointCoordinates(membersElementsMasterUpper(i));

//        std::string resultFileName2 = resultDir + "/Coordinates" + std::to_string(membersElementsMasterUpper(i)) + ".dat";
//        coordinates.WriteToFile(resultFileName2, "  ");
//    }

    myIntegrationScheme.SetMinTimeStepPlot(1.);
    myIntegrationScheme.SetLastTimePlot(0.);

    myIntegrationScheme.SetToleranceForce(1.e-10);
    myIntegrationScheme.SetMaxNumIterations(50);
    myIntegrationScheme.SetTimeStep(timeStep);
    myIntegrationScheme.SetPerformLineSearch(false);
    myIntegrationScheme.Solve(simulationTime);

#ifdef ENABLE_VISUALIZE
    //set result directory
//    if (boost::filesystem::exists(resultDir))
//    {
//        if (boost::filesystem::is_directory(resultDir))
//        {
//            boost::filesystem::remove_all(resultDir);
//        }
//    }

//    // create result directory
//    boost::filesystem::create_directory(resultDir);

    myStructure->AddVisualizationComponent(groupElementsIGAlayer, NuTo::eVisualizeWhat::DISPLACEMENTS);
    myStructure->ElementGroupExportVtkDataFile(groupElementsIGAlayer, resultDir+"/ElementsLayer.vtu", true);
    // master + slave fe only
    myStructure->AddVisualizationComponent(groupElementsSlave, NuTo::eVisualizeWhat::DISPLACEMENTS);
    myStructure->AddVisualizationComponent(groupElementsSlave, NuTo::eVisualizeWhat::ENGINEERING_STRAIN);
    myStructure->AddVisualizationComponent(groupElementsSlave, NuTo::eVisualizeWhat::ENGINEERING_STRESS);
    myStructure->ElementGroupExportVtkDataFile(groupElementsSlave, resultDir+"/ElementsSlave.vtu", true);
    // master + slave fe only
    int visualizationGroup = myStructure->GroupUnion(groupElementsSlave, groupElementsMaster);
    myStructure->AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::DISPLACEMENTS);
    myStructure->AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ENGINEERING_STRAIN);
    myStructure->AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ENGINEERING_STRESS);
    myStructure->ElementGroupExportVtkDataFile(visualizationGroup, resultDir+"/ElementsSlaveMaster.vtu", true);
    // master + layer
    visualizationGroup = myStructure->GroupUnion(groupElementsIGAlayer, groupElementsMaster);
    myStructure->AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::DISPLACEMENTS);
    myStructure->ElementGroupExportVtkDataFile(visualizationGroup, resultDir+"/ElementsMaster.vtu", true);
    // slave contact
    myStructure->AddVisualizationComponent(groupElementsSlaveLower, NuTo::eVisualizeWhat::DISPLACEMENTS);
    myStructure->ElementGroupExportVtkDataFile(groupElementsSlaveLower, resultDir+"/ElementsSlaveContact.vtu", true);
#endif
}

void dynamicSolve(NuTo::Structure *myStructure, const std::string &resultDir, double offset, int groupElementsIGAlayer, int groupElementsSlave, int groupElementsMaster, int groupElementsSlaveLower)
{
    int slaveNodesGroupId = myStructure->GroupCreate(NuTo::eGroupId::Nodes);
    myStructure->GroupAddNodesFromElements(slaveNodesGroupId, groupElementsSlave);

    NuTo::FullVector<double,Eigen::Dynamic>  velocities(2);
    velocities(0) =   0.;
    velocities(1) = -10.;
    myStructure->NodeGroupSetDisplacements(slaveNodesGroupId, 1, velocities);

    double disp = -offset;

    //  initial values
    myStructure->GroupAddNodesFromElements(slaveNodesGroupId, groupElementsSlave);
    NuTo::FullVector<double,Eigen::Dynamic>  dispVec(2);
    dispVec(0) =   0.;
    dispVec(1) = disp;
    myStructure->NodeGroupSetDisplacements(slaveNodesGroupId, 0, dispVec);


    myStructure->CalculateMaximumIndependentSets();
    myStructure->NodeBuildGlobalDofs();


    NuTo::RungeKutta4 myIntegrationScheme(myStructure);

    double simulationTime(0.0005);

//    myIntegrationScheme.SetMinTimeStep(0.0001 * myIntegrationScheme.GetMaxTimeStep());
    double h = 1.e-6;
    myIntegrationScheme.SetTimeStep(h);
    int numOutputs = 700;
    myIntegrationScheme.SetMinTimeStepPlot(simulationTime/numOutputs);

    //set output during the simulation to false
    myStructure->SetShowTime(false);

    myIntegrationScheme.AddResultTime("Time");

    //set result directory
    bool deleteResultDirectoryFirst(true);
    myIntegrationScheme.SetResultDirectory(resultDir, deleteResultDirectoryFirst);

#ifdef ENABLE_VISUALIZE
    myStructure->AddVisualizationComponent(groupElementsIGAlayer, NuTo::eVisualizeWhat::DISPLACEMENTS);
    // slave contact elements
    myStructure->AddVisualizationComponent(groupElementsSlaveLower, NuTo::eVisualizeWhat::DISPLACEMENTS);
    myStructure->AddVisualizationComponent(groupElementsSlaveLower, NuTo::eVisualizeWhat::ENGINEERING_STRAIN);
    myStructure->AddVisualizationComponent(groupElementsSlaveLower, NuTo::eVisualizeWhat::ENGINEERING_STRESS);
    // master + slave fe only
    int visualizationGroup = myStructure->GroupUnion(groupElementsSlave, groupElementsMaster);
    myStructure->AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::DISPLACEMENTS);
    myStructure->AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ENGINEERING_STRAIN);
    myStructure->AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ENGINEERING_STRESS);
    // master
    myStructure->AddVisualizationComponent(groupElementsMaster, NuTo::eVisualizeWhat::DISPLACEMENTS);
#endif

    //solve (perform Newton raphson iteration
    myIntegrationScheme.Solve(simulationTime);
}

void testInterpolationTypes(const std::string &path,
                            const std::string &resultDir,
                            const std::string &fileNameMaster)
{
    NuTo::Structure *myStructure = new NuTo::Structure(2);

    NuTo::FullMatrix<int, Eigen::Dynamic, Eigen::Dynamic> groupIndicesMaster = myStructure->ImportFromGmsh(path + fileNameMaster,
                                                                                                          NuTo::ElementData::eElementDataType::CONSTITUTIVELAWIP,
                                                                                                          NuTo::IpData::eIpDataType::NOIPDATA);

    int masterInterpolationType = myStructure->InterpolationTypeCreate(NuTo::Interpolation::eShapeType::QUAD2D);

    myStructure->InterpolationTypeAdd(masterInterpolationType, NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    myStructure->InterpolationTypeAdd(masterInterpolationType, NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    int masterElementsGroupId = groupIndicesMaster.GetValue(0, 0);
    myStructure->ElementGroupSetInterpolationType(masterElementsGroupId, masterInterpolationType);
    myStructure->ElementConvertToInterpolationType(masterElementsGroupId, 1.e-6, 1.e-3);
    //myStructure->InterpolationTypeSetIntegrationType(masterInterpolationType,NuTo::eIntegrationType::IntegrationType2D4NLobatto25Ip,NuTo::IpData::eIpDataType::STATICDATA);

    double Thickness = 1.;
    int section = myStructure->SectionCreate("PLANE_STRESS");
    myStructure->SectionSetThickness(section, Thickness);

    int constitutiveLaw = myStructure->ConstitutiveLawCreate("Linear_Elastic_Engineering_Stress");
    double E = 1.e5;
    double nue = 0.3;
    double rho = 1.;
    myStructure->ConstitutiveLawSetParameterDouble(constitutiveLaw, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, E);
    myStructure->ConstitutiveLawSetParameterDouble(constitutiveLaw, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, nue);
    myStructure->ConstitutiveLawSetParameterDouble(constitutiveLaw, NuTo::Constitutive::eConstitutiveParameter::DENSITY, rho);
    myStructure->ElementTotalSetSection(section);
    myStructure->ElementTotalSetConstitutiveLaw(constitutiveLaw);

    int masterNodesGroupId = myStructure->GroupCreate(NuTo::eGroupId::Nodes);
    myStructure->GroupAddNodesFromElements(masterNodesGroupId, masterElementsGroupId);

    auto LambdaGetMasterNodesLower = [](NuTo::NodeBase* rNodePtr) -> bool
    {
        double Tol = 1.e-6;
        if (rNodePtr->GetNum(NuTo::Node::eDof::COORDINATES)>0 && rNodePtr->IsDof(NuTo::Node::eDof::DISPLACEMENTS))
        {
            double y = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[1];
            if ((y >= -2. - Tol && y <= -2. + Tol))
            {
                return true;
            }
        }
        return false;
    };  // LambdaGetMasterNodesLower

    // DBC for the bottom of the master
    int groupNodesMasterLower = myStructure->GroupCreate(NuTo::eGroupId::Nodes);
    myStructure->GroupAddNodeFunction(groupNodesMasterLower, LambdaGetMasterNodesLower);

    NuTo::FullVector<int,Eigen::Dynamic> rMembers;
    myStructure->NodeGroupGetMembers(groupNodesMasterLower, rMembers);

    int countDBC = 0;

    NuTo::FullVector<double,Eigen::Dynamic> direction(2);
    direction << 0, 1;
    for(int i = 0; i < rMembers.rows(); i++)
    {
        countDBC = myStructure->ConstraintLinearSetDisplacementNode(rMembers(i), direction, 0.0);
    }
    direction << 1, 0;
    countDBC = myStructure->ConstraintLinearSetDisplacementNode(rMembers(0), direction, 0.0);


    auto LambdaGetMasterNodesUpper = [](NuTo::NodeBase* rNodePtr) -> bool
    {
        double Tol = 1.e-6;
        if (rNodePtr->GetNum(NuTo::Node::eDof::COORDINATES)>0 && rNodePtr->IsDof(NuTo::Node::eDof::DISPLACEMENTS))
        {
            double y = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[1];
            if (y >= -1. - Tol && y <= -1. + Tol)
            {
                return true;
            }
        }
        return false;
    };  // LambdaGetMasterNodesUpper

    int groupNodesMasterUpper = myStructure->GroupCreate(NuTo::eGroupId::Nodes);

    myStructure->GroupAddNodeFunction(groupNodesMasterUpper, LambdaGetMasterNodesUpper);

    myStructure->NodeGroupGetMembers(groupNodesMasterUpper, rMembers);

    int groupElementsMasterUpper = myStructure->GroupCreate(NuTo::eGroupId::Elements);
    myStructure->GroupAddElementsFromNodes(groupElementsMasterUpper, groupNodesMasterUpper, false);

    myStructure->ElementGroupGetMembers(groupElementsMasterUpper, rMembers);

    // Pressure
    myStructure->SetNumLoadCases(1);
    double Stress = 10.;
    myStructure->LoadSurfacePressureCreate2D(0, groupElementsMasterUpper, groupNodesMasterUpper, Stress);

    myStructure->Info();

    myStructure->CalculateMaximumIndependentSets();
    myStructure->NodeBuildGlobalDofs();

#ifdef ENABLE_VISUALIZE

    //set result directory
    if (boost::filesystem::exists(resultDir))
    {
        if (boost::filesystem::is_directory(resultDir))
        {
            boost::filesystem::remove_all(resultDir);
        }
    }

    // create result directory
    boost::filesystem::create_directory(resultDir);

    myStructure->AddVisualizationComponent(masterElementsGroupId, NuTo::eVisualizeWhat::DISPLACEMENTS);
    myStructure->AddVisualizationComponent(masterElementsGroupId, NuTo::eVisualizeWhat::ENGINEERING_STRAIN);
    myStructure->AddVisualizationComponent(masterElementsGroupId, NuTo::eVisualizeWhat::ENGINEERING_STRESS);
    myStructure->ElementGroupExportVtkDataFile(masterElementsGroupId, resultDir+"/ElementsMaster.vtu", true);

    myStructure->ExportVtkDataFileElements(resultDir+"/Elements.vtu", true);
#endif

    myStructure->SolveGlobalSystemStaticElastic();

//    NuTo::NewmarkDirect myIntegrationScheme(myStructure);
//    double timeStep = 1.;
//    double simulationTime = 1.;

//    myIntegrationScheme.SetResultDirectory(resultDir, true);

//    myIntegrationScheme.SetMinTimeStepPlot(1.);
//    myIntegrationScheme.SetLastTimePlot(0.);

//    myIntegrationScheme.SetToleranceForce(1.e-10);
//    myIntegrationScheme.SetMaxNumIterations(50);
//    myIntegrationScheme.SetTimeStep(timeStep);
//    myIntegrationScheme.SetPerformLineSearch(false);
//    myIntegrationScheme.Solve(simulationTime);

}

int main(int argc, char* argv[])
{
    std::string path           = "/home/potto1/mechanics/otto/gmsh/";
    std::string fileNameSlave  = "";
    std::string fileNameMaster = "";

    std::string resultDir = "";
    double offset = 0.;//1.e-1;
    double penalty;
    int    contactAlgo; //mortar = 0, non-mortar = 1
    int groupElementsIGAlayer(0), groupElementsSlave(0), groupElementsMaster(0), groupElementsSlaveUpper(0), groupNodesSlaveUpper(0), groupElementsSlaveLower(0), groupElementsMasterUpper(0);

    NuTo::Structure *myStructure = NULL;

    /*********************************************static - penalty parameter study******************************************************/
    fileNameSlave  = "CircleUnstructuredQudrilat.msh";
    fileNameMaster = "masterHalfspaceCombined.msh";

    resultDir      = "./ResultsStatic";
    testInterpolationTypes(path, resultDir, fileNameMaster);

    return 0;

    contactAlgo = 1;
    penalty = 1.e10;
    resultDir      = "./ResultsContactStatic_e10";

    myStructure = new NuTo::Structure(2);

    myStructure->SetNumTimeDerivatives(0);
    ContactHertzQuarterCircle(path,
                              fileNameSlave,
                              fileNameMaster,
                              contactAlgo,
                              myStructure,
                              penalty,
                              offset,
                              NuTo::eIntegrationType::IntegrationType1D2NGauss12Ip,
                              groupElementsIGAlayer, groupElementsSlave, groupElementsMaster, groupElementsSlaveUpper, groupNodesSlaveUpper, groupElementsSlaveLower, groupElementsMasterUpper);

    myStructure->SetNumTimeDerivatives(0);
    staticSolve(myStructure, resultDir, offset, groupElementsIGAlayer, groupElementsSlave, groupElementsMaster, groupElementsSlaveUpper, groupNodesSlaveUpper, groupElementsSlaveLower, groupElementsMasterUpper);

    delete myStructure;

        return 0;

//------------------------------------------------------------------------------------------------------------------------


    penalty = 5.e10;
    resultDir      = "./ResultsContactStatic_5_e10";

    myStructure = new NuTo::Structure(2);

    myStructure->SetNumTimeDerivatives(0);
    ContactHertzQuarterCircle(path,
                              fileNameSlave,
                              fileNameMaster,
                              contactAlgo,
                              myStructure,
                              penalty,
                              offset,
                              NuTo::eIntegrationType::IntegrationType1D2NGauss12Ip,
                              groupElementsIGAlayer, groupElementsSlave, groupElementsMaster, groupElementsSlaveUpper, groupNodesSlaveUpper, groupElementsSlaveLower, groupElementsMasterUpper);

    myStructure->SetNumTimeDerivatives(0);
    staticSolve(myStructure, resultDir, offset,
                groupElementsIGAlayer, groupElementsSlave, groupElementsMaster, groupElementsSlaveUpper, groupNodesSlaveUpper, groupElementsSlaveLower, groupElementsMasterUpper);
    delete myStructure;

//------------------------------------------------------------------------------------------------------------------------


    penalty = 1.e11;
    resultDir      = "./ResultsContactStatic_e11";

    myStructure = new NuTo::Structure(2);

    myStructure->SetNumTimeDerivatives(0);
    ContactHertzQuarterCircle(path,
                              fileNameSlave,
                              fileNameMaster,
                              contactAlgo,
                              myStructure,
                              penalty,
                              offset,
                              NuTo::eIntegrationType::IntegrationType1D2NGauss12Ip,
                              groupElementsIGAlayer, groupElementsSlave, groupElementsMaster, groupElementsSlaveUpper, groupNodesSlaveUpper, groupElementsSlaveLower, groupElementsMasterUpper);

    myStructure->SetNumTimeDerivatives(0);
    staticSolve(myStructure, resultDir, offset,
                groupElementsIGAlayer, groupElementsSlave, groupElementsMaster, groupElementsSlaveUpper, groupNodesSlaveUpper, groupElementsSlaveLower, groupElementsMasterUpper);

    delete myStructure;

//------------------------------------------------------------------------------------------------------------------------


    penalty = 5.e11;
    resultDir      = "./ResultsContactStatic_5_e11";

    myStructure = new NuTo::Structure(2);

    myStructure->SetNumTimeDerivatives(0);
    ContactHertzQuarterCircle(path,
                              fileNameSlave,
                              fileNameMaster,
                              contactAlgo,
                              myStructure,
                              penalty,
                              offset,
                              NuTo::eIntegrationType::IntegrationType1D2NGauss12Ip,
                              groupElementsIGAlayer, groupElementsSlave, groupElementsMaster, groupElementsSlaveUpper, groupNodesSlaveUpper, groupElementsSlaveLower, groupElementsMasterUpper);

    myStructure->SetNumTimeDerivatives(0);
    staticSolve(myStructure, resultDir, offset,
                groupElementsIGAlayer, groupElementsSlave, groupElementsMaster, groupElementsSlaveUpper, groupNodesSlaveUpper, groupElementsSlaveLower, groupElementsMasterUpper);

    delete myStructure;

    /*********************************************static - mesh study******************************************************/
    penalty = 1.e11;
    resultDir      = "./ResultsContactStatic_e11_Fine";
    fileNameSlave  = "CircleUnstructuredQudrilat.msh";
    fileNameMaster = "masterHalfspaceCombined.msh";

    myStructure = new NuTo::Structure(2);

    myStructure->SetNumTimeDerivatives(0);
    ContactHertzQuarterCircle(path,
                              fileNameSlave,
                              fileNameMaster,
                              contactAlgo,
                              myStructure,
                              penalty,
                              offset,
                              NuTo::eIntegrationType::IntegrationType1D2NGauss12Ip,
                              groupElementsIGAlayer, groupElementsSlave, groupElementsMaster, groupElementsSlaveUpper, groupNodesSlaveUpper, groupElementsSlaveLower, groupElementsMasterUpper);

    myStructure->SetNumTimeDerivatives(0);
    staticSolve(myStructure, resultDir, offset,
                groupElementsIGAlayer, groupElementsSlave, groupElementsMaster, groupElementsSlaveUpper, groupNodesSlaveUpper, groupElementsSlaveLower, groupElementsMasterUpper);
    delete myStructure;


    return 0;

    /*********************************************dynamic*****************************************************/

    myStructure->SetNumTimeDerivatives(2);
    dynamicSolve(myStructure, resultDir, offset, groupElementsIGAlayer, groupElementsSlave, groupElementsMaster, groupElementsSlaveLower);

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
//    dynamicExplicitContactTest(resultDir, path, fileNameSlave, fileNameMaster, penalty, NuTo::eIntegrationType::IntegrationType1D2NLobatto5Ip, contactAlgo, 4, 4, 4, 4);

    /*********************************************dynamic*****************************************************/

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
//    dynamicExplicitContactTest(resultDir, path, fileNameSlave, fileNameMaster, penalty, NuTo::eIntegrationType::IntegrationType1D2NLobatto5Ip, contactAlgo, 20, 20, 20, 20);
}
