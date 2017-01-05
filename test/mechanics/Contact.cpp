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
#include "nuto/mechanics/elements/ElementBase.h"
#include "nuto/mechanics/elements/ElementDataBase.h"

#include "nuto/mechanics/elements/IpDataEnum.h"
#include "nuto/mechanics/MechanicsException.h"

#define _USE_MATH_DEFINES
#include <cmath>

#ifdef ENABLE_VISUALIZE
#include "nuto/visualize/VisualizeEnum.h"
#endif

#define PRINTRESULT true

void ContactHertzQuarterCircle(const std::string &path,
                               const std::string &fileNameSlave,
                               const std::string &fileNameMaster,
                               int rContactAlgo,
                               NuTo::Structure *myStructure,
                               double rPenalty,
                               double rGapMax,
                               NuTo::eIntegrationType rIntegrationType,
                               NuTo::Interpolation::eTypeOrder rElementTypeIdentDisps,
                               NuTo::eIntegrationType rIntegrationTypeDisps,
                               int &groupElementsIGAlayer,
                               int &slaveElementsGroupId,
                               int &masterElementsGroupId,
                               int &groupElementsSlaveUpper,
                               int &groupNodesSlaveUpper,
                               int &groupElementsSlaveLower,
                               int &groupElementsMasterUpper,
                               int &slaveContactElementsGroupId,
                               int &masterContactElementsGroupId)
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
    myStructure->InterpolationTypeAdd(slaveInterpolationType, NuTo::Node::eDof::DISPLACEMENTS, rElementTypeIdentDisps);
    slaveElementsGroupId = groupIndicesSlave.GetValue(0, 0);
    myStructure->ElementGroupSetInterpolationType(slaveElementsGroupId, slaveInterpolationType);
    myStructure->ElementConvertToInterpolationType(slaveElementsGroupId);
    myStructure->InterpolationTypeSetIntegrationType(slaveInterpolationType,rIntegrationTypeDisps,NuTo::IpData::eIpDataType::STATICDATA);

    int slaveNodesGroupId = myStructure->GroupCreate(NuTo::eGroupId::Nodes);
    myStructure->GroupAddNodesFromElements(slaveNodesGroupId, slaveElementsGroupId,  NuTo::Node::eDof::DISPLACEMENTS);

    // ===> build contact elements slave elements
    auto LambdaGetSlaveNodesLower = [](NuTo::NodeBase* rNodePtr) -> bool
    {
        double Tol = 1.e-4;
        //double angleincr = (15. * M_PI)/180.0;
        //double anglemax = 3.*M_PI/2. + angleincr;
        double radius = 1.;
        if ( rNodePtr->GetNum(NuTo::Node::eDof::DISPLACEMENTS)>0 )
        {
            double x = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[0];
            double y = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[1];
            double r = sqrt(x*x+y*y);
            //double angle = std::atan2(y,x);
            if (r <= radius + Tol && r >= radius - Tol && x >= 0.0 && x <= 0.125)//angle <= anglemax)
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
        if (rNodePtr->GetNum(NuTo::Node::eDof::DISPLACEMENTS)>0)
        {
            double y = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[1];
            if ((y >= 0. - Tol && y <= 0. + Tol))
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
    myStructure->InterpolationTypeAdd(masterInterpolationType, NuTo::Node::eDof::DISPLACEMENTS, rElementTypeIdentDisps);
    masterElementsGroupId = groupIndicesMaster.GetValue(0, 0);
    myStructure->ElementGroupSetInterpolationType(masterElementsGroupId, masterInterpolationType);
    myStructure->ElementConvertToInterpolationType(masterElementsGroupId);
    myStructure->InterpolationTypeSetIntegrationType(masterInterpolationType, rIntegrationTypeDisps, NuTo::IpData::eIpDataType::STATICDATA);

    int masterNodesGroupId = myStructure->GroupCreate(NuTo::eGroupId::Nodes);
    myStructure->GroupAddNodesFromElements(masterNodesGroupId, masterElementsGroupId,  NuTo::Node::eDof::DISPLACEMENTS);

    auto LambdaGetMasterNodesLower = [](NuTo::NodeBase* rNodePtr) -> bool
    {
        double Tol = 1.e-6;
        if (rNodePtr->GetNum(NuTo::Node::eDof::DISPLACEMENTS)>0)
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
        if (rNodePtr->GetNum(NuTo::Node::eDof::DISPLACEMENTS) > 0)
        {
            double x = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[0];
            if (x >= -Tol && x <= Tol)
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
        if (rNodePtr->GetNum(NuTo::Node::eDof::DISPLACEMENTS)>0)
        {
            double x = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[0];
            double y = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[1];
            if (y >= -1. - Tol && y <= -1. + Tol && x >= 0. - Tol && x <= 0.125 + Tol)
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
//    std::cout << coordinatesAndIDs << std::endl;

    // create IGA curve
    int degree = 3;

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

    int constitutiveLawPC = myStructure->ConstitutiveLawCreate("Contact_Constitutive_Law");
    std::function<double(double)> constitutiveContactLaw  =
    [rPenalty](double rGap) -> double
    {
        if(rGap<0)
            return rPenalty*rGap;
        else
            return 0.;
    };

    std::function<double(double)> constitutiveContactLawDerivative =
    [rPenalty](double rGap) -> double
    {
        if(rGap<=0)
            return rPenalty;
        else
            return 0.;
    };

    double a = rPenalty/(2*rGapMax);
    std::function<double(double)> constitutiveContactLawNLin =
    [a](double rGap) -> double
    {
        if(rGap<=0)
            return a*rGap*rGap;
        else
            return 0.;
    };

    std::function<double(double)> constitutiveContactLawDerivativeNLin =
    [a](double rGap) -> double
    {
        if(rGap<=0)
            return 2*a*rGap;
        else
            return 0.;
    };


    myStructure->ConstitutiveLawSetParameterFunction(constitutiveLawPC, NuTo::Constitutive::eConstitutiveParameter::CONSTITUTIVE_LAW_FUNCTION, constitutiveContactLaw);
    myStructure->ConstitutiveLawSetParameterFunction(constitutiveLawPC, NuTo::Constitutive::eConstitutiveParameter::CONSTITUTIVE_LAW_DERIVATIVE_FUNCTION, constitutiveContactLawDerivative);

    slaveContactElementsGroupId = myStructure->NuTo::Structure::ContactElementsCreate<2,1>(groupElementsSlaveLower,
                                                                                           groupNodesSlaveLower,
                                                                                           elementsMaster,
                                                                                           rIntegrationType,
                                                                                           rContactAlgo,
                                                                                           constitutiveLawPC);

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

    // create boundary elements for visualizing the results
    masterContactElementsGroupId = myStructure->BoundaryElementsCreate(groupElementsMasterUpper, groupNodesMasterUpper, NULL);
    NuTo::FullVector<int,Eigen::Dynamic> rMembersMasterContactBoundaryVisualize;
    myStructure->ElementGroupGetMembers(masterContactElementsGroupId, rMembersMasterContactBoundaryVisualize);
    for(int i = 0; i < rMembersMasterContactBoundaryVisualize.rows(); i++)
    {
        NuTo::ElementBase* elementPtr = myStructure->ElementGetElementPtr(rMembersMasterContactBoundaryVisualize(i));
        NuTo::IpData::eIpDataType ipDataType = elementPtr->GetIpDataType(0);
        elementPtr->SetIntegrationType(myStructure->GetPtrIntegrationType(rIntegrationType), ipDataType);
    }
}


void staticSolve(NuTo::Structure *myStructure,
                 const std::string &resultDir,
                 NuTo::eIntegrationType rIntegrationType,
                 int groupElementsIGAlayer,
                 int groupElementsSlave,
                 int groupElementsMaster,
                 int groupElementsSlaveUpper,
                 int groupNodesSlaveUpper,
                 int groupElementsSlaveLower,
                 int groupElementsMasterUpper,
                 int slaveContactElementsGroupId,
                 int masterContactElementsGroupId)
{
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

    // ===> Neumann BCs
    myStructure->SetNumLoadCases(1);
    double Stress = 10.;
    myStructure->LoadSurfacePressureCreate2D(0, groupElementsSlaveUpper, groupNodesSlaveUpper, Stress);

    double disp = 0.;

    // ===> DBCs
//    double disp = - 0.0008;

//    NuTo::FullVector<int,Eigen::Dynamic> rMembersSlaveDBC;
//    myStructure->NodeGroupGetMembers(groupNodesSlaveUpper, rMembersSlaveDBC);

//    NuTo::FullVector<double,Eigen::Dynamic> direction(2);
//    direction << 0, 1;
//    for(int i = 0; i < rMembersSlaveDBC.rows(); i++)
//    {
//        myStructure->ConstraintLinearSetDisplacementNode(rMembersSlaveDBC(i), direction, disp);
//    }

    // ===> initial values
    disp = - 0.0005;
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
    std::vector<int> IPIdsMaster;
    std::vector<int> IPIdsSlave;

    switch(rIntegrationType)
    {
    case NuTo::eIntegrationType::IntegrationType1D2NGauss1Ip:
    {
        IPIdsMaster.clear();
        IPIdsMaster.push_back(0);

        IPIdsSlave.clear();
        IPIdsSlave.push_back(0);
        break;
    }
    case NuTo::eIntegrationType::IntegrationType1D2NGauss2Ip:
    {
        IPIdsMaster.clear();
        IPIdsMaster.push_back(1);
        IPIdsMaster.push_back(0);

        IPIdsSlave.clear();
        IPIdsSlave.push_back(0);
        IPIdsSlave.push_back(1);
        break;
    }
    case NuTo::eIntegrationType::IntegrationType1D2NGauss3Ip:
    {
        IPIdsMaster.clear();
        IPIdsMaster.push_back(2);
        IPIdsMaster.push_back(1);
        IPIdsMaster.push_back(0);

        IPIdsSlave.clear();
        IPIdsSlave.push_back(0);
        IPIdsSlave.push_back(1);
        IPIdsSlave.push_back(2);
        break;
    }
    case NuTo::eIntegrationType::IntegrationType1D2NGauss4Ip:
    {
        IPIdsMaster.clear();
        IPIdsMaster.push_back(3);
        IPIdsMaster.push_back(2);
        IPIdsMaster.push_back(1);
        IPIdsMaster.push_back(0);

        IPIdsSlave.clear();
        IPIdsSlave.push_back(0);
        IPIdsSlave.push_back(1);
        IPIdsSlave.push_back(2);
        IPIdsSlave.push_back(3);
        break;
    }
    case NuTo::eIntegrationType::IntegrationType1D2NGauss5Ip:
    {
        IPIdsMaster.clear();
        IPIdsMaster.push_back(4);
        IPIdsMaster.push_back(3);
        IPIdsMaster.push_back(2);
        IPIdsMaster.push_back(1);
        IPIdsMaster.push_back(0);

        IPIdsSlave.clear();
        IPIdsSlave.push_back(0);
        IPIdsSlave.push_back(1);
        IPIdsSlave.push_back(2);
        IPIdsSlave.push_back(3);
        IPIdsSlave.push_back(4);
        break;
    }
    case NuTo::eIntegrationType::IntegrationType1D2NGauss12Ip:
    {
        IPIdsMaster.clear();
        IPIdsMaster.push_back(11);
        IPIdsMaster.push_back(10);
        IPIdsMaster.push_back(9);
        IPIdsMaster.push_back(8);
        IPIdsMaster.push_back(7);
        IPIdsMaster.push_back(6);
        IPIdsMaster.push_back(5);
        IPIdsMaster.push_back(4);
        IPIdsMaster.push_back(3);
        IPIdsMaster.push_back(2);
        IPIdsMaster.push_back(1);
        IPIdsMaster.push_back(0);

        IPIdsSlave.clear();
        IPIdsSlave.push_back(0);
        IPIdsSlave.push_back(1);
        IPIdsSlave.push_back(2);
        IPIdsSlave.push_back(3);
        IPIdsSlave.push_back(4);
        IPIdsSlave.push_back(5);
        IPIdsSlave.push_back(6);
        IPIdsSlave.push_back(7);
        IPIdsSlave.push_back(8);
        IPIdsSlave.push_back(9);
        IPIdsSlave.push_back(10);
        IPIdsSlave.push_back(11);
        break;
    }
    default:
        throw NuTo::MechanicsException("[NuTo::Test::Contact] No Integration Type Defined.");
    }

    myIntegrationScheme.AddResultElementGroupIpData("ContactStressMaster1",  masterContactElementsGroupId, 1, IPIdsMaster, NuTo::IpData::eIpStaticDataType::ENGINEERING_STRESS);
    myIntegrationScheme.AddResultElementGroupIpData("ContactStressSlave1",  slaveContactElementsGroupId, 1, IPIdsSlave, NuTo::IpData::eIpStaticDataType::ENGINEERING_STRESS);

    myIntegrationScheme.SetMinTimeStepPlot(1.);
    myIntegrationScheme.SetLastTimePlot(0.);

    myIntegrationScheme.SetToleranceForce(1.e-10);
    myIntegrationScheme.SetMaxNumIterations(50);
    myIntegrationScheme.SetTimeStep(timeStep);
    myIntegrationScheme.SetPerformLineSearch(false);
    myIntegrationScheme.Solve(simulationTime);

#ifdef ENABLE_VISUALIZE
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
    myStructure->AddVisualizationComponent(groupElementsSlaveLower, NuTo::eVisualizeWhat::ENGINEERING_STRAIN);
    myStructure->AddVisualizationComponent(groupElementsSlaveLower, NuTo::eVisualizeWhat::ENGINEERING_STRESS);
    myStructure->ElementGroupExportVtkDataFile(groupElementsSlaveLower, resultDir+"/ElementsSlaveContact.vtu", true);

    myStructure->ExportVtkDataFileNodes(resultDir+"/Nodes.vtu", true);
#endif
}

void run(const std::string &resultDir, const std::string &path, const std::string &fileNameSlave, const std::string &fileNameMaster, double penalty, double gap, int contactAlgo, NuTo::eIntegrationType rIntegrationType,
         NuTo::Interpolation::eTypeOrder rElementTypeIdentDisps,
         NuTo::eIntegrationType rIntegrationTypeDisps)
{
    int groupElementsIGAlayer(0),
        groupElementsSlave(0),
        groupElementsMaster(0),
        groupElementsSlaveUpper(0),
        groupNodesSlaveUpper(0),
        groupElementsSlaveLower(0),
        groupElementsMasterUpper(0),
        slaveContactElementsGroupId(0),
        masterContactElementsGroupId(0);

    NuTo::Structure *myStructure = new NuTo::Structure(2);
    myStructure->SetNumTimeDerivatives(0);
    ContactHertzQuarterCircle(path,
                              fileNameSlave,
                              fileNameMaster,
                              contactAlgo,
                              myStructure,
                              penalty,
                              gap,
                              rIntegrationType,
                              rElementTypeIdentDisps,
                              rIntegrationTypeDisps,
                              groupElementsIGAlayer,
                              groupElementsSlave,
                              groupElementsMaster,
                              groupElementsSlaveUpper,
                              groupNodesSlaveUpper,
                              groupElementsSlaveLower,
                              groupElementsMasterUpper,
                              slaveContactElementsGroupId,
                              masterContactElementsGroupId);

    myStructure->SetNumTimeDerivatives(0);
    staticSolve(myStructure, resultDir, rIntegrationType,
                groupElementsIGAlayer,
                groupElementsSlave,
                groupElementsMaster,
                groupElementsSlaveUpper,
                groupNodesSlaveUpper,
                groupElementsSlaveLower,
                groupElementsMasterUpper,
                slaveContactElementsGroupId,
                masterContactElementsGroupId);
}


int main(int argc, char* argv[])
{
    std::string path = "./";

    std::string fileNameSlave  = "";
    std::string fileNameMaster = "";
    std::string resultDir = "";

    double penalty;
    double gap = -0.013;
    int    contactAlgo = 1; //mortar = 0, non-mortar = 1

    path = "/home/potto1/mechanics/otto/gmsh/";
    fileNameSlave  = "CircleUnstructuredQudrilat.msh";
    fileNameMaster = "masterHalfspaceCombined.msh";

    // ===> penalty parameter study <=== //
    resultDir = "./ResultsContactStatic_1_10Lobatto2";
    penalty = 1.e10;
    run(resultDir, path, fileNameSlave, fileNameMaster, penalty, gap, contactAlgo,  NuTo::eIntegrationType::IntegrationType1D2NGauss12Ip, NuTo::Interpolation::eTypeOrder::LOBATTO2, NuTo::eIntegrationType::IntegrationType2D4NLobatto9Ip);

    return 0;

    resultDir = "./ResultsContactStatic_1_10";
    penalty = 1.e10;
    run(resultDir, path, fileNameSlave, fileNameMaster, penalty, gap, contactAlgo,  NuTo::eIntegrationType::IntegrationType1D2NGauss12Ip, NuTo::Interpolation::eTypeOrder::EQUIDISTANT2, NuTo::eIntegrationType::IntegrationType2D4NGauss4Ip);

    contactAlgo = 0;

    resultDir = "./ResultsContactStatic_1_10Lobatto2";
    penalty = 1.e10;
    run(resultDir, path, fileNameSlave, fileNameMaster, penalty, gap, contactAlgo,  NuTo::eIntegrationType::IntegrationType1D2NGauss12Ip, NuTo::Interpolation::eTypeOrder::LOBATTO2, NuTo::eIntegrationType::IntegrationType2D4NLobatto9Ip);

    resultDir = "./ResultsContactStatic_1_10";
    penalty = 1.e10;
    run(resultDir, path, fileNameSlave, fileNameMaster, penalty, gap, contactAlgo,  NuTo::eIntegrationType::IntegrationType1D2NGauss12Ip, NuTo::Interpolation::eTypeOrder::EQUIDISTANT2, NuTo::eIntegrationType::IntegrationType2D4NGauss4Ip);

    return 0;

    resultDir = "./ResultsContactStatic_1_9";
    penalty = 1.e9;
    run(resultDir, path, fileNameSlave, fileNameMaster, penalty, gap, contactAlgo,  NuTo::eIntegrationType::IntegrationType1D2NGauss12Ip, NuTo::Interpolation::eTypeOrder::EQUIDISTANT2, NuTo::eIntegrationType::IntegrationType2D4NGauss4Ip);



    return 0;

    resultDir = "./ResultsContactStatic_1_10_Order3";
    penalty = 1.e10;
    run(resultDir, path, fileNameSlave, fileNameMaster, penalty, gap, contactAlgo,  NuTo::eIntegrationType::IntegrationType1D2NGauss12Ip, NuTo::Interpolation::eTypeOrder::LOBATTO3, NuTo::eIntegrationType::IntegrationType2D4NLobatto16Ip);

    resultDir = "./ResultsContactStatic_1_10_Order4";
    penalty = 1.e10;
    run(resultDir, path, fileNameSlave, fileNameMaster, penalty, gap, contactAlgo,  NuTo::eIntegrationType::IntegrationType1D2NGauss12Ip, NuTo::Interpolation::eTypeOrder::LOBATTO4, NuTo::eIntegrationType::IntegrationType2D4NLobatto25Ip);


    resultDir = "./ResultsContactStatic_5_10";
    penalty = 5.e10;
    run(resultDir, path, fileNameSlave, fileNameMaster, penalty, gap, contactAlgo,  NuTo::eIntegrationType::IntegrationType1D2NGauss12Ip, NuTo::Interpolation::eTypeOrder::LOBATTO2, NuTo::eIntegrationType::IntegrationType2D4NLobatto9Ip);

    resultDir = "./ResultsContactStatic_1_11";
    penalty = 1.e11;
    run(resultDir, path, fileNameSlave, fileNameMaster, penalty, gap, contactAlgo,  NuTo::eIntegrationType::IntegrationType1D2NGauss12Ip, NuTo::Interpolation::eTypeOrder::LOBATTO2, NuTo::eIntegrationType::IntegrationType2D4NLobatto9Ip);

    resultDir = "./ResultsContactStatic_5_11";
    penalty = 5.e11;
    run(resultDir, path, fileNameSlave, fileNameMaster, penalty, gap, contactAlgo,  NuTo::eIntegrationType::IntegrationType1D2NGauss12Ip, NuTo::Interpolation::eTypeOrder::LOBATTO2, NuTo::eIntegrationType::IntegrationType2D4NLobatto9Ip);

    // ===> integration order <=== //


    return 0;
}
