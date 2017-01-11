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

#include "nuto/mechanics/timeIntegration/NewmarkDirect.h"

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
    std::cout << "Node: " << node << ", Coordinates: " << nodeCoordinates.Trans() << std::endl;
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
                std::cout << "Node: " << node << ", Coordinates: " << nodeCoordinates.Trans() << std::endl;
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
                    std::cout << "Node: " << node << ", Coordinates: " << nodeCoordinates.Trans() << std::endl;
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

            if (rElementTypeIdent == NuTo::Interpolation::eTypeOrder::EQUIDISTANT1)
            {
                int temp = elementIncidence(2);
                elementIncidence(2) = elementIncidence(3);
                elementIncidence(3) = temp;

            }

            int myElement = myStructure->ElementCreate(interpolationType, elementIncidence);
            std::cout << "Element: " << myElement << ", Nodes: " << elementIncidence.Trans() << std::endl;
            myStructure->GroupAddElement(rElementGroup, myElement);
        }
    }
    return node;
}

void ContactTest(const std::string &resultDir,
                 NuTo::Interpolation::eTypeOrder rElementTypeIdent,
                 int rNumNodesPerElementInOneDir,
                 NuTo::FullVector<double, Eigen::Dynamic>& nodeCoordinatesFirstElement,
                 double rPenalty,
                 NuTo::eIntegrationType rIntegrationType,
                 int rContactAlgo,
                 int numElXSlave, int numElYSlave,
                 int numElXMaster, int numElYMaster)
{
    NuTo::Structure myStructure(2);
    myStructure.SetNumTimeDerivatives(0);

    #ifdef _OPENMP
    int numThreads = 4;
    myStructure.SetNumProcessors(numThreads);
    #endif
    /////////////////////////////////////////////////////////////////////
    // ====> create SLAVE mesh from gmsh (smth. impacting a rectangle) //
    /////////////////////////////////////////////////////////////////////

    double startySlave = 0.0;

    std::set<NuTo::Node::eDof> setOfDOFSSlave;
    setOfDOFSSlave.insert(NuTo::Node::eDof::COORDINATES);
    setOfDOFSSlave.insert(NuTo::Node::eDof::DISPLACEMENTS);

    int groupNodesSlave = myStructure.GroupCreate(NuTo::eGroupId::Nodes);
    int groupElementsSlave = myStructure.GroupCreate(NuTo::eGroupId::Elements);
    int node = buildStructure2D(rElementTypeIdent, rNumNodesPerElementInOneDir, nodeCoordinatesFirstElement, numElXSlave, numElYSlave, 1., 1., 0., 0., 0, &myStructure, groupNodesSlave, groupElementsSlave, setOfDOFSSlave);

    // ===> build contact elements slave elements
    auto LambdaGetSlaveNodesLower = [startySlave](NuTo::NodeBase* rNodePtr) -> bool
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
    };  // LambdaGetSlaveNodesLower

    int groupNodesSlaveLower = myStructure.GroupCreate(NuTo::eGroupId::Nodes);
    myStructure.GroupAddNodeFunction(groupNodesSlaveLower, groupNodesSlave, LambdaGetSlaveNodesLower);

    int groupElementsSlaveLower = myStructure.GroupCreate(NuTo::eGroupId::Elements);
    myStructure.GroupAddElementsFromNodes(groupElementsSlaveLower, groupNodesSlaveLower, false);

    //////////////////////////////
    // ====> create MASTER mesh //
    //////////////////////////////

    std::set<NuTo::Node::eDof> setOfDOFSMaster;
    setOfDOFSMaster.insert(NuTo::Node::eDof::COORDINATES);
    setOfDOFSMaster.insert(NuTo::Node::eDof::DISPLACEMENTS);

    int groupNodesMaster = myStructure.GroupCreate(NuTo::eGroupId::Nodes);
    int groupElementsMaster = myStructure.GroupCreate(NuTo::eGroupId::Elements);
    buildStructure2D(rElementTypeIdent, rNumNodesPerElementInOneDir, nodeCoordinatesFirstElement, numElXMaster, numElYMaster, 1., 1., 0., -1., node, &myStructure, groupNodesMaster, groupElementsMaster, setOfDOFSMaster);

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
    myStructure.GroupAddNodeFunction(groupNodesMasterUpper, groupNodesMaster, LambdaGetMasterNodesUpper);
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
    int degree = 1;

    Eigen::MatrixXd A;

    NuTo::NURBSCurve curve(degree, coordinatesAndIDs.block(0, 0, coordinates.rows(), coordinates.cols()), A);

//    for(int i = 0; i < 0; i++)
//    {
//        curve.DuplicateKnots();
//    }

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
    double rho = 0.;

    int constitutiveLaw = myStructure.ConstitutiveLawCreate("Linear_Elastic_Engineering_Stress");
    myStructure.ConstitutiveLawSetParameterDouble(constitutiveLaw, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, E);
    myStructure.ConstitutiveLawSetParameterDouble(constitutiveLaw, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, nue);
    myStructure.ConstitutiveLawSetParameterDouble(constitutiveLaw, NuTo::Constitutive::eConstitutiveParameter::DENSITY, rho);
    myStructure.ElementTotalSetSection(section);
    myStructure.ElementTotalSetConstitutiveLaw(constitutiveLaw);

    ////////////////////////////
    // ===> CONTACT ELEMENTS  //
    ////////////////////////////

    int constitutiveLawPC = myStructure.ConstitutiveLawCreate("Contact_Constitutive_Law");
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

    myStructure.ConstitutiveLawSetParameterFunction(constitutiveLawPC, NuTo::Constitutive::eConstitutiveParameter::CONSTITUTIVE_LAW_FUNCTION, constitutiveContactLaw);
    myStructure.ConstitutiveLawSetParameterFunction(constitutiveLawPC, NuTo::Constitutive::eConstitutiveParameter::CONSTITUTIVE_LAW_DERIVATIVE_FUNCTION, constitutiveContactLawDerivative);

    int slaveContactElementsGroupId = myStructure.NuTo::Structure::ContactElementsCreate<2,1>(groupElementsSlaveLower,
                                                                                              groupNodesSlaveLower,
                                                                                              elementsMaster,
                                                                                              rIntegrationType,
                                                                                              rContactAlgo,
                                                                                              constitutiveLawPC);

    ////////////////////////////
    // ===> IGA coupling      //
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

    int countDBC = 0;
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

    // ===> DBCs
    auto LambdaGetMasterNodesLower = [](NuTo::NodeBase* rNodePtr) -> bool
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
    };  // LambdaGetMasterNodesLower

    // DBC for the bottom of the master
    int groupNodesMasterLower = myStructure.GroupCreate(NuTo::eGroupId::Nodes);
    myStructure.GroupAddNodeFunction(groupNodesMasterLower, groupNodesMaster, LambdaGetMasterNodesLower);

    NuTo::FullVector<int,Eigen::Dynamic> rMembersMasterDBC;
    myStructure.NodeGroupGetMembers(groupNodesMasterLower, rMembersMasterDBC);

    NuTo::FullVector<double,Eigen::Dynamic> direction(2);
    direction << 0, 1;
    for(int i = 0; i < rMembersMasterDBC.rows(); i++)
    {
        myStructure.ConstraintLinearSetDisplacementNode(rMembersMasterDBC(i), direction, 0.0);
    }

    // slave upper
    auto LambdaGetSlaveNodesUpper = [](NuTo::NodeBase* rNodePtr) -> bool
    {
        double Tol = 1.e-6;
        if (rNodePtr->GetNum(NuTo::Node::eDof::COORDINATES)>0)
        {
            double y = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[1];
            if ((y >= 1. - Tol && y <= 1. + Tol))
            {
                return true;
            }
        }
        return false;
    };  // LambdaGetSlaveNodesUpper

    int groupNodesSlaveUpper = myStructure.GroupCreate(NuTo::eGroupId::Nodes);
    myStructure.GroupAddNodeFunction(groupNodesSlaveUpper, groupNodesSlave, LambdaGetSlaveNodesUpper);
    NuTo::FullVector<int, Eigen::Dynamic> membersDBCSlave = myStructure.GroupGetMemberIds(groupNodesSlaveUpper);

    double disp = - 0.00008;

    NuTo::FullVector<int,Eigen::Dynamic> rMembersSlaveDBC;
    myStructure.NodeGroupGetMembers(groupNodesSlaveUpper, rMembersSlaveDBC);

    direction << 0, 1;
    for(int i = 0; i < membersDBCSlave.rows(); i++)
    {
        myStructure.ConstraintLinearSetDisplacementNode(membersDBCSlave(i), direction, disp);
    }

    // slave master left
    auto LambdaGetSlaveNodesLeft = [](NuTo::NodeBase* rNodePtr) -> bool
    {
        double Tol = 1.e-6;
        if (rNodePtr->GetNum(NuTo::Node::eDof::COORDINATES)>0)
        {
            double x = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[0];
            if ((x >= 0. - Tol && x <= 0. + Tol))
            {
                return true;
            }
        }
        return false;
    };  // LambdaGetSlaveNodesLeft

    int groupNodesLeft = myStructure.GroupCreate(NuTo::eGroupId::Nodes);
    myStructure.GroupAddNodeFunction(groupNodesLeft, groupNodesSlave, LambdaGetSlaveNodesLeft);
    NuTo::FullVector<int, Eigen::Dynamic> rMembersLeftDBC = myStructure.GroupGetMemberIds(groupNodesLeft);

    direction << 1, 0;
    for(int i = 0; i < rMembersLeftDBC.rows(); i++)
    {
        myStructure.ConstraintLinearSetDisplacementNode(rMembersLeftDBC(i), direction, 0.);
    }

    // slave master left right
    auto LambdaGetSlaveNodesRight = [](NuTo::NodeBase* rNodePtr) -> bool
    {
        double Tol = 1.e-6;
        if (rNodePtr->GetNum(NuTo::Node::eDof::COORDINATES)>0)
        {
            double x = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[0];
            if ((x >= 1. - Tol && x <= 1. + Tol))
            {
                return true;
            }
        }
        return false;
    };  // LambdaGetSlaveNodesLeft

    int groupNodesRight = myStructure.GroupCreate(NuTo::eGroupId::Nodes);
    myStructure.GroupAddNodeFunction(groupNodesRight, groupNodesSlave, LambdaGetSlaveNodesRight);
    NuTo::FullVector<int, Eigen::Dynamic> rMembersRightDBC = myStructure.GroupGetMemberIds(groupNodesRight);


    direction << 1, 0;
    for(int i = 0; i < rMembersRightDBC.rows(); i++)
    {
        myStructure.ConstraintLinearSetDisplacementNode(rMembersRightDBC(i), direction, 0.);
    }


//    // ===> initial values
//    disp = - 0.0001;
//    int slaveNodesGroupId = myStructure->GroupCreate(NuTo::eGroupId::Nodes);
//    myStructure->GroupAddNodesFromElements(slaveNodesGroupId, groupElementsSlave, NuTo::Node::eDof::DISPLACEMENTS);
//    NuTo::FullVector<double,Eigen::Dynamic>  dispVec(2);
//    dispVec(0) =   0.;
//    dispVec(1) = disp;
//    myStructure->NodeGroupSetDisplacements(slaveNodesGroupId, 0, dispVec);

    myStructure.CalculateMaximumIndependentSets();
    myStructure.NodeBuildGlobalDofs();

    NuTo::NewmarkDirect myIntegrationScheme(&myStructure);
    double timeStep = 1.;
    double simulationTime = 1.;

    myIntegrationScheme.SetResultDirectory(resultDir, true);

    std::vector<int> IPIdsSlave;
    std::vector<int> IPIdsMaster;

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
    case NuTo::eIntegrationType::IntegrationType1D2NLobatto6Ip:
    {
        IPIdsMaster.clear();

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
        break;
    }
    case NuTo::eIntegrationType::IntegrationType1D2NLobatto4Ip:
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

    // create boundary elements for visualizing the results
//    int groupElementsMasterUpper = myStructure.GroupCreate(NuTo::eGroupId::Elements);
//    myStructure.GroupAddElementsFromNodes(groupElementsMasterUpper, groupNodesMasterUpper, false);

//    int masterContactElementsGroupId = myStructure.BoundaryElementsCreate(groupElementsMasterUpper, groupNodesMasterUpper, NULL);
//    NuTo::FullVector<int,Eigen::Dynamic> rMembersMasterContactBoundaryVisualize;
//    myStructure.ElementGroupGetMembers(masterContactElementsGroupId, rMembersMasterContactBoundaryVisualize);
//    for(int i = 0; i < rMembersMasterContactBoundaryVisualize.rows(); i++)
//    {
//        NuTo::ElementBase* elementPtr = myStructure.ElementGetElementPtr(rMembersMasterContactBoundaryVisualize(i));
//        NuTo::IpData::eIpDataType ipDataType = elementPtr->GetIpDataType(0);
//        elementPtr->SetIntegrationType(myStructure.GetPtrIntegrationType(rIntegrationType), ipDataType);
//    }


//    myIntegrationScheme.AddResultElementGroupIpData("ContactStressMaster1",  masterContactElementsGroupId, 1, IPIdsMaster, NuTo::IpData::eIpStaticDataType::ENGINEERING_STRESS);

    myIntegrationScheme.AddResultElementGroupIpData("ContactStressSlave1",  slaveContactElementsGroupId, 1, IPIdsSlave, NuTo::IpData::eIpStaticDataType::ENGINEERING_STRESS);

    myIntegrationScheme.SetMinTimeStepPlot(1.);
    myIntegrationScheme.SetLastTimePlot(0.);

    myIntegrationScheme.SetToleranceForce(1.e-10);
    myIntegrationScheme.SetMaxNumIterations(50);
    myIntegrationScheme.SetTimeStep(timeStep);
    myIntegrationScheme.SetPerformLineSearch(false);
    myIntegrationScheme.Solve(simulationTime);

#ifdef ENABLE_VISUALIZE
    int visualizationGroup = myStructure.GroupUnion(groupElementsSlave, groupElementsMaster);
    myStructure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::DISPLACEMENTS);
    myStructure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ENGINEERING_STRAIN);
    myStructure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ENGINEERING_STRESS);
    myStructure.ElementGroupExportVtkDataFile(visualizationGroup, resultDir+"/Elements.vtu", true);
#endif
}

int main()
{

    NuTo::FullVector<double, Eigen::Dynamic> nodeCoordinates(2);
    nodeCoordinates(0) = 0;
    nodeCoordinates(1) = 2;

    std::string resultDir = "./ResultsStaticElementElementGauss1";
    ContactTest(resultDir,
                NuTo::Interpolation::eTypeOrder::EQUIDISTANT1,
                2,
                nodeCoordinates,
                1.e9,
                NuTo::eIntegrationType::IntegrationType1D2NGauss12Ip, 1, 1, 1, 1, 1);

    return 0;

    nodeCoordinates.resize(3);
    nodeCoordinates(0) = 0;
    nodeCoordinates(1) = 1;
    nodeCoordinates(2) = 2;

    resultDir = "./ResultsStaticElementElementLobatto2";
    ContactTest(resultDir,
                NuTo::Interpolation::eTypeOrder::LOBATTO2,
                3,
                nodeCoordinates,
                1.e9,
                NuTo::eIntegrationType::IntegrationType1D2NGauss12Ip, 1, 5, 5, 1, 1);

    nodeCoordinates.resize(4);
    nodeCoordinates.fill(0.);
    NuTo::FullVector<double, Eigen::Dynamic> ones(4); ones.fill(1);
    NuTo::IntegrationType1D2NLobatto4Ip Lobatto1D2N4Ip;
    for (int i = 0; i < 4; i++) Lobatto1D2N4Ip.GetLocalIntegrationPointCoordinates1D(i, nodeCoordinates(i));
    nodeCoordinates+=ones;


    resultDir = "./ResultsStaticElementElementLobatto3";
    ContactTest(resultDir,
                NuTo::Interpolation::eTypeOrder::LOBATTO3,
                4,
                nodeCoordinates,
                1.e9,
                NuTo::eIntegrationType::IntegrationType1D2NGauss12Ip, 1, 5, 5, 1, 1);
    return 0;
}
