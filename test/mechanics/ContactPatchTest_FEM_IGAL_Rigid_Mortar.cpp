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
#include "nuto/mechanics/integrationtypes/IntegrationType1D2NGaussNIp.h"

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

#include "nuto/mechanics/elements/ContinuumContactElement.h"

#include "nuto/mechanics/interpolationtypes/InterpolationTypeEnum.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeBase.h"

#include "nuto/mechanics/interpolationtypes/InterpolationType.h"
#include "nuto/mechanics/interpolationtypes/InterpolationBase.h"


#include <math.h>

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
//    std::cout << "Node: " << node << ", Coordinates: " << nodeCoordinates.Trans() << std::endl;
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
//                std::cout << "Node: " << node << ", Coordinates: " << nodeCoordinates.Trans() << std::endl;
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
//                    std::cout << "Node: " << node << ", Coordinates: " << nodeCoordinates.Trans() << std::endl;
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

//            NuTo::ElementBase* el = myStructure->ElementGetElementPtr(myElement);
//            Eigen::VectorXd nodalCurrent;
//            nodalCurrent.setZero(0);
//            nodalCurrent = el->ExtractNodeValues(0, NuTo::Node::eDof::DISPLACEMENTS);
//            nodalCurrent = el->ExtractNodeValues(0, NuTo::Node::eDof::COORDINATES);

//            std::cout << "Element: " << myElement << ", Nodes: " << elementIncidence.Trans() << std::endl;
            myStructure->GroupAddElement(rElementGroup, myElement);
        }
    }
    myStructure->ElementConvertToInterpolationType(rElementGroup, 1.e-6, 10);
    return node;
}

void SetDBCPatchTestRigidIGA(NuTo::Structure &myStructure, int groupNodesSlave, int groupNodesIGAlayer, int groupNodesMaster, int &countDBC)
{
    Eigen::Vector2d direction(0,0);
    double xMin(0.), xMax(0.);
    double yMin(0.), yMax(0.);
    auto LambdaNodes = [&](NuTo::NodeBase* rNodePtr) -> bool
    {
        double Tol = 1.e-6;
        if (rNodePtr->GetNum(NuTo::Node::eDof::COORDINATES)>0)
        {
            double x = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[0];
            double y = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[1];
            bool xR = false;
            bool yR = false;

            if (x >= xMin - Tol && x <= xMax + Tol) xR = true;
            if (y >= yMin - Tol && y <= yMax + Tol) yR = true;

            if (xR == true && yR == true) return true;
        }
        return false;
    };

    NuTo::FullMatrix<int, Eigen::Dynamic> members;
    // ===> master bottom <=== //
    xMin = 0.; xMax = std::numeric_limits<double>::infinity();
    yMin = yMax = 0.;
    int groupNodesMasterLower = myStructure.GroupCreate(NuTo::eGroupId::Nodes);
    myStructure.GroupAddNodeFunction(groupNodesMasterLower, groupNodesMaster, LambdaNodes);

    members = myStructure.GroupGetMemberIds(groupNodesMasterLower);

    direction << 0.,1.;
    countDBC = myStructure.ConstraintLinearSetDisplacementNodeGroup(groupNodesMasterLower, direction, 0.0);
    direction << 1.,0.;
    countDBC = myStructure.ConstraintLinearSetDisplacementNodeGroup(groupNodesMasterLower, direction, 0.0);

    // ===> initial values <=== //
    NuTo::FullVector<double,Eigen::Dynamic>  dispVec(2);
    dispVec(0) =   0.;
    dispVec(1) = -0.00001;
    myStructure.NodeGroupSetDisplacements(groupNodesSlave, 0, dispVec);
    if(groupNodesIGAlayer >= 0)
        myStructure.NodeGroupSetDisplacements(groupNodesIGAlayer, 0, dispVec);

    // ===> slave left <=== //
    xMin = 1.; xMax = 1.;
    yMin = 1.; yMax = 1.;
    int groupNodesLeft = myStructure.GroupCreate(NuTo::eGroupId::Nodes);
    myStructure.GroupAddNodeFunction(groupNodesLeft, groupNodesSlave, LambdaNodes);

    members = myStructure.GroupGetMemberIds(groupNodesLeft);

    direction << 1, 0;
    countDBC = myStructure.ConstraintLinearSetDisplacementNodeGroup(groupNodesLeft, direction, 0.0);

    // ===> PATCH TEST BOUNDARY <=== //
    double Stress = 10.;
    myStructure.SetNumLoadCases(1);
    // ===> slave top <=== //
    xMin = 1.; xMax = 2.;
    yMin = 1.; yMax = 1.;
    int groupNodesSlaveUpper = myStructure.GroupCreate(NuTo::eGroupId::Nodes);
    myStructure.GroupAddNodeFunction(groupNodesSlaveUpper, LambdaNodes);

    members = myStructure.GroupGetMemberIds(groupNodesSlaveUpper);

    int groupElementsSlaveUpper = myStructure.GroupCreate(NuTo::eGroupId::Elements);
    myStructure.GroupAddElementsFromNodes(groupElementsSlaveUpper, groupNodesSlaveUpper, false);

    members = myStructure.GroupGetMemberIds(groupElementsSlaveUpper);

    myStructure.LoadSurfacePressureCreate2D(0, groupElementsSlaveUpper, groupNodesSlaveUpper, Stress);

    // ===> slave bottom <=== //
//    xMin = 1.; xMax = 2.;
//    yMin = 0.; yMax = 0.;
//    int groupNodesSlaveBottom = myStructure.GroupCreate(NuTo::eGroupId::Nodes);
//    myStructure.GroupAddNodeFunction(groupNodesSlaveBottom, groupNodesSlave, LambdaNodes);
//    direction << 0, 1;
//    countDBC = myStructure.ConstraintLinearSetDisplacementNodeGroup(groupNodesSlaveBottom, direction, 0.0);


//    NuTo::FullVector<int,Eigen::Dynamic> rMembersSlaveBottom;
//    myStructure.NodeGroupGetMembers(groupNodesSlaveBottom, rMembersSlaveBottom);
//    direction << 1.,0.;
//    countDBC = myStructure.ConstraintLinearSetDisplacementNode(rMembersSlaveBottom(0), direction, 0.0);

    countDBC++;
}

void AddIGALayerMortar(NuTo::Structure *myStructure,
                       const std::function<bool(NuTo::NodeBase *)> &rFunction,
                       int rNodesGroupId,
                       int rDegree,
                       int &groupElementsIGAlayer,
                       int &groupNodesIGAlayer,
                       Eigen::Matrix<std::pair<int, int>, Eigen::Dynamic, Eigen::Dynamic> &rElements)
{
    // Nodes on the part to interpolate
    int groupFENodes = myStructure->GroupCreate(NuTo::eGroupId::Nodes);
    myStructure->GroupAddNodeFunction(groupFENodes, rNodesGroupId, rFunction);
    NuTo::FullVector<int,Eigen::Dynamic> idsFENodes;
    myStructure->NodeGroupGetMembers(groupFENodes, idsFENodes);

    // Matrix containing the ids and coordinates of the FE nodes => 'coordinatesAndIDs'
    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> coordinates;
    myStructure->NodeGroupGetCoordinates(groupFENodes, coordinates);

    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> coordinatesAndIDs;
    coordinatesAndIDs.resize(coordinates.rows(), coordinates.cols() + 1);
    coordinatesAndIDs.block(0, 0, coordinates.rows(), coordinates.cols()) = coordinates;

    for(int i = 0; i < coordinates.rows(); i++)
        coordinatesAndIDs(i, coordinates.cols()) = idsFENodes(i);

    coordinatesAndIDs.SortRow(0);

    Eigen::MatrixXd controlPointsCoordinates(coordinatesAndIDs.rows(), 2);
    for(int i = 0; i < coordinatesAndIDs.rows(); i++)
    {
        controlPointsCoordinates(i,0) = coordinatesAndIDs(i,0);
        controlPointsCoordinates(i,1) = coordinatesAndIDs(i,1);
    }

    int numKnots = controlPointsCoordinates.rows() + rDegree + 1;

    Eigen::VectorXd weights(controlPointsCoordinates.rows());
    weights.fill(1.);

    Eigen::VectorXd knots(numKnots);
    int numElements  = controlPointsCoordinates.rows() - rDegree;

    for (int i = 0; i <= rDegree; i++) knots(i) = 0.;
    for (int i = rDegree + 1; i <= rDegree + numElements - 1 ; i++) knots(i) = knots(i-1) + 1./numElements;
    for (int i = rDegree + numElements ; i < numKnots; i++) knots(i) = 1;

    NuTo::NURBSCurve curve(knots, controlPointsCoordinates, weights, rDegree);

    for(int i = 0; i < 0; i++) curve.DuplicateKnots();

    // Build the IGA layer
    std::set<NuTo::Node::eDof> setOfDOFS;
    setOfDOFS.insert(NuTo::Node::eDof::COORDINATES);
    setOfDOFS.insert(NuTo::Node::eDof::DISPLACEMENTS);

    groupNodesIGAlayer = myStructure->GroupCreate("Nodes");
    groupElementsIGAlayer = myStructure->GroupCreate("Elements");

    Eigen::VectorXi nodeIDs;
    rElements = curve.buildIGAStructure(*myStructure, setOfDOFS, groupElementsIGAlayer, groupNodesIGAlayer, "IGA1DLAYER", nodeIDs);
}

void AddIGALayerMortarOther(NuTo::Structure *myStructure,
                            double startCoordinate,
                            double length,
                            int    numRefinements,
                            int    rDegree,
                            int    &groupElementsIGAlayer,
                            int    &groupNodesIGAlayer,
                            Eigen::Matrix<std::pair<int, int>, Eigen::Dynamic, Eigen::Dynamic> &rElements)
{
    int numControlPoints = rDegree+1;
    Eigen::MatrixXd controlPointsCoordinates(numControlPoints, 2);
    for(int i = 0; i < numControlPoints; i++)
    {
        controlPointsCoordinates(i,0) = startCoordinate + i*length/(numControlPoints - 1);
        controlPointsCoordinates(i,1) = 0.;
    }

    int numKnots = controlPointsCoordinates.rows() + rDegree + 1;

    Eigen::VectorXd weights(controlPointsCoordinates.rows());
    weights.fill(1.);

    Eigen::VectorXd knots(numKnots);
    int numElements  = controlPointsCoordinates.rows() - rDegree;

    for (int i = 0; i <= rDegree; i++) knots(i) = 0.;
    for (int i = rDegree + 1; i <= rDegree + numElements - 1 ; i++) knots(i) = knots(i-1) + 1./numElements;
    for (int i = rDegree + numElements ; i < numKnots; i++) knots(i) = 1;

    NuTo::NURBSCurve curve(knots, controlPointsCoordinates, weights, rDegree);

    for(int i = 0; i < numRefinements; i++) curve.DuplicateKnots();

    // Build the IGA layer
    std::set<NuTo::Node::eDof> setOfDOFS;
    setOfDOFS.insert(NuTo::Node::eDof::COORDINATES);
    setOfDOFS.insert(NuTo::Node::eDof::DISPLACEMENTS);

    groupNodesIGAlayer = myStructure->GroupCreate("Nodes");
    groupElementsIGAlayer = myStructure->GroupCreate("Elements");

    Eigen::VectorXi nodeIDs;
    rElements = curve.buildIGAStructure(*myStructure, setOfDOFS, groupElementsIGAlayer, groupNodesIGAlayer, "IGA1DLAYER", nodeIDs);
}

void AddIGALayerMortarOther2(NuTo::Structure *myStructure,
                             double startCoordinate,
                             double length,
                             int    numRefinements,
                             int    rDegree,
                             int    &groupElementsIGAlayer,
                             int    &groupNodesIGAlayer,
                             Eigen::Matrix<std::pair<int, int>, Eigen::Dynamic, Eigen::Dynamic> &rElements)
{
    int numControlPoints = rDegree+2;
    Eigen::MatrixXd controlPointsCoordinates;
    controlPointsCoordinates.setZero(numControlPoints, 2);

    controlPointsCoordinates(0,0)  = 1.;
    controlPointsCoordinates(1,0)  = 1.25;
    controlPointsCoordinates(2,0)  = 1.75;
    controlPointsCoordinates(3,0)  = 2.0;

    int numKnots = controlPointsCoordinates.rows() + rDegree + 1;

    Eigen::VectorXd weights(controlPointsCoordinates.rows());
    weights.fill(1.);

    Eigen::VectorXd knots(numKnots);
    int numElements  = controlPointsCoordinates.rows() - rDegree;

    for (int i = 0; i <= rDegree; i++) knots(i) = 0.;
    for (int i = rDegree + 1; i <= rDegree + numElements - 1 ; i++) knots(i) = knots(i-1) + 1./numElements;
    for (int i = rDegree + numElements ; i < numKnots; i++) knots(i) = 1;

    knots(3) = 0.4;

    NuTo::NURBSCurve curve(knots, controlPointsCoordinates, weights, rDegree);

    //for(int i = 0; i < numRefinements; i++) curve.DuplicateKnots();

    // Build the IGA layer
    std::set<NuTo::Node::eDof> setOfDOFS;
    setOfDOFS.insert(NuTo::Node::eDof::COORDINATES);
    setOfDOFS.insert(NuTo::Node::eDof::DISPLACEMENTS);

    groupNodesIGAlayer = myStructure->GroupCreate("Nodes");
    groupElementsIGAlayer = myStructure->GroupCreate("Elements");

    Eigen::VectorXi nodeIDs;
    rElements = curve.buildIGAStructure(*myStructure, setOfDOFS, groupElementsIGAlayer, groupNodesIGAlayer, "IGA1DLAYER", nodeIDs);
}

void AddIGALayerMortarInterp(NuTo::Structure *myStructure,
                             const std::function<bool(NuTo::NodeBase *)> &rFunction,
                             int rNodesGroupId,
                             int rDegree,
                             int &groupElementsIGAlayer,
                             int &groupNodesIGAlayer,
                             Eigen::Matrix<std::pair<int, int>, Eigen::Dynamic, Eigen::Dynamic> &rElements)
{
    // Nodes on the part to interpolate
    int groupFENodes = myStructure->GroupCreate(NuTo::eGroupId::Nodes);
    myStructure->GroupAddNodeFunction(groupFENodes, rNodesGroupId, rFunction);
    NuTo::FullVector<int,Eigen::Dynamic> idsFENodes;
    myStructure->NodeGroupGetMembers(groupFENodes, idsFENodes);

    // Matrix containing the ids and coordinates of the FE nodes => 'coordinatesAndIDs'
    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> coordinates;
    myStructure->NodeGroupGetCoordinates(groupFENodes, coordinates);

    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> coordinatesAndIDs;
    coordinatesAndIDs.resize(coordinates.rows(), coordinates.cols() + 1);
    coordinatesAndIDs.block(0, 0, coordinates.rows(), coordinates.cols()) = coordinates;

    for(int i = 0; i < coordinates.rows(); i++)
        coordinatesAndIDs(i, coordinates.cols()) = idsFENodes(i);

    coordinatesAndIDs.SortRow(0);

    Eigen::MatrixXd A;
    NuTo::NURBSCurve curve(rDegree, coordinatesAndIDs.block(0, 0, coordinates.rows(), coordinates.cols()), A);

    for(int i = 0; i < 0; i++)
    {
        curve.DuplicateKnots();
    }

    std::set<NuTo::Node::eDof> setOfDOFS;
    setOfDOFS.insert(NuTo::Node::eDof::COORDINATES);
    setOfDOFS.insert(NuTo::Node::eDof::DISPLACEMENTS);

    groupNodesIGAlayer = myStructure->GroupCreate("Nodes");
    groupElementsIGAlayer = myStructure->GroupCreate("Elements");

    Eigen::VectorXi nodeIDs;
    rElements = curve.buildIGAStructure(*myStructure,
                                        setOfDOFS,
                                        groupElementsIGAlayer,
                                        groupNodesIGAlayer,
                                        "IGA1DLAYER",
                                        nodeIDs);
}

void ContactTestOneElementLayerSlave(const std::string &resultDir,
                                     NuTo::Interpolation::eTypeOrder rElementTypeIdent,
                                     int rNumNodesPerElementInOneDir,
                                     NuTo::FullVector<double, Eigen::Dynamic>& nodeCoordinatesFirstElement,
                                     int rDegree,
                                     double rPenalty,
                                     NuTo::eIntegrationType rIntegrationType,
                                     int rContactAlgo,
                                     int numElXSlave,
                                     int numElYSlave,
                                     int numRefinementsLayer)
{
    NuTo::Structure myStructure(2);
    myStructure.SetNumTimeDerivatives(0);

#ifdef _OPENMP
    int numThreads = 4;
    myStructure.SetNumProcessors(numThreads);
#endif

    // *** slave *** //
    std::set<NuTo::Node::eDof> setOfDOFSSlave;
    setOfDOFSSlave.insert(NuTo::Node::eDof::COORDINATES);
    setOfDOFSSlave.insert(NuTo::Node::eDof::DISPLACEMENTS);

    int groupNodesSlave = myStructure.GroupCreate(NuTo::eGroupId::Nodes);
    int groupElementsSlave = myStructure.GroupCreate(NuTo::eGroupId::Elements);
    double startX = 1.;
    double length = 1.;
    buildStructure2D(rElementTypeIdent, rNumNodesPerElementInOneDir, nodeCoordinatesFirstElement,
                     numElXSlave, numElYSlave, 1., length, startX, 0., 0, &myStructure, groupNodesSlave, groupElementsSlave, setOfDOFSSlave);

    // ===> build contact elements slave elements
    auto LambdaGetSlaveNodesLower = [](NuTo::NodeBase* rNodePtr) -> bool
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

    // *** master - rigid *** //
    int numPoints = 10;
    double LengthMaster = 3.;
    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> points(numPoints,2);
    for(int i = 0; i < numPoints; i++)
    {
        points(i,0) = i*LengthMaster/(numPoints-1);
        points(i,1) = 0.;
    }

    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> AInv;
    NuTo::NURBSCurve curve(rDegree, points, AInv);

    for(int i = 0; i < 1; i++) curve.DuplicateKnots();

    std::set<NuTo::Node::eDof> setOfDOFS;
    setOfDOFS.insert(NuTo::Node::eDof::COORDINATES);
    setOfDOFS.insert(NuTo::Node::eDof::DISPLACEMENTS);

    int groupNodesIGA    = myStructure.GroupCreate("Nodes");
    int groupElementsIGA = myStructure.GroupCreate("Elements");

    Eigen::VectorXi nodeIDs;
    // master side fixed
    Eigen::Matrix<std::pair<int, int>, Eigen::Dynamic, Eigen::Dynamic> elementsMaster =
    curve.buildIGAStructure(myStructure, setOfDOFS, groupElementsIGA, groupNodesIGA, "IGA1DLAYER", nodeIDs);

    int countDBC = 0;

    // contact
    int groupElementsIGAlayer, groupNodesIGAlayer;
    Eigen::Matrix<std::pair<int, int>, Eigen::Dynamic, Eigen::Dynamic> rElements;
    AddIGALayerMortarInterp(&myStructure,
                            LambdaGetSlaveNodesLower,
                            groupNodesSlave,
                            rDegree,
                            groupElementsIGAlayer,
                            groupNodesIGAlayer,
                            rElements);

//    AddIGALayerMortarOther(&myStructure,
//                           startX,
//                           length,
//                           numRefinementsLayer,
//                           2, // degree
//                           groupElementsIGAlayer,
//                           groupNodesIGAlayer,
//                           rElements);

//    AddIGALayerMortarOther2(&myStructure,
//                            startX,
//                            length,
//                            numRefinementsLayer,
//                            2, // degree
//                            groupElementsIGAlayer,
//                            groupNodesIGAlayer,
//                            rElements);

    SetDBCPatchTestRigidIGA(myStructure, groupNodesSlave, groupNodesIGAlayer, groupNodesIGA,  countDBC);

    ///////////////////
    // ===> material //
    ///////////////////

    double Thickness = 1.;
    int section = myStructure.SectionCreate("PLANE_STRESS");
    myStructure.SectionSetThickness(section, Thickness);

    double E = rPenalty;
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
        if(rGap < 0.)
            return rPenalty*rGap;
        else
            return 0.;
    };

    std::function<double(double)> constitutiveContactLawDerivative =
            [rPenalty](double rGap) -> double
    {
        if(rGap < 0.)
            return rPenalty;
        else
            return 0.;
    };

    myStructure.ConstitutiveLawSetParameterFunction(constitutiveLawPC, NuTo::Constitutive::eConstitutiveParameter::CONSTITUTIVE_LAW_FUNCTION, constitutiveContactLaw);
    myStructure.ConstitutiveLawSetParameterFunction(constitutiveLawPC, NuTo::Constitutive::eConstitutiveParameter::CONSTITUTIVE_LAW_DERIVATIVE_FUNCTION, constitutiveContactLawDerivative);

    // contact
    myStructure.NuTo::Structure::ContactElementsCreate<1,1>(groupElementsIGAlayer, groupNodesIGAlayer, elementsMaster, rIntegrationType, rContactAlgo, constitutiveLawPC);

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

    myStructure.CalculateMaximumIndependentSets();
    myStructure.NodeBuildGlobalDofs();

    // layer tying
    int cceID = myStructure.NuTo::Structure::ContactElementsCreate<2,1>(groupElementsSlaveLower, groupNodesSlaveLower,
                                                                        rElements,
                                                                        rIntegrationType, rContactAlgo, constitutiveLawPC);

    NuTo::ElementBase* elementBase = myStructure.ElementGetElementPtr(cceID);
    NuTo::ContinuumContactElement<2,1> &contactElement = elementBase->AsContinuumContactElement21();

    Eigen::MatrixXd D, M;
    std::unordered_map<int, int> mappingGlobal2LocalSlaveNode, mappingGlobal2LocalMasterNode;
    contactElement.ComputeMeshTyingMatrix(D, M, mappingGlobal2LocalSlaveNode, mappingGlobal2LocalMasterNode);

//    std::cout << D << std::endl << std::endl;
//    std::cout << M << std::endl;

    if(1)
    {
        for(auto &itSlaveRow : mappingGlobal2LocalSlaveNode) // global
        {
            int localSlaveNodeIdRow  = itSlaveRow.second;

            for(int dof = 0; dof < myStructure.GetDimension(); dof++)
            {
                std::unordered_map<int, int>::iterator itSlaveCol = mappingGlobal2LocalSlaveNode.begin();
                int localSlaveNodeIdCol = itSlaveCol->second;
                int globalSlaveNodeId = itSlaveCol->first;
                // define equation
                myStructure.ConstraintLinearEquationCreate(countDBC, globalSlaveNodeId, NuTo::Node::eDof::DISPLACEMENTS, dof, D(localSlaveNodeIdRow, localSlaveNodeIdCol), 0.);
//                std::cout << D(localSlaveNodeIdRow, localSlaveNodeIdCol) << " ";
                itSlaveCol++;

                for(;itSlaveCol != mappingGlobal2LocalSlaveNode.end() ;itSlaveCol++)
                {
                    localSlaveNodeIdCol  = itSlaveCol->second;
                    globalSlaveNodeId = itSlaveCol->first;
                    // add slave ingredients
                    myStructure.ConstraintLinearEquationAddTerm(countDBC, globalSlaveNodeId, NuTo::Node::eDof::DISPLACEMENTS, dof, D(localSlaveNodeIdRow, localSlaveNodeIdCol));
//                    std::cout << D(localSlaveNodeIdRow, localSlaveNodeIdCol) << " ";
                    // master side
                }

//                std::cout << " || ";
                for(auto &itMaster : mappingGlobal2LocalMasterNode)
                {
                    int globalMasterNodeId = itMaster.first;
                    int localMasterNodeIdRow  = itMaster.second;
                    // add master ingredients
                    myStructure.ConstraintLinearEquationAddTerm(countDBC, globalMasterNodeId, NuTo::Node::eDof::DISPLACEMENTS, dof, M(localSlaveNodeIdRow, localMasterNodeIdRow));

//                    std::cout << M(localSlaveNodeIdRow, localMasterNodeIdRow) << " ";
                }
//                std::cout << std::endl;
                countDBC++;
            }
        }
    }
    else
    {

        Eigen::MatrixXd D_small, M_small;
        D_small.resize(mappingGlobal2LocalSlaveNode.size(), mappingGlobal2LocalSlaveNode.size());
        M_small.resize(mappingGlobal2LocalSlaveNode.size(), mappingGlobal2LocalMasterNode.size());
        int i = 0, jSlave = 0, jMaster  = 0;
        for(auto &itSlaveRow : mappingGlobal2LocalSlaveNode)
        {
            int row = itSlaveRow.second;
            jSlave = 0;
            for(auto &itSlaveCol : mappingGlobal2LocalSlaveNode)
            {
                int col = itSlaveCol.second;
                D_small(i,jSlave) = D(row,col);
                jSlave++;
            }
            jMaster = 0;
            for(auto &itMasterCol : mappingGlobal2LocalMasterNode)
            {
                int colMaster = itMasterCol.second;
                M_small(i,jMaster) = M(row,colMaster);
                jMaster++;
            }
            i++;
        }

    //    std::cout << std::endl << D_small << std::endl << std::endl;
    //    std::cout << M_small << std::endl << std::endl;
        Eigen::MatrixXd P = (D_small.inverse())*M_small;
    //    std::cout << D_small.inverse()*(D_small) << std::endl << std::endl;
    //    std::cout << P << std::endl << std::endl;

        i = 0;
        for(auto &itSlaveRow : mappingGlobal2LocalSlaveNode)
        {
            int globalSlaveNodeId = itSlaveRow.first;
            for(int dof = 0; dof < myStructure.GetDimension(); dof++)
            {
                myStructure.ConstraintLinearEquationCreate(countDBC, globalSlaveNodeId, NuTo::Node::eDof::DISPLACEMENTS, dof, 1., 0.);
                int j = 0;
                for(auto &itMasterCol : mappingGlobal2LocalMasterNode)
                {
                    int globalMasterNodeId = itMasterCol.first;
                    myStructure.ConstraintLinearEquationAddTerm(countDBC, globalMasterNodeId, NuTo::Node::eDof::DISPLACEMENTS, dof, P(i, j));
                    j++;
                }
                countDBC++;
            }
            i++;
        }
    }

    myStructure.ElementDelete(cceID);

    // end- layer tying

    NuTo::NewmarkDirect myIntegrationScheme(&myStructure);
    double timeStep = 1.;
    double simulationTime = 1.;

    myIntegrationScheme.SetResultDirectory(resultDir, true);

    std::vector<int> IPIdsSlave;
    std::vector<int> IPIdsMaster;

    NuTo::IntegrationTypeBase* integrPtr = myStructure.GetPtrIntegrationType(rIntegrationType);

    IPIdsSlave.clear();
    IPIdsMaster.clear();

    int  numIPs = integrPtr->GetNumIntegrationPoints();
    for(int i = 0; i < numIPs; i++)
    {
        IPIdsSlave.push_back(i);
        IPIdsMaster.push_back(numIPs - i - 1);
    }

    // create boundary elements for visualizing the results slave
    int slaveContactElementsGroupId = myStructure.BoundaryElementsCreate(groupElementsSlaveLower, groupNodesSlaveLower, NULL);
    NuTo::FullVector<int,Eigen::Dynamic> rMembersSlaveContactBoundaryVisualize;
    myStructure.ElementGroupGetMembers(slaveContactElementsGroupId, rMembersSlaveContactBoundaryVisualize);
    for(int i = 0; i < rMembersSlaveContactBoundaryVisualize.rows(); i++)
    {
        NuTo::ElementBase* elementPtr = myStructure.ElementGetElementPtr(rMembersSlaveContactBoundaryVisualize(i));
        NuTo::IpData::eIpDataType ipDataType = elementPtr->GetIpDataType(0);
        elementPtr->SetIntegrationType(myStructure.GetPtrIntegrationType(rIntegrationType), ipDataType);
    }

    // plot the stresses on boundary elements
    myIntegrationScheme.AddResultElementGroupIpData("ContactStressSlave1",  slaveContactElementsGroupId, 1, IPIdsSlave, NuTo::IpData::eIpStaticDataType::ENGINEERING_STRESS);
    // plot the stresses overall
    NuTo::FullVector<int,Eigen::Dynamic> rMembersSlaveVisualize;
    myStructure.ElementGroupGetMembers(groupElementsSlave, rMembersSlaveVisualize);
    NuTo::ElementBase* elementPtrSlave = myStructure.ElementGetElementPtr(rMembersSlaveVisualize(0));
    const NuTo::InterpolationBase& interpolationTypeDispSlave = elementPtrSlave->GetInterpolationType()->Get(NuTo::Node::eDof::DISPLACEMENTS);
    NuTo::eIntegrationType ittype = interpolationTypeDispSlave.GetStandardIntegrationType();
    integrPtr = myStructure.GetPtrIntegrationType(ittype);
    numIPs = integrPtr->GetNumIntegrationPoints();
    std::vector<int> IPIds;
    for(int i = 0; i < numIPs; i++) IPIds.push_back(i);
    myIntegrationScheme.AddResultElementGroupIpData("StressSlave", groupElementsSlave, 1, IPIds, NuTo::IpData::eIpStaticDataType::ENGINEERING_STRESS);

    myIntegrationScheme.SetMinTimeStepPlot(1.);
    myIntegrationScheme.SetLastTimePlot(0.);

    myIntegrationScheme.SetToleranceForce(1.e-10);
    myIntegrationScheme.SetMaxNumIterations(50);
    myIntegrationScheme.SetTimeStep(timeStep);
    myIntegrationScheme.SetPerformLineSearch(false);
    myIntegrationScheme.Solve(simulationTime);

#ifdef ENABLE_VISUALIZE
    int visualizationGroup = groupElementsSlave;
    myStructure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::DISPLACEMENTS);
    myStructure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ENGINEERING_STRAIN);
    myStructure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ENGINEERING_STRESS);
    myStructure.ElementGroupExportVtkDataFile(visualizationGroup, resultDir+"/Elements.vtu", true);

    myStructure.AddVisualizationComponent(groupElementsIGAlayer, NuTo::eVisualizeWhat::DISPLACEMENTS);
    myStructure.ElementGroupExportVtkDataFile(groupElementsIGAlayer, resultDir+"/ElementsLayerSlave.vtu", true);

    myStructure.AddVisualizationComponent(groupElementsIGA, NuTo::eVisualizeWhat::DISPLACEMENTS);
    myStructure.ElementGroupExportVtkDataFile(groupElementsIGA, resultDir+"/ElementsLayer.vtu", true);
#endif
}


#define EPS 2e-16

void gauleg(double x1, double x2, double x[], double w[], int n)
{
    int m,j,i;
    double z1,z,xm,xl,pp,p3,p2,p1;

    m=(n+1)/2;
    xm=0.5*(x2+x1);
    xl=0.5*(x2-x1);
    for (i=1;i<=m;i++) {
        z=cos(M_PI*(i-0.25)/(n+0.5));
        do {
            p1=1.0;
            p2=0.0;
            for (j=1;j<=n;j++) {
                p3=p2;
                p2=p1;
                p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
            }
            pp=n*(z*p1-p2)/(z*z-1.0);
            z1=z;
            z=z1-p1/pp;
        } while (fabs(z-z1) > EPS);
        x[i-1]=xm-xl*z;
        x[n+1-i-1]=xm+xl*z;
        w[i-1]=2.0*xl/((1.0-z*z)*pp*pp);
        w[n+1-i-1]=w[i-1];
    }
}
#undef EPS

int main()
{
    int n = 40;
    double *x, *w;
    x = new double[n];
    w = new double[n];
    gauleg(-1, 1, x, w, n);
    std::cout << "Absc: \n";
    for(int  i = 0; i < n; i++)
    {
        std::cout << std::setprecision(16) << x[i] << " ";
    }
    std::cout << "\nWeights: \n";
    for(int  i = 0; i < n; i++)
    {
        std::cout << std::setprecision(16) << w[i] << " ";
    }
    std::cout << std::flush;
    delete[] x;
    delete[] w;


    std::string resultDir = "";
    int degree = 0;

    // LOBATTO 2
    NuTo::FullVector<double, Eigen::Dynamic> nodeCoordinates(3);

    int contactAlgorithm = 0;
    degree = 2;

    nodeCoordinates.resize(3);
    nodeCoordinates(0) = 0;
    nodeCoordinates(1) = 1;
    nodeCoordinates(2) = 2;
    contactAlgorithm = 0;

    resultDir = "./ResultsStaticFEM_IGAL_Rigid_Mortar_2E_5IP";

    ContactTestOneElementLayerSlave(resultDir,
                                    NuTo::Interpolation::eTypeOrder::LOBATTO2,
                                    3,
                                    nodeCoordinates,
                                    degree,
                                    1.e6,
                                    NuTo::eIntegrationType::IntegrationType1D2NGauss5Ip, contactAlgorithm, 2, 1,
                                    3);

    resultDir = "./ResultsStaticFEM_IGAL_Rigid_Mortar_2E_20IP";

    ContactTestOneElementLayerSlave(resultDir,
                                    NuTo::Interpolation::eTypeOrder::LOBATTO2,
                                    3,
                                    nodeCoordinates,
                                    degree,
                                    1.e6,
                                    NuTo::eIntegrationType::IntegrationType1D2NGauss20Ip, contactAlgorithm, 2, 1,
                                    3);



//    resultDir = "./ResultsStaticFEM_IGAL_Rigid_Mortar_10E_5IP";

//    ContactTestOneElementLayerSlave(resultDir,
//                                    NuTo::Interpolation::eTypeOrder::LOBATTO2,
//                                    3,
//                                    nodeCoordinates,
//                                    degree,
//                                    1.e6,
//                                    NuTo::eIntegrationType::IntegrationType1D2NGauss5Ip, contactAlgorithm, 10, 1,
//                                    3);

//    resultDir = "./ResultsStaticFEM_IGAL_Rigid_Mortar_10E_12IP";

//    ContactTestOneElementLayerSlave(resultDir,
//                                    NuTo::Interpolation::eTypeOrder::LOBATTO2,
//                                    3,
//                                    nodeCoordinates,
//                                    degree,
//                                    1.e6,
//                                    NuTo::eIntegrationType::IntegrationType1D2NGauss12Ip, contactAlgorithm, 10, 1,
//                                    3);

//    resultDir = "./ResultsStaticFEM_IGAL_Rigid_Mortar_20E_5IP";

//    ContactTestOneElementLayerSlave(resultDir,
//                                    NuTo::Interpolation::eTypeOrder::LOBATTO2,
//                                    3,
//                                    nodeCoordinates,
//                                    degree,
//                                    1.e6,
//                                    NuTo::eIntegrationType::IntegrationType1D2NGauss5Ip, contactAlgorithm, 20, 1,
//                                    4);

//    resultDir = "./ResultsStaticFEM_IGAL_Rigid_Mortar_20E_12IP";

//    ContactTestOneElementLayerSlave(resultDir,
//                                    NuTo::Interpolation::eTypeOrder::LOBATTO2,
//                                    3,
//                                    nodeCoordinates,
//                                    degree,
//                                    1.e6,
//                                    NuTo::eIntegrationType::IntegrationType1D2NGauss12Ip, contactAlgorithm, 20, 1,
//                                    3);


    return 0;
}
