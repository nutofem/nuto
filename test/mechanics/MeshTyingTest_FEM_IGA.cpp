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

//            NuTo::ElementBase* el = myStructure->ElementGetElementPtr(myElement);
//            Eigen::VectorXd nodalCurrent;
//            nodalCurrent.setZero(0);
//            nodalCurrent = el->ExtractNodeValues(0, NuTo::Node::eDof::DISPLACEMENTS);
//            nodalCurrent = el->ExtractNodeValues(0, NuTo::Node::eDof::COORDINATES);

//            std::cout << "Element: " << myElement << ", Nodes: " << elementIncidence.Trans() << std::endl;
            myStructure->GroupAddElement(rElementGroup, myElement);
        }
    }
    myStructure->ElementConvertToInterpolationType(rElementGroup, 1.e-6,10);
    return node;
}


NuTo::NURBSSurface buildRect2DIGA(int rDegree, double x0, double y0, double Height, double Length,
                                  int &groupElements, int &groupNodes, int refinementX, int refinementY,
                                  NuTo::Structure *myStructure)
{
    /** Knots and control points **/
    int numElementsX = 1;
    int numElementsY = 1;

    Eigen::Vector2i degree(rDegree,rDegree);

    int numKnotsX = 2*(degree(0)+1) + numElementsX - 1;
    int numKnotsY = 2*(degree(1)+1) + numElementsY - 1;

    Eigen::VectorXd knotsX(numKnotsX);
    Eigen::VectorXd knotsY(numKnotsY);

    for (int i = 0; i <= degree(0); i++) knotsX(i) = 0.;
    for (int i = degree(0) + 1; i <= degree(0) + numElementsX - 1 ; i++) knotsX(i) = knotsX(i-1) + 1./numElementsX;
    for (int i = degree(0) + numElementsX ; i < numKnotsX; i++) knotsX(i) = 1;

    for (int i = 0; i <= degree(1); i++) knotsY(i) = 0.;
    for (int i = degree(1) + 1; i <= degree(1) + numElementsY - 1 ; i++) knotsY(i) = knotsY(i-1) + 1./numElementsY;
    for (int i = degree(1) + numElementsY ; i < numKnotsY; i++) knotsY(i) = 1;

    int numControlPointsX = (numKnotsX -1) - degree(0);
    int numControlPointsY = (numKnotsY -1) - degree(1);

    Eigen::Matrix<Eigen::VectorXd, Eigen::Dynamic, Eigen::Dynamic> controlPoints(numControlPointsY, numControlPointsX);

    double incrx = (Length/(numControlPointsX-1));
    double incry = (Height/(numControlPointsY-1));
    for(int i = 0; i < numControlPointsY; i++)
    {
        for(int j = 0; j < numControlPointsX; j++)
        {
            controlPoints(i,j) = Eigen::Vector2d(x0 + incrx*j, y0 + incry*i);
        }
    }

    Eigen::MatrixXd weights(numControlPointsY, numControlPointsX);
    weights.setOnes(numControlPointsY, numControlPointsX);

    NuTo::NURBSSurface surface(degree, knotsX, knotsY, controlPoints, weights);

    for(int i = 0; i < refinementX-1; i++) surface.DuplicateKnots(0);
    for(int i = 0; i < refinementY-1; i++) surface.DuplicateKnots(1);

    std::cout << "Control point:\t" << std::endl;
    for(int i = 0; i < surface.GetNumControlPoints(1); i++)
    {
        std::cout << std::endl;
        for(int j = 0; j < surface.GetNumControlPoints(0); j++)
            std::cout << "(" << surface.GetControlPoint(i,j).transpose() << ")\t";
    }

    std::cout << std::flush;

    std::set<NuTo::Node::eDof> setOfDOFS;
    setOfDOFS.insert(NuTo::Node::eDof::COORDINATES);
    setOfDOFS.insert(NuTo::Node::eDof::DISPLACEMENTS);

    surface.buildIGAStructure(*myStructure, setOfDOFS, groupElements, groupNodes);

    return surface;
}

void SetDBCPatchTestIGAIGA(NuTo::Structure &myStructure, int groupNodesSlave, int groupNodesMaster, int &countDBC)
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

            if (x >= xMin - Tol && x <= xMax + Tol)
            {
                xR = true;
            }
            if (y >= yMin - Tol && y <= yMax + Tol)
            {
                yR = true;
            }

            if (xR == true && yR == true) return true;
        }
        return false;
    };

    // ===> master bottom <=== //
    xMin = 0.; xMax = std::numeric_limits<double>::infinity();
    yMin = yMax = -1.;
    int groupNodesMasterLower = myStructure.GroupCreate(NuTo::eGroupId::Nodes);
    myStructure.GroupAddNodeFunction(groupNodesMasterLower, LambdaNodes);
    direction << 0.,1.;
    countDBC = myStructure.ConstraintLinearSetDisplacementNodeGroup(groupNodesMasterLower, direction, 0.0);

    NuTo::FullVector<int,Eigen::Dynamic> rMembersMasterBottom;
    myStructure.NodeGroupGetMembers(groupNodesMasterLower, rMembersMasterBottom);
    direction << 1.,0.;
    countDBC = myStructure.ConstraintLinearSetDisplacementNode(rMembersMasterBottom(0), direction, 0.0);

    NuTo::FullVector<int, Eigen::Dynamic> members;

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

    countDBC++;
}

void MeshTyingIGA_IGA_Patch(const std::string &resultDir,
                            int rDegree,
                            double rPenalty,
                            NuTo::Interpolation::eTypeOrder rElementTypeIdent,
                            int rNumNodesPerElementInOneDir,
                            NuTo::FullVector<double, Eigen::Dynamic>& nodeCoordinatesFirstElement,
                            NuTo::eIntegrationType rIntegrationType,
                            int rContactAlgo,
                            int numElXSlave, int numElYSlave,
                            int refinementElXMaster, int refinementElYMaster)
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
    int groupNodesSlave = myStructure.GroupCreate(NuTo::eGroupId::Nodes);
    int groupElementsSlave = myStructure.GroupCreate(NuTo::eGroupId::Elements);
    double startX = 1.; double startY = 0.;
    double length = 1.; double height = 1.;

    std::set<NuTo::Node::eDof> setOfDOFSSlave;
    setOfDOFSSlave.insert(NuTo::Node::eDof::COORDINATES);
    setOfDOFSSlave.insert(NuTo::Node::eDof::DISPLACEMENTS);

    buildStructure2D(rElementTypeIdent, rNumNodesPerElementInOneDir, nodeCoordinatesFirstElement,
                     numElXSlave, numElYSlave, height, length, startX, startY, 0, &myStructure, groupNodesSlave, groupElementsSlave, setOfDOFSSlave);

    // ===> build contact elements slave elements
    auto LambdaGetSlaveNodesLower = [](NuTo::NodeBase* rNodePtr) -> bool
    {
        double Tol = 1.e-6;
        if (rNodePtr->GetNum(NuTo::Node::eDof::COORDINATES)>0)
        {
            double y = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[1];
            if (y >= 0. - Tol && y <= 0. + Tol)
            {
                return true;
            }
        }
        return false;
    };  // LambdaGetSlaveNodesLower

    int groupNodesSlaveLower = myStructure.GroupCreate(NuTo::eGroupId::Nodes);
    myStructure.GroupAddNodeFunction(groupNodesSlaveLower, groupNodesSlave, LambdaGetSlaveNodesLower);

    NuTo::FullMatrix<int, Eigen::Dynamic> members = myStructure.GroupGetMemberIds(groupNodesSlaveLower);
    std::cout << "Members Nodes Slave lower: \n" << members << std::endl;

    int groupElementsSlaveLower = myStructure.GroupCreate(NuTo::eGroupId::Elements);
    myStructure.GroupAddElementsFromNodes(groupElementsSlaveLower, groupNodesSlaveLower, false);
    members = myStructure.GroupGetMemberIds(groupElementsSlaveLower);
    std::cout << "Members Elements Slave lower: \n" << members << std::endl;

    //////////////////////////////
    // ====> create MASTER mesh //
    //////////////////////////////
    int groupNodesMaster = myStructure.GroupCreate(NuTo::eGroupId::Nodes);
    int groupElementsMaster = myStructure.GroupCreate(NuTo::eGroupId::Elements);
    startX =  1.; startY = -1.;
    length =  1.; height =  1.;
    buildRect2DIGA(rDegree, startX, startY, height, length,
                   groupElementsMaster, groupNodesMaster,  refinementElXMaster,
                   refinementElYMaster, &myStructure);

    int groupNodesMasterLower = myStructure.GroupCreate(NuTo::eGroupId::Nodes);
    myStructure.GroupAddNodeFunction(groupNodesMasterLower, groupNodesMaster, LambdaGetSlaveNodesLower);

    members = myStructure.GroupGetMemberIds(groupNodesMasterLower);
    std::cout << "Members Nodes Master lower: \n" << members << std::endl;

    int groupElementsMasterLower = myStructure.GroupCreate(NuTo::eGroupId::Elements);
    myStructure.GroupAddElementsFromNodes(groupElementsMasterLower, groupNodesMasterLower, false);
    members = myStructure.GroupGetMemberIds(groupElementsMasterLower);
    std::cout << "Members Elements Master lower: \n" << members << std::endl;

    int countDBC;
    Eigen::VectorXd min1Coords(2,1);
    min1Coords << 0,0;
    Eigen::VectorXd max2Coords(2,1);
    max2Coords << 4,0;
    SetDBCPatchTestIGAIGA(myStructure, groupNodesSlave, groupNodesMaster, countDBC);

    ///////////////////
    // ===> material //
    ///////////////////

    double Thickness = 1.;
    int section = myStructure.SectionCreate("PLANE_STRESS");
    myStructure.SectionSetThickness(section, Thickness);

    double E   = 1.e5;
    double nue = 0.3;

    int constitutiveLaw = myStructure.ConstitutiveLawCreate("Linear_Elastic_Engineering_Stress");
    myStructure.ConstitutiveLawSetParameterDouble(constitutiveLaw, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, E);
    myStructure.ConstitutiveLawSetParameterDouble(constitutiveLaw, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, nue);
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

    myStructure.ConstitutiveLawSetParameterFunction(constitutiveLawPC,
                                                    NuTo::Constitutive::eConstitutiveParameter::CONSTITUTIVE_LAW_FUNCTION,
                                                    constitutiveContactLaw);
    myStructure.ConstitutiveLawSetParameterFunction(constitutiveLawPC,
                                                    NuTo::Constitutive::eConstitutiveParameter::CONSTITUTIVE_LAW_DERIVATIVE_FUNCTION,
                                                    constitutiveContactLawDerivative);

    // *** mesh tying *** //

    int cceID = myStructure.ContactElementsCreate<2,2>(groupElementsSlaveLower, groupNodesSlaveLower,
                                                       groupElementsMasterLower, groupNodesMasterLower,
                                                       rIntegrationType, rContactAlgo, constitutiveLawPC);

    NuTo::ElementBase* elementBase = myStructure.ElementGetElementPtr(cceID);
    NuTo::ContinuumContactElement<2,2> &contactElement = elementBase->AsContinuumContactElement22();

    Eigen::MatrixXd D, M;
    std::unordered_map<int, int> mappingGlobal2LocalSlaveNode, mappingGlobal2LocalMasterNode;
    contactElement.ComputeMeshTyingMatrix(D, M, mappingGlobal2LocalSlaveNode, mappingGlobal2LocalMasterNode);

    std::cout << D << std::endl;

    std::cout << "===============================" << std::endl;

    std::cout << M << std::endl;

    for(int dof = 0; dof < myStructure.GetDimension(); dof++)
    {
        for(auto &itSlaveRow : mappingGlobal2LocalSlaveNode) // global
        {
            int localSlaveNodeIdRow  = itSlaveRow.second;

            std::unordered_map<int, int>::iterator itSlaveCol = mappingGlobal2LocalSlaveNode.begin();
            int localSlaveNodeIdCol = itSlaveCol->second;
            int globalSlaveNodeId = itSlaveCol->first;
            // define equation
            myStructure.ConstraintLinearEquationCreate(countDBC, globalSlaveNodeId, NuTo::Node::eDof::DISPLACEMENTS, dof, D(localSlaveNodeIdRow, localSlaveNodeIdCol), 0.);
            std::cout << D(localSlaveNodeIdRow, localSlaveNodeIdCol) << " ";
            itSlaveCol++;

            for(;itSlaveCol != mappingGlobal2LocalSlaveNode.end() ;itSlaveCol++)
            {
                localSlaveNodeIdCol  = itSlaveCol->second;
                globalSlaveNodeId = itSlaveCol->first;
                // add slave ingredients
                myStructure.ConstraintLinearEquationAddTerm(countDBC, globalSlaveNodeId, NuTo::Node::eDof::DISPLACEMENTS, dof, D(localSlaveNodeIdRow, localSlaveNodeIdCol));
                std::cout << D(localSlaveNodeIdRow, localSlaveNodeIdCol) << " ";
                // master side
            }

            std::cout << " || ";
            for(auto &itMaster : mappingGlobal2LocalMasterNode)
            {
                int globalMasterNodeId = itMaster.first;
                int localMasterNodeIdRow  = itMaster.second;
                // add master ingredients
                myStructure.ConstraintLinearEquationAddTerm(countDBC, globalMasterNodeId, NuTo::Node::eDof::DISPLACEMENTS, dof, M(localSlaveNodeIdRow, localMasterNodeIdRow));

                std::cout << M(localSlaveNodeIdRow, localMasterNodeIdRow) << " ";
            }
            std::cout << std::endl;

            countDBC++;
        }
    }

    myStructure.ElementDelete(cceID);

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

    NuTo::NewmarkDirect myIntegrationScheme(&myStructure);
    double timeStep = 1.;
    double simulationTime = 1.;

    myIntegrationScheme.SetResultDirectory(resultDir, true);

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
    std::string resultDir = "./ResultsStatic_Tying_FEM_IGA_Patch_Mortar";
    int contactAlgorithm = 0;
    int degree = 2;
    double penalty = 1.e9;

    // LOBATTO 2
    NuTo::FullVector<double, Eigen::Dynamic> nodeCoordinates(3);

    nodeCoordinates.resize(3);
    nodeCoordinates(0) = 0;
    nodeCoordinates(1) = 1;
    nodeCoordinates(2) = 2;

    MeshTyingIGA_IGA_Patch(resultDir,
                           degree,
                           penalty,
                           NuTo::Interpolation::eTypeOrder::LOBATTO2,
                           3,
                           nodeCoordinates,
                           NuTo::eIntegrationType::IntegrationType1D2NGauss5Ip, contactAlgorithm, 5, 5, 3, 3);

//    resultDir = "./ResultsStatic_Tying_IGA_IGA_Patch_NonMortar";
//    contactAlgorithm = 1;
//    MeshTyingIGA_IGA_Patch(resultDir,
//                             degree,
//                             penalty,
//                             NuTo::eIntegrationType::IntegrationType1D2NGauss12Ip, contactAlgorithm, 1, 1, 2, 1);

    return 0;
}
