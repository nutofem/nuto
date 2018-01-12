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

//    std::cout << "Control point:\t" << std::endl;
//    for(int i = 0; i < surface.GetNumControlPoints(1); i++)
//    {
//        std::cout << std::endl;
//        for(int j = 0; j < surface.GetNumControlPoints(0); j++)
//            std::cout << "(" << surface.GetControlPoint(i,j).transpose() << ")\t";
//    }

//    std::cout << std::flush;

    std::set<NuTo::Node::eDof> setOfDOFS;
    setOfDOFS.insert(NuTo::Node::eDof::COORDINATES);
    setOfDOFS.insert(NuTo::Node::eDof::DISPLACEMENTS);

    surface.buildIGAStructure(*myStructure, setOfDOFS, groupElements, groupNodes);

    return surface;
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
    dispVec(1) = -1.e-5;
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

void ContactTestRigidIGA_L_IGA(const std::string &resultDir,
                               int rDegree,
                               double rPenalty,
                               NuTo::eIntegrationType rIntegrationType,
                               int rContactAlgo,
                               int numElXSlave, int numElYSlave)
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
    buildRect2DIGA(rDegree, 1., 0., 1., 1., groupElementsSlave, groupNodesSlave,  numElXSlave, numElYSlave, &myStructure);

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
//    std::cout << "Members Nodes Slave lower: \n" << members << std::endl;

    int groupElementsSlaveLower = myStructure.GroupCreate(NuTo::eGroupId::Elements);
    myStructure.GroupAddElementsFromNodes(groupElementsSlaveLower, groupNodesSlaveLower, false);
    members = myStructure.GroupGetMemberIds(groupElementsSlaveLower);
//    std::cout << "Members Elements Slave lower: \n" << members << std::endl;

    //////////////////////////////
    // ====> create MASTER mesh //
    //////////////////////////////
    // create IGA curve
    int numPoints = 10;
    double Length = 3.;
    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> points(numPoints,2);
    for(int i = 0; i < numPoints; i++)
    {
        points(i,0) = i*Length/(numPoints-1);
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
    Eigen::Matrix<std::pair<int, int>,Eigen::Dynamic,Eigen::Dynamic> elementsMaster = curve.buildIGAStructure(myStructure, setOfDOFS, groupElementsIGA, groupNodesIGA, "IGA1DLAYER", nodeIDs);

    int countDBC;
//    SetDBC(myStructure, groupNodesSlave, groupNodesIGA,  countDBC);
//    SetDBCPatchTest(myStructure, groupNodesSlave, groupNodesIGA,  countDBC);
    SetDBCPatchTestRigidIGA(myStructure, groupNodesSlave, -1, groupNodesIGA,  countDBC);

    ///////////////////
    // ===> material //
    ///////////////////

    double Thickness = 1.;
    int section = myStructure.SectionCreate("PLANE_STRESS");
    myStructure.SectionSetThickness(section, Thickness);

    double E = 1.e5;
    double nue = 0.;
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

    myStructure.NuTo::Structure::ContactElementsCreate<2,1>(groupElementsSlaveLower, groupNodesSlaveLower, elementsMaster, rIntegrationType, rContactAlgo, constitutiveLawPC);

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
    int visualizationGroup = groupElementsSlave;
    myStructure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::DISPLACEMENTS);
    myStructure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ENGINEERING_STRAIN);
    myStructure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ENGINEERING_STRESS);
    myStructure.ElementGroupExportVtkDataFile(visualizationGroup, resultDir+"/Elements.vtu", true);

    myStructure.AddVisualizationComponent(groupElementsIGA, NuTo::eVisualizeWhat::DISPLACEMENTS);
    myStructure.ElementGroupExportVtkDataFile(groupElementsIGA, resultDir+"/ElementsLayer.vtu", true);
#endif
}


int main()
{
    std::string resultDir = "";
    resultDir = "./ResultsStaticIGA_IGA_L_Rigid";

    int factor = 3;
    int degree = 2;
    int contactAlgorithm = 0;

    ContactTestRigidIGA_L_IGA(resultDir,
                              degree,
                              1.e6,
                              NuTo::eIntegrationType::IntegrationType1D2NGauss12Ip, contactAlgorithm, 1*factor, 1*factor);


    return 0;
}
