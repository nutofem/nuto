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

void SetDBCPatchTestIGAIGA(NuTo::Structure &myStructure,
                           int groupNodesSlave, int groupNodesMaster,
                           const Eigen::VectorXd &min1Coords,  const Eigen::VectorXd &max1Coords,  const Eigen::VectorXd &min2Coords,  const Eigen::VectorXd &max2Coords,
                           int &countDBC)
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
    // ===> slave master left <=== //
    xMin = 1.; xMax = 1.;
    yMin = 1.; yMax = 1.;
    int groupNodesLeft = myStructure.GroupCreate(NuTo::eGroupId::Nodes);
    myStructure.GroupAddNodeFunction(groupNodesLeft, groupNodesSlave, LambdaNodes);
    members = myStructure.GroupGetMemberIds(groupNodesLeft);
    direction << 1, 0;
    countDBC = myStructure.ConstraintLinearSetDisplacementNodeGroup(groupNodesLeft, direction, 0.0);


    // ===> PATCH TEST BOUNDARY <=== //
    double Stress = 10.;
    myStructure.SetNumLoadCases(3);
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

    // ===> master top <=== //
    // LEFT
    xMin = 0.; xMax = max1Coords(0); // depends on the refinement ....
    yMin = 0.; yMax = 0.;
    int groupNodesMasterUpperLeft = myStructure.GroupCreate(NuTo::eGroupId::Nodes);
    myStructure.GroupAddNodeFunction(groupNodesMasterUpperLeft, groupNodesMaster, LambdaNodes);
    members = myStructure.GroupGetMemberIds(groupNodesMasterUpperLeft);

    int groupElementsMasterUpperLeft = myStructure.GroupCreate(NuTo::eGroupId::Elements);
    myStructure.GroupAddElementsFromNodes(groupElementsMasterUpperLeft, groupNodesMasterUpperLeft, false);

    myStructure.LoadSurfacePressureCreate2D(1, groupElementsMasterUpperLeft, groupNodesMasterUpperLeft, Stress);
    members = myStructure.GroupGetMemberIds(groupElementsMasterUpperLeft);
    // RIGHT
    xMin = min2Coords(0); xMax = 4.;
    yMin = 0.; yMax = 0.;
    int groupNodesMasterUpperRight = myStructure.GroupCreate(NuTo::eGroupId::Nodes);
    myStructure.GroupAddNodeFunction(groupNodesMasterUpperRight, groupNodesMaster, LambdaNodes);
    members = myStructure.GroupGetMemberIds(groupNodesMasterUpperRight);

    int groupElementsMasterUpperRight = myStructure.GroupCreate(NuTo::eGroupId::Elements);
    myStructure.GroupAddElementsFromNodes(groupElementsMasterUpperRight, groupNodesMasterUpperRight, false);
    myStructure.LoadSurfacePressureCreate2D(2, groupElementsMasterUpperRight, groupNodesMasterUpperRight, Stress);
    members = myStructure.GroupGetMemberIds(groupElementsMasterUpperRight);

    // ===> initial values <=== //
    NuTo::FullVector<double,Eigen::Dynamic>  dispVec(2);
    dispVec(0) =   0.;
    dispVec(1) = -1.e-6;
    myStructure.NodeGroupSetDisplacements(groupNodesSlave, 0, dispVec);
    countDBC++;
}

void ContactTestIGA_IGA_Patch(const std::string &resultDir,
                              int rDegree,
                              double rPenalty,
                              NuTo::eIntegrationType rIntegrationType,
                              int rContactAlgo,
                              int refinementXSlave, int refinementYSlave,
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
    buildRect2DIGA(rDegree, startX, startY, height, length,
                   groupElementsSlave, groupNodesSlave,  refinementXSlave, refinementYSlave,
                   &myStructure);

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
    startX =  0.; startY = -1.;
    length =  4.; height =  1.;
    NuTo::NURBSSurface surfMaster = buildRect2DIGA(rDegree, startX, startY, height, length,
                                                   groupElementsMaster, groupNodesMaster,  refinementElXMaster,
                                                   refinementElYMaster, &myStructure);

    Eigen::VectorXd max1Coords = surfMaster.GetControlPoint(0, rDegree + std::pow(2,refinementElYMaster-1) - 1);
    Eigen::VectorXd min2Coords = surfMaster.GetControlPoint(0, rDegree + 2*std::pow(2,refinementElYMaster-1)-2);

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
    SetDBCPatchTestIGAIGA(myStructure, groupNodesSlave, groupNodesMaster,  min1Coords, max1Coords, min2Coords, max2Coords, countDBC);

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

    myStructure.ConstitutiveLawSetParameterFunction(constitutiveLawPC,
                                                    NuTo::Constitutive::eConstitutiveParameter::CONSTITUTIVE_LAW_FUNCTION,
                                                    constitutiveContactLaw);
    myStructure.ConstitutiveLawSetParameterFunction(constitutiveLawPC,
                                                    NuTo::Constitutive::eConstitutiveParameter::CONSTITUTIVE_LAW_DERIVATIVE_FUNCTION,
                                                    constitutiveContactLawDerivative);

    myStructure.NuTo::Structure::ContactElementsCreate<2,2>(groupElementsSlaveLower, groupNodesSlaveLower, groupElementsMasterLower,
                                                            groupNodesMasterLower, rIntegrationType, rContactAlgo, constitutiveLawPC);

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
    std::string resultDir = "./ResultsStatic_IGA_IGA_Patch_Mortar";
    int contactAlgorithm = 0;
    int degree = 2;
    int factor = 2;
    double penalty = 1.e9;

    ContactTestIGA_IGA_Patch(resultDir,
                             degree,
                             penalty,
                             NuTo::eIntegrationType::IntegrationType1D2NGauss12Ip, contactAlgorithm, 2*factor, 1*factor, 1*(factor+2), 1*factor);

    resultDir = "./ResultsStatic_IGA_IGA_Patch_NonMortar";
    contactAlgorithm = 1;
    ContactTestIGA_IGA_Patch(resultDir,
                             degree,
                             penalty,
                             NuTo::eIntegrationType::IntegrationType1D2NGauss12Ip, contactAlgorithm, 2*factor, 1*factor, 1*(factor+2), 1*factor);

    return 0;
}
