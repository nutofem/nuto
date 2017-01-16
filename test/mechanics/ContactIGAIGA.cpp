#include <iostream>
#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"

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

#include "nuto/mechanics/constitutive/ConstitutiveEnum.h"
#include "nuto/mechanics/interpolationtypes/InterpolationTypeEnum.h"
#include "nuto/mechanics/elements/ElementDataEnum.h"
#include "nuto/mechanics/elements/IpDataEnum.h"

#include "nuto/mechanics/IGA/NURBSCurve.h"
#include "nuto/mechanics/IGA/NURBSSurface.h"

#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/nodes/NodeEnum.h"

#include "nuto/mechanics/groups/GroupEnum.h"

#include <boost/filesystem.hpp>

#include "nuto/mechanics/timeIntegration/NewmarkDirect.h"

#ifdef ENABLE_VISUALIZE
#include "nuto/visualize/VisualizeEnum.h"
#endif

NuTo::NURBSSurface buildRect2D(double x0, double y0, double Height, double Length, int &groupElements, int &groupNodes)
{
    /** Knots and control points **/
    int numElementsX = 1;
    int numElementsY = 1;

    Eigen::Vector2i degree(2,2);

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

    for(int i = 0; i < 2; i++)
    {
        surface.DuplicateKnots(0);
        surface.DuplicateKnots(1);
    }

    std::set<NuTo::Node::eDof> setOfDOFS;
    setOfDOFS.insert(NuTo::Node::eDof::COORDINATES);
    setOfDOFS.insert(NuTo::Node::eDof::DISPLACEMENTS);

    surface.buildIGAStructure(*myStructure, setOfDOFS, groupElements, groupNodes);
}


void SetDBCPatchTest(NuTo::Structure &myStructure, int groupNodesSlave, int groupNodesMaster, int &countDBC)
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

    // ===> initial values <=== //
    NuTo::FullVector<double,Eigen::Dynamic>  dispVec(2);
    dispVec(0) =   0.;
    dispVec(1) = -0.000001;
    myStructure.NodeGroupSetDisplacements(groupNodesSlave, 0, dispVec);

    // ===> slave master left <=== //
    xMin = 1.; xMax = 1.;
    yMin = 1.; yMax = 1.;
    int groupNodesLeft = myStructure.GroupCreate(NuTo::eGroupId::Nodes);
    myStructure.GroupAddNodeFunction(groupNodesLeft, groupNodesSlave, LambdaNodes);
    direction << 1, 0;
    countDBC = myStructure.ConstraintLinearSetDisplacementNodeGroup(groupNodesLeft, direction, 0.0);

    // ===> slave master right <=== //
//    xMin = 2.; xMax = 2.;
//    yMin = 1.; yMax = 1.;
//    int groupNodesRight = myStructure.GroupCreate(NuTo::eGroupId::Nodes);
//    myStructure.GroupAddNodeFunction(groupNodesRight, groupNodesSlave, LambdaNodes);
//    direction << 1, 0;
//    countDBC = myStructure.ConstraintLinearSetDisplacementNodeGroup(groupNodesRight, direction, 0.0);

    // ===> PATCH TEST BOUNDARY <=== //
    double Stress = 10.;
    myStructure.SetNumLoadCases(3);
    // ===> slave top <=== //
    xMin = 1.; xMax = 2.;
    yMin = 1.; yMax = 1.;
    int groupNodesSlaveUpper = myStructure.GroupCreate(NuTo::eGroupId::Nodes);
    myStructure.GroupAddNodeFunction(groupNodesSlaveUpper, LambdaNodes);
    int groupElementsSlaveUpper = myStructure.GroupCreate(NuTo::eGroupId::Elements);
    myStructure.GroupAddElementsFromNodes(groupElementsSlaveUpper, groupNodesSlaveUpper, false);
    myStructure.LoadSurfacePressureCreate2D(0, groupElementsSlaveUpper, groupNodesSlaveUpper, Stress);

    // ===> master top <=== //
    NuTo::FullVector<int, Eigen::Dynamic> members;
    // LEFT
    xMin = 0.; xMax = 1.;
    yMin = 0.; yMax = 0.;
    int groupNodesMasterUpperLeft = myStructure.GroupCreate(NuTo::eGroupId::Nodes);
    myStructure.GroupAddNodeFunction(groupNodesMasterUpperLeft, groupNodesMaster, LambdaNodes);
    members = myStructure.GroupGetMemberIds(groupNodesMasterUpperLeft);
    int groupElementsMasterUpperLeft = myStructure.GroupCreate(NuTo::eGroupId::Elements);
    myStructure.GroupAddElementsFromNodes(groupElementsMasterUpperLeft, groupNodesMasterUpperLeft, false);
    myStructure.LoadSurfacePressureCreate2D(1, groupElementsMasterUpperLeft, groupNodesMasterUpperLeft, Stress);
    members = myStructure.GroupGetMemberIds(groupElementsMasterUpperLeft);
    // RIGHT
    xMin = 2.0; xMax = 3.0;
    yMin = 0.; yMax = 0.;
    int groupNodesMasterUpperRight = myStructure.GroupCreate(NuTo::eGroupId::Nodes);
    myStructure.GroupAddNodeFunction(groupNodesMasterUpperRight, groupNodesMaster, LambdaNodes);
    members = myStructure.GroupGetMemberIds(groupNodesMasterUpperRight);
    int groupElementsMasterUpperRight = myStructure.GroupCreate(NuTo::eGroupId::Elements);
    myStructure.GroupAddElementsFromNodes(groupElementsMasterUpperRight, groupNodesMasterUpperRight, false);
    myStructure.LoadSurfacePressureCreate2D(2, groupElementsMasterUpperRight, groupNodesMasterUpperRight, Stress);
    members = myStructure.GroupGetMemberIds(groupElementsMasterUpperRight);

    countDBC++;
}

void SetDBC(NuTo::Structure &myStructure, int groupNodesSlave, int groupNodesMaster, int &countDBC)
{
    Eigen::Vector2d direction(0,0);
    double xValue = 0.;
    double yValue = 0.;
    auto LambdaNodes = [&](NuTo::NodeBase* rNodePtr) -> bool
    {
        double Tol = 1.e-6;
        if (rNodePtr->GetNum(NuTo::Node::eDof::COORDINATES)>0)
        {
            double x = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[0];
            double y = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[1];
            bool xR = false;
            bool yR = false;

            if      (xValue == std::numeric_limits<double>::infinity()) xR = true;
            else if (x >= xValue - Tol && x <= xValue + Tol) xR = true;

            if      (yValue == std::numeric_limits<double>::infinity()) yR = true;
            else if (y >= yValue - Tol && y <= yValue + Tol) yR = true;

            if (xR == true && yR == true) return true;
        }
        return false;
    };

    // ===> master bottom <=== //
    xValue = std::numeric_limits<double>::infinity();
    yValue = -1.;
    int groupNodesMasterLower = myStructure.GroupCreate(NuTo::eGroupId::Nodes);
    myStructure.GroupAddNodeFunction(groupNodesMasterLower, LambdaNodes);
    direction << 0.,1.;
    countDBC = myStructure.ConstraintLinearSetDisplacementNodeGroup(groupNodesMasterLower, direction, 0.0);

    NuTo::FullVector<int,Eigen::Dynamic> rMembersMasterBottom;
    myStructure.NodeGroupGetMembers(groupNodesMasterLower, rMembersMasterBottom);
    direction << 1.,0.;
    countDBC = myStructure.ConstraintLinearSetDisplacementNode(rMembersMasterBottom(0), direction, 0.0);

    // ===> slave top <=== //
    xValue = std::numeric_limits<double>::infinity();
    yValue = 1.;
    int groupNodesSlaveUpper = myStructure.GroupCreate(NuTo::eGroupId::Nodes);
    myStructure.GroupAddNodeFunction(groupNodesSlaveUpper, LambdaNodes);
    double disp = - 0.00008;
    direction << 0.,1.;
    countDBC = myStructure.ConstraintLinearSetDisplacementNodeGroup(groupNodesSlaveUpper, direction, disp);

    // ===> slave master left <=== //
    xValue = 1.;
    yValue = 1.;
    int groupNodesLeft = myStructure.GroupCreate(NuTo::eGroupId::Nodes);
    myStructure.GroupAddNodeFunction(groupNodesLeft, groupNodesSlave, LambdaNodes);
    direction << 1, 0;
    countDBC = myStructure.ConstraintLinearSetDisplacementNodeGroup(groupNodesLeft, direction, 0.0);

//    // ===> slave master right <=== //
//    xValue = 1.;
//    yValue = std::numeric_limits<double>::infinity();
//    int groupNodesRight = myStructure.GroupCreate(NuTo::eGroupId::Nodes);
//    myStructure.GroupAddNodeFunction(groupNodesRight, LambdaNodes);
//    direction << 1, 0;
//    countDBC = myStructure.ConstraintLinearSetDisplacementNodeGroup(groupNodesRight, direction, 0.0);

    countDBC++;

}

void AddIGALayer(NuTo::Structure *myStructure,
                const std::function<bool(NuTo::NodeBase *)> &rFunction,
                int rNodesGroupId,
                int &groupElementsIGAlayer,
                int &groupNodesIGAlayer,
                Eigen::Matrix<std::pair<int, int>, Eigen::Dynamic, Eigen::Dynamic> &rElements,
                int &countDBC)
{
    // Nodes on the part to interpolate
    int groupFENodes = myStructure->GroupCreate(NuTo::eGroupId::Nodes);
    myStructure->GroupAddNodeFunction(groupFENodes, rNodesGroupId, rFunction);
    NuTo::FullVector<int,Eigen::Dynamic> idsFENodes;
    myStructure->NodeGroupGetMembers(groupFENodes, idsFENodes);

    // Corresponding elements
    int groupFE = myStructure->GroupCreate(NuTo::eGroupId::Elements);
    myStructure->GroupAddElementsFromNodes(groupFE, groupFENodes, false);

    // Matrix containing the ids and coordinates of the FE nodes => 'coordinatesAndIDs'
    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> coordinates;
    myStructure->NodeGroupGetCoordinates(groupFENodes, coordinates);

    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> coordinatesAndIDs;
    coordinatesAndIDs.resize(coordinates.rows(), coordinates.cols() + 1);
    coordinatesAndIDs.block(0,0,coordinates.rows(), coordinates.cols()) = coordinates;

    for(int i = 0; i < coordinates.rows(); i++)
        coordinatesAndIDs(i,coordinates.cols()) = idsFENodes(i);

    coordinatesAndIDs.SortRow(0);

    // Interpolations
    NuTo::NURBSCurve curve(rDegree, coordinatesAndIDs.block(0, 0, coordinates.rows(), coordinates.cols()), A);

    // Build the IGA layer
    std::set<NuTo::Node::eDof> setOfDOFS;
    setOfDOFS.insert(NuTo::Node::eDof::COORDINATES);
    setOfDOFS.insert(NuTo::Node::eDof::DISPLACEMENTS);

    groupNodesIGAlayer = myStructure->GroupCreate("Nodes");
    groupElementsIGAlayer  = myStructure->GroupCreate("Elements");

    rElements = curve.buildIGAStructure(*myStructure, setOfDOFS, groupElementsIGAlayer, groupNodesIGAlayer, "IGA1DLAYER");

    // Matrix containing the ids and coordinates of the IGA layer control points (aka. nodes) => 'coordinatesAndIDsLayer'
    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> coordinatesLayer;
    myStructure->NodeGroupGetCoordinates(groupNodesIGAlayer, coordinatesLayer);

    NuTo::FullVector<int,Eigen::Dynamic> rMembersMasterContactBoundaryLayer;
    myStructure->NodeGroupGetMembers(groupNodesIGAlayer, rMembersMasterContactBoundaryLayer);

    // all together
    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> coordinatesAndIDsLayer;
    coordinatesAndIDsLayer.resize(coordinatesLayer.rows(), coordinatesLayer.cols() + 1);
    coordinatesAndIDsLayer.block(0,0,coordinatesLayer.rows(), coordinatesLayer.cols()) = coordinatesLayer;

    for(int i = 0; i < coordinatesLayer.rows(); i++)
        coordinatesAndIDsLayer(i,coordinatesLayer.cols()) = rMembersMasterContactBoundaryLayer(i);

    coordinatesAndIDsLayer.SortRow(0);

    // build the FE - IGA coupling
    int dim = coordinatesAndIDs.cols() - 1;
    for(int node = 0; node < coordinatesAndIDs.rows(); node++)
    {
        for(int dof = 0; dof < dim; dof++)
        {
            myStructure->ConstraintLinearEquationCreate (countDBC, coordinatesAndIDs(node, dim), NuTo::Node::eDof::DISPLACEMENTS, dof,  1., 0.);
            for(int controlPoint = 0; controlPoint < coordinatesAndIDsLayer.rows(); controlPoint++)
            {
                myStructure->ConstraintLinearEquationAddTerm(countDBC, coordinatesAndIDsLayer(controlPoint, dim), NuTo::Node::eDof::DISPLACEMENTS, dof, -A(node, controlPoint));
            }
            countDBC++;
        }
    }
}

void constantStress(const std::string &resultDir)
{
    /** Structure 2D **/
    NuTo::Structure* myStructure = new NuTo::Structure(2);

#ifdef _OPENMP
    int numThreads = 1;
    myStructure->SetNumProcessors(numThreads);
#endif

    // ===> Slave <=== //
    int groupNodesSlave = myStructure.GroupCreate(NuTo::eGroupId::Nodes);
    int groupElementsSlave = myStructure.GroupCreate(NuTo::eGroupId::Elements);
    buildRect2D(1, 0, 1, 1, groupNodesSlave, groupElementsSlave);

    // ===> Master <=== //
    int groupNodesMaster = myStructure.GroupCreate(NuTo::eGroupId::Nodes);
    int groupElementsMaster = myStructure.GroupCreate(NuTo::eGroupId::Elements);
    buildRect2D(0, -1, 1, 3, groupNodesMaster, groupElementsMaster);

    // ===> Boundary <=== //
    int countDBC;
    SetDBC(myStructure, groupNodesSlave, groupNodesMaster,  countDBC);

    // ===> Slave side <=== //
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

    int groupNodesSlaveLower = myStructure->GroupCreate(NuTo::eGroupId::Nodes);
    myStructure->GroupAddNodeFunction(groupNodesSlaveLower, groupNodesSlave, LambdaGetSlaveNodesLower);

    int groupElementsSlaveLower = myStructure->GroupCreate(NuTo::eGroupId::Elements);
    myStructure->GroupAddElementsFromNodes(groupElementsSlaveLower, groupNodesSlaveLower, false);

    // ===> Master Side <=== //

    auto LambdaGetMasterNodesUpper = [](NuTo::NodeBase* rNodePtr) -> bool
    {
        double Tol = 1.e-6;
        if (rNodePtr->GetNum(NuTo::Node::eDof::DISPLACEMENTS)>0)
        {
            double x = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[0];
            double y = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[1];
            if (y >= -Tol && y <= Tol && x >= 0.5 - Tol && x <= 2.5 + Tol)
            {
                return true;
            }
        }
        return false;
    };

    int groupNodesMasterLower = myStructure->GroupCreate(NuTo::eGroupId::Nodes);
    myStructure->GroupAddNodeFunction(groupNodesMasterLower, groupNodesMaster, LambdaGetMasterNodesUpper);

    int groupElementsMasterLower = myStructure->GroupCreate(NuTo::eGroupId::Elements);
    myStructure->GroupAddElementsFromNodes(groupElementsMasterLower, groupNodesMasterLower, false);



    // ===> Material <=== //
    double  YoungsModulus = 20000.;
    double  PoissonRatio = 0.3;
    // there should be no dependency, because the pressure at the boundary is predefined
    double  Thickness = 1.;
    int mySection = myStructure->SectionCreate("Plane_Stress");
    myStructure->SectionSetThickness(mySection,Thickness);
    int myMatLin = myStructure->ConstitutiveLawCreate("LINEAR_ELASTIC_ENGINEERING_STRESS");
    myStructure->ConstitutiveLawSetParameterDouble(myMatLin,NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS,YoungsModulus);
    myStructure->ConstitutiveLawSetParameterDouble(myMatLin,NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO,PoissonRatio);
    myStructure->ElementTotalSetConstitutiveLaw(myMatLin);
    myStructure->ElementTotalSetSection(mySection);

    // ===> Contact element <=== //
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

    myStructure->ConstitutiveLawSetParameterFunction(constitutiveLawPC, NuTo::Constitutive::eConstitutiveParameter::CONSTITUTIVE_LAW_FUNCTION, constitutiveContactLaw);
    myStructure->ConstitutiveLawSetParameterFunction(constitutiveLawPC, NuTo::Constitutive::eConstitutiveParameter::CONSTITUTIVE_LAW_DERIVATIVE_FUNCTION, constitutiveContactLawDerivative);

    myStructure->NuTo::Structure::ContactElementsCreate<2,2>(groupElementsSlaveLower, groupNodesSlaveLower, elementsMaster, rIntegrationType, rContactAlgo, constitutiveLawPC);

    // ===> Solution <=== //

    if (boost::filesystem::exists(resultDir))
    {
        if (boost::filesystem::is_directory(resultDir))
        {
            boost::filesystem::remove_all(resultDir);
        }
    }

    // create result directory
    boost::filesystem::create_directory(resultDir);

    myStructure->CalculateMaximumIndependentSets();
    myStructure->NodeBuildGlobalDofs();

    NuTo::NewmarkDirect myIntegrationScheme(myStructure);
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
    int visualizationGroup = groupElements;
    myStructure->AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::DISPLACEMENTS);
    myStructure->AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ENGINEERING_STRAIN);
    myStructure->AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ENGINEERING_STRESS);
    myStructure->ElementGroupExportVtkDataFile(visualizationGroup, resultDir+"/Elements.vtu", true);
#endif

}

int main()
{
    NuTo::Structure* myStructure = NULL;

    std::string resultDir = "./ResultsIGA";

    if (boost::filesystem::exists(resultDir))
    {
        if (boost::filesystem::is_directory(resultDir))
        {
            boost::filesystem::remove_all(resultDir);
        }
    }

    // create result directory
    boost::filesystem::create_directory(resultDir);

    double solution = 0.;

    // 2D constant stress
    constantStress(solution, 1, resultDir);
}
