#include "BoostUnitTest.h"
#include "PlateWithHoleAnalytic.h"
#include <iostream>
#include <boost/filesystem.hpp>

#include "mechanics/structures/StructureBase.h"
#include "mechanics/structures/unstructured/Structure.h"

#include "math/SparseMatrixCSRGeneral.h"
#include "math/SparseMatrixCSRSymmetric.h"
#include "math/SparseDirectSolverMUMPS.h"
#include "math/SparseDirectSolverPardiso.h"

#include "mechanics/MechanicsEnums.h"
#include "visualize/VisualizeEnum.h"

#include "mechanics/iga/BSplineCurve.h"
#include "mechanics/iga/BSplineSurface.h"

#include "mechanics/sections/SectionPlane.h"

#include "mechanics/constraints/ConstraintCompanion.h"

#include "mechanics/groups/Group.h"


/*  ||>*----*----*----*----*
       |    |    |    |    | -->
       |    |    |    |    |
    ||>*----*----*----*----*     Sigma
       |    |    |    |    | -->
       |    |    |    |    |
    ||>*----*----*----*----*
       ^
*/
NuTo::BSplineSurface buildRect2D(double x0, double y0, double Height, double Length)
{
    /** Knots and control points **/
    int numElementsX = 1;
    int numElementsY = 1;

    Eigen::Vector2i degree(2, 2);

    int numKnotsX = 2 * (degree(0) + 1) + numElementsX - 1;
    int numKnotsY = 2 * (degree(1) + 1) + numElementsY - 1;

    Eigen::VectorXd knotsX(numKnotsX);
    Eigen::VectorXd knotsY(numKnotsY);

    for (int i = 0; i <= degree(0); i++)
        knotsX(i) = 0.;
    for (int i = degree(0) + 1; i <= degree(0) + numElementsX - 1; i++)
        knotsX(i) = knotsX(i - 1) + 1. / numElementsX;
    for (int i = degree(0) + numElementsX; i < numKnotsX; i++)
        knotsX(i) = 1;

    for (int i = 0; i <= degree(1); i++)
        knotsY(i) = 0.;
    for (int i = degree(1) + 1; i <= degree(1) + numElementsY - 1; i++)
        knotsY(i) = knotsY(i - 1) + 1. / numElementsY;
    for (int i = degree(1) + numElementsY; i < numKnotsY; i++)
        knotsY(i) = 1;

    int numControlPointsX = (numKnotsX - 1) - degree(0);
    int numControlPointsY = (numKnotsY - 1) - degree(1);

    Eigen::Matrix<Eigen::VectorXd, Eigen::Dynamic, Eigen::Dynamic> controlPoints(numControlPointsY, numControlPointsX);

    double incrx = (Length / (numControlPointsX - 1));
    double incry = (Height / (numControlPointsY - 1));
    for (int i = 0; i < numControlPointsY; i++)
    {
        for (int j = 0; j < numControlPointsX; j++)
        {
            controlPoints(i, j) = Eigen::Vector2d(x0 + incrx * j, y0 + incry * i);
        }
    }

    Eigen::MatrixXd weights(numControlPointsY, numControlPointsX);
    weights.setOnes(numControlPointsY, numControlPointsX);

    return NuTo::BSplineSurface(degree, knotsX, knotsY, controlPoints, weights);
}

void SolveAndVisualize(NuTo::Structure* s, std::string name)
{
    s->SolveGlobalSystemStaticElastic();

    int visualizationGroup = s->GroupGetElementsTotal();
    s->AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::DISPLACEMENTS);
    s->AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ENGINEERING_STRAIN);
    s->AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ENGINEERING_STRESS);

    std::string resultDir = "./ResultsIGA";
    boost::filesystem::create_directory(resultDir);
    s->ExportVtkDataFileElements(resultDir + "/Elements" + name + ".vtu");
    s->ExportVtkDataFileNodes(resultDir + "/Nodes" + name + ".vtu");
}

BOOST_AUTO_TEST_CASE(IGA_ConstantStress)
{
    /** parameters **/
    double YoungsModulus = 20000.;
    double PoissonRatio = 0.3;
    double Height = 5.;
    // there should be no dependency, because the pressure at the boundary is predefined
    double thickness = 2.123548;
    double Length = 10;
    double Stress = 10.;


    NuTo::BSplineSurface surface = buildRect2D(0, 0, Height, Length);

    /** Structure 2D **/
    NuTo::Structure* s = new NuTo::Structure(2);

    /** create section **/
    auto mySection = NuTo::SectionPlane::Create(thickness, false);

    /** create constitutive law **/
    int myMatLin = s->ConstitutiveLawCreate("LINEAR_ELASTIC_ENGINEERING_STRESS");
    s->ConstitutiveLawSetParameterDouble(myMatLin, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS,
                                         YoungsModulus);
    s->ConstitutiveLawSetParameterDouble(myMatLin, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO,
                                         PoissonRatio);

    std::set<NuTo::Node::eDof> setOfDOFS;
    setOfDOFS.insert(NuTo::Node::eDof::COORDINATES);
    setOfDOFS.insert(NuTo::Node::eDof::DISPLACEMENTS);

    int groupNodes = s->GroupCreate("Nodes");
    int groupElements = s->GroupCreate("Elements");

    surface.buildIGAStructure(*s, setOfDOFS, groupElements, groupNodes);

    s->Info();

    /** assign constitutive law **/
    s->ElementTotalSetConstitutiveLaw(myMatLin);
    s->ElementTotalSetSection(mySection);

    /** Boundary condition **/
    for (int i = 0; i < surface.GetNumControlPoints(1); i++)
    {
        const auto& node = *s->NodeGetNodePtr(i * surface.GetNumControlPoints(0));
        s->Constraints().Add(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Constraint::Component(node, {NuTo::eDirection::X}));
    }

    for (int i = 0; i < surface.GetNumControlPoints(0); i++)
    {
        const auto& node = *s->NodeGetNodePtr(i);
        s->Constraints().Add(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Constraint::Component(node, {NuTo::eDirection::Y}));
    }

    // right boundary
    int groupNumberNodesLeft = s->GroupCreate("NODES");
    for (int i = 1; i <= surface.GetNumControlPoints(1); i++)
    {
        s->GroupAddNode(groupNumberNodesLeft, i * surface.GetNumControlPoints(0) - 1);
    }

    int groupNumberElementsLeft = s->GroupCreate("ELEMENTS");
    for (int i = 1; i <= surface.GetNumIGAElements(1); i++)
    {
        s->GroupAddElement(groupNumberElementsLeft, i * surface.GetNumIGAElements(0) - 1);
    }

    s->LoadSurfacePressureCreate2D(groupNumberElementsLeft, groupNumberNodesLeft, -Stress);

    s->CalculateMaximumIndependentSets();
    s->NodeBuildGlobalDofs();

    s->NodeInfo(10);

    SolveAndVisualize(s, "Rectangle1");

    double displacementCorrect = (Stress * Length) / YoungsModulus;
    for (int id : s->GroupGetMemberIds(groupNumberNodesLeft))
    {
        auto displ = s->NodeGetNodePtr(id)->Get(NuTo::Node::eDof::DISPLACEMENTS);
        BOOST_CHECK_CLOSE(displ[0], displacementCorrect, 1.e-6);
    }
}


BOOST_AUTO_TEST_CASE(IGA_PlateWithHoleNeumann)
{
    NuTo::Structure* s = new NuTo::Structure(2);

    /*********/
    // Mesh  //
    /*********/

    int refine = 5;
    int noPtsX = 3;
    int noPtsY = 4;

    Eigen::Matrix<Eigen::VectorXd, Eigen::Dynamic, Eigen::Dynamic> controlPts(noPtsX, noPtsY);

    controlPts(0, 0) = Eigen::Vector2d(-1, 0);
    controlPts(0, 1) = Eigen::Vector2d(-1, 0.4142135623730951);
    controlPts(0, 2) = Eigen::Vector2d(-0.4142135623730951, 1);
    controlPts(0, 3) = Eigen::Vector2d(0, 1);

    controlPts(1, 0) = Eigen::Vector2d(-2.5, 0);
    controlPts(1, 1) = Eigen::Vector2d(-2.5, 0.75);
    controlPts(1, 2) = Eigen::Vector2d(-0.75, 2.5);
    controlPts(1, 3) = Eigen::Vector2d(0, 2.5);

    controlPts(2, 0) = Eigen::Vector2d(-4, 0);
    controlPts(2, 1) = Eigen::Vector2d(-4, 4);
    controlPts(2, 2) = Eigen::Vector2d(-4, 4);
    controlPts(2, 3) = Eigen::Vector2d(0, 4);

    Eigen::MatrixXd weights(noPtsX, noPtsY);
    weights.setOnes(noPtsX, noPtsY);
    weights(0, 1) = (1. + 1. / sqrt(2)) / 2.;
    weights(0, 2) = (1. + 1. / sqrt(2)) / 2.;


    Eigen::VectorXd knotsX(7);
    knotsX << 0, 0, 0, 0.5, 1, 1, 1;
    Eigen::VectorXd knotsY(6);
    knotsY << 0, 0, 0, 1, 1, 1;

    Eigen::Vector2i degree(2, 2);

    NuTo::BSplineSurface surface(degree, knotsX, knotsY, controlPts, weights);

    for (int i = 0; i < refine; i++)
    {
        surface.DuplicateKnots(0);
        surface.DuplicateKnots(1);
    }

    /**************/
    // Structure  //
    /**************/

    std::set<NuTo::Node::eDof> setOfDOFS;
    setOfDOFS.insert(NuTo::Node::eDof::COORDINATES);
    setOfDOFS.insert(NuTo::Node::eDof::DISPLACEMENTS);

    int groupNodes = s->GroupCreate("Nodes");
    int groupElements = s->GroupCreate("Elements");

    surface.buildIGAStructure(*s, setOfDOFS, groupElements, groupNodes);

    /** create section **/
    double thickness = 1.;
    auto mySection = NuTo::SectionPlane::Create(thickness, false);

    /** create constitutive law **/
    double YoungsModulus = 1.e5;
    double PoissonRatio = 0.3;
    int myMatLin = s->ConstitutiveLawCreate("LINEAR_ELASTIC_ENGINEERING_STRESS");
    s->ConstitutiveLawSetParameterDouble(myMatLin, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS,
                                         YoungsModulus);
    s->ConstitutiveLawSetParameterDouble(myMatLin, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO,
                                         PoissonRatio);

    s->ElementTotalSetConstitutiveLaw(myMatLin);
    s->ElementTotalSetSection(mySection);

    /**********************/
    // Boundary condition //
    /**********************/

    int groupElementsLeft = s->GroupCreate("ELEMENTS");
    int groupElementsUpper = s->GroupCreate("ELEMENTS");

    int start = (surface.GetNumIGAElements(1) - 1) * surface.GetNumIGAElements(0);
    for (int i = start; i < start + surface.GetNumIGAElements(0) / 2; i++)
    {
        s->GroupAddElement(groupElementsLeft, i);
    }

    start += surface.GetNumIGAElements(0) / 2;
    for (int i = start; i < surface.GetNumIGAElements(); i++)
    {
        s->GroupAddElement(groupElementsUpper, i);
    }

    auto& groupLeft = s->GroupGetNodesAtCoordinate(NuTo::eDirection::X, -4);
    auto& groupRight = s->GroupGetNodesAtCoordinate(NuTo::eDirection::X, 0);
    auto& groupBottom = s->GroupGetNodesAtCoordinate(NuTo::eDirection::Y, 0);
    auto& groupTop = s->GroupGetNodesAtCoordinate(NuTo::eDirection::Y, 4);

    s->Constraints().Add(NuTo::Node::eDof::DISPLACEMENTS,
                         NuTo::Constraint::Component(groupBottom, {NuTo::eDirection::Y}));
    s->Constraints().Add(NuTo::Node::eDof::DISPLACEMENTS,
                         NuTo::Constraint::Component(groupRight, {NuTo::eDirection::X}));

    s->LoadSurfacePressureFunctionCreate2D(groupElementsLeft, s->GroupGetId(&groupLeft),
                                           NuTo::Test::PlateWithHoleAnalytical::PressureLeft);
    s->LoadSurfacePressureFunctionCreate2D(groupElementsUpper, s->GroupGetId(&groupTop),
                                           NuTo::Test::PlateWithHoleAnalytical::PressureTop);

    BOOST_CHECK_NO_THROW(SolveAndVisualize(s, "Hole5"));
}
