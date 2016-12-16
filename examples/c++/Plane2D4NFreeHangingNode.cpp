#include "nuto/math/MathException.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/math/SparseDirectSolverMUMPS.h"
#include "nuto/mechanics/dofSubMatrixStorage/BlockFullMatrix.h"
#include "nuto/mechanics/dofSubMatrixStorage/BlockScalar.h"
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/structures/StructureOutputBlockVector.h"
#include "nuto/mechanics/MechanicsEnums.h"
#include "nuto/visualize/VisualizeEnum.h"

int main()
{
    // create structure
    NuTo::Structure myStructure(2);
    int myInterpolationType = myStructure.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::QUAD2D);
    myStructure.InterpolationTypeAdd(
            myInterpolationType, NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    myStructure.InterpolationTypeAdd(
            myInterpolationType, NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    myStructure.Info();

    // create nodes
    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> Coordinates(2, 11);
    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> Displacements(2, 11);

    Coordinates(0, 0) = 0;
    Coordinates(1, 0) = 0;
    Coordinates(0, 1) = 1;
    Coordinates(1, 1) = 0;
    Coordinates(0, 2) = 2;
    Coordinates(1, 2) = 0;
    Coordinates(0, 3) = 0;
    Coordinates(1, 3) = 1;
    Coordinates(0, 4) = 1;
    Coordinates(1, 4) = 1;
    Coordinates(0, 5) = 2;
    Coordinates(1, 5) = 1;
    Coordinates(0, 6) = 0;
    Coordinates(1, 6) = 2;
    Coordinates(0, 7) = 1;
    Coordinates(1, 7) = 2;
    Coordinates(0, 8) = 2;
    Coordinates(1, 8) = 2;

    Coordinates(0, 9) = 4;
    Coordinates(1, 9) = 0;
    Coordinates(0, 10) = 4;
    Coordinates(1, 10) = 2;
    Coordinates.Info();

    NuTo::FullVector<int, Eigen::Dynamic> Nodes = myStructure.NodesCreate(Coordinates);

    // create elements
    std::vector<int> Incidences(4);
    NuTo::FullVector<int, Eigen::Dynamic> Elements(5);

    // element1
    Incidences[0] = Nodes(0);
    Incidences[1] = Nodes(1);
    Incidences[2] = Nodes(4);
    Incidences[3] = Nodes(3);
    Elements[0] = myStructure.ElementCreate(myInterpolationType, Incidences);

    // element2
    Incidences[0] = Nodes(1);
    Incidences[1] = Nodes(2);
    Incidences[2] = Nodes(5);
    Incidences[3] = Nodes(4);
    Elements[1] = myStructure.ElementCreate(myInterpolationType, Incidences);

    // element3
    Incidences[0] = Nodes(3);
    Incidences[1] = Nodes(4);
    Incidences[2] = Nodes(7);
    Incidences[3] = Nodes(6);
    Elements[2] = myStructure.ElementCreate(myInterpolationType, Incidences);

    // element4
    Incidences[0] = Nodes(4);
    Incidences[1] = Nodes(5);
    Incidences[2] = Nodes(8);
    Incidences[3] = Nodes(7);
    Elements[3] = myStructure.ElementCreate(myInterpolationType, Incidences);

    // element5
    Incidences[0] = Nodes(2, 0);
    Incidences[1] = Nodes(9, 0);
    Incidences[2] = Nodes(10, 0);
    Incidences[3] = Nodes(8, 0);
    Elements[4] = myStructure.ElementCreate(myInterpolationType, Incidences);

    Coordinates.Info();

    // create constitutive law
    int myMatLin =  myStructure.ConstitutiveLawCreate("LINEAR_ELASTIC_ENGINEERING_STRESS");
    myStructure.ConstitutiveLawSetParameterDouble(
            myMatLin, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, 10);
    myStructure.ConstitutiveLawSetParameterDouble(
            myMatLin, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, 0.1);

    // create section
    int mySection1 = myStructure.SectionCreate("PLANE_STRAIN");
    myStructure.SectionSetThickness(mySection1, 0.01);

    // assign material, section and integration type
    myStructure.ElementTotalConvertToInterpolationType(1.e-6, 3);
    myStructure.ElementTotalSetConstitutiveLaw(myMatLin);
    myStructure.ElementTotalSetSection(mySection1);

    // build global stiffness matrix and equivalent load vector which correspond to prescribed boundary values
    auto hessianElement0 = myStructure.ElementBuildHessian0(0);
    auto hessianElement4 = myStructure.ElementBuildHessian0(4);

    hessianElement0.Export().WriteToFile("stiffness", "   ");
    hessianElement4.Export().WriteToFile("stiffnessCoarse", "   ");

    // boundary conditions
    NuTo::FullVector<double, Eigen::Dynamic> direction(2);
    direction(0) = 1;
    direction(1) = 0;
    myStructure.ConstraintLinearSetDisplacementNode(0, direction, 0.0);
    myStructure.ConstraintLinearSetDisplacementNode(3, direction, 0.0);
    myStructure.ConstraintLinearSetDisplacementNode(6, direction, 0.0);
    direction(0) = 0;
    direction(1) = 1;
    for (int i = 0; i < 11; ++i)
    {
        myStructure.ConstraintLinearSetDisplacementNode(i, direction, 0.0);
    }

    std::cout << "Displacement control" << std::endl;
    // boundary displacments
    double BoundaryDisplacement = 1;
    direction(0) = 1;
    direction(1) = 0;
    myStructure.ConstraintLinearSetDisplacementNode(9, direction, BoundaryDisplacement);
    myStructure.ConstraintLinearSetDisplacementNode(10, direction, BoundaryDisplacement);

    // start analysis
    myStructure.SolveGlobalSystemStaticElastic(1);
    auto residual = myStructure.BuildGlobalInternalGradient() - myStructure.BuildGlobalExternalLoadVector(1);

    std::cout << "residual: " << residual.J.CalculateNormL2() << std::endl;

    // visualize element
    int visualizationGroup = myStructure.GroupCreate(NuTo::eGroupId::Elements);
    myStructure.GroupAddElementsTotal(visualizationGroup);

    myStructure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::DISPLACEMENTS);
    myStructure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ENGINEERING_STRAIN);
    myStructure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ENGINEERING_STRESS);
    myStructure.ExportVtkDataFileElements("Plane2D4N.vtk");
    return 0;
}
