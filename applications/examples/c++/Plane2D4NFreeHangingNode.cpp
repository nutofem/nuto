#include "base/serializeStream/SerializeStreamOut.h"
#include "math/MathException.h"
#include "math/SparseDirectSolverMUMPS.h"
#include "mechanics/dofSubMatrixStorage/BlockFullMatrix.h"
#include "mechanics/dofSubMatrixStorage/BlockScalar.h"
#include "mechanics/MechanicsException.h"
#include "mechanics/sections/SectionPlane.h"
#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/structures/StructureOutputBlockVector.h"
#include "mechanics/MechanicsEnums.h"
#include "visualize/VisualizeEnum.h"

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
    // create nodes
    Eigen::MatrixXd coordinates(2, 11);

    coordinates(0, 0) = 0;
    coordinates(1, 0) = 0;
    coordinates(0, 1) = 1;
    coordinates(1, 1) = 0;
    coordinates(0, 2) = 2;
    coordinates(1, 2) = 0;
    coordinates(0, 3) = 0;
    coordinates(1, 3) = 1;
    coordinates(0, 4) = 1;
    coordinates(1, 4) = 1;
    coordinates(0, 5) = 2;
    coordinates(1, 5) = 1;
    coordinates(0, 6) = 0;
    coordinates(1, 6) = 2;
    coordinates(0, 7) = 1;
    coordinates(1, 7) = 2;
    coordinates(0, 8) = 2;
    coordinates(1, 8) = 2;

    coordinates(0, 9) = 4;
    coordinates(1, 9) = 0;
    coordinates(0, 10) = 4;
    coordinates(1, 10) = 2;

    std::cout << coordinates << std::endl;

    std::vector<int> nodeIds = myStructure.NodesCreate(coordinates);

    std::vector<int> Incidences(4);

    // element1
    Incidences[0] = nodeIds[0];
    Incidences[1] = nodeIds[1];
    Incidences[2] = nodeIds[4];
    Incidences[3] = nodeIds[3];
    myStructure.ElementCreate(myInterpolationType, Incidences);

    // element2
    Incidences[0] = nodeIds[1];
    Incidences[1] = nodeIds[2];
    Incidences[2] = nodeIds[5];
    Incidences[3] = nodeIds[4];
    myStructure.ElementCreate(myInterpolationType, Incidences);

    // element3
    Incidences[0] = nodeIds[3];
    Incidences[1] = nodeIds[4];
    Incidences[2] = nodeIds[7];
    Incidences[3] = nodeIds[6];
    myStructure.ElementCreate(myInterpolationType, Incidences);

    // element4
    Incidences[0] = nodeIds[4];
    Incidences[1] = nodeIds[5];
    Incidences[2] = nodeIds[8];
    Incidences[3] = nodeIds[7];
    myStructure.ElementCreate(myInterpolationType, Incidences);

    // element5
    Incidences[0] = nodeIds[2 ];
    Incidences[1] = nodeIds[9 ];
    Incidences[2] = nodeIds[10];
    Incidences[3] = nodeIds[8 ];
    myStructure.ElementCreate(myInterpolationType, Incidences);

    // create constitutive law
    int myMatLin =  myStructure.ConstitutiveLawCreate("LINEAR_ELASTIC_ENGINEERING_STRESS");
    myStructure.ConstitutiveLawSetParameterDouble(
            myMatLin, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, 10);
    myStructure.ConstitutiveLawSetParameterDouble(
            myMatLin, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, 0.1);

    // create section
    auto mySection1 = NuTo::SectionPlane::Create(0.01, true);

    // assign material, section and integration type
    myStructure.ElementTotalConvertToInterpolationType();
    myStructure.ElementTotalSetConstitutiveLaw(myMatLin);
    myStructure.ElementTotalSetSection(mySection1);

    // build global stiffness matrix and equivalent load vector which correspond to prescribed boundary values
    auto hessianElement0 = myStructure.ElementBuildHessian0(0);
    auto hessianElement3 = myStructure.ElementBuildHessian0(4);

    NuTo::SerializeStreamOut stiffnessFile("stiffness", false);
    stiffnessFile.SaveMatrix(hessianElement0.Export());
    NuTo::SerializeStreamOut stiffnessCoarseFile("stiffnessCoarse", false);
    stiffnessCoarseFile.SaveMatrix(hessianElement3.Export());

    // boundary conditions
    Eigen::VectorXd direction(2);
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
