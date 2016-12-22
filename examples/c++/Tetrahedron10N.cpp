#include "math/FullMatrix.h"
#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/dofSubMatrixStorage/BlockFullMatrix.h"
#include "mechanics/MechanicsEnums.h"
#include "visualize/VisualizeEnum.h"

int main()
{
    // create structure
    NuTo::Structure myStructure(3);
    int myInterpolationType = myStructure.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::TETRAHEDRON3D);
    myStructure.InterpolationTypeAdd(
            myInterpolationType, NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);
    myStructure.InterpolationTypeAdd(
            myInterpolationType, NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);

    // create nodes
    NuTo::FullVector<double, Eigen::Dynamic> Coordinates(3);
    std::vector<int> Incidence(10);

    // create nodes
    Coordinates(0) = 0.0;
    Coordinates(1) = 0.0;
    Coordinates(2) = 0.0;
    Incidence[0] = myStructure.NodeCreate(Coordinates);

    Coordinates(0) = 1.0;
    Coordinates(1) = 0.0;
    Coordinates(2) = 0.0;
    Incidence[1] = myStructure.NodeCreate(Coordinates);

    Coordinates(0) = 0.0;
    Coordinates(1) = 1.0;
    Coordinates(2) = 0.0;
    Incidence[2] = myStructure.NodeCreate(Coordinates);

    Coordinates(0) = 0.0;
    Coordinates(1) = 0.0;
    Coordinates(2) = 1.0;
    Incidence[3] = myStructure.NodeCreate(Coordinates);

    Coordinates(0) = 0.5;
    Coordinates(1) = 0.0;
    Coordinates(2) = 0.0;
    Incidence[4] = myStructure.NodeCreate(Coordinates);

    Coordinates(0) = 0.5;
    Coordinates(1) = 0.5;
    Coordinates(2) = 0.0;
    Incidence[5] = myStructure.NodeCreate(Coordinates);

    Coordinates(0) = 0.0;
    Coordinates(1) = 0.5;
    Coordinates(2) = 0.0;
    Incidence[6] = myStructure.NodeCreate(Coordinates);

    Coordinates(0) = 0.0;
    Coordinates(1) = 0.0;
    Coordinates(2) = 0.5;
    Incidence[7] = myStructure.NodeCreate(Coordinates);

    Coordinates(0) = 0.0;
    Coordinates(1) = 0.5;
    Coordinates(2) = 0.5;
    Incidence[8] = myStructure.NodeCreate(Coordinates);

    Coordinates(0) = 0.5;
    Coordinates(1) = 0.0;
    Coordinates(2) = 0.5;
    Incidence[9] = myStructure.NodeCreate(Coordinates);

    // create element
    int myElement1 = myStructure.ElementCreate(myInterpolationType, Incidence);
    myStructure.ElementTotalConvertToInterpolationType(1e-6, 10);

    // create constitutive law
    int myMatLin = myStructure.ConstitutiveLawCreate("Linear_Elastic_Engineering_Stress");
    myStructure.ConstitutiveLawSetParameterDouble(
            myMatLin, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, 1);
    myStructure.ConstitutiveLawSetParameterDouble(
            myMatLin, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, 0);

    // create section
    int mySection = myStructure.SectionCreate("Volume");

    // assign constitutive law
    myStructure.ElementSetConstitutiveLaw(myElement1, myMatLin);
    // assign section
    myStructure.ElementSetSection(myElement1, mySection);

    auto Ke = myStructure.ElementBuildHessian0(myElement1);

    // visualize results
    int visualizationGroup = myStructure.GroupCreate(NuTo::eGroupId::Elements);
    myStructure.GroupAddElementsTotal(visualizationGroup);

    myStructure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::DISPLACEMENTS);
    myStructure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ENGINEERING_STRAIN);
    myStructure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ENGINEERING_STRESS);
    myStructure.ExportVtkDataFileElements("TetrahedronElements10N.vtu", true);
    myStructure.ExportVtkDataFileNodes("TetrahedronNodes10N.vtu", true);

    return 0;
}
