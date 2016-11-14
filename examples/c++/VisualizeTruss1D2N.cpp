#include "nuto/math/MathException.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/MechanicsEnums.h"
#include "nuto/visualize/VisualizeEnum.h"

int main()
{
    // create structure
    NuTo::Structure myStructure(1);

    // create nodes
    NuTo::FullVector<double, Eigen::Dynamic> Coordinates(1);
    NuTo::FullVector<double, Eigen::Dynamic> Displacements(1);

    Coordinates(0) = 1;
    myStructure.NodeCreate(0, Coordinates);
    Coordinates(0) = 6;
    myStructure.NodeCreate(1, Coordinates);
    Coordinates(0) = 10;
    myStructure.NodeCreate(2, Coordinates);

    int interpolationType = myStructure.InterpolationTypeCreate("TRUSS1D");
    myStructure.InterpolationTypeAdd(
            interpolationType, NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    myStructure.InterpolationTypeAdd(
            interpolationType, NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);

    // create element
    NuTo::FullVector<int, Eigen::Dynamic> Incidences(2);

    Incidences(0) = 0;
    Incidences(1) = 1;
    int myElement1 = myStructure.ElementCreate(interpolationType, Incidences);

    Incidences(0) = 1;
    Incidences(1) = 2;
    int myElement2 = myStructure.ElementCreate(interpolationType, Incidences);

    myStructure.ElementTotalConvertToInterpolationType(1e-6, 3);

    Displacements(0) = 0;
    myStructure.NodeSetDisplacements(0, Displacements);
    Displacements(0) = 0.3;
    myStructure.NodeSetDisplacements(1, Displacements);
    Displacements(0) = 0.6;
    myStructure.NodeSetDisplacements(2, Displacements);

    // create constitutive law
    int myMatLin = myStructure.ConstitutiveLawCreate("Linear_Elastic_Engineering_Stress");
    myStructure.ConstitutiveLawSetParameterDouble(
            myMatLin, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, 10);
    myStructure.ConstitutiveLawSetParameterDouble(
            myMatLin, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, 0.1);

    // create section
    int mySection1 = myStructure.SectionCreate("Truss");
    myStructure.SectionSetArea(mySection1, 0.01);

    // assign material, section and integration type
    myStructure.InterpolationTypeSetIntegrationType(
            interpolationType, NuTo::eIntegrationType::IntegrationType1D2NGauss2Ip);

    myStructure.ElementSetConstitutiveLaw(myElement1, myMatLin);
    myStructure.ElementSetSection(myElement1, mySection1);
    myStructure.ElementSetConstitutiveLaw(myElement2, myMatLin);
    myStructure.ElementSetSection(myElement2, mySection1);

    // visualize element
    int visualizationGroup = myStructure.GroupCreate(NuTo::eGroupId::Elements);
    myStructure.GroupAddElementsTotal(visualizationGroup);

    myStructure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::DISPLACEMENTS);
    myStructure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ENGINEERING_STRAIN);
    myStructure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ENGINEERING_STRESS);
    myStructure.ExportVtkDataFileElements("Truss1D2N.vtk");

    return 0;
}
