#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/MechanicsException.h"

#define PRINTRESULT false

void CheckTriangle(NuTo::IntegrationType::eIntegrationType rIntegrationType)
{
    //create structure
    NuTo::Structure myStructure(2);
    //create nodes
    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> nodeCoordinates(2, 4);
    nodeCoordinates <<
            0, 10,  10,  0,
            0,  0,  10, 10;
    myStructure.NodesCreate(nodeCoordinates);

    int interpolationType = myStructure.InterpolationTypeCreate(NuTo::Interpolation::TRIANGLE2D);
    myStructure.InterpolationTypeAdd(interpolationType, NuTo::Node::COORDINATES, NuTo::Interpolation::EQUIDISTANT1);
    myStructure.InterpolationTypeAdd(interpolationType, NuTo::Node::DISPLACEMENTS, NuTo::Interpolation::EQUIDISTANT4);
    myStructure.InterpolationTypeSetIntegrationType(interpolationType, rIntegrationType, NuTo::IpData::NOIPDATA);


    //create element
    NuTo::FullMatrix<int, Eigen::Dynamic, Eigen::Dynamic> nodeNumbers(3, 2);
    nodeNumbers << 0, 2,
                   1, 3,
                   2, 0;

    myStructure.ElementsCreate(interpolationType, nodeNumbers);
    myStructure.ElementTotalConvertToInterpolationType();


    //Calculate maximum independent sets for parallelization (openmp)
    myStructure.CalculateMaximumIndependentSets();

    //create constitutive law
    int myMatLin = myStructure.ConstitutiveLawCreate("LinearElasticEngineeringStress");
    myStructure.ConstitutiveLawSetParameterDouble(myMatLin, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, 10);
    myStructure.ConstitutiveLawSetParameterDouble(myMatLin, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, 0.2);

    //create section
    int mySection = myStructure.SectionCreate("Plane_Strain");
    myStructure.SectionSetThickness(mySection, 1);

    //assign constitutive law
    myStructure.ElementTotalSetConstitutiveLaw(myMatLin);
    myStructure.ElementTotalSetSection(mySection);

    NuTo::ElementBase* element = myStructure.ElementGetElementPtr(0);
    for (int i = 0; i < element->GetNumNodes(NuTo::Node::DISPLACEMENTS); ++i)
    {
        NuTo::NodeBase* node = element->GetNode(i, NuTo::Node::DISPLACEMENTS);
        Eigen::Vector2d coord = node->GetCoordinates2D();
        Eigen::Vector2d displ(coord.x()*coord.x(), 0);
        node->SetDisplacements2D(displ);
    }

    element = myStructure.ElementGetElementPtr(1);
    for (int i = 0; i < element->GetNumNodes(NuTo::Node::DISPLACEMENTS); ++i)
    {
        NuTo::NodeBase* node = element->GetNode(i, NuTo::Node::DISPLACEMENTS);
        Eigen::Vector2d coord = node->GetCoordinates2D();
        Eigen::Vector2d displ(coord.x()*coord.x(), 0);
        node->SetDisplacements2D(displ);

    }

    int visualizationGroup = myStructure.GroupCreate(NuTo::Groups::eGroupId::Elements);
    myStructure.GroupAddElementsTotal(visualizationGroup);

    myStructure.AddVisualizationComponent(visualizationGroup, NuTo::VisualizeBase::DISPLACEMENTS);
    myStructure.AddVisualizationComponent(visualizationGroup, NuTo::VisualizeBase::ENGINEERING_STRAIN);

    myStructure.ExportVtkDataFileElements("./voronoi.vtk");

}

int main()
{
    try
    {
        CheckTriangle(NuTo::IntegrationType::IntegrationType2D3NGauss12IpDetail);
    }
    catch (NuTo::Exception& e)
    {
        std::cout << e.ErrorMessage() << std::endl;
        return -1;
    }

    return 0;
}
