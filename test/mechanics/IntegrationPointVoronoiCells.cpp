#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/constitutive/ConstitutiveEnum.h"
#include "nuto/mechanics/elements/ElementBase.h"
#include "nuto/mechanics/elements/IpDataEnum.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeEnum.h"
#include "nuto/mechanics/interpolationtypes/InterpolationTypeEnum.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/nodes/NodeEnum.h"
#include "nuto/mechanics/groups/GroupEnum.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/MechanicsException.h"


#include "nuto/visualize/VisualizeEnum.h"


#define PRINTRESULT false

void CheckTriangle(NuTo::eIntegrationType rIntegrationType)
{
    //create structure
    NuTo::Structure myStructure(2);
    //create nodes
    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> nodeCoordinates(2, 4);
    nodeCoordinates <<
            0, 10,  10,  0,
            0,  0,  10, 10;
    myStructure.NodesCreate(nodeCoordinates);

    int interpolationType = myStructure.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::TRIANGLE2D);
    myStructure.InterpolationTypeAdd(interpolationType, NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    myStructure.InterpolationTypeAdd(interpolationType, NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT4);
    myStructure.InterpolationTypeSetIntegrationType(interpolationType, rIntegrationType, NuTo::IpData::eIpDataType::NOIPDATA);


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
    int myMatLin = myStructure.ConstitutiveLawCreate("Linear_Elastic_Engineering_Stress");
    myStructure.ConstitutiveLawSetParameterDouble(myMatLin, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, 10);
    myStructure.ConstitutiveLawSetParameterDouble(myMatLin, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, 0.2);

    //create section
    int mySection = myStructure.SectionCreate("Plane_Strain");
    myStructure.SectionSetThickness(mySection, 1);

    //assign constitutive law
    myStructure.ElementTotalSetConstitutiveLaw(myMatLin);
    myStructure.ElementTotalSetSection(mySection);

    NuTo::ElementBase* element = myStructure.ElementGetElementPtr(0);
    for (int i = 0; i < element->GetNumNodes(NuTo::Node::eDof::DISPLACEMENTS); ++i)
    {
        NuTo::NodeBase* node = element->GetNode(i, NuTo::Node::eDof::DISPLACEMENTS);
        Eigen::Vector2d coord = node->Get(NuTo::Node::eDof::COORDINATES);
        Eigen::Vector2d displ(coord.x()*coord.x(), 0);
        node->Set(NuTo::Node::eDof::DISPLACEMENTS, displ);
    }

    element = myStructure.ElementGetElementPtr(1);
    for (int i = 0; i < element->GetNumNodes(NuTo::Node::eDof::DISPLACEMENTS); ++i)
    {
        NuTo::NodeBase* node = element->GetNode(i, NuTo::Node::eDof::DISPLACEMENTS);
        Eigen::Vector2d coord = node->Get(NuTo::Node::eDof::COORDINATES);
        Eigen::Vector2d displ(coord.x()*coord.x(), 0);
        node->Set(NuTo::Node::eDof::DISPLACEMENTS, displ);

    }

    int visualizationGroup = myStructure.GroupCreate(NuTo::eGroupId::Elements);
    myStructure.GroupAddElementsTotal(visualizationGroup);

    myStructure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::DISPLACEMENTS);
    myStructure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ENGINEERING_STRAIN);

    myStructure.ExportVtkDataFileElements("./voronoi.vtk");

}

int main()
{
    try
    {
        CheckTriangle(NuTo::eIntegrationType::IntegrationType2D3NGauss12IpDetail);
    }
    catch (NuTo::Exception& e)
    {
        std::cout << e.ErrorMessage() << std::endl;
        return -1;
    }

    return 0;
}
