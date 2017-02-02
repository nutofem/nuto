#include "math/MathException.h"
#include "mechanics/MechanicsException.h"
#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/MechanicsEnums.h"
#include "visualize/VisualizeEnum.h"

int main()
{
    try
    {
        // create structure
        NuTo::Structure myStructure(2);
        myStructure.Info();

        // create nodes
        Eigen::MatrixXd coordinates(2, 9);

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

        std::cout << coordinates << std::endl;

        std::vector<int> nodeIds = myStructure.NodesCreate(coordinates);

        // create elements
        Eigen::MatrixXi elementIds(4, 4);

        // element1
        elementIds(0, 0) = nodeIds[0];
        elementIds(1, 0) = nodeIds[1];
        elementIds(2, 0) = nodeIds[4];
        elementIds(3, 0) = nodeIds[3];
        // element2
        elementIds(0, 1) = nodeIds[1];
        elementIds(1, 1) = nodeIds[2];
        elementIds(2, 1) = nodeIds[5];
        elementIds(3, 1) = nodeIds[4];
        // element3
        elementIds(0, 2) = nodeIds[3];
        elementIds(1, 2) = nodeIds[4];
        elementIds(2, 2) = nodeIds[7];
        elementIds(3, 2) = nodeIds[6];
        // element4
        elementIds(0, 3) = nodeIds[4];
        elementIds(1, 3) = nodeIds[5];
        elementIds(2, 3) = nodeIds[8];
        elementIds(3, 3) = nodeIds[7];

        std::cout << elementIds << std::endl;

        int interpolationType = myStructure.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::QUAD2D);
        myStructure.InterpolationTypeAdd(interpolationType, NuTo::Node::eDof::COORDINATES,
                NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
        myStructure.InterpolationTypeAdd(interpolationType, NuTo::Node::eDof::DISPLACEMENTS,
                NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
        myStructure.InterpolationTypeSetIntegrationType(interpolationType,
                NuTo::eIntegrationType::IntegrationType2D4NGauss4Ip);

        myStructure.ElementsCreate(interpolationType, elementIds);

        myStructure.ElementTotalConvertToInterpolationType(1.e-6, 10);

        // create constitutive law
        int myMatLin = myStructure.ConstitutiveLawCreate("LINEAR_ELASTIC_ENGINEERING_STRESS");
        myStructure.ConstitutiveLawSetParameterDouble(myMatLin,
                NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, 10);
        myStructure.ConstitutiveLawSetParameterDouble(myMatLin,
                NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, 0.1);

        // create section
        int mySection1 = myStructure.SectionCreate("PLANE_STRAIN");
        myStructure.SectionSetThickness(mySection1, 0.01);

        // assign material, section and integration type
        myStructure.ElementTotalSetConstitutiveLaw(myMatLin);
        myStructure.ElementTotalSetSection(mySection1);

        // visualize element
        int visualizationGroup = myStructure.GroupCreate(NuTo::eGroupId::Elements);
        myStructure.GroupAddElementsTotal(visualizationGroup);

        myStructure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::DISPLACEMENTS);
        myStructure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ENGINEERING_STRAIN);
        myStructure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ENGINEERING_STRESS);
        myStructure.ExportVtkDataFileElements("Plane2D4N.vtk");
    }
    catch (NuTo::MathException& e)
    {
        std::cout << e.ErrorMessage() << std::endl;
    }
    catch (NuTo::Exception& e)
    {
        std::cout << e.ErrorMessage() << std::endl;
    }
    catch (...)
    {
        std::cout << "Unexpected" << std::endl;
    }

    return 0;
}
