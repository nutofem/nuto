// $Id: VisualizeTruss1D2N.cpp 147 2009-12-02 17:12:54Z eckardt4 $

#include "nuto/math/MathException.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"

#include "nuto/base/Debug.h"

int main()
{
    try
    {
        // create structure
        NuTo::Structure myStructure(2);
        myStructure.Info();

        // create nodes
        NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> Coordinates(2, 9);
        NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> Displacements(2, 9);

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

        DBG_POSITION_INFO("Coordinates Matrix")
        Coordinates.Info();

        NuTo::FullVector<int, Eigen::Dynamic> Nodes = myStructure.NodesCreate(Coordinates);

        // create elements
        NuTo::FullMatrix<int, Eigen::Dynamic, Eigen::Dynamic> Incidences(4, 4);

        // element1
        Incidences(0, 0) = Nodes(0, 0);
        Incidences(1, 0) = Nodes(1, 0);
        Incidences(2, 0) = Nodes(4, 0);
        Incidences(3, 0) = Nodes(3, 0);
        // element2
        Incidences(0, 1) = Nodes(1, 0);
        Incidences(1, 1) = Nodes(2, 0);
        Incidences(2, 1) = Nodes(5, 0);
        Incidences(3, 1) = Nodes(4, 0);
        // element3
        Incidences(0, 2) = Nodes(3, 0);
        Incidences(1, 2) = Nodes(4, 0);
        Incidences(2, 2) = Nodes(7, 0);
        Incidences(3, 2) = Nodes(6, 0);
        // element4
        Incidences(0, 3) = Nodes(4, 0);
        Incidences(1, 3) = Nodes(5, 0);
        Incidences(2, 3) = Nodes(8, 0);
        Incidences(3, 3) = Nodes(7, 0);

        DBG_POSITION_INFO("Incidence Matrix")
        Coordinates.Info();

        int interpolationType = myStructure.InterpolationTypeCreate(NuTo::Interpolation::QUAD2D);
        myStructure.InterpolationTypeAdd(interpolationType, NuTo::Node::COORDINATES, NuTo::Interpolation::EQUIDISTANT1);
        myStructure.InterpolationTypeAdd(interpolationType, NuTo::Node::DISPLACEMENTS, NuTo::Interpolation::EQUIDISTANT1);
        myStructure.InterpolationTypeSetIntegrationType(interpolationType, NuTo::IntegrationType::IntegrationType2D4NGauss4Ip, NuTo::IpData::NOIPDATA);

        NuTo::FullVector<int, Eigen::Dynamic> Elements = myStructure.ElementsCreate(interpolationType, Incidences);

        myStructure.ElementTotalConvertToInterpolationType(1.e-6, 10);

        // create constitutive law
        int myMatLin = myStructure.ConstitutiveLawCreate("LINEARELASTICENGINEERINGSTRESS");
        myStructure.ConstitutiveLawSetParameterDouble(myMatLin,NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, 10);
        myStructure.ConstitutiveLawSetParameterDouble(myMatLin,NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, 0.1);

        // create section
        int mySection1 = myStructure.SectionCreate("PLANE_STRAIN");
        myStructure.SectionSetThickness(mySection1, 0.01);

        // assign material, section and integration type
        myStructure.ElementTotalSetConstitutiveLaw(myMatLin);
        myStructure.ElementTotalSetSection(mySection1);

#ifdef ENABLE_VISUALIZE
        // visualize element
        int visualizationGroup = myStructure.GroupCreate(NuTo::Groups::eGroupId::Elements);
        myStructure.GroupAddElementsTotal(visualizationGroup);

        myStructure.AddVisualizationComponent(visualizationGroup, NuTo::VisualizeBase::DISPLACEMENTS);
        myStructure.AddVisualizationComponent(visualizationGroup, NuTo::VisualizeBase::ENGINEERING_STRAIN);
        myStructure.AddVisualizationComponent(visualizationGroup, NuTo::VisualizeBase::ENGINEERING_STRESS);
        myStructure.ExportVtkDataFileElements("Plane2D4N.vtk");
#endif
    } catch (NuTo::MathException& e)
    {
        std::cout << e.ErrorMessage() << std::endl;
    } catch (NuTo::Exception& e)
    {
        std::cout << e.ErrorMessage() << std::endl;
    } catch (...)
    {
        std::cout << "Unexpected" << std::endl;
    }

    return 0;
}
