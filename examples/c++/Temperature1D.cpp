#include <iostream>
#include "nuto/math/FullMatrix.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"
#include "nuto/math/SparseDirectSolverMUMPS.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"

int main()
{
	// geometry/mesh
	double area = 1.0;
	double length = 1.0;
	int num_elements = 10;

    // boundaries
    double boundary_temperature = 20.0;
    double boundary_flux = 10.0;

    // material
    double conductivity = 1.0;

	// create one-dimensional structure
	NuTo::Structure myStructure(1);

	// create section
	auto truss = myStructure.SectionCreate("Truss");
	myStructure.SectionSetArea(truss, area);

	// create material law
    auto material = myStructure.ConstitutiveLawCreate("Heat_Conduction");
    myStructure.ConstitutiveLawSetParameterDouble(material, 
        NuTo::Constitutive::eConstitutiveParameter::THERMAL_CONDUCTIVITY, conductivity);

	// create nodes
	NuTo::FullVector<double,Eigen::Dynamic> nodeCoordinates(1);
	for(int node = 0; node < num_elements + 1; node++)
	{
		nodeCoordinates(0) = node * length/num_elements;
        myStructure.NodeCreate(node, nodeCoordinates);
	}

    auto InterpolationType = myStructure.InterpolationTypeCreate("Truss1D");
    myStructure.InterpolationTypeAdd(InterpolationType, NuTo::Node::COORDINATES,
        NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    myStructure.InterpolationTypeAdd(InterpolationType, NuTo::Node::TEMPERATURE,
        NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    myStructure.InterpolationTypeSetIntegrationType(InterpolationType, "1D2NGauss2Ip", "noipdata");

	// create elements
	NuTo::FullVector<int,Eigen::Dynamic> elementIncidence(2);
	for(int element = 0; element < num_elements; element++)
	{
		elementIncidence(0) = element;
		elementIncidence(1) = element + 1;
        myStructure.ElementCreate(InterpolationType, elementIncidence);
		myStructure.ElementSetSection(element, truss);
		myStructure.ElementSetConstitutiveLaw(element, material);
	}

    myStructure.ElementTotalConvertToInterpolationType();

	// set boundary conditions and loads
	NuTo::FullVector<double,Eigen::Dynamic> direction(1);
	direction(0) = 1.0;
	myStructure.ConstraintLinearSetTemperatureNode(0, boundary_temperature);
	myStructure.SetNumLoadCases(1);
    myStructure.LoadCreateNodeHeatFlux(0, num_elements, direction, boundary_flux);

	// start analysis
    myStructure.SolveGlobalSystemStaticElastic(0);
    auto residual = myStructure.BuildGlobalInternalGradient() - myStructure.BuildGlobalExternalLoadVector(0);
    std::cout << "residual: " << residual.J.CalculateNormL2() << std::endl;

	// visualize results
	int visualizationGroup = myStructure.GroupCreate(NuTo::Groups::eGroupId::Elements);
    myStructure.GroupAddElementsTotal(visualizationGroup);

    myStructure.AddVisualizationComponent(visualizationGroup, NuTo::VisualizeBase::TEMPERATURE);
    myStructure.ExportVtkDataFileElements("Temperature1D.vtk");

	return 0;
}
