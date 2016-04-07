// $Id$

#include <iostream>
#include "nuto/math/FullMatrix.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"
#include "nuto/math/SparseDirectSolverMUMPS.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"


int main()
{
	// definitions
	double YoungsModulus = 20000.;
	double Area = 100. * 100.;
	double Length = 1000.;
	int NumElements = 10;
	double Force = 1.;
	bool EnableDisplacementControl = false;
	double BoundaryDisplacement = 0.1;

	// create one-dimensional structure
	NuTo::Structure myStructure(1);

	// create section
	int Section1 = myStructure.SectionCreate("Truss");
	myStructure.SectionSetArea(Section1, Area);

	// create material law
    int Material1 = myStructure.ConstitutiveLawCreate("Linear_Elastic_Engineering_Stress");
    myStructure.ConstitutiveLawSetParameterDouble(Material1,NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, YoungsModulus);

	// create nodes
	NuTo::FullVector<double,Eigen::Dynamic> nodeCoordinates(1);
	for(int node = 0; node < NumElements + 1; node++)
	{
		std::cout << "create node: " << node << " coordinates: " << node * Length/NumElements << std::endl;
		nodeCoordinates(0) = node * Length/NumElements;
        myStructure.NodeCreate(node, nodeCoordinates);
	}

    int InterpolationType = myStructure.InterpolationTypeCreate("Truss1D");
    myStructure.InterpolationTypeAdd(InterpolationType, NuTo::Node::COORDINATES,    NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    myStructure.InterpolationTypeAdd(InterpolationType, NuTo::Node::DISPLACEMENTS,  NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);

	// create elements
	NuTo::FullVector<int,Eigen::Dynamic> elementIncidence(2);
	for(int element = 0; element < NumElements; element++)
	{
		std::cout <<  "create element: " << element << " nodes: " << element << "," << element+1 << std::endl;
		elementIncidence(0) = element;
		elementIncidence(1) = element + 1;
        myStructure.ElementCreate(InterpolationType, elementIncidence);
		myStructure.ElementSetSection(element,Section1);
		myStructure.ElementSetConstitutiveLaw(element,Material1);
	}

    myStructure.ElementTotalConvertToInterpolationType();

	// set boundary conditions and loads
	NuTo::FullVector<double,Eigen::Dynamic> direction(1);
	direction(0) = 1;
	myStructure.ConstraintLinearSetDisplacementNode(0, direction, 0.0);
	myStructure.SetNumLoadCases(1);
	if(EnableDisplacementControl)
	{
		std::cout << "Displacement control" << std::endl;
		myStructure.ConstraintLinearSetDisplacementNode(NumElements, direction, BoundaryDisplacement);
	}
	else
	{
		std::cout << "Load control" << std::endl;
		myStructure.LoadCreateNodeForce(0,NumElements, direction, Force);
	}

	// start analysis
    myStructure.SolveGlobalSystemStaticElastic(1);
    auto residual = myStructure.BuildGlobalInternalGradient() - myStructure.BuildGlobalExternalLoadVector(1);

    std::cout << "residual: " << residual.J.CalculateNormL2() << std::endl;

	// visualize results
	int visualizationGroup = myStructure.GroupCreate(NuTo::Groups::eGroupId::Elements);
    myStructure.GroupAddElementsTotal(visualizationGroup);

    myStructure.AddVisualizationComponent(visualizationGroup, NuTo::VisualizeBase::DISPLACEMENTS);
    myStructure.AddVisualizationComponent(visualizationGroup, NuTo::VisualizeBase::ENGINEERING_STRAIN);
    myStructure.AddVisualizationComponent(visualizationGroup, NuTo::VisualizeBase::ENGINEERING_STRESS);
    myStructure.ExportVtkDataFileElements("Truss1D2N.vtk");

	return 0;
}
