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
	int Material1 = myStructure.ConstitutiveLawCreate("LinearElastic");
	myStructure.ConstitutiveLawSetYoungsModulus(Material1, YoungsModulus);

	// create nodes
	NuTo::FullVector<double,Eigen::Dynamic> nodeCoordinates(1);
	for(int node = 0; node < NumElements + 1; node++)
	{
		std::cout << "create node: " << node << " coordinates: " << node * Length/NumElements << std::endl;
		nodeCoordinates(0) = node * Length/NumElements;
		myStructure.NodeCreate(node, "displacements", nodeCoordinates);
	}

	// create elements
	NuTo::FullVector<int,Eigen::Dynamic> elementIncidence(2);
	for(int element = 0; element < NumElements; element++)
	{
		std::cout <<  "create element: " << element << " nodes: " << element << "," << element+1 << std::endl;
		elementIncidence(0) = element;
		elementIncidence(1) = element + 1;
		myStructure.ElementCreate(element, "Truss1D2N", elementIncidence);
		myStructure.ElementSetSection(element,Section1);
		myStructure.ElementSetConstitutiveLaw(element,Material1);
	}

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
	// build global dof numbering
	myStructure.NodeBuildGlobalDofs();

	// build global stiffness matrix and equivalent load vector which correspond to prescribed boundary values
	NuTo::SparseMatrixCSRGeneral<double> stiffnessMatrix;
	NuTo::FullVector<double,Eigen::Dynamic> dispForceVector;
	myStructure.BuildGlobalCoefficientMatrix0(stiffnessMatrix, dispForceVector);

	// build global external load vector
	NuTo::FullVector<double,Eigen::Dynamic> extForceVector;
	myStructure.BuildGlobalExternalLoadVector(0,extForceVector);

	// calculate right hand side
	NuTo::FullVector<double,Eigen::Dynamic> rhsVector = dispForceVector + extForceVector;

	// solve
	NuTo::SparseDirectSolverMUMPS mySolver;
	NuTo::FullVector<double,Eigen::Dynamic> displacementVector;
	stiffnessMatrix.SetOneBasedIndexing();
#ifdef HAVE_MUMPS
	mySolver.Solve(stiffnessMatrix, rhsVector, displacementVector);

	// write displacements to node
	myStructure.NodeMergeActiveDofValues(displacementVector);

	// calculate residual
	NuTo::FullVector<double,Eigen::Dynamic> intForceVector;
	myStructure.BuildGlobalGradientInternalPotentialVector(intForceVector);
	NuTo::FullVector<double,Eigen::Dynamic> residualVector = extForceVector - intForceVector;
	std::cout << "residual: " << residualVector.Norm() << std::endl;

#ifdef ENABLE_VISUALIZE
	// visualize results
    myStructure.AddVisualizationComponentDisplacements();
    myStructure.AddVisualizationComponentEngineeringStrain();
    myStructure.AddVisualizationComponentEngineeringStress();
	myStructure.ExportVtkDataFile("Truss1D2N.vtk");
#endif

#else
    std::cout << "MUMPS not available - can't solve system of equations " << std::endl;
#endif // HAVE_MUMPS
	return 0;
}
