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
	myStructure.SectionCreate("Section1","Truss");
	myStructure.SectionSetArea("Section1", Area);

	// create material law
	int Material1 = myStructure.ConstitutiveLawCreate("LinearElastic");
	myStructure.ConstitutiveLawSetYoungsModulus(Material1, YoungsModulus);

	// create nodes
	NuTo::FullMatrix<double> nodeCoordinates(1,1);
	for(int node = 0; node < NumElements + 1; node++)
	{
		std::cout << "create node: " << node << " coordinates: " << node * Length/NumElements << std::endl;
		nodeCoordinates(0, 0) = node * Length/NumElements;
		myStructure.NodeCreate(node, "displacements", nodeCoordinates);
	}

	// create elements
	NuTo::FullMatrix<int> elementIncidence(2,1);
	for(int element = 0; element < NumElements; element++)
	{
		std::cout <<  "create element: " << element << " nodes: " << element << "," << element+1 << std::endl;
		elementIncidence(0, 0) = element;
		elementIncidence(1, 0) = element + 1;
		myStructure.ElementCreate(element, "Truss1D2N", elementIncidence);
		myStructure.ElementSetSection(element,"Section1");
		myStructure.ElementSetConstitutiveLaw(element,Material1);
	}

	// set boundary conditions and loads
	NuTo::FullMatrix<double> direction(1,1);
	direction(0,0) = 1;
	myStructure.ConstraintSetDisplacementNode(0, direction, 0.0);
	if(EnableDisplacementControl)
	{
		std::cout << "Displacement control" << std::endl;
		myStructure.ConstraintSetDisplacementNode(NumElements, direction, BoundaryDisplacement);
	}
	else
	{
		std::cout << "Load control" << std::endl;
		myStructure.LoadCreateNodeForce(NumElements, direction, Force);
	}

	// start analysis
	// build global dof numbering
	myStructure.NodeBuildGlobalDofs();

	// build global stiffness matrix and equivalent load vector which correspond to prescribed boundary values
	NuTo::SparseMatrixCSRGeneral<double> stiffnessMatrix;
	NuTo::FullMatrix<double> dispForceVector;
	myStructure.BuildGlobalCoefficientMatrix0(stiffnessMatrix, dispForceVector);

	// build global external load vector
	NuTo::FullMatrix<double> extForceVector;
	myStructure.BuildGlobalExternalLoadVector(extForceVector);

	// calculate right hand side
	NuTo::FullMatrix<double> rhsVector = dispForceVector + extForceVector;

	// solve
	NuTo::SparseDirectSolverMUMPS mySolver;
	NuTo::FullMatrix<double> displacementVector;
	stiffnessMatrix.SetOneBasedIndexing();
#ifdef HAVE_MUMPS
	mySolver.Solve(stiffnessMatrix, rhsVector, displacementVector);

	// write displacements to node
	myStructure.NodeMergeActiveDofValues(displacementVector);

	// calculate residual
	NuTo::FullMatrix<double> intForceVector;
	myStructure.BuildGlobalGradientInternalPotentialVector(intForceVector);
	NuTo::FullMatrix<double> residualVector = extForceVector - intForceVector;
	std::cout << "residual: " << residualVector.Norm() << std::endl;

	// visualize results
	myStructure.ExportVtkDataFile("Truss1D2N.vtk","displacements engineering_strain engineering_stress");

#else
    std::cout << "MUMPS not available - can't solve system of equations " << std::endl;
#endif // HAVE_MUMPS
	return 0;
}
