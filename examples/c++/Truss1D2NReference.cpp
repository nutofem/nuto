// $Id$


//Author: Andrea keßler
//changed form Truss1D2N
#include <iostream>
#include "nuto/math/FullMatrix.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"
#include "nuto/math/SparseDirectSolverMUMPS.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"

int main()
{
	// definitions
	double YoungsModulus = 10000.;
	double YoungsModulusDiff = 20000.; //difference from node 0 to lats node, +/-
	double Area = 10. * 10.;
    //double Length = 6000.;
	int NumElements = 4;
	double Force = 1.;
	bool EnableDisplacementControl = true;
	double BoundaryDisplacement = 1;
	std::ofstream output;

	// create one-dimensional structure
	NuTo::Structure myStructure(1);

	// create section
	int Section1 = myStructure.SectionCreate("Truss");
	myStructure.SectionSetArea(Section1, Area);

	// create material law
	int Material1 = myStructure.ConstitutiveLawCreate("LINEARELASTICENGINEERINGSTRESS");
	int Material2 = myStructure.ConstitutiveLawCreate("LINEARELASTICENGINEERINGSTRESS");
    myStructure.ConstitutiveLawSetParameterDouble(Material1,NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, YoungsModulus);
    myStructure.ConstitutiveLawSetParameterDouble(Material2,NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, YoungsModulus+YoungsModulusDiff);
//    myStructure.ConstitutiveLawSetPoissonsRatio(Material1,0.2);
//    myStructure.ConstitutiveLawSetPoissonsRatio(Material2,0.2);

	std::string filename="Truss1D2NReference-nodes";
	output.open(filename.c_str());
	if(output)
	{
		output<<"#"<<filename<<"\n";
		output<<"nodes\n";
	}
	else
		std::cout<<"Truss1DMultiphase] error output file.\n";
	// create nodes
	NuTo::FullVector<double,Eigen::Dynamic> nodeCoordinates(1);
	nodeCoordinates(0) =  0.;
	output<<nodeCoordinates(0)<<"\n";
    myStructure.NodeCreate(0, nodeCoordinates);
	nodeCoordinates(0) =  2000.;
	output<<nodeCoordinates(0)<<"\n";
    myStructure.NodeCreate(1, nodeCoordinates);
	nodeCoordinates(0) =  3000.0;
	output<<nodeCoordinates(0)<<"\n";
    myStructure.NodeCreate(2, nodeCoordinates);
	nodeCoordinates(0) =  4000.0;
	output<<nodeCoordinates(0)<<"\n";
    myStructure.NodeCreate(3, nodeCoordinates);
	nodeCoordinates(0) =  6000.0;
	output<<nodeCoordinates(0)<<"\n";
    myStructure.NodeCreate(4, nodeCoordinates);
	output.close();


    int interpolationType = myStructure.InterpolationTypeCreate("TRUSS1D");
    myStructure.InterpolationTypeAdd(interpolationType, NuTo::Node::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    myStructure.InterpolationTypeAdd(interpolationType, NuTo::Node::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);

//	for(int node = 0; node < NumElements + 1; node++)
//	{
//		std::cout << "create node: " << node << " coordinates: " << node * Length/NumElements << std::endl;
//		nodeCoordinates(0) = node * Length/NumElements;
//		myStructure.NodeCreate(node, "displacements", nodeCoordinates);
//	}

	filename="Truss1D2NReference-E";
	output.open(filename.c_str());
	if(output)
	{
		output<<"#"<<filename<<"\n";
		output<<"modul\n";
	}
	else
		std::cout<<__LINE__<<" Truss1D2NReference] error output file.\n";
	// create elements
	NuTo::FullVector<int,Eigen::Dynamic> elementIncidence(2);
	for(int element = 0; element < NumElements; element++)
	{
		std::cout <<  "create element: " << element << " nodes: " << element << "," << element+1 << std::endl;
		elementIncidence(0) = element;
		elementIncidence(1) = element + 1;
        myStructure.ElementCreate(interpolationType, elementIncidence);
        myStructure.ElementSetSection(element,Section1);
//			myStructure.ElementSetConstitutiveLaw(element,Material1);
		if(element<2)
		{
			myStructure.ElementSetConstitutiveLaw(element,Material1);
            output<<myStructure.ConstitutiveLawGetParameterDouble(Material1,NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS)<<"\n";
		}
		else
		{
			myStructure.ElementSetConstitutiveLaw(element,Material2);
            output<<myStructure.ConstitutiveLawGetParameterDouble(Material2,NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS)<<"\n";
		}

	}
	output.close();

    myStructure.ElementTotalConvertToInterpolationType(1e-6,3);

	filename="Truss1D2NReference-ipCoor";
	output.open(filename.c_str());
	if(output)
	{
		output<<"#"<<filename<<"\n";
		output<<"ipCoor\n";
		output<<"1000\n2500\n3500\n5000\n";
		output.close();
	}
	else
		std::cout<<__LINE__<<" Truss1D2NReference] error output file.\n";

	// set boundary conditions and loads
	NuTo::FullVector<double,Eigen::Dynamic> direction(1);
	direction(0) = 1;
	myStructure.ConstraintLinearSetDisplacementNode(0, direction, 0.0);
	if(EnableDisplacementControl)
	{
		std::cout << "Displacement control" << std::endl;
		myStructure.ConstraintLinearSetDisplacementNode(NumElements, direction, BoundaryDisplacement);
	}
	else
	{
		std::cout << "Load control" << std::endl;
		myStructure.LoadCreateNodeForce(1,NumElements, direction, Force);
	}
	// start analysis
	// build global dof numbering
	myStructure.NodeBuildGlobalDofs();

	// build global stiffness matrix and equivalent load vector which correspond to prescribed boundary values
	NuTo::SparseMatrixCSRVector2General<double> stiffnessMatrixVec;
	NuTo::FullVector<double,Eigen::Dynamic> dispForceVector;
	myStructure.CalculateMaximumIndependentSets();
	myStructure.BuildGlobalCoefficientMatrix0(stiffnessMatrixVec, dispForceVector);
	NuTo::SparseMatrixCSRGeneral<double> stiffnessMatrix(stiffnessMatrixVec);

	// build global external load vector
	NuTo::FullVector<double,Eigen::Dynamic> extForceVector;
	myStructure.BuildGlobalExternalLoadVector(1,extForceVector);

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

	NuTo::FullVector<double,Eigen::Dynamic> rDisplacements;
	filename="Truss1D2NReference-disp";
	output.open(filename.c_str());
	if(output)
	{
		output<<"#"<<filename<<"\n";
		output<<"disp\n 0.\n";
		if(EnableDisplacementControl)
		{
			for(int element = 1; element < NumElements; element++)
			{
				myStructure.NodeGetDisplacements(element,rDisplacements);
				output<<rDisplacements<<"\n";
			}
			output<<"1.\n";
		}
		else
			output<<displacementVector;
		output.close();
	}


	// calculate residual
	NuTo::FullVector<double,Eigen::Dynamic> intForceVector;
	myStructure.BuildGlobalGradientInternalPotentialVector(intForceVector);
	NuTo::FullVector<double,Eigen::Dynamic> residualVector = extForceVector - intForceVector;
	std::cout << "residual: " << residualVector.Norm() << std::endl;

	// calculate Stress and strain per element
	NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> rEngineeringStress0;
	NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> rEngineeringStrain0;


	myStructure.ElementGetEngineeringStrain(0,rEngineeringStrain0);
	myStructure.ElementGetEngineeringStress(0,rEngineeringStress0);
	std::cout << "element 0: strain\n" << rEngineeringStrain0 << std::endl;
	std::cout << "           stress\n" << rEngineeringStress0 << std::endl;

	myStructure.ElementGetEngineeringStrain(1,rEngineeringStrain0);
	myStructure.ElementGetEngineeringStress(1,rEngineeringStress0);
	std::cout << "element 1: strain\n" << rEngineeringStrain0 << std::endl;
	std::cout << "           stress\n" << rEngineeringStress0 << std::endl;

	myStructure.ElementGetEngineeringStrain(2,rEngineeringStrain0);
	myStructure.ElementGetEngineeringStress(2,rEngineeringStress0);
	std::cout << "element 2: strain\n" << rEngineeringStrain0 << std::endl;
	std::cout << "           stress\n" << rEngineeringStress0 << std::endl;

	myStructure.ElementGetEngineeringStrain(3,rEngineeringStrain0);
	myStructure.ElementGetEngineeringStress(3,rEngineeringStress0);
	std::cout << "element 3: strain\n" << rEngineeringStrain0 << std::endl;
	std::cout << "           stress\n" << rEngineeringStress0 << std::endl;


	filename="Truss1D2NReference-strain00";
	output.open(filename.c_str());
	if(output)
	{
		output<<"#"<<filename<<"\n";
		output<<"strain00\n";
		for(int element = 0; element < NumElements; element++)
		{
			myStructure.ElementGetEngineeringStrain(element,rEngineeringStrain0);
			output<<rEngineeringStrain0.GetRow(0)<<"\n";
		}
		output.close();
	}
	else
		std::cout<<__LINE__<<" Truss1D2NReference] error output file.\n";


	filename="Truss1D2NReference-stress00";
	output.open(filename.c_str());
	if(output)
	{
		output<<"#"<<filename<<"\n";
		output<<"stress00\n";
		for(int element = 0; element < NumElements; element++)
		{
			myStructure.ElementGetEngineeringStress(element,rEngineeringStress0);
			output<<rEngineeringStress0.GetRow(0)<<"\n";

		}
		output.close();
	}
	else
		std::cout<<__LINE__<<" Truss1D2NReference] error output file.\n";

#ifdef ENABLE_VISUALIZE
	// visualize results
    myStructure.AddVisualizationComponentDisplacements();
    myStructure.AddVisualizationComponentEngineeringStrain();
    myStructure.AddVisualizationComponentEngineeringStress();
    myStructure.AddVisualizationComponentConstitutive();
	myStructure.ExportVtkDataFileElements("Truss1D2NReference.vtk");
#endif

#else
    std::cout << "MUMPS not available - can't solve system of equations " << std::endl;
#endif // HAVE_MUMPS
	return 0;
}
