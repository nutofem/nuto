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
	int Material1 = myStructure.ConstitutiveLawCreate("LINEAR_ELASTIC_ENGINEERING_STRESS");
	int Material2 = myStructure.ConstitutiveLawCreate("LINEAR_ELASTIC_ENGINEERING_STRESS");
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
	myStructure.CalculateMaximumIndependentSets();
    myStructure.SolveGlobalSystemStaticElastic(1);

    auto displacementVector = myStructure.NodeExtractDofValues(0);

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
    auto residual = myStructure.BuildGlobalInternalGradient() - myStructure.BuildGlobalExternalLoadVector(1);

    std::cout << "residual: " << residual.J.CalculateNormL2() << std::endl;

	std::cout << "element 0: strain\n" << myStructure.ElementGetEngineeringStrain(0) << std::endl;
    std::cout << "           stress\n" << myStructure.ElementGetEngineeringStress(0) << std::endl;

    std::cout << "element 0: strain\n" << myStructure.ElementGetEngineeringStrain(1) << std::endl;
    std::cout << "           stress\n" << myStructure.ElementGetEngineeringStress(1) << std::endl;

    std::cout << "element 0: strain\n" << myStructure.ElementGetEngineeringStrain(2) << std::endl;
    std::cout << "           stress\n" << myStructure.ElementGetEngineeringStress(2) << std::endl;

    std::cout << "element 0: strain\n" << myStructure.ElementGetEngineeringStrain(3) << std::endl;
    std::cout << "           stress\n" << myStructure.ElementGetEngineeringStress(3) << std::endl;

	filename="Truss1D2NReference-strain00";
	output.open(filename.c_str());
	if(output)
	{
		output<<"#"<<filename<<"\n";
		output<<"strain00\n";
		for(int element = 0; element < NumElements; element++)
		{
			output<<myStructure.ElementGetEngineeringStrain(element).GetRow(0)<<"\n";
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
            output<<myStructure.ElementGetEngineeringStress(element).GetRow(0)<<"\n";
		}
		output.close();
	}
	else
		std::cout<<__LINE__<<" Truss1D2NReference] error output file.\n";

	// visualize results
    int visualizationGroup = myStructure.GroupCreate(NuTo::Groups::eGroupId::Elements);
    myStructure.GroupAddElementsTotal(visualizationGroup);

    myStructure.AddVisualizationComponent(visualizationGroup, NuTo::VisualizeBase::DISPLACEMENTS);
    myStructure.AddVisualizationComponent(visualizationGroup, NuTo::VisualizeBase::ENGINEERING_STRAIN);
    myStructure.AddVisualizationComponent(visualizationGroup, NuTo::VisualizeBase::ENGINEERING_STRESS);
    myStructure.AddVisualizationComponent(visualizationGroup, NuTo::VisualizeBase::SECTION);
    myStructure.AddVisualizationComponent(visualizationGroup, NuTo::VisualizeBase::CONSTITUTIVE);

	myStructure.ExportVtkDataFileElements("Truss1D2NReference.vtk");

	return 0;
}
