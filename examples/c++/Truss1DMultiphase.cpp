// $Id$

//Author: Andrea Keßler
#include <iostream>
#include "nuto/math/FullMatrix.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"
#include "nuto/math/SparseDirectSolverMUMPS.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/elements/Truss1D3N.h"

// arguments:
// 		0 - 2 for const to quad,
// 		nbrNodes of multiphase element (2 or 3)
//		nbrIps of multiphase element (1 to 4)
int main(int argc, char *argv[])
{
	// material function in multiphase element is:
	// quadratic

	std::cout<<"[Truss1DMultiphase] start.\n";
//	std::string modulFunction="Quadratic";
	std::string modulFunction="Linear";
//	std::string modulFunction="Constant";

	std::string elemType="Truss1D3N";
	int nbrNodes=3;
	std::string ips="2Ip";
	int nbrIps=2;


	if(argc==4) //arguments (progamm, ..)
	{
		if(atoi(argv[1])==0)
			modulFunction="Constant";
		else if(atoi(argv[1])==1)
			modulFunction="Linear";
		else if(atoi(argv[1])==2)
			modulFunction="Quadratic";
		else
			std::cout<<"[Truss1DMultiphase] error input argument.\n";

		elemType="Truss1D";
		elemType+=argv[2];
		elemType+="N";
		nbrNodes=atoi(argv[2]);
		ips=argv[3];
		ips+="Ip";
		nbrIps=atoi(argv[3]);
	}
	std::cout<<"[Truss1DMultiphase] "<<elemType+ips+modulFunction<<"\n";
	std::cout<<"[Truss1DMultiphase] nbrNodes: "<<nbrNodes<<" nbrIps: "<<nbrIps<<"\n";

	std::string intType="1D2NGauss"+ips;
	std::ofstream output;

	// definitions
	double YoungsModulus = 10000.; //basis modulus at node 0
	double YoungsModulusDiff = 20000.; //difference from node 0 to last node, +/-
	int NumElements = 7;
	double Area = 10. * 10.;
	double Length = 6000.;
	double Force = 1.;
	bool EnableDisplacementControl = true;
	double BoundaryDisplacement = 1;

	// create one-dimensional structure
	NuTo::Structure myStructure(1);
	myStructure.SetVerboseLevel(0);
	// create section
	int Section1 = myStructure.SectionCreate("Truss");
	myStructure.SectionSetArea(Section1, Area);

	// create material law
	int Material1 = myStructure.ConstitutiveLawCreate("LINEARELASTICENGINEERINGSTRESS");
	int Material2 = myStructure.ConstitutiveLawCreate("LINEARELASTICENGINEERINGSTRESS");
	myStructure.ConstitutiveLawSetYoungsModulus(Material1, YoungsModulus);
//    myStructure.ConstitutiveLawSetPoissonsRatio(Material1,0.2);
//    myStructure.ConstitutiveLawSetPoissonsRatio(Material2,0.2);
	myStructure.ConstitutiveLawSetYoungsModulus(Material2, YoungsModulus+YoungsModulusDiff);

	// create nodes
	NuTo::FullVector<double,Eigen::Dynamic> nodeCoordinates(1);
	std::string filename=elemType+ips+"nodes";
	output.open(filename.c_str());
	if(output)
	{
		output<<"#"<<filename<<"\n";
		output<<"nodes\n";
	}
	else
		std::cout<<"[Truss1DMultiphase] error output file.\n";

	// +1 for one element with 3 nodes
	int numNodes=NumElements+nbrNodes-1;
	for(int node = 0; node < numNodes; node++)
	{
		if(nbrNodes>2)
		{
//		std::cout<<__LINE__<<" TEST\n";
			if(node<(int)(NumElements + 1)/2)
			{
				nodeCoordinates(0) = node * Length/NumElements;
				std::cout << "create node: " << node << " coordinates: " << nodeCoordinates(0) << std::endl;
				myStructure.NodeCreate(node, "displacements", nodeCoordinates);
				output<<nodeCoordinates(0)<<"\n";

			}
			else if(node==(int)(NumElements + 1)/2)
			{
				// create middle node
				if(nbrNodes==3)
				{
					nodeCoordinates(0) = Length/2;
					myStructure.NodeCreate(node, "displacements", nodeCoordinates);
					std::cout << "create node: " << node << " coordinates: " << nodeCoordinates(0) << std::endl;
					output<<nodeCoordinates(0)<<"\n";
				}
				else
					std::cout<<"Truss1DMultiphase] error only 3 nodes implemented.\n";

			}
			else
			{
				nodeCoordinates(0) = (node-1) * Length/NumElements;
				std::cout << "create node: " << node << " coordinates: " << nodeCoordinates(0) << std::endl;
				myStructure.NodeCreate(node, "displacements", nodeCoordinates);
				output<<nodeCoordinates(0)<<"\n";
			}
		}
		else
		{
//		std::cout<<__LINE__<<" TEST\n";
			if(nbrNodes<2)
				std::cout<<"[Truss1DMultiphase] error number of nodes is to small.\n";
			//nbrNodes==2
			nodeCoordinates(0) = node * Length/NumElements;
			std::cout << "create node: " << node << " coordinates: " << nodeCoordinates(0) << std::endl;
			myStructure.NodeCreate(node, "displacements", nodeCoordinates);
			output<<nodeCoordinates(0)<<"\n";
		}
	}
	output.close();

	filename=elemType+ips+"E";
	output.open(filename.c_str());
	if(output)
	{
		output<<"#"<<filename<<"\n";
		output<<"modul\n";
	}
	else
		std::cout<<__LINE__<<" Truss1D3NMultiphase] error output file.\n";

	// create elements
	NuTo::FullVector<int,Eigen::Dynamic> elementIncidence(2);
	NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>ipCoor;
	for(int element = 0; element < NumElements; element++)
	{
		if(element<(NumElements-1)/2)
		{
			std::cout <<  "create element: " << element << " nodes: " << element << "," << element+1 << std::endl;
			elementIncidence(0) = element;
			elementIncidence(1) = element + 1;
			myStructure.ElementCreate(element, "Truss1D2N", elementIncidence);
			myStructure.ElementSetSection(element,Section1);
			myStructure.ElementSetConstitutiveLaw(element,Material1);
			output<<myStructure.ConstitutiveLawGetYoungsModulus(Material1)<<"\n";
		}
		else if (element==(NumElements-1)/2)
		{
			NuTo::FullVector<int,Eigen::Dynamic> elementIncidence2(nbrNodes);
			std::cout <<  "create element: " << element << " nodes: " ;
			for(int node=0;node<nbrNodes;node++)
			{
				elementIncidence2(node) = element+node;
				std::cout<< elementIncidence2(node) << "," ;
			}
			std::cout<< std::endl;
			myStructure.ElementCreate(element, elemType, elementIncidence2,"VariableConstitutiveLawIp","StaticData");
			myStructure.ElementSetSection(element,Section1);
			myStructure.ElementSetIntegrationType(element,intType,"STATICDATA");
			myStructure.ElementGetIntegrationPointCoordinates(element,ipCoor);
			std::cout<<"[Truss1DMultiphase] create material data for mp.\n";


			std::vector<int> mpMat;
//			std::vector<double> youngsMod;
			double youngsMod;
			for(int thisIp=0;thisIp<nbrIps;thisIp++)
			{
				//ip Coor global coordinates
				double thisIpCoor=ipCoor.GetValue(0,thisIp);
				// local coordinates
				double ipCoorLocal=(thisIpCoor-0.5*Length)/(Length/NumElements);
				std::cout<<__LINE__<<"[Truss1D3NMultiphase] thisIpCoor "<<thisIpCoor<<" "<<ipCoorLocal<<"\n";

				// linear variation of modulus
				if(modulFunction=="Linear")
					//NEW: fct for modulus: x\in (-1,1) ->f(x)=\Delta E/2 *(x+1) + E_1
					youngsMod=(YoungsModulus+(1+ipCoorLocal)*YoungsModulusDiff/2);
				else if (modulFunction=="Quadratic")
				//NEW: fct for modulus: x\in (-1,1) ->f(x)=-1/4 * \Delta E x^3+ 3/4 * \Delta E x+ E_1 + \Delta E/2
					youngsMod=(-0.25 * YoungsModulusDiff * ipCoorLocal *ipCoorLocal * ipCoorLocal
						+  0.75 * YoungsModulusDiff * ipCoorLocal
						+ YoungsModulus + 0.5 * YoungsModulusDiff);
				else if (modulFunction=="Constant")
					youngsMod=(YoungsModulus + 0.5 * YoungsModulusDiff);
				else
					std::cout<<"[Truss1D3NMultiphase] ERROR wrong type of modulus function.\n";

				mpMat.push_back(myStructure.ConstitutiveLawCreate("LINEARELASTICENGINEERINGSTRESS"));
				myStructure.ConstitutiveLawSetYoungsModulus(mpMat[thisIp],youngsMod);
				myStructure.ElementSetConstitutiveLaw(element,thisIp,mpMat[thisIp]);
				output<<myStructure.ConstitutiveLawGetYoungsModulus(mpMat[thisIp])<<"\n";
				std::cout<<"mat["<<thisIp<<"]="<<myStructure.ConstitutiveLawGetYoungsModulus(mpMat[thisIp])<<"\n";
//				std::cout<<"mat["<<thisIp<<"] fac="<<(((ipCoor.GetValue(0,thisIp))-2000)/2000)<<"\n";

			}
//			myStructure.ElementSetConstitutiveLaw(element,0,Material1);
//			myStructure.ElementSetConstitutiveLaw(element,1,Material2);
		}
		else
		{
			std::cout <<  "create element: " << element <<
					" nodes: " << element+nbrNodes-2 << "," << element+nbrNodes-1 << std::endl;
			elementIncidence(0) = element+nbrNodes-2;
			elementIncidence(1) = element+nbrNodes-1;
			myStructure.ElementCreate(element, "Truss1D2N", elementIncidence);
			myStructure.ElementSetSection(element,Section1);
			myStructure.ElementSetConstitutiveLaw(element,Material2);
			output<<myStructure.ConstitutiveLawGetYoungsModulus(Material2)<<"\n";
		}
	}
	output.close();
	//write ipCoor to file
	filename=elemType+ips+"ipCoor";
	output.open(filename.c_str());
	if(output)
	{
		output<<"#"<<filename<<"\n";
		output<<"ipCoor\n";
		for(int element=0;element<NumElements;element++)
		{
			myStructure.ElementGetIntegrationPointCoordinates(element,ipCoor);
			if(element==(NumElements-1)/2)
			{
				for(int ip=0;ip<nbrIps;ip++)
					output<<ipCoor.GetValue(0,ip)<<"\n";
			}
			else
				output<<ipCoor.GetValue(0,0)<<"\n";
		}
		output.close();
	}
	else
		std::cout<<__LINE__<<" Truss1D3NMultiphase] error output file.\n";


	// set boundary conditions and loads
	NuTo::FullVector<double,Eigen::Dynamic> direction(1);
	direction(0) = 1;
	myStructure.ConstraintLinearSetDisplacementNode(0, direction, 0.0);
	if(EnableDisplacementControl)
	{
		std::cout << "Displacement control" << std::endl;
		myStructure.ConstraintLinearSetDisplacementNode(numNodes-1, direction, BoundaryDisplacement);
	}
	else
	{
		std::cout << "Load control" << std::endl;
		myStructure.LoadCreateNodeForce(1,numNodes-1, direction, Force);
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

//	std::cout<<__LINE__<<" [Truss1D3NMultiphase]   \n";

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
	filename=elemType+ips+"disp";
	output.open(filename.c_str());
	if(output)
	{
		output<<"#"<<filename<<"\n";
		output<<"disp \n";
		std::cout<<"disp \n";
		if(EnableDisplacementControl)
		{
			for(int node = 0; node < numNodes; node++)
			{
				myStructure.NodeGetDisplacements(node,rDisplacements);

				std::cout<<rDisplacements<<"\n";
				output<<rDisplacements<<"\n";
			}
//			output<<"1.\n";
		}
		else
		{
			output<<"0. \n";
			output<<displacementVector;
			std::cout<<"disp \n0. \n";
			std::cout<<displacementVector;
		}
		std::cout<<"\n";
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

//	std::cout<<__LINE__<<" Test\n";
	for(int element=0;element<NumElements;element++)
	{
		myStructure.ElementGetEngineeringStrain(element,rEngineeringStrain0);
		myStructure.ElementGetEngineeringStress(element,rEngineeringStress0);
		std::cout << "element "<<element<<": strain\n" << rEngineeringStrain0 << std::endl;
		std::cout << "           stress\n" << rEngineeringStress0 << std::endl;
	}

	filename=elemType+ips+"strain00";
	output.open(filename.c_str());
	if(output)
	{
		output<<"#"<<filename<<"\n";
		output<<"strain00\n";
		for(int element = 0; element < NumElements; element++)
		{
			myStructure.ElementGetEngineeringStrain(element,rEngineeringStrain0);

			if(element==(NumElements-1)/2)
			{
				for(int ip=0;ip<nbrIps;ip++)
					output<<rEngineeringStrain0.GetValue(0,ip)<<"\n";
			}
			else
				output<<rEngineeringStrain0.GetRow(0)<<"\n";

		}
		output.close();
	}
	else
		std::cout<<__LINE__<<" Truss1D3NMultiphase] error output file.\n";


	filename=elemType+ips+"stress00";
	output.open(filename.c_str());
	if(output)
	{
		output<<"#"<<filename<<"\n";
		output<<"stress00\n";
		for(int element = 0; element < NumElements; element++)
		{
			myStructure.ElementGetEngineeringStress(element,rEngineeringStress0);
			if(element==(NumElements-1)/2)
			{
				for(int ip=0;ip<nbrIps;ip++)
					output<<rEngineeringStress0.GetValue(0,ip)<<"\n";
			}
			else
				output<<rEngineeringStress0.GetRow(0)<<"\n";

		}
		output.close();
	}
	else
		std::cout<<__LINE__<<" Truss1D3NMultiphase] error output file.\n";

#ifdef ENABLE_VISUALIZE
	// visualize results
    myStructure.AddVisualizationComponentDisplacements();
    myStructure.AddVisualizationComponentEngineeringStrain();
    myStructure.AddVisualizationComponentEngineeringStress();
    myStructure.AddVisualizationComponentConstitutive();

	myStructure.ExportVtkDataFileElements(elemType+ips+"Multiphase.vtk");
#endif

#else
    std::cout << "MUMPS not available - can't solve system of equations " << std::endl;
#endif // HAVE_MUMPS
	return 0;
}
