#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"

#include "nuto/math/SparseMatrixCSRGeneral.h"
#include <eigen3/Eigen/Core>


int main()
{
    try
    {
	//create structure
	NuTo::Structure myStructure(2);

	//3x3 nodes 2x2 element grid
	//create nodes
    NuTo::FullMatrix<double> Coordinates(2,1);
	Coordinates(0,0) = 0.0;
	Coordinates(1,0) = 0.0;
	int node1 = myStructure.NodeCreate("displacements",Coordinates);

	Coordinates(0,0) = 1.0;
	Coordinates(1,0) = 0.0;
	int node2 = myStructure.NodeCreate("displacements",Coordinates);

	Coordinates(0,0) = 2.0;
	Coordinates(1,0) = 0.0;
	int node3 = myStructure.NodeCreate("displacements",Coordinates);

	Coordinates(0,0) = 1.0;
	Coordinates(1,0) = 1.0;
	int node5 = myStructure.NodeCreate("displacements",Coordinates);

	Coordinates(0,0) = 0.0;
	Coordinates(1,0) = 1.0;
	int node4 = myStructure.NodeCreate("displacements",Coordinates);

	Coordinates(0,0) = 2.0;
	Coordinates(1,0) = 1.0;
	int node6 = myStructure.NodeCreate("displacements",Coordinates);

	Coordinates(0,0) = 0.0;
    Coordinates(1,0) = 2.0;
	int node7 = myStructure.NodeCreate("displacements",Coordinates);

	Coordinates(0,0) = 1.0;
	Coordinates(1,0) = 2.0;
	int node8 = myStructure.NodeCreate("displacements",Coordinates);

	Coordinates(0,0) = 2.0;
	Coordinates(1,0) = 2.0;
	int node9 = myStructure.NodeCreate("displacements",Coordinates);

	//create elements
    NuTo::FullMatrix<int> Incidence(4,1);
	Incidence(0,0) = node1;
	Incidence(1,0) = node2;
	Incidence(2,0) = node5;
	Incidence(3,0) = node4;
    int myElement1 = myStructure.ElementCreate("PLANE2D4N",Incidence,"ConstitutiveLawIpNonlocal","StaticDataNonlocal");
    myStructure.ElementSetIntegrationType(myElement1,"2D4NGauss1Ip","StaticDataNonlocal");
		
	Incidence(0,0) = node2;
	Incidence(1,0) = node3;
	Incidence(2,0) = node6;
	Incidence(3,0) = node5;
    int myElement2 = myStructure.ElementCreate("PLANE2D4N",Incidence,"ConstitutiveLawIpNonlocal","StaticDataNonlocal");
    myStructure.ElementSetIntegrationType(myElement2,"2D4NGauss4Ip","StaticDataNonlocal");
		
	Incidence(0,0) = node4;
	Incidence(1,0) = node5;
	Incidence(2,0) = node8;
	Incidence(3,0) = node7;
    int myElement3 = myStructure.ElementCreate("PLANE2D4N",Incidence,"ConstitutiveLawIpNonlocal","StaticDataNonlocal");
    myStructure.ElementSetIntegrationType(myElement3,"2D4NGauss1Ip","StaticDataNonlocal");
		
/*	Incidence(0,0) = node5;
	Incidence(1,0) = node6;
	Incidence(2,0) = node9;
	Incidence(3,0) = node8;
    int myElement4 = myStructure.ElementCreate("PLANE2D4N",Incidence,"ConstitutiveLawIpNonlocal","StaticDataNonlocal");
    myStructure.ElementSetIntegrationType(myElement4,"2D4NGauss4Ip","StaticDataNonlocal");
*/
    NuTo::FullMatrix<int> Incidence3(3,1);

    Incidence3(0,0) = node5;
   	Incidence3(1,0) = node6;
   	Incidence3(2,0) = node9;
   	int myElement4 = myStructure.ElementCreate("PLANE2D3N",Incidence3,"ConstitutiveLawIpNonlocal","StaticDataNonlocal");
    myStructure.ElementSetIntegrationType(myElement4,"2D3NGauss1Ip","StaticDataNonlocal");

    Incidence3(0,0) = node5;
   	Incidence3(1,0) = node9;
   	Incidence3(2,0) = node8;
   	int myElement5 = myStructure.ElementCreate("PLANE2D3N",Incidence3,"ConstitutiveLawIpNonlocal","StaticDataNonlocal");
    myStructure.ElementSetIntegrationType(myElement5,"2D3NGauss1Ip","StaticDataNonlocal");

	//create constitutive law
	int myMatDamage = myStructure.ConstitutiveLawCreate("NonlocalDamagePlasticity");
	myStructure.ConstitutiveLawSetYoungsModulus(myMatDamage,9);
	myStructure.ConstitutiveLawSetPoissonsRatio(myMatDamage,0.25);
	myStructure.ConstitutiveLawSetNonlocalRadius(myMatDamage,2.);
	myStructure.ConstitutiveLawSetTensileStrength(myMatDamage,2);
	myStructure.ConstitutiveLawSetCompressiveStrength(myMatDamage,20);
	myStructure.ConstitutiveLawSetBiaxialCompressiveStrength(myMatDamage,25);
	myStructure.ConstitutiveLawSetFractureEnergy(myMatDamage,0.2);

	int myMatLin = myStructure.ConstitutiveLawCreate("LinearElastic");
	myStructure.ConstitutiveLawSetYoungsModulus(myMatLin,10);
	myStructure.ConstitutiveLawSetPoissonsRatio(myMatLin,0.25);

	//create section
	int mySection = myStructure.SectionCreate("Plane_Strain");
	myStructure.SectionSetThickness(mySection,5);

	//assign constitutive law 
	myStructure.ElementTotalSetSection(mySection);
	myStructure.ElementTotalSetConstitutiveLaw(myMatDamage);

	//Build nonlocal elements
	myStructure.BuildNonlocalData(myMatDamage);

#ifdef ENABLE_VISUALIZE
	// visualize results
	myStructure.AddVisualizationComponentNonlocalWeights(myElement1,0);

	myStructure.AddVisualizationComponentNonlocalWeights(myElement2,0);
	myStructure.AddVisualizationComponentNonlocalWeights(myElement2,1);
	myStructure.AddVisualizationComponentNonlocalWeights(myElement2,2);
	myStructure.AddVisualizationComponentNonlocalWeights(myElement2,3);
#endif

	//build maximum independent sets for openmp parallel assembly
	myStructure.CalculateMaximumIndependentSets();

    //calculate linear elastic matrix
	NuTo::SparseMatrixCSRVector2General<double> stiffnessMatrix;
	NuTo::FullMatrix<double> dispForceVector;

	myStructure.ElementTotalUpdateTmpStaticData();
	myStructure.BuildGlobalCoefficientMatrix0(stiffnessMatrix, dispForceVector);
	stiffnessMatrix.RemoveZeroEntries(0,1e-14);

	NuTo::FullMatrix<double> fullStiffnessMatrixElastic(stiffnessMatrix);
	std::cout<<"stiffnessMatrix elastic"<<std::endl;
	fullStiffnessMatrixElastic.Info();
	NuTo::FullMatrix<double> displacements;
	NuTo::FullMatrix<double> dependentDofs;
	NuTo::FullMatrix<double> intForce;
	NuTo::FullMatrix<double> intForce2;

	//check the stiffness twice, once in the initial deformed state
	//and once after the update (should be linear elastic)
	//loadstep 0 : uniform plastic loading
	//loadstep 1 : unloading to zero
	//loadstep 2 : nonuniform loading, some elements unloading
	for (int theLoadStep=0; theLoadStep<3; theLoadStep++)
	{
		//apply displacements
		double rightDisp;
		switch (theLoadStep)
		{
		case 0:
		    rightDisp = 0.5;
		break;
		case 1:
		    rightDisp = 0.0;
		break;
		case 2:
			rightDisp = 0.6;
		}

		NuTo::FullMatrix<double>matrixRightDisp(2,1);
		matrixRightDisp.SetValue(0,0,rightDisp);
		matrixRightDisp.SetValue(1,0,0.);

		myStructure.NodeSetDisplacements(node3,matrixRightDisp);
		myStructure.NodeSetDisplacements(node6,matrixRightDisp);
		myStructure.NodeSetDisplacements(node9,matrixRightDisp);

		NuTo::FullMatrix<double>matrixCenterDisp(2,1);
		if (theLoadStep!=2)
		    matrixCenterDisp.SetValue(0,0,0.5*rightDisp);
		else
			matrixCenterDisp.SetValue(0,0,0.4*rightDisp);
		matrixCenterDisp.SetValue(1,0,0.);

		myStructure.NodeSetDisplacements(node2,matrixCenterDisp);
		myStructure.NodeSetDisplacements(node5,matrixCenterDisp);
		myStructure.NodeSetDisplacements(node8,matrixCenterDisp);

		NuTo::FullMatrix<double>matrixLeftDisp(2,1);
		matrixLeftDisp.SetValue(0,0,0.0);
		matrixLeftDisp.SetValue(1,0,0.);

		myStructure.NodeSetDisplacements(node1,matrixLeftDisp);
		myStructure.NodeSetDisplacements(node4,matrixLeftDisp);
		myStructure.NodeSetDisplacements(node7,matrixLeftDisp);

	    //NuTo::FullMatrix<double>matrixLeftDispNode1(2,1);
	    //matrixLeftDispNode1.SetValue(0,0,0.0);
	    //matrixLeftDispNode1.SetValue(1,0,0.);
	    //myStructure.NodeSetDisplacements(node1,matrixLeftDispNode1);


		myStructure.ElementTotalUpdateTmpStaticData();

		myStructure.BuildGlobalCoefficientMatrix0(stiffnessMatrix, dispForceVector);
		stiffnessMatrix.RemoveZeroEntries(0,1e-14);

		NuTo::FullMatrix<double> fullStiffnessMatrix(stiffnessMatrix);

		std::cout<<"stiffnessMatrix analytic"<<std::endl;
		fullStiffnessMatrix.Info();

		myStructure.NodeExtractDofValues(displacements,dependentDofs);
		myStructure.BuildGlobalGradientInternalPotentialVector(intForce);

		double delta(1e-8);
		NuTo::FullMatrix<double> stiffnessMatrixCD(displacements.GetNumRows(),displacements.GetNumRows());

		//check with central differences
		for (int count2=0; count2<displacements.GetNumRows(); count2++)
		{
			displacements(count2,0) = displacements(count2,0) + delta;
			myStructure.NodeMergeActiveDofValues(displacements);
			myStructure.ElementTotalUpdateTmpStaticData();
			myStructure.BuildGlobalGradientInternalPotentialVector(intForce2);
			//std::cout<<"intForce delta "<< count << std::endl;
			//intForce2.Info();
			stiffnessMatrixCD.SetColumn(count2,(intForce2-intForce)*(1/delta));
			displacements(count2,0) = displacements(count2,0) - delta;
			myStructure.NodeMergeActiveDofValues(displacements);
			myStructure.ElementTotalUpdateTmpStaticData();
		}
		std::cout << "stiffnessMatrixCD" << std::endl;
		stiffnessMatrixCD.Info();
		int row,col;
		double maxerror((fullStiffnessMatrix-stiffnessMatrixCD).Abs().Max(row,col));
	    switch(theLoadStep)
	    {
	    case 0:
	        std::cout << "max difference in stiffness matrix for uniform plastic loading " << maxerror << " row " << row << " col " << col << std::endl;
	    break;
	    case 1:
	    {
	        std::cout << "max difference in stiffness matrix for unloading " << maxerror << " row " << row << " col " << col << std::endl;
	    	double omega(fullStiffnessMatrix(0,0)/fullStiffnessMatrixElastic(0,0));
	    	double maxerror2((fullStiffnessMatrixElastic*omega-stiffnessMatrix).Abs().Max(row,col));
			std::cout<<"stiffnessMatrix elastic*omega"<<std::endl;
			(fullStiffnessMatrixElastic*omega-stiffnessMatrix).Info();
    	    std::cout << "max difference in stiffness matrix for unloading and scaled elastic matrix " << maxerror2 << " row " << row << " col " << col << std::endl;
	    }
    	break;
	    case 2:
	        std::cout << "max difference in stiffness matrix for nonuniform plastic loading/unloading " << maxerror << " row " << row << " col " << col << std::endl;
	    break;
	    }
		//update the structure
		myStructure.ElementTotalUpdateStaticData();
	}

#ifdef ENABLE_VISUALIZE
	myStructure.AddVisualizationComponentDisplacements();
	myStructure.AddVisualizationComponentEngineeringStrain();
	myStructure.AddVisualizationComponentEngineeringStress();
	myStructure.AddVisualizationComponentDamage();
	myStructure.AddVisualizationComponentEngineeringPlasticStrain();
	myStructure.ExportVtkDataFile("NonlocalDamagePlasticityModel.vtk");
#endif

    NuTo::FullMatrix<double> shift(2,1);
    shift(0,0) = 3+myStructure.ConstitutiveLawGetNonlocalRadius(myMatDamage);
    shift(1,0) = 0;
    myStructure.CopyAndTranslate(shift);

#ifdef ENABLE_VISUALIZE
	myStructure.ExportVtkDataFile("NonlocalDamagePlasticityModelCopied.vtk");
#endif

    }
    catch (NuTo::Exception& e)
    {
        std::cout << e.ErrorMessage() << std::endl;
    }
    return 0;
}
