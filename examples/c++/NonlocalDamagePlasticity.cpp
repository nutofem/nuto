#include "nuto/math/FullMatrix.h"
#include "nuto/math/FullVector.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"

//just for test
#include "nuto/mechanics/constitutive/mechanics/DeformationGradient2D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain2D.h"
#include "nuto/mechanics/constitutive/mechanics/NonlocalDamagePlasticity.h"
#include <eigen2/Eigen/Core>


int main()
{
    try
    {
	//create structure
	NuTo::Structure myStructure(2);

	//3x3 nodes 2x2 element grid
	//create nodes
	NuTo::FullVector<double,Eigen::Dynamic> Coordinates(2);
	Coordinates(0) = 0.0;
	Coordinates(1) = 0.0;
	int node1 = myStructure.NodeCreate("displacements",Coordinates);

	Coordinates(0) = 1.0;
	Coordinates(1) = 0.0;
	int node2 = myStructure.NodeCreate("displacements",Coordinates);

	Coordinates(0) = 2.0;
	Coordinates(1) = 0.0;
	int node3 = myStructure.NodeCreate("displacements",Coordinates);

	Coordinates(0) = 0.0;
	Coordinates(1) = 1.0;
	int node4 = myStructure.NodeCreate("displacements",Coordinates);

	Coordinates(0) = 1.0;
	Coordinates(1) = 1.0;
	int node5 = myStructure.NodeCreate("displacements",Coordinates);

	Coordinates(0) = 2.0;
	Coordinates(1) = 1.0;
	int node6 = myStructure.NodeCreate("displacements",Coordinates);

	Coordinates(0) = 0.0;
	Coordinates(1,0) = 2.0;
	int node7 = myStructure.NodeCreate("displacements",Coordinates);

	Coordinates(0) = 1.0;
	Coordinates(1) = 2.0;
	int node8 = myStructure.NodeCreate("displacements",Coordinates);

	Coordinates(0) = 2.0;
	Coordinates(1) = 2.0;
	int node9 = myStructure.NodeCreate("displacements",Coordinates);

	//create elements
	NuTo::FullVector<int,Eigen::Dynamic> Incidence(4);
	Incidence(0) = node1;
	Incidence(1) = node2;
	Incidence(2) = node5;
	Incidence(3) = node4;
    int myElement1 = myStructure.ElementCreate("PLANE2D4N",Incidence,"ConstitutiveLawIpNonlocal","StaticDataNonlocal");
//    myStructure.ElementSetIntegrationType(myElement1,"2D4NGauss1Ip","StaticDataNonlocal");

	Incidence(0) = node2;
	Incidence(1) = node3;
	Incidence(2) = node6;
	Incidence(3) = node5;
    int myElement2 = myStructure.ElementCreate("PLANE2D4N",Incidence,"ConstitutiveLawIpNonlocal","StaticDataNonlocal");
		
	Incidence(0) = node4;
	Incidence(1) = node5;
	Incidence(2) = node8;
	Incidence(3) = node7;
    int myElement3 = myStructure.ElementCreate("PLANE2D4N",Incidence,"ConstitutiveLawIpNonlocal","StaticDataNonlocal");
		
	Incidence(0) = node5;
	Incidence(1) = node6;
	Incidence(2) = node9;
	Incidence(3) = node8;
    int myElement4 = myStructure.ElementCreate("PLANE2D4N",Incidence,"ConstitutiveLawIpNonlocal","StaticDataNonlocal");

	//create constitutive law
	int myMatLin = myStructure.ConstitutiveLawCreate("NonlocalDamagePlasticity");
	myStructure.ConstitutiveLawSetYoungsModulus(myMatLin,10);
	myStructure.ConstitutiveLawSetPoissonsRatio(myMatLin,0.25);
	myStructure.ConstitutiveLawSetNonlocalRadius(myMatLin,1);
	myStructure.ConstitutiveLawSetTensileStrength(myMatLin,2);
	myStructure.ConstitutiveLawSetCompressiveStrength(myMatLin,20);
	myStructure.ConstitutiveLawSetBiaxialCompressiveStrength(myMatLin,25);
	myStructure.ConstitutiveLawSetFractureEnergy(myMatLin,0.2);

	//create section
	myStructure.SectionCreate("mySection","Plane_Strain");
	myStructure.SectionSetThickness("mySection",1);

	//assign constitutive law 
	myStructure.ElementTotalSetSection("mySection");
	myStructure.ElementTotalSetConstitutiveLaw(myMatLin);

	//Build nonlocal elements
	myStructure.BuildNonlocalData(myMatLin);

	// visualize results
	myStructure.AddVisualizationComponentNonlocalWeights(myElement1,0);
	//myStructure.AddVisualizationComponentNonlocalWeights(myElement1,1);
	//myStructure.AddVisualizationComponentNonlocalWeights(myElement1,2);
	//myStructure.AddVisualizationComponentNonlocalWeights(myElement1,3);
	myStructure.AddVisualizationComponentNonlocalWeights(myElement2,0);
	myStructure.AddVisualizationComponentNonlocalWeights(myElement2,1);
	myStructure.AddVisualizationComponentNonlocalWeights(myElement2,2);
	myStructure.AddVisualizationComponentNonlocalWeights(myElement2,3);
	myStructure.AddVisualizationComponentNonlocalWeights(myElement3,0);
	myStructure.AddVisualizationComponentNonlocalWeights(myElement3,1);
	myStructure.AddVisualizationComponentNonlocalWeights(myElement3,2);
	myStructure.AddVisualizationComponentNonlocalWeights(myElement3,3);
	myStructure.AddVisualizationComponentNonlocalWeights(myElement4,0);
	myStructure.AddVisualizationComponentNonlocalWeights(myElement4,1);
	myStructure.AddVisualizationComponentNonlocalWeights(myElement4,2);
	myStructure.AddVisualizationComponentNonlocalWeights(myElement4,3);
	myStructure.ExportVtkDataFile("PlaneNonlocalWeights.vtk");

	//apply displacements
	double rightDisp(1);
	NuTo::FullVector<double,Eigen::Dynamic>matrixRightDisp(2);
	matrixRightDisp.SetValue(0,rightDisp);
	matrixRightDisp.SetValue(1,0.);

	myStructure.NodeSetDisplacements(node3,matrixRightDisp);
	myStructure.NodeSetDisplacements(node6,matrixRightDisp);
	myStructure.NodeSetDisplacements(node9,matrixRightDisp);

	NuTo::FullVector<double,Eigen::Dynamic> matrixCenterDisp(2);
	matrixCenterDisp.SetValue(0,0.5*rightDisp);
	matrixCenterDisp.SetValue(1,0.);

	myStructure.NodeSetDisplacements(node2,matrixCenterDisp);
	myStructure.NodeSetDisplacements(node5,matrixCenterDisp);
	myStructure.NodeSetDisplacements(node8,matrixCenterDisp);

	NuTo::FullVector<double,Eigen::Dynamic> matrixLeftDisp(2);
	matrixLeftDisp.SetValue(0,0.0);
	matrixLeftDisp.SetValue(1,0.);

	myStructure.NodeSetDisplacements(node2,matrixLeftDisp);
	myStructure.NodeSetDisplacements(node5,matrixLeftDisp);
	myStructure.NodeSetDisplacements(node8,matrixLeftDisp);

	myStructure.ElementTotalUpdateTmpStaticData();
	myStructure.ElementTotalUpdateTmpStaticData();

	//calculate the stiffness matrix
	NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> Ke;
        NuTo::FullMatrix<int,Eigen::Dynamic,Eigen::Dynamic> rowIndex;
    	NuTo::FullMatrix<int,Eigen::Dynamic,Eigen::Dynamic> colIndex;
    	myStructure.ElementStiffness(myElement1,Ke,rowIndex,colIndex);

    }
    catch (NuTo::Exception& e)
    {
        std::cout << e.ErrorMessage() << std::endl;
    }
    return 0;
}
