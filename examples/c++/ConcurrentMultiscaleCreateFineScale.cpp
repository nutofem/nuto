#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/structures/unstructured/StructureIp.h"
#include "nuto/mechanics/MechanicsException.h"

#include "nuto/math/SparseDirectSolverMUMPS.h"
#include "nuto/math/SparseMatrixCSRVector2General.h"

#define MAXNUMNEWTONITERATIONS 20
//#define MAXNORMRHS 100
#define PRINTRESULT true
#define MIN_DELTA_STRAIN_FACTOR 1e-7

// there is still an error at the very end of the calculation
// uncomment at the end of void NuTo::NonlocalDamagePlasticity::YieldSurfaceRankine2DRoundedDerivatives
// the check for the derivatives of the yield surface to see the problem

int main()
{
try
{
    //
    bool square(true);

    //create structure
    NuTo::StructureIp myStructureFineScale(2);
    myStructureFineScale.SetShowTime(true);


    NuTo::FullMatrix<int> createdGroupIds;
    if (square)
    {
        //myStructureFineScale.ImportFromGmsh("/home/unger3/develop/nuto/examples/c++/ConcurrentMultiscale.msh","displacements", "ConstitutiveLawIpNonlocal", "StaticDataNonlocal",createdGroupIds);
        //myStructureFineScale.ImportFromGmsh("/home/unger3/develop/nuto/examples/c++/ConcurrentMultiscaleFineScaleCrack.msh","displacements", "ConstitutiveLawIpNonlocal", "StaticDataNonlocal",createdGroupIds);
        myStructureFineScale.ImportFromGmsh("/home/unger3/develop/nuto/examples/c++/ConcurrentMultiscaleFineScaleSquare.msh","displacements", "ConstitutiveLawIpNonlocal", "StaticDataNonlocal",createdGroupIds);
        //myStructureFineScale.ImportFromGmsh("/home/unger3/develop/nuto/examples/c++/ConcurrentMultiscaleFineScale.msh","displacements", "ConstitutiveLawIpNonlocal", "StaticDataNonlocal",createdGroupIds);
    }
    else
    {
        //myStructureFineScale.ImportFromGmsh("/home/unger3/develop/nuto/examples/c++/ConcurrentMultiscaleRoundCrack.msh","displacements", "ConstitutiveLawIpNonlocal", "StaticDataNonlocal",createdGroupIds);
        myStructureFineScale.ImportFromGmsh("/home/unger3/develop/nuto/examples/c++/ConcurrentMultiscaleRound.msh","displacements", "ConstitutiveLawIpNonlocal", "StaticDataNonlocal",createdGroupIds);
        //myStructureFineScale.ImportFromGmsh("/home/unger3/develop/nuto/examples/c++/ConcurrentMultiscaleFineScaleSquare.msh","displacements", "ConstitutiveLawIpNonlocal", "StaticDataNonlocal",createdGroupIds);
    }

    //create constitutive law nonlocal damage
    int myMatDamage = myStructureFineScale.ConstitutiveLawCreate("NonlocalDamagePlasticity");
    double YoungsModulusDamage(20000);
    myStructureFineScale.ConstitutiveLawSetYoungsModulus(myMatDamage,YoungsModulusDamage);
    myStructureFineScale.ConstitutiveLawSetPoissonsRatio(myMatDamage,0.0);
    double nonlocalRadius(15);
    myStructureFineScale.ConstitutiveLawSetNonlocalRadius(myMatDamage,nonlocalRadius);
    double fct(2);
    myStructureFineScale.ConstitutiveLawSetTensileStrength(myMatDamage,fct);
    myStructureFineScale.ConstitutiveLawSetCompressiveStrength(myMatDamage,fct*10);
    myStructureFineScale.ConstitutiveLawSetBiaxialCompressiveStrength(myMatDamage,fct*12.5);
    myStructureFineScale.ConstitutiveLawSetFractureEnergy(myMatDamage,.2);

    //create constitutive law linear elastic (finally not used, since the elements are deleted)
    int myMatLinear = myStructureFineScale.ConstitutiveLawCreate("LinearElastic");
    double YoungsModulusLE(20000);
    myStructureFineScale.ConstitutiveLawSetYoungsModulus(myMatLinear,YoungsModulusLE);
    myStructureFineScale.ConstitutiveLawSetPoissonsRatio(myMatLinear,0.0);

    //create section
    double thickness(1);
    int mySectionParticle = myStructureFineScale.SectionCreate("Plane_Strain");
    myStructureFineScale.SectionSetThickness(mySectionParticle,thickness);

    int mySectionMatrix = myStructureFineScale.SectionCreate("Plane_Strain");
    myStructureFineScale.SectionSetThickness(mySectionMatrix,thickness);

    //assign constitutive law
//    myStructureFineScale.ElementGroupSetSection(101,mySectionParticle);
//    myStructureFineScale.ElementGroupSetSection(102,mySectionMatrix);
//    myStructureFineScale.ElementGroupSetConstitutiveLaw(101,myMatLinear);
//    myStructureFineScale.ElementGroupSetConstitutiveLaw(102,myMatDamage);
//    bool deleteNodes=true;
//    myStructureFineScale.ElementGroupDelete(101,deleteNodes);

    myStructureFineScale.ElementTotalSetSection(mySectionMatrix);
//change here and nonlocal elements afterwards
    myStructureFineScale.ElementTotalSetConstitutiveLaw(myMatDamage);
    //myStructureFineScale.ElementTotalSetConstitutiveLaw(myMatLinear);

    //Build nonlocal elements
    myStructureFineScale.BuildNonlocalData(myMatDamage);

    int GrpNodes_Boundary;
    if (square)
    {
        //Create groups to apply the periodic boundary conditions
        double lX(100);
        double lY(100);
        //left boundary
        int GrpNodes_Left = myStructureFineScale.GroupCreate("Nodes");
        int direction = 0; //either 0,1,2
        double min(0.);
        double max(0.);
        myStructureFineScale.GroupAddNodeCoordinateRange(GrpNodes_Left,direction,min,max);

        //right boundary
        int GrpNodes_Right = myStructureFineScale.GroupCreate("Nodes");
        direction = 0; //either 0,1,2
        min=lX;
        max=lX;
        myStructureFineScale.GroupAddNodeCoordinateRange(GrpNodes_Right,direction,min,max);

        //top boundary
        int GrpNodes_Top = myStructureFineScale.GroupCreate("Nodes");
        direction=1;
        min=lY;
        max=lY;
        myStructureFineScale.GroupAddNodeCoordinateRange(GrpNodes_Top,direction,min,max);

        //bottom boundary
        int GrpNodes_Bottom = myStructureFineScale.GroupCreate("Nodes");
        direction=1;
        min=0;
        max=0;
        myStructureFineScale.GroupAddNodeCoordinateRange(GrpNodes_Bottom,direction,min,max);

        //top right node
        int GrpNodes_BottomTop = myStructureFineScale.GroupUnion(GrpNodes_Bottom,GrpNodes_Top);
        int GrpNodes_LeftRight = myStructureFineScale.GroupUnion(GrpNodes_Left,GrpNodes_Right);
        //GrpNodes_Boundary = myStructureFineScale.GroupUnion(GrpNodes_BottomTop,GrpNodes_LeftRight);

        std::cout << "all nodes are boundary nodes!!!!!!!!!!" << std::endl;
        direction=1;
        min=-1;
        max=101;
        GrpNodes_Boundary = myStructureFineScale.GroupCreate("Nodes");
        myStructureFineScale.GroupAddNodeCoordinateRange(GrpNodes_Boundary,direction,min,max);
    }
    else
    {
        GrpNodes_Boundary = myStructureFineScale.GroupCreate("Nodes");
        NuTo::FullMatrix<double> center(2,1);
        center(0,0) = 0.;
        center(1,0) = 0.;
        double rmin=49.999;
        double rmax=50.001;
        myStructureFineScale.GroupAddNodeRadiusRange(GrpNodes_Boundary,center,rmin,rmax);

    }

    myStructureFineScale.SetGroupBoundaryNodes(GrpNodes_Boundary);
    myStructureFineScale.GroupInfo(10);
	//update conre mat
	myStructureFineScale.NodeBuildGlobalDofs();

	myStructureFineScale.AddVisualizationComponentSection();
	myStructureFineScale.AddVisualizationComponentConstitutive();
	myStructureFineScale.AddVisualizationComponentDisplacements();
	myStructureFineScale.AddVisualizationComponentEngineeringStrain();
	myStructureFineScale.AddVisualizationComponentEngineeringStress();
	myStructureFineScale.AddVisualizationComponentDamage();
	myStructureFineScale.AddVisualizationComponentEngineeringPlasticStrain();
	myStructureFineScale.AddVisualizationComponentPrincipalEngineeringStress();
	myStructureFineScale.Info();

#ifdef ENABLE_SERIALIZATION
    myStructureFineScale.Save("myStructureFineScale.xml","xml");
    myStructureFineScale.Save("myStructureFineScale.bin","binary");
    myStructureFineScale.Restore("myStructureFineScale.bin","binary");
#endif
}
catch (NuTo::Exception& e)
{
    std::cout << e.ErrorMessage() << std::endl;
}
return 0;
}
