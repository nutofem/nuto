// $Id$

#include "nuto/math/MathException.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/math/FullVector.h"
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"

#include "nuto/base/Debug.h"

int main()
{
	try
	{
		// create structure
		NuTo::Structure myStructure(2);
		myStructure.Info();

		// create nodes
		NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> Coordinates(2,9);
		NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> Displacements(2,9);

		Coordinates(0,0) = 0; Coordinates(1,0) = 0;
		Coordinates(0,1) = 1; Coordinates(1,1) = 0;
		Coordinates(0,2) = 2; Coordinates(1,2) = 0;
		Coordinates(0,3) = 0; Coordinates(1,3) = 1;
		Coordinates(0,4) = 1; Coordinates(1,4) = 1;
		Coordinates(0,5) = 2; Coordinates(1,5) = 1;
		Coordinates(0,6) = 0; Coordinates(1,6) = 2;
		Coordinates(0,7) = 1; Coordinates(1,7) = 2;
		Coordinates(0,8) = 2; Coordinates(1,8) = 2;

		DBG_POSITION_INFO("Coordinates Matrix")
		Coordinates.Info();

		NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> Coordinates2(2,9, std::vector< double >{0,0 , 1,0 , 2,0 , 0,1 , 1,1 , 2,1 , 0,2 , 1,2 , 2,2});
		DBG_POSITION_INFO("Coordinates2 Matrix")
		Coordinates2.Info();
		
		NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> DiffCoords = Coordinates - Coordinates2;
		DBG_POSITION_INFO("DiffCoords Matrix")
		DiffCoords.Info();
		
		NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> InitCrackCoords1(2,3, std::vector< double >{0.0,1.5, 1,1.75 , 2,0.5 });
		DBG_POSITION_INFO("InitCrackCoords1 Matrix")
		InitCrackCoords1.Info();
		NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> InitCrackCoords5(2,4, std::vector< double >{-1.0,-1.0, 0.5,1 , 1,1.5 , 1.5,0});
		DBG_POSITION_INFO("InitCrackCoords5 Matrix")
		InitCrackCoords5.Info();
		
		NuTo::FullMatrix<int,Eigen::Dynamic,Eigen::Dynamic> Nodes = myStructure.NodesCreate("displacements", Coordinates);
		
		NuTo::FullVector<int,Eigen::Dynamic> CrackNodes1 = myStructure.NodesCreate("coordinates", InitCrackCoords1);
		NuTo::FullVector<int,Eigen::Dynamic> CrackNodes5 = myStructure.NodesCreate("coordinates", InitCrackCoords5);
		DBG_POSITION_INFO("CrackNodes1")
		CrackNodes1.Info(5);
		DBG_POSITION_INFO("CrackNodes5")
		CrackNodes5.Info(5);

		// create elements
		NuTo::FullMatrix<int,Eigen::Dynamic,Eigen::Dynamic> Incidences(4,4);

		// element1
		Incidences(0,0) = Nodes(0,0);
		Incidences(1,0) = Nodes(1,0);
		Incidences(2,0) = Nodes(4,0);
		Incidences(3,0) = Nodes(3,0);
		// element2
		Incidences(0,1) = Nodes(1,0);
		Incidences(1,1) = Nodes(2,0);
		Incidences(2,1) = Nodes(5,0);
		Incidences(3,1) = Nodes(4,0);
		// element3
		Incidences(0,2) = Nodes(3,0);
		Incidences(1,2) = Nodes(4,0);
		Incidences(2,2) = Nodes(7,0);
		Incidences(3,2) = Nodes(6,0);
		// element4
		Incidences(0,3) = Nodes(4,0);
		Incidences(1,3) = Nodes(5,0);
		Incidences(2,3) = Nodes(8,0);
		Incidences(3,3) = Nodes(7,0);

		DBG_POSITION_INFO("Incidence Matrix")
		Incidences.Info();

		NuTo::FullVector<int,Eigen::Dynamic> Elements = myStructure.ElementsCreate("Plane2D4N", Incidences, "CONSTITUTIVELAWIPCRACK" , "NOIPDATA");

	    // create constitutive law
	    int myMatLin = myStructure.ConstitutiveLawCreate("LinearElastic");
	    myStructure.ConstitutiveLawSetYoungsModulus(myMatLin,10);
	    myStructure.ConstitutiveLawSetPoissonsRatio(myMatLin,0.1);

	    // create section
	    //~ int mySection1 = myStructure.SectionCreate("PLANE_STRAIN");
	    int mySection1 = myStructure.SectionCreate("PLANE_STRESS");
	    myStructure.SectionSetThickness(mySection1,0.01);

	    // assign material, section and integration type
DBG_POSITION_INFO("before")
	    myStructure.ElementTotalSetIntegrationType("2D4NGauss4Ip","NOIPDATA");
DBG_POSITION_INFO("after")
	    myStructure.ElementTotalSetConstitutiveLaw(myMatLin);
	    myStructure.ElementTotalSetSection(mySection1);

	    // Crack initiation
	    int crack1= myStructure.CrackCreate();
	    int crack2= myStructure.CrackCreate();
	    int crack3= myStructure.CrackCreate();
	    int crack4= myStructure.CrackCreate();
DBG_PRINT_VAL(crack1)
DBG_PRINT_VAL(crack2)
DBG_PRINT_VAL(crack3)
DBG_PRINT_VAL(crack4)
	    myStructure.CrackInfo(5);
	    int crack5= myStructure.CrackCreate();
DBG_PRINT_VAL(crack5)
	    myStructure.CrackInfo(5);
	    for (int i=0; i < CrackNodes1.GetNumRows(); i++ )
	    {
			myStructure.CrackPushBack(crack1, myStructure.NodeGetNodePtr(CrackNodes1(i,0)));
		}
	    for (int i=0; i < CrackNodes5.GetNumRows(); i++ )
	    {
			myStructure.CrackPushBack(crack5, myStructure.NodeGetNodePtr(CrackNodes5(i,0)));
		}
		//myStructure.CrackPushFront(crack5, myStructure.NodeGetNodePtr(CrackNodes(InitCrackCoords.GetNumColumns()-1,0)));
	    myStructure.CrackDelete(crack2);
	    myStructure.CrackDelete(crack3);
	    myStructure.CrackDelete(crack4);
	    myStructure.CrackInfo(5);

	    // merge crack into elements
	    myStructure.InitiateCracks();


        // visualize element
#ifdef ENABLE_VISUALIZE
        myStructure.AddVisualizationComponentDisplacements();
        myStructure.AddVisualizationComponentElement();
        myStructure.AddVisualizationComponentSection();
        myStructure.AddVisualizationComponentCracks();
        myStructure.AddVisualizationComponentEngineeringStrain();
        myStructure.AddVisualizationComponentEngineeringStress();
        myStructure.ExportVtkDataFile("CrackPlane2D4N.vtk");
#endif
	}
	catch (NuTo::MathException& e)
	{
		std::cout << e.ErrorMessage() << std::endl;
	}
	catch (NuTo::Exception& e)
	{
		std::cout << e.ErrorMessage() << std::endl;
	}
	catch (...)
	{
		std::cout << "Unexpected" << std::endl;
	}

    return 0;
}
