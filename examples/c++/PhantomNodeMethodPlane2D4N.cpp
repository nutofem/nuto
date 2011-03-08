// $Id$

#include "nuto/math/MathException.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/math/SparseDirectSolverMUMPS.h"
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"

#include "nuto/base/Debug.h"

int main()
{
//	int numIp=16;
//	std::cout << (int)sqrt((double)numIp) << std::endl;
//	numIp=((int)sqrt((double)numIp))*((int)sqrt((double)numIp));
//	std::cout << numIp << std::endl;
//	for(int i=0; i<numIp; ++i)
//	{
//		if(i%(int)sqrt((double)numIp)==0) std::cout << std::endl << i/(int)sqrt((double)numIp) << ":";
//		std::cout << "\t" << i%(int)sqrt((double)numIp);
//	}
//	std::cout << std::endl;
//
//	for(int i=0; i<numIp; ++i)
//	{
//		std::cout << "coordinates: [" 	<< (0.5+i%(int)sqrt((double)numIp))/(sqrt((double)numIp))*2-1 <<
//									";" << (0.5+i/(int)sqrt((double)numIp))/(sqrt((double)numIp))*2-1 << "]" << std::endl;
//	}
//
//	unsigned int mNumIp1D=(unsigned int)sqrt((double)numIp);
//	for(int i=0; i<(mNumIp1D+1)*(mNumIp1D+1); ++i)
//	{
//		std::cout << "VisulizeCoordinates: [" 	<< (i%(mNumIp1D+1))/((double)mNumIp1D)*2-1 <<
//											";" << (i/(mNumIp1D+1))/((double)mNumIp1D)*2-1 << "]" << std::endl;
//	}

	try
	{
		// create structure
		NuTo::Structure myStructure(2);
		myStructure.Info();

		// create nodes
		NuTo::FullMatrix<double> Displacements(2,12);

		NuTo::FullMatrix<double> Coordinates(2,12, std::vector< double >{0,0 ,
																		1,0 ,
																		2,0 ,
																		0,1 ,
																		1,1 ,
																		2,1 ,
																		0,2 ,
																		1,2 ,
																		2,2 ,
																		0,3 ,
																		1,3 ,
																		2,3});
/*		NuTo::FullMatrix<double> Coordinates(2,12, std::vector< double >{0,0 ,
																		1,0 ,
																		2,0 ,
																		0,1 ,
																		1.9,1 ,
																		2,1 ,
																		0,2 ,
																		1,2 ,
																		2,2 ,
																		0,3 ,
																		1,3 ,
																		2,3});
*/		DBG_POSITION_INFO("Coordinates Matrix")
		Coordinates.Info();

/*		NuTo::FullMatrix<double> InitCrackCoords(2,4, std::vector< double >{
																		-1.0,-1.0 ,
																		 0.5, 1.0 ,
																		 1.0, 1.5 ,
																		 1.5, 0.0});*/
//NuTo::FullMatrix<double> InitCrackCoords(2,3, std::vector< double >{
//																-1.0, 2.0 ,
//																-1.0, 1.5 ,
//																 1.5, 1.5 });
NuTo::FullMatrix<double> InitCrackCoords(2,2, std::vector< double >{
																-1.0, 1.5 ,
																 0.5, 1.5 });
//NuTo::FullMatrix<double> InitCrackCoords(2,2, std::vector< double >{
//																-1.0, 1.5 ,
//																 0.5, 1.6 });
//NuTo::FullMatrix<double> InitCrackCoords(2,2, std::vector< double >{
//																 1.5, 0.5 ,
//																 2.5, 0.5 });
		DBG_POSITION_INFO("InitCrackCoords Matrix")
		InitCrackCoords.Info();

		NuTo::FullMatrix<int> Nodes = myStructure.NodesCreate("displacements", Coordinates);

		NuTo::FullMatrix<int> CrackNodes = myStructure.NodesCreate("coordinates", InitCrackCoords);
		DBG_POSITION_INFO("CrackNodes")
		CrackNodes.Info(5);

		// create elements
		NuTo::FullMatrix<int> Incidences(4,6, std::vector< int >{
			Nodes(0,0), Nodes(1,0), Nodes(4,0) , Nodes(3,0) ,
			Nodes(1,0), Nodes(2,0), Nodes(5,0) , Nodes(4,0) ,
			Nodes(3,0), Nodes(4,0), Nodes(7,0) , Nodes(6,0) ,
			Nodes(4,0), Nodes(5,0), Nodes(8,0) , Nodes(7,0) ,
			Nodes(6,0), Nodes(7,0), Nodes(10,0), Nodes(9,0) ,
			Nodes(7,0), Nodes(8,0), Nodes(11,0), Nodes(10,0)});
		DBG_POSITION_INFO("Incidence Matrix")
		Incidences.Info();
	    NuTo::FullMatrix<int> Elements = myStructure.ElementsCreate("Plane2D4N", Incidences, "CONSTITUTIVELAWIPCRACK" , "NOIPDATA");

	    // create constitutive law
	    int myMatLin = myStructure.ConstitutiveLawCreate("LinearElastic");
	    myStructure.ConstitutiveLawSetYoungsModulus(myMatLin,10);
	    myStructure.ConstitutiveLawSetPoissonsRatio(myMatLin,0.1);

	    // create section
	    int mySection1 = myStructure.SectionCreate("PLANE_STRAIN");
	    myStructure.SectionSetThickness(mySection1,0.01);

	    // assign material, section and integration type
//	    myStructure.ElementTotalSetIntegrationType("2D4Nmod10Ip","NOIPDATA");
//	    myStructure.ElementTotalSetIntegrationType("2D4Nconst10Ip","NOIPDATA");
	    myStructure.ElementTotalSetIntegrationType("2D4NGauss4Ip","NOIPDATA");
	    // Crack initiation
	    int crack1= myStructure.CrackCreate();
	    for (int i=0; i < CrackNodes.GetNumRows(); i++ )
			myStructure.CrackPushBack(crack1, myStructure.NodeGetNodePtr(CrackNodes(i,0)));
	    myStructure.CrackInfo(5);

	    // merge crack into elements
	    //! @todo [DA] muss intern ablaufen: sowas wie InitiatePhantomNode oder so
	    NuTo::Structure::elementBasePtrSet_t crackedElems;
	    myStructure.InitiateCrack(crack1, crackedElems);
	    BOOST_FOREACH(NuTo::Structure::elementBasePtr_t thisCrackedElem, crackedElems)
	    {
	    	double globCoors[]={10.0,1.0};
	    	double locCoords[2];

	    	if(thisCrackedElem->GetLocalPointCoordinates(globCoors , locCoords))
	    	{
	    		DBG_PRINT_VEC(locCoords)
	    	}else{
	    		DBG_POSITION_INFO("point is not in Element")
	    	}

	    	//! @todo intersect element's corners with the crack

	    }

	    BOOST_FOREACH(NuTo::Structure::elementBasePtr_t thisCrackedElem, crackedElems)
	    	    myStructure.ElementSetIntegrationType(thisCrackedElem->ElementGetId(),"2D4Nconst100Ip","NOIPDATA");

	    myStructure.InitiatePhantomNodeMethod(crackedElems);

	    myStructure.NodeInfo(5);
	    myStructure.ElementInfo(5);

	    myStructure.IntegrationTypeInfo(1);

	    myStructure.ElementTotalSetConstitutiveLaw(myMatLin);
	    myStructure.ElementTotalSetSection(mySection1);

        // visualize element
        myStructure.AddVisualizationComponentDisplacements();
        myStructure.AddVisualizationComponentElement();
        myStructure.AddVisualizationComponentSection();
        myStructure.AddVisualizationComponentCracks();
        myStructure.AddVisualizationComponentEngineeringStrain();
        myStructure.AddVisualizationComponentEngineeringStress();
        myStructure.ExportVtkDataFile("PhantomNodeMethodPlane2D4N.vtk");

        // boundary conditions
        // x-direction
        NuTo::FullMatrix<double> direction(2,1);
        direction(0,0)= 1;
        direction(1,0)= 0;
        myStructure.ConstraintLinearSetDisplacementNode(0, direction, 0.0);
        // y-direction bottom
        direction(0,0)= 0;
        direction(1,0)= 1;
        myStructure.ConstraintLinearSetDisplacementNode(0, direction, 0.0);
        myStructure.ConstraintLinearSetDisplacementNode(1, direction, 0.0);
        myStructure.ConstraintLinearSetDisplacementNode(2, direction, 0.0);
        // y-direction top
        direction(0,0)= 0;
        direction(1,0)= 1;
        myStructure.ConstraintLinearSetDisplacementNode( 9, direction, 0.5);
        myStructure.ConstraintLinearSetDisplacementNode(10, direction, 0.5);
        myStructure.ConstraintLinearSetDisplacementNode(11, direction, 0.5);

        // start analysis
        // build global dof numbering
        myStructure.NodeBuildGlobalDofs();
        // build global stiffness matrix and equivalent load vector which correspond to prescribed boundary values
        NuTo::SparseMatrixCSRGeneral<double> stiffnessMatrix;
        NuTo::FullMatrix<double> dispForceVector;
        myStructure.BuildGlobalCoefficientMatrix0(stiffnessMatrix, dispForceVector);
        stiffnessMatrix.RemoveZeroEntries(0,1e-14);
        
        //~ NuTo::FullMatrix<double> stiffnessFullMatrix;
        //~ stiffnessMatrix.WriteEntriesToFullMatrix(stiffnessFullMatrix);
        //~ stiffnessFullMatrix.WriteToFile("stiffnessMatrix.txt"," ");

        // build global external load vector
        NuTo::FullMatrix<double> extForceVector;
        myStructure.BuildGlobalExternalLoadVector(extForceVector);
        //extForceVector.Info();

        // calculate right hand side
        NuTo::FullMatrix<double> rhsVector = dispForceVector + extForceVector;
        rhsVector.WriteToFile("rhsVector.txt"," ");

        // solve
        NuTo::SparseDirectSolverMUMPS mySolver;
        NuTo::FullMatrix<double> displacementVector;
        stiffnessMatrix.SetOneBasedIndexing();
        mySolver.Solve(stiffnessMatrix, rhsVector, displacementVector);
        displacementVector.WriteToFile("displacementVector.txt"," ");

        // write displacements to node
        myStructure.NodeMergeActiveDofValues(displacementVector);

        // calculate residual
        NuTo::FullMatrix<double> intForceVector;
        myStructure.BuildGlobalGradientInternalPotentialVector(intForceVector);
        NuTo::FullMatrix<double> residualVector = extForceVector - intForceVector;
        std::cout << "residual: " << residualVector.Norm() << std::endl;

        // visualize results
        myStructure.ExportVtkDataFile("PhantomNodeMethodPlane2D4N-deformed.vtk");
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
