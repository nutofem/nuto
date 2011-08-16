// $Id$

#include "nuto/math/MathException.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/math/SparseDirectSolverMUMPS.h"
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"
#include "nuto/base/Debug.h"

int main()
{
	// definitions
	const double Width = 5.;
	const double Height = 10.;
	const size_t  NumElementsX = 13;
	const size_t  NumElementsY = 25;

	try
	{
		// create structure
		NuTo::Structure myStructure(2);
		myStructure.Info();
#ifdef SHOW_TIME
		myStructure.SetShowTime(true);
#endif //SHOW_TIME
		// create nodes
		NuTo::FullMatrix<double> nodeCoordinates(2,1);
		size_t node = 0;
		for(size_t yCount = 0; yCount < NumElementsY + 1; yCount++)
		{
			nodeCoordinates(1,0) = (double)yCount * Height/(double)NumElementsY;
			for(size_t xCount = 0; xCount < NumElementsX + 1; xCount++)
			{
				nodeCoordinates(0,0) = (double)xCount * Width/(double)NumElementsX;
				//std::cout << "node: " << node << " coordinates: " << nodeCoordinates.GetValue(0,0) << "," << nodeCoordinates.GetValue(1,0) << "," << nodeCoordinates.GetValue(2,0) << std::endl;
				myStructure.NodeCreate(node, "displacements", nodeCoordinates);
				node ++;
			}
		}

		// create elements
		NuTo::FullMatrix<int> elementIncidence(4,1);
		size_t element = 0;
		for(size_t yCount = 0; yCount < NumElementsY; yCount++)
		{
			for(size_t xCount = 0; xCount < NumElementsX; xCount++)
			{
				size_t node1 = yCount * (NumElementsX + 1) + xCount;
				elementIncidence(0,0) = node1;
				elementIncidence(1,0) = node1 + 1;
				elementIncidence(2,0) = node1 + NumElementsX + 2;
				elementIncidence(3,0) = node1 + NumElementsX + 1;
				//std::cout << "element: " << element << " incidence: " << std::endl;
				//elementIncidence.Info();
				myStructure.ElementCreate(element, "plane2d4n", elementIncidence, "CONSTITUTIVELAWIPCRACK" , "NOIPDATA");
				element ++;
			}
		}
		//~ NuTo::FullMatrix<double> InitCrackCoords(2,4, std::vector< double >{
																 //~ 1.6, 5 ,
																 //~ 2.5, 5 ,
																 //~ 2.6, 7.5 ,
																 //~ 3.4, 7.5 });
		//~ NuTo::FullMatrix<double> InitCrackCoords(2,2, std::vector< double >{
																 //~ 1.6, 5 ,
																 //~ 3.4, 5 });
		NuTo::FullMatrix<double> InitCrackCoords(2,2, std::vector< double >{
																 1.6, 5 ,
																 6, 5 });
		DBG_POSITION_INFO("InitCrackCoords Matrix")
		InitCrackCoords.Info();

		NuTo::FullMatrix<int> CrackNodes = myStructure.NodesCreate("coordinates", InitCrackCoords);
		DBG_POSITION_INFO("CrackNodes")
		CrackNodes.Info(5);

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

		myStructure.Info();
	    // merge crack into elements
	    /*
	    NuTo::Structure::elementBasePtrSet_t crackedElems;
	    myStructure.InitiateCrack(crack1, crackedElems);
	    
	    BOOST_FOREACH(NuTo::Structure::elementBasePtr_t thisCrackedElem, crackedElems)
	    	    myStructure.ElementSetIntegrationType(thisCrackedElem->ElementGetId(),"2D4Nconst100Ip","NOIPDATA");

	    myStructure.InitiatePhantomNodeMethod(crackedElems);
		*/
		
		myStructure.InitiatePhantomNodeMethod(1000);

	    //~ myStructure.NodeInfo(5);
	    //~ myStructure.ElementInfo(5);
//~ 
	    //~ myStructure.IntegrationTypeInfo(1);

	    myStructure.ElementTotalSetConstitutiveLaw(myMatLin);
	    myStructure.ElementTotalSetSection(mySection1);

#ifdef ENABLE_VISUALIZE
        // visualize element
        myStructure.AddVisualizationComponentDisplacements();
        myStructure.AddVisualizationComponentElement();
        myStructure.AddVisualizationComponentSection();
        myStructure.AddVisualizationComponentCracks();
        myStructure.AddVisualizationComponentEngineeringStrain();
        myStructure.AddVisualizationComponentEngineeringStress();
        myStructure.ExportVtkDataFile("PhantomNodeMethodPlane2D4N.vtk");
#endif

        // boundary conditions
        // x-direction
        NuTo::FullMatrix<double> direction(2,1);
        direction(0,0)= 1;
        direction(1,0)= 0;
        myStructure.ConstraintLinearSetDisplacementNode(0, direction, 0.0);
        // y-direction bottom
        direction(0,0)= 0;
        direction(1,0)= 1;
		for(size_t xCount = 0; xCount < NumElementsX + 1; xCount++)
		{
			size_t node = xCount;
			//std::cout << "node: " << node << std::endl;
			myStructure.ConstraintLinearSetDisplacementNode(node, direction, 0.0);
		}
        // y-direction top
        direction(0,0)= 0;
        direction(1,0)= 1;
		for(size_t xCount = 0; xCount < NumElementsX + 1; xCount++)
		{
			size_t node = xCount + (NumElementsX + 1) * NumElementsY;
			//std::cout << "node: " << node << std::endl;
			myStructure.ConstraintLinearSetDisplacementNode(node, direction, 0.5);
		}

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

#ifdef ENABLE_VISUALIZE
        // visualize results
        myStructure.ExportVtkDataFile("PhantomNodeMethodPlane2D4N-deformed.vtk");
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
