// $Id$

#include "nuto/math/MathException.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"

#include "nuto/math/SparseDirectSolverMUMPS.h"


#include "nuto/base/Debug.h"

int main()
{
	try
	{
		// create structure
		NuTo::Structure myStructure(2);
		myStructure.Info();

		// create nodes
		NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> Coordinates(2,11);
		NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> Displacements(2,11);

		Coordinates(0,0) = 0; Coordinates(1,0) = 0;
		Coordinates(0,1) = 1; Coordinates(1,1) = 0;
		Coordinates(0,2) = 2; Coordinates(1,2) = 0;
		Coordinates(0,3) = 0; Coordinates(1,3) = 1;
		Coordinates(0,4) = 1; Coordinates(1,4) = 1;
		Coordinates(0,5) = 2; Coordinates(1,5) = 1;
		Coordinates(0,6) = 0; Coordinates(1,6) = 2;
		Coordinates(0,7) = 1; Coordinates(1,7) = 2;
		Coordinates(0,8) = 2; Coordinates(1,8) = 2;

		Coordinates(0,9) = 4; Coordinates(1,9) = 0;
		Coordinates(0,10) = 4; Coordinates(1,10) = 2;
		DBG_POSITION_INFO("Coordinates Matrix")
		Coordinates.Info();

		NuTo::FullVector<int,Eigen::Dynamic> Nodes = myStructure.NodesCreate("displacements", Coordinates);

		// create elements
		NuTo::FullMatrix<int,Eigen::Dynamic,Eigen::Dynamic> Incidences(4,5);

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

		// element5
		Incidences(0,4) = Nodes(2,0);
		Incidences(1,4) = Nodes(9,0);
		Incidences(2,4) = Nodes(10,0);
		Incidences(3,4) = Nodes(8,0);

		DBG_POSITION_INFO("Incidence Matrix")
		Coordinates.Info();

		NuTo::FullVector<int,Eigen::Dynamic> Elements = myStructure.ElementsCreate("Plane2D4N", Incidences);

	    // create constitutive law
	    int myMatLin = myStructure.ConstitutiveLawCreate("LINEARELASTICENGINEERINGSTRESS");
	    myStructure.ConstitutiveLawSetYoungsModulus(myMatLin,10);
	    myStructure.ConstitutiveLawSetPoissonsRatio(myMatLin,0.1);

	    // create section
	    int mySection1 = myStructure.SectionCreate("PLANE_STRAIN");
	    myStructure.SectionSetThickness(mySection1,0.01);

	    // assign material, section and integration type
	    myStructure.ElementTotalSetIntegrationType("2D4NGauss4Ip","NOIPDATA");
	    myStructure.ElementTotalSetConstitutiveLaw(myMatLin);
	    myStructure.ElementTotalSetSection(mySection1);


		// build global stiffness matrix and equivalent load vector which correspond to prescribed boundary values
		NuTo::FullVector<int,Eigen::Dynamic> rows;
		NuTo::FullVector<int,Eigen::Dynamic> coluums;
		NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> stiffnessMatrix;
		myStructure.ElementStiffness(0,stiffnessMatrix,rows,coluums );
//		if(mVerboseLevel>3)
		stiffnessMatrix.WriteToFile("stiffness","   ");
		myStructure.ElementStiffness(4,stiffnessMatrix,rows,coluums );
//		if(mVerboseLevel>3)
		stiffnessMatrix.WriteToFile("stiffnessCoarse","   ");


	   // boundary conditions
		NuTo::FullVector<double,Eigen::Dynamic> direction(2);
		direction(0)= 1;
		direction(1)= 0;
		myStructure.ConstraintLinearSetDisplacementNode(0, direction, 0.0);
		myStructure.ConstraintLinearSetDisplacementNode(3, direction, 0.0);
		myStructure.ConstraintLinearSetDisplacementNode(6, direction, 0.0);
		direction(0)= 0;
		direction(1)= 1;
		for(int i=0;i<11;++i)
		{
			myStructure.ConstraintLinearSetDisplacementNode(i, direction, 0.0);

		}
//		myStructure.ConstraintLinearSetDisplacementNode(0, direction, 0.0);
//		myStructure.ConstraintLinearSetDisplacementNode(1, direction, 0.0);
//		myStructure.ConstraintLinearSetDisplacementNode(2, direction, 0.0);
//		myStructure.ConstraintLinearSetDisplacementNode(9, direction, 0.0);

		std::cout << "Displacement control" << std::endl;
		// boundary displacments
		double BoundaryDisplacement=1;
		direction(0)= 1;
		direction(1)= 0;
		myStructure.ConstraintLinearSetDisplacementNode(9, direction, BoundaryDisplacement);
		myStructure.ConstraintLinearSetDisplacementNode(10, direction, BoundaryDisplacement);

	   // start analysis
		// build global dof numbering
		myStructure.NodeBuildGlobalDofs();

		// build global stiffness matrix and equivalent load vector which correspond to prescribed boundary values
		NuTo::SparseMatrixCSRVector2General<double> stiffnessMatrixVector2;
		NuTo::FullVector<double,Eigen::Dynamic> dispForceVector;
		myStructure.CalculateMaximumIndependentSets();
		myStructure.BuildGlobalCoefficientMatrix0(stiffnessMatrixVector2, dispForceVector);
		stiffnessMatrixVector2.RemoveZeroEntries(0,1e-14);
	        //NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> A(stiffnessMatrix);
	        //A.WriteToFile("stiffnessMatrix.txt"," ");
	        //stiffnessMatrix.Info();
	        //dispForceVector.Info();

		// build global external load vector
		NuTo::FullVector<double,Eigen::Dynamic> extForceVector;
		myStructure.BuildGlobalExternalLoadVector(0,extForceVector);
		//extForceVector.Info();

		// calculate right hand side
		NuTo::FullVector<double,Eigen::Dynamic> rhsVector = dispForceVector + extForceVector;
		rhsVector.WriteToFile("rhsVector.txt"," ");

	        // solve
		NuTo::SparseDirectSolverMUMPS mySolver;
		NuTo::FullVector<double,Eigen::Dynamic> displacementVector;
		NuTo::SparseMatrixCSRGeneral<double> stiffnessMatrixGlob(stiffnessMatrixVector2);
		stiffnessMatrixGlob.SetOneBasedIndexing();
		mySolver.Solve(stiffnessMatrixGlob, rhsVector, displacementVector);
		displacementVector.WriteToFile("displacementVector.txt"," ");

		// write displacements to node
		myStructure.NodeMergeActiveDofValues(displacementVector);

		// calculate residual
		NuTo::FullVector<double,Eigen::Dynamic> intForceVector;
		myStructure.BuildGlobalGradientInternalPotentialVector(intForceVector);
		NuTo::FullVector<double,Eigen::Dynamic> residualVector = extForceVector - intForceVector;
		std::cout << "residual: " << residualVector.Norm() << std::endl;


#ifdef ENABLE_VISUALIZE
            // visualize element
        myStructure.AddVisualizationComponentDisplacements();
        myStructure.AddVisualizationComponentEngineeringStrain();
        myStructure.AddVisualizationComponentEngineeringStress();
        myStructure.ExportVtkDataFileElements("Plane2D4N.vtk");
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
