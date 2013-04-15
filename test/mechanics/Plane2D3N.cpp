#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/MechanicsException.h"

#include "nuto/math/SparseDirectSolverMUMPS.h"
#include "nuto/math/SparseMatrixCSRVector2General.h"

#define MAXNUMNEWTONITERATIONS 20
#define PRINTRESULT true
#define MIN_DELTA_STRAIN_FACTOR 1e-7

int main()
{
try
{
    //create structure
    NuTo::Structure myStructure(2);

    //create nodes
    myStructure.NodesCreate("displacements", nuto.DoubleFullMatrix(2,8,(	 
    0 ,  0 ,
    10 ,  0 ,
    2 ,  2 ,
    8 ,  3 ,
    4 ,  7 ,
    8 ,  7 ,
    0 , 10 ,
    10 , 10	)));

    //create element
    elementIncidence = nuto.IntFullMatrix(3,10,(	
    0,1,3,
    0,2,6,
    0,3,2,
    1,7,3,
    2,4,6,
    2,3,4,
    3,5,4,
    3,7,5,
    4,5,6,
    5,7,6 ));
    myStructure.ElementsCreate("PLANE2D3N", elementIncidence);

    //Calculate maximum independent sets for parallelization (openmp)
    myStructure.CalculateMaximumIndependentSets();

    //create constitutive law
    int myMatLin = myStructure.ConstitutiveLawCreate("LinearElasticEngineeringStress");
    myStructure.ConstitutiveLawSetYoungsModulus(myMatLin,10);
    myStructure.ConstitutiveLawSetPoissonsRatio(myMatLin,0.2);

    //create section
    int mySection = myStructure.SectionCreate("Plane_Strain");
    myStructure.SectionSetThickness(mySection,1);

    //assign constitutive law 
    myStructure.ElementTotalSetConstitutiveLaw(myMatLin);
    myStructure.ElementTotalSetSection(mySection);

    //set displacements of right node
    myStructure.ConstraintLinearSetDisplacementNode(0, nuto.DoubleFullMatrix(2,1,(1,0)), 0.0);
    myStructure.ConstraintLinearSetDisplacementNode(0, nuto.DoubleFullMatrix(2,1,(0,1)), 0.0);
    myStructure.ConstraintLinearSetDisplacementNode(6, nuto.DoubleFullMatrix(2,1,(1,0)), 0.0);
    myStructure.ConstraintLinearSetDisplacementNode(1, nuto.DoubleFullMatrix(2,1,(1,0)), 1.0);
    myStructure.ConstraintLinearSetDisplacementNode(7, nuto.DoubleFullMatrix(2,1,(1,0)), 1.0);

    myStructure.SetVerboseLevel(10);
    myStructure.Info();
    //start analysis
    //build global dof numbering
    myStructure.NodeBuildGlobalDofs();

    // build global stiffness matrix and equivalent load vector which correspond to prescribed boundary values
    NuTo::SparseMatrixCSRVector2General<double> stiffnessMatrixCSRVector2(0,0);
    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> dispForceVector;
    myStructure.BuildGlobalCoefficientMatrix0(stiffnessMatrixCSRVector2, dispForceVector);
    NuTo::SparseMatrixCSRGeneral<double> stiffnessMatrix (stiffnessMatrixCSRVector2);

    // build global external load vector
    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> extForceVector;
    myStructure.BuildGlobalExternalLoadVector(extForceVector);

    // calculate right hand side
    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> rhsVector = dispForceVector + extForceVector;

    // solve
    NuTo::SparseDirectSolverMUMPS mySolver;
    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> displacementVector;
    stiffnessMatrix.SetOneBasedIndexing();
    mySolver.Solve(stiffnessMatrix, rhsVector, displacementVector);

    // write displacements to node
    myStructure.NodeMergeActiveDofValues(displacementVector);
    
    // calculate residual
    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> intForceVector = nuto.DoubleFullMatrix();
    myStructure.BuildGlobalGradientInternalPotentialVector(intForceVector);
    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> residualVector (extForceVector - intForceVector);
    if ((residualVector).Abs().Max()[0]>1e-8)
    {
	std::cout << "[Plane2D3N] : residual force vector is not zero." << std::endl;
	error = true;
    }

    //calculate engineering strain of myelement1 at all integration points
    //the size the matrix is not important and reallocated within the procedure
    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> EngineeringStrain(6,1);
    //correct strain
    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> EngineeringStrainCorrect;
    EngineeringStrainCorrect(0,0) = 0.1;
    EngineeringStrainCorrect(1,0) = -0.025;

    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> EngineeringStress(6,1);
    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> EngineeringStressCorrect(6,1);
    EngineeringStressCorrect(0,0) = 1.0416666666666667;
    EngineeringStressCorrect(2,0) = 0.2083333333333333;
 
    for (int element=0; element<8; element++)
    {
	myStructure.ElementGetEngineeringStrain(element, EngineeringStrain);

	if (printResult)
	{
	    std::cout << "EngineeringStrainCorrect" << std::endl;
	    EngineeringStrainCorrect.Info()
	    std::cout << "EngineeringStrain" << std::end;
	    EngineeringStrain.Info()
	}

	if ((EngineeringStrain-EngineeringStrainCorrect).Abs().Max()[0]>1e-8)
	{
	    if (!printResult)
	    {
		std::cout << "EngineeringStrainCorrect" << std::endl;
		EngineeringStrainCorrect.Info()
		std::cout << "EngineeringStrain" << std::end;
		EngineeringStrain.Info()
	    }
	    std::cout << "[Plane2D3N] : strain is not correct." << std::endl;
	    error = true;
	}

	//calculate engineering strain at all integration points
	myStructure.ElementGetEngineeringStress(element, EngineeringStress)
	//correct stress

	if (printResult)
	{
	    std::cout << "EngineeringStressCorrect" << std::endl;
	    EngineeringStressCorrect.Info()
	    std::cout << "EngineeringStress" << std::end;
	    EngineeringStress.Info()
	}

	if ((EngineeringStress-EngineeringStressCorrect).Abs().Max()[0]>1e-8)
	{
	    if (!printResult)
	    {
		std::cout << "EngineeringStressCorrect" << std::endl;
		EngineeringStressCorrect.Info()
		std::cout << "EngineeringStress" << std::end;
		EngineeringStress.Info()
	    }
	    std::cout << "[Plane2D3N] : stress is not correct." << std::endl;
	    error = true;
	}


    // visualize results
    myStructure.AddVisualizationComponentDisplacements();
    myStructure.AddVisualizationComponentEngineeringStrain();
    myStructure.AddVisualizationComponentEngineeringStress();
    myStructure.ExportVtkDataFile( "Plane2D3N.vtk");
}
catch (NuTo::Exception& e)
{
    std::cout << e.ErrorMessage() << std::endl;
    return -1;
}
return 0;
}
