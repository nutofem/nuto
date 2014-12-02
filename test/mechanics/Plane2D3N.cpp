#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/MechanicsException.h"

#include "nuto/math/SparseDirectSolverMUMPS.h"
#include "nuto/math/SparseMatrixCSRVector2General.h"

#define PRINTRESULT false

int main()
{
    try
    {
        //create structure
        NuTo::Structure myStructure(2);

        //create nodes
        NuTo::FullMatrix<double, Eigen::Dynamic,Eigen::Dynamic> nodeCoordinates(2,8);
        nodeCoordinates <<
                0, 10, 2, 8, 4, 8,  0, 10,
                0,  0, 2, 3, 7, 7, 10, 10;

        myStructure.NodesCreate("displacements",nodeCoordinates);

        //create element
        NuTo::FullMatrix<int, Eigen::Dynamic, Eigen::Dynamic> nodeNumbers(3,10);
        nodeNumbers <<
                0, 0, 0, 1, 2, 2, 3, 3, 4, 5,
                1, 2, 3, 7, 4, 3, 5, 7, 5, 7,
                3, 6, 2, 3, 6, 4, 4, 5, 6, 6;

        myStructure.ElementsCreate("PLANE2D3N", nodeNumbers);

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

        myStructure.NodeBuildGlobalDofs();
        myStructure.CheckStiffness();

        //set displacements of right node
        NuTo::FullVector<double, 2> dirX({1,0});
        NuTo::FullVector<double, 2> dirY({0,1});

        myStructure.ConstraintLinearSetDisplacementNode(0, dirX, 0.0);
        myStructure.ConstraintLinearSetDisplacementNode(0, dirY, 0.0);
        myStructure.ConstraintLinearSetDisplacementNode(6, dirX, 0.0);
        myStructure.ConstraintLinearSetDisplacementNode(1, dirX, 1.0);
        myStructure.ConstraintLinearSetDisplacementNode(7, dirX, 1.0);

        myStructure.SetVerboseLevel(10);
        myStructure.Info();
        //start analysis
        //build global dof numbering
        myStructure.NodeBuildGlobalDofs();

        // build global stiffness matrix and equivalent load vector which correspond to prescribed boundary values
        NuTo::SparseMatrixCSRVector2General<double> stiffnessMatrixCSRVector2(0,0);
        NuTo::FullVector<double,Eigen::Dynamic> dispForceVector;
        myStructure.BuildGlobalCoefficientMatrix0(stiffnessMatrixCSRVector2, dispForceVector);
        NuTo::SparseMatrixCSRGeneral<double> stiffnessMatrix (stiffnessMatrixCSRVector2);

        // build global external load vector
        NuTo::FullVector<double,Eigen::Dynamic> extForceVector;
        myStructure.BuildGlobalExternalLoadVector(0,extForceVector);

        // calculate right hand side
        NuTo::FullVector<double,Eigen::Dynamic> rhsVector = dispForceVector + extForceVector;

        // solve
        NuTo::SparseDirectSolverMUMPS mySolver;
        NuTo::FullVector<double,Eigen::Dynamic> displacementVector;
        stiffnessMatrix.SetOneBasedIndexing();
        mySolver.Solve(stiffnessMatrix, rhsVector, displacementVector);

        // write displacements to node
        myStructure.NodeMergeActiveDofValues(displacementVector);

        // calculate residual
        NuTo::FullVector<double,Eigen::Dynamic> intForceVector;
        myStructure.BuildGlobalGradientInternalPotentialVector(intForceVector);
        NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> residualVector (extForceVector - intForceVector);
        if (residualVector.Abs().Max() > 1e-8)
        {
            std::cout << "[Plane2D3N] : residual force vector is not zero." << std::endl;
        }

        //calculate engineering strain of myelement1 at all integration points
        //the size the matrix is not important and reallocated within the procedure
        NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> EngineeringStrain(6,1);
        //correct strain
        NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> EngineeringStrainCorrect(6,1);
        EngineeringStrainCorrect(0,0) = 0.1;
        EngineeringStrainCorrect(1,0) = -0.025;

        NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> EngineeringStress(6,1);
        NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> EngineeringStressCorrect(6,1);
        EngineeringStressCorrect(0,0) = 1.0416666666666667;
        EngineeringStressCorrect(2,0) = 0.2083333333333333;

        for (int element=0; element<8; element++)
        {
            myStructure.ElementGetEngineeringStrain(element, EngineeringStrain);

            if (PRINTRESULT)
            {
                std::cout << "EngineeringStrainCorrect" << std::endl;
                EngineeringStrainCorrect.Info();
                std::cout << "EngineeringStrain" << std::endl;
                EngineeringStrain.Info();
            }

            if ((EngineeringStrain-EngineeringStrainCorrect).cwiseAbs().maxCoeff() > 1e-8)
            {
                if (!PRINTRESULT)
                {
                    std::cout << "EngineeringStrainCorrect" << std::endl;
                    EngineeringStrainCorrect.Info();
                    std::cout << "EngineeringStrain" << std::endl;
                    EngineeringStrain.Info();
                }
                throw NuTo::Exception("[Plane2D3N] : strain is not correct.");
            }

            //calculate engineering strain at all integration points
            myStructure.ElementGetEngineeringStress(element, EngineeringStress);
            //correct stress

            if (PRINTRESULT)
            {
                std::cout << "EngineeringStressCorrect" << std::endl;
                EngineeringStressCorrect.Info();
                std::cout << "EngineeringStress" << std::endl;
                EngineeringStress.Info();
            }

            if ((EngineeringStress-EngineeringStressCorrect).cwiseAbs().maxCoeff() > 1e-8)
            {
                if (!PRINTRESULT)
                {
                std::cout << "EngineeringStressCorrect" << std::endl;
                EngineeringStressCorrect.Info();
                std::cout << "EngineeringStress" << std::endl;
                EngineeringStress.Info();
                }
                throw NuTo::Exception("[Plane2D3N] : stress is not correct.");
            }
        }


        // visualize results
        myStructure.AddVisualizationComponentDisplacements();
        myStructure.AddVisualizationComponentEngineeringStrain();
        myStructure.AddVisualizationComponentEngineeringStress();
        myStructure.ExportVtkDataFileElements( "Plane2D3N.vtk");
    }
    catch (NuTo::Exception& e)
    {
        std::cout << e.ErrorMessage() << std::endl;
        return -1;
    }

    return 0;
}
