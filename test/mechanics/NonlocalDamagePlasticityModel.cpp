#include "nuto/math/FullMatrix.h"
#include "nuto/math/FullVector.h"
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
    NuTo::FullVector<double,Eigen::Dynamic> Coordinates(2);
    Coordinates(0) = 0.0;
    Coordinates(1) = 0.0;
    int node1 = myStructure.NodeCreate(Coordinates);

    Coordinates(0) = 1.0;
    Coordinates(1) = 0.0;
    int node2 = myStructure.NodeCreate(Coordinates);

    Coordinates(0) = 2.0;
    Coordinates(1) = 0.0;
    int node3 = myStructure.NodeCreate(Coordinates);

    Coordinates(0) = 1.0;
    Coordinates(1) = 1.0;
    int node5 = myStructure.NodeCreate(Coordinates);

    Coordinates(0) = 0.0;
    Coordinates(1) = 1.0;
    int node4 = myStructure.NodeCreate(Coordinates);

    Coordinates(0) = 2.0;
    Coordinates(1) = 1.0;
    int node6 = myStructure.NodeCreate(Coordinates);

    Coordinates(0) = 0.0;
    Coordinates(1) = 2.0;
    int node7 = myStructure.NodeCreate(Coordinates);

    Coordinates(0) = 1.0;
    Coordinates(1) = 2.0;
    int node8 = myStructure.NodeCreate(Coordinates);

    Coordinates(0) = 2.0;
    Coordinates(1) = 2.0;
    int node9 = myStructure.NodeCreate(Coordinates);

    int interpolationType = myStructure.InterpolationTypeCreate("Quad2D");
    myStructure.InterpolationTypeAdd(interpolationType, NuTo::Node::COORDINATES, NuTo::Interpolation::EQUIDISTANT1);
    myStructure.InterpolationTypeAdd(interpolationType, NuTo::Node::DISPLACEMENTS, NuTo::Interpolation::EQUIDISTANT1);

    int interpolationType3 = myStructure.InterpolationTypeCreate("Triangle2D");
    myStructure.InterpolationTypeAdd(interpolationType3, NuTo::Node::COORDINATES, NuTo::Interpolation::EQUIDISTANT1);
    myStructure.InterpolationTypeAdd(interpolationType3, NuTo::Node::DISPLACEMENTS, NuTo::Interpolation::EQUIDISTANT1);

    //create elements
    NuTo::FullVector<int,Eigen::Dynamic> Incidence(4);
    Incidence(0) = node1;
    Incidence(1) = node2;
    Incidence(2) = node5;
    Incidence(3) = node4;
    int myElement1 = myStructure.ElementCreate(interpolationType,Incidence,"ConstitutiveLawIpNonlocal","StaticDataNonlocal");

    Incidence(0) = node2;
    Incidence(1) = node3;
    Incidence(2) = node6;
    Incidence(3) = node5;
    int myElement2 = myStructure.ElementCreate(interpolationType,Incidence,"ConstitutiveLawIpNonlocal","StaticDataNonlocal");

    Incidence(0) = node4;
    Incidence(1) = node5;
    Incidence(2) = node8;
    Incidence(3) = node7;
    int myElement3 = myStructure.ElementCreate(interpolationType,Incidence,"ConstitutiveLawIpNonlocal","StaticDataNonlocal");

/*  Incidence(0,0) = node5;
    Incidence(1,0) = node6;
    Incidence(2,0) = node9;
    Incidence(3,0) = node8;
    int myElement4 = myStructure.ElementCreate(interpolationType,Incidence,"ConstitutiveLawIpNonlocal","StaticDataNonlocal");
    myStructure.ElementSetIntegrationType(myElement4,"2D4NGauss4Ip","StaticDataNonlocal");
*/
    NuTo::FullVector<int,Eigen::Dynamic> Incidence3(3);

    Incidence3(0) = node5;
    Incidence3(1) = node6;
    Incidence3(2) = node9;
    int myElement4 = myStructure.ElementCreate(interpolationType3,Incidence3,"ConstitutiveLawIpNonlocal","StaticDataNonlocal");

    Incidence3(0) = node5;
    Incidence3(1) = node9;
    Incidence3(2) = node8;
    int myElement5 = myStructure.ElementCreate(interpolationType3,Incidence3,"ConstitutiveLawIpNonlocal","StaticDataNonlocal");

    myStructure.ElementTotalConvertToInterpolationType();

    //create constitutive law
    int myMatDamage = myStructure.ConstitutiveLawCreate("NonlocalDamagePlasticityEngineeringStress");
    myStructure.ConstitutiveLawSetYoungsModulus(myMatDamage,9);
    myStructure.ConstitutiveLawSetPoissonsRatio(myMatDamage,0.25);
    myStructure.ConstitutiveLawSetNonlocalRadius(myMatDamage,2.);
    myStructure.ConstitutiveLawSetTensileStrength(myMatDamage,2);
    myStructure.ConstitutiveLawSetCompressiveStrength(myMatDamage,20);
    myStructure.ConstitutiveLawSetBiaxialCompressiveStrength(myMatDamage,25);
    myStructure.ConstitutiveLawSetFractureEnergy(myMatDamage,0.2);

    int myMatLin = myStructure.ConstitutiveLawCreate("LinearElasticEngineeringStress");
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
    NuTo::FullVector<double,Eigen::Dynamic> dispForceVector;

    myStructure.ElementTotalUpdateTmpStaticData();
    myStructure.BuildGlobalCoefficientMatrix0(stiffnessMatrix, dispForceVector);
    stiffnessMatrix.RemoveZeroEntries(0,1e-14);

    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> fullStiffnessMatrixElastic(stiffnessMatrix);
    std::cout<<"stiffnessMatrix elastic"<<std::endl;
    fullStiffnessMatrixElastic.Info();
    NuTo::FullVector<double,Eigen::Dynamic> displacements;
    NuTo::FullVector<double,Eigen::Dynamic> dependentDofs;
    NuTo::FullVector<double,Eigen::Dynamic> intForce;
    NuTo::FullVector<double,Eigen::Dynamic> intForce2;

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
        break;
        }

        NuTo::FullVector<double,Eigen::Dynamic> matrixRightDisp(2);
        matrixRightDisp.SetValue(0,0,rightDisp);
        matrixRightDisp.SetValue(1,0,0.);

        myStructure.NodeSetDisplacements(node3,matrixRightDisp);
        myStructure.NodeSetDisplacements(node6,matrixRightDisp);
        myStructure.NodeSetDisplacements(node9,matrixRightDisp);

        NuTo::FullVector<double,Eigen::Dynamic> matrixCenterDisp(2);
        if (theLoadStep!=2)
            matrixCenterDisp.SetValue(0,0,0.5*rightDisp);
        else
            matrixCenterDisp.SetValue(0,0,0.4*rightDisp);
        matrixCenterDisp.SetValue(1,0,0.);

        myStructure.NodeSetDisplacements(node2,matrixCenterDisp);
        myStructure.NodeSetDisplacements(node5,matrixCenterDisp);
        myStructure.NodeSetDisplacements(node8,matrixCenterDisp);

        NuTo::FullVector<double,Eigen::Dynamic> matrixLeftDisp(2);
        matrixLeftDisp.SetValue(0,0,0.0);
        matrixLeftDisp.SetValue(1,0,0.);

        myStructure.NodeSetDisplacements(node1,matrixLeftDisp);
        myStructure.NodeSetDisplacements(node4,matrixLeftDisp);
        myStructure.NodeSetDisplacements(node7,matrixLeftDisp);

        //NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>matrixLeftDispNode1(2,1);
        //matrixLeftDispNode1.SetValue(0,0,0.0);
        //matrixLeftDispNode1.SetValue(1,0,0.);
        //myStructure.NodeSetDisplacements(node1,matrixLeftDispNode1);


        myStructure.ElementTotalUpdateTmpStaticData();

        myStructure.BuildGlobalCoefficientMatrix0(stiffnessMatrix, dispForceVector);
        stiffnessMatrix.RemoveZeroEntries(0,1e-14);

        NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> fullStiffnessMatrix(stiffnessMatrix);

        std::cout<<"stiffnessMatrix analytic"<<std::endl;
        fullStiffnessMatrix.Info();

        myStructure.NodeExtractDofValues(displacements,dependentDofs);
        myStructure.BuildGlobalGradientInternalPotentialVector(intForce);

        double delta(1e-8);
        NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> stiffnessMatrixCD(displacements.GetNumRows(),displacements.GetNumRows());

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
        double maxerror((fullStiffnessMatrix-stiffnessMatrixCD).cwiseAbs().maxCoeff(&row,&col));
        switch(theLoadStep)
        {
        case 0:
            std::cout << "max difference in stiffness matrix for uniform plastic loading " << maxerror << " row " << row << " col " << col << std::endl;
        break;
        case 1:
        {
            std::cout << "max difference in stiffness matrix for unloading " << maxerror << " row " << row << " col " << col << std::endl;
            double omega(fullStiffnessMatrix(0,0)/fullStiffnessMatrixElastic(0,0));
            NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> diffMatrix(fullStiffnessMatrixElastic*omega-NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>(stiffnessMatrix));
            double maxerror2(diffMatrix.cwiseAbs().maxCoeff(&row,&col));
            std::cout<<"stiffnessMatrix elastic*omega"<<std::endl;
            diffMatrix.Info();
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
    myStructure.ExportVtkDataFileElements("NonlocalDamagePlasticityModel.vtk");
#endif

    NuTo::FullVector<double,Eigen::Dynamic> shift(2);
    shift(0) = 3.+myStructure.ConstitutiveLawGetNonlocalRadius(myMatDamage);
    shift(1) = 0;
    myStructure.CopyAndTranslate(shift);

#ifdef ENABLE_VISUALIZE
    myStructure.ExportVtkDataFileElements("NonlocalDamagePlasticityModelCopied.vtk");
#endif

    }
    catch (NuTo::Exception& e)
    {
        std::cout << e.ErrorMessage() << std::endl;
    }
    return 0;
}
