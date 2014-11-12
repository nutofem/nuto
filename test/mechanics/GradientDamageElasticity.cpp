#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain1D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain3D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress1D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress3D.h"
#include "nuto/mechanics/constitutive/mechanics/GradientDamageEngineeringStress.h"
#include "nuto/mechanics/MechanicsException.h"

#include "nuto/math/SparseDirectSolverMUMPS.h"
#include "nuto/math/SparseMatrixCSRVector2General.h"
#include "nuto/mechanics/timeIntegration/NewmarkDirect.h"
#include "nuto/base/Logger.h"

#define PRINTRESULT
#include <eigen3/Eigen/Eigenvalues>




int main()
{
try {









    // 1D structure
    NuTo::Structure myStructure(1);


    int numElements = 3;
    double length = 10;
    double area1D = 1;

    double displacementBCLeft = -5.e-4;
    double displacementBCRight = -displacementBCLeft;

    NuTo::FullVector<double, Eigen::Dynamic> myDamageLaw(3);
    myDamageLaw(0) = NuTo::Constitutive::eDamageLawType::ISOTROPIC_EXPONENTIAL_SOFTENING;
    myDamageLaw(1) = 1.e-4;
    myDamageLaw(2) = 0.005;

    int myNumberConstitutiveLaw = myStructure.ConstitutiveLawCreate("GradientDamageEngineeringStress");

    myStructure.ConstitutiveLawSetDensity       (myNumberConstitutiveLaw, 1.0);
    myStructure.ConstitutiveLawSetYoungsModulus (myNumberConstitutiveLaw, 20000.);
    myStructure.ConstitutiveLawSetPoissonsRatio (myNumberConstitutiveLaw, 0.3);
    myStructure.ConstitutiveLawSetNonlocalRadius(myNumberConstitutiveLaw, 1.0);
    myStructure.ConstitutiveLawSetDamageLaw(myNumberConstitutiveLaw, myDamageLaw);



    // create nodes
    int numNodes= numElements * 2 + 1;



    double lengthElement=length / numElements;

    // create nodes
    NuTo::FullVector<double,Eigen::Dynamic> nodeCoordinates(1);
    for(int node = 0; node < numNodes; node++)
    {
        //std::cout << "create node: " << node << " coordinates: " << node * l_e << std::endl;
        nodeCoordinates(0) = node * lengthElement;
        if (node % 2 == 0)
            myStructure.NodeCreate(node, "displacements nonlocalEqStrain", nodeCoordinates, 0);
        else
            myStructure.NodeCreate(node, "displacements", nodeCoordinates, 0);
    }

    // Apply displacement BC via linear constraints

    int groupNodesLeft = myStructure.GroupCreate("Nodes");
    int groupNodesRight = myStructure.GroupCreate("Nodes");
    myStructure.GroupAddNodeCoordinateRange(groupNodesLeft,0,0,0);
    myStructure.GroupAddNodeCoordinateRange(groupNodesRight,0,length,length);

    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> directionConstraint(1,1);
    directionConstraint(0,0) = 1;
    myStructure.ConstraintLinearSetDisplacementNodeGroup(groupNodesLeft , directionConstraint, displacementBCLeft);
    myStructure.ConstraintLinearSetDisplacementNodeGroup(groupNodesRight, directionConstraint, displacementBCRight);
    myStructure.NodeBuildGlobalDofs();

    // check node Getter/Setter

    NuTo::FullVector<double, Eigen::Dynamic> actDofValues, depDofValues, actDofValuesStructure, depDofValuesStructure;
    myStructure.NodeExtractDofValues(actDofValues, depDofValues);
    actDofValues += 13.37;
    depDofValues += 61.74;
    myStructure.NodeMergeDofValues(actDofValues, depDofValues);
    myStructure.NodeExtractDofValues(actDofValuesStructure, depDofValuesStructure);
    double normDof = (actDofValues - actDofValuesStructure).norm() + (depDofValues - depDofValuesStructure).norm();
    if (normDof > 1.e-15)
        throw NuTo::Exception("[GradientDamageElasticity] Error while merging or extracting Dof values");

    // create elements

    // TODO: create section and weakened section

    for (int element = 0; element < numElements; element ++)
    {
        NuTo::FullVector<int, 3> nodeIds;
        nodeIds(0) = 2*element;
        nodeIds(1) = 2*element+1;
        nodeIds(2) = 2*element+2;
        myStructure.ElementCreate("Truss1D3N", nodeIds);
    }

    myStructure.Info();






} catch (NuTo::Exception& e) {
    std::cout << "Error executing GradientDamageElasticity "<< std::endl;
    std::cout << e.ErrorMessage() << std::endl;
    return -1;
}
    std::cout << "GradientDamageElasticity terminated normally."<< std::endl;
    return 0;
}
