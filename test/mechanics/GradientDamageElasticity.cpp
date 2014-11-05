#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain1D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain3D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress1D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress3D.h"
#include "nuto/mechanics/constitutive/mechanics/GradientDamagePlasticityEngineeringStress.h"
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



    double l=100;
    int numNodes = 51;
//    int numElements = numNodes-1;
    double l_e=l/(numNodes-1);
//    double area1D = 1.;

    // create nodes
    NuTo::FullVector<double,Eigen::Dynamic> nodeCoordinates(1);
    for(int node = 0; node < numNodes ; node++)
    {
        //std::cout << "create node: " << node << " coordinates: " << node * l_e << std::endl;
        nodeCoordinates(0) = node *l_e;
//        myStructure.NodeCreate(node, "displacements nonlocalEqStrain", nodeCoordinates);
    }
//    int nodeLeft = 0;
//    int nodeRight = numNodes-1;





} catch (NuTo::Exception& e) {
    std::cout << "Error executing GradientDamageElasticity "<< std::endl;
    std::cout << e.ErrorMessage() << std::endl;
    return -1;
}
    std::cout << "GradientDamageElasticity terminated normally."<< std::endl;
    return 0;
}
