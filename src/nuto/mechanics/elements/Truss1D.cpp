// $Id$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION

#include "nuto/mechanics/constitutive/mechanics/DeformationGradient1D.h"
#include "nuto/mechanics/elements/Truss1D.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include <assert.h>


NuTo::Truss1D::Truss1D(NuTo::StructureBase* rStructure, ElementData::eElementDataType rElementDataType,
		IntegrationType::eIntegrationType rIntegrationType, IpData::eIpDataType rIpDataType) :
        Truss(rStructure, rElementDataType, rIntegrationType, rIpDataType)
{}

//! @brief calculates the local coordinates of the nodes
//! @param localCoordinates vector with already correct size allocated
//! this can be checked with an assertation
void NuTo::Truss1D::CalculateLocalCoordinates(std::vector<double>& rLocalCoordinates)const
{
    assert((int)rLocalCoordinates.size()==GetNumNodes());
    for (int theNode=0; theNode<GetNumNodes(); theNode++)
    {
        GetNode(theNode)->GetCoordinates1D(&(rLocalCoordinates[theNode]));
    }
}

//! @brief calculates the local displacements of the nodes
//! @param time derivative (0 temperature, 1 temperature rate, 2 second time derivative of temperature)
//! @param localDisplacements vector with already correct size allocated
//! this can be checked with an assertation
void NuTo::Truss1D::CalculateLocalDisplacements(int rTimeDerivative, std::vector<double>& rLocalDisplacements)const
{
    assert(rTimeDerivative>=0);
    assert(rTimeDerivative<3);
    assert((int)rLocalDisplacements.size()==GetNumNodes());
    for (int theNode=0; theNode<GetNumNodes(); theNode++)
    {
        GetNode(theNode)->GetDisplacements1D(rTimeDerivative, &(rLocalDisplacements[theNode]));
    }
}

// interpolate geometry
void NuTo::Truss1D::InterpolateCoordinatesFrom1D(double rLocalCoordinates, double rGlobalCoordinates[3]) const
{
    // calculate shape functions
    std::vector<double> ShapeFunctions(this->GetNumNodes());
    this->CalculateShapeFunctions(rLocalCoordinates, ShapeFunctions);

    // start interpolation
    rGlobalCoordinates[0] = 0.0;
    rGlobalCoordinates[1] = 0.0;
    rGlobalCoordinates[2] = 0.0;
    for (int theNode = 0; theNode < this->GetNumNodes(); theNode++)
    {
        // get node coordinate
        double NodeCoordinate;
        GetNode(theNode)->GetCoordinates1D(&NodeCoordinate);

        // add node contribution
        rGlobalCoordinates[0] += ShapeFunctions[theNode] *  NodeCoordinate;
    }
}

// interpolate displacements
void NuTo::Truss1D::InterpolateDisplacementsFrom1D(double rLocalCoordinates, double rGlobalDisplacements[3]) const
{
    // calculate shape functions
    std::vector<double> ShapeFunctions(this->GetNumNodes());
    this->CalculateShapeFunctions(rLocalCoordinates, ShapeFunctions);

    // start interpolation
    rGlobalDisplacements[0] = 0.0;
    rGlobalDisplacements[1] = 0.0;
    rGlobalDisplacements[2] = 0.0;
    for (int theNode = 0; theNode < this->GetNumNodes(); theNode++)
    {
        // get node displacements
        double NodeDisplacement;
        GetNode(theNode)->GetDisplacements1D(&NodeDisplacement);

        // add node contribution
        rGlobalDisplacements[0] += ShapeFunctions[theNode] *  NodeDisplacement;
    }
}

// build global row dofs
void NuTo::Truss1D::CalculateGlobalRowDofs(std::vector<int>& rGlobalRowDofs, int rNumDispDofs, int rNumTempDofs,
		int rNumNonlocalEqPlasticStrainDofs, int rNumNonlocalTotalStrainDofs) const
{
    rGlobalRowDofs.resize(rNumDispDofs+rNumTempDofs+rNumNonlocalEqPlasticStrainDofs+rNumNonlocalTotalStrainDofs);
    int indexDisp(0);
    int indexNonlocalTotalStrain(0);
    int numNodes(this->GetNumNodes());
    for (int nodeCount = 0; nodeCount < numNodes; nodeCount++)
    {
        const NodeBase * nodePtr(GetNode(nodeCount));
        if (rNumDispDofs>0)
        {
            for (int countDisp=0; countDisp<nodePtr->GetNumDisplacements(); countDisp++)
            {
        	    rGlobalRowDofs[indexDisp] = nodePtr->GetDofDisplacement(countDisp);
        	    indexDisp++;
            }
        }
        if (nodePtr->GetNumTemperatures()>0 && rNumTempDofs>0)
        {
            rGlobalRowDofs[rNumDispDofs + nodeCount ] = nodePtr->GetDofTemperature();
        }
        if (nodePtr->GetNumNonlocalEqPlasticStrain()>0 && rNumNonlocalEqPlasticStrainDofs>0)
        {
            for (int countNonlocalEqPlasticStrain=0; countNonlocalEqPlasticStrain<nodePtr->GetNumNonlocalEqPlasticStrain(); countNonlocalEqPlasticStrain++)
            {
				rGlobalRowDofs[rNumDispDofs + rNumTempDofs + nodeCount + numNodes*countNonlocalEqPlasticStrain] = nodePtr->GetDofNonlocalEqPlasticStrain(countNonlocalEqPlasticStrain);
            }
        }
        if (nodePtr->GetNumNonlocalTotalStrain()>0 && rNumNonlocalTotalStrainDofs>0)
        {
            for (int countNonlocalTotalStrain=0; countNonlocalTotalStrain<nodePtr->GetNumNonlocalTotalStrain(); countNonlocalTotalStrain++)
            {
				rGlobalRowDofs[rNumDispDofs + rNumTempDofs + rNumNonlocalEqPlasticStrainDofs + indexNonlocalTotalStrain] = nodePtr->GetDofNonlocalTotalStrain(countNonlocalTotalStrain);
				indexNonlocalTotalStrain++;
            }
        }
    }

}

// build global column dof
void NuTo::Truss1D::CalculateGlobalColumnDofs(std::vector<int>& rGlobalColumnDofs, int rNumDispDofs, int rNumTempDofs, int rNumNonlocalEqPlasticStrainDofs, int rNumNonlocalTotalStrainDofs) const
{
    this->CalculateGlobalRowDofs(rGlobalColumnDofs,rNumDispDofs,rNumTempDofs,rNumNonlocalEqPlasticStrainDofs,rNumNonlocalTotalStrainDofs);
}

// check element definition
void NuTo::Truss1D::CheckElement()
{
    // check nodes
    for (int nodeCount = 0; nodeCount < this->GetNumNodes(); nodeCount++)
    {
        int numCoordinates(GetNode(nodeCount)->GetNumCoordinates());
    	if (numCoordinates<1 || numCoordinates>3)
        {
            throw MechanicsException("[NuTo::Truss1D::CheckElement] invalid node type (check node definition for coordinates).");
        }
    }

    // check node ordering (element length must be positive) and for changing sign in jacobian determinant
    //calculate local coordinates
    std::vector<double> localNodeCoord(this->GetNumLocalDofs());
    this->CalculateLocalCoordinates(localNodeCoord);

    // check number of integration points
    if (this->GetNumIntegrationPoints() < 1)
    {
        throw MechanicsException("[NuTo::Truss1D::CheckElement] invalid integration type.");
    }

    // check sign of the jacobian determinant of the first integration point
    double localIPCoord;
    this->GetLocalIntegrationPointCoordinates(0, localIPCoord);

    std::vector<double> derivativeShapeFunctions(this->GetLocalDimension()*this->GetNumShapeFunctions());
    this->CalculateDerivativeShapeFunctions(localIPCoord, derivativeShapeFunctions);

    double detJacobian = this->DetJacobian(derivativeShapeFunctions,localNodeCoord);

    // reorder nodes if determinant is negative
    if (detJacobian < 0.0)
    {
        this->ReorderNodes();
        // recalculate node coordinates after reordering
        this->CalculateLocalCoordinates(localNodeCoord);
    }

    // check jacobian determinant for all integration points for positive sign and calculate element length
    double length = 0;
    for (int ipCount = 0; ipCount < this->GetNumIntegrationPoints(); ipCount++)
    {
        // calculate jacobian determinant
        this->GetLocalIntegrationPointCoordinates(ipCount, localIPCoord);
        this->CalculateDerivativeShapeFunctions(localIPCoord, derivativeShapeFunctions);
        detJacobian = this->DetJacobian(derivativeShapeFunctions,localNodeCoord);
        if (detJacobian <= 0)
        {
            throw MechanicsException("[NuTo::Truss1D::CheckElement] element is not properly defined by this nodes (zero or negative jacobian determinant).");
        }
        length += this->GetIntegrationPointWeight(ipCount) * detJacobian;
    }

    // check element length
    if (length < 1e-14)
    {
        throw MechanicsException("[NuTo::Truss1D::CheckElement] element with zero length (check nodes).");
    }
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::Truss1D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::Truss1D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::Truss1D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::Truss1D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::Truss1D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::Truss1D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::Truss1D::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize Truss1D" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Truss);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize Truss1D" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::Truss1D)
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::Truss1D)
#endif // ENABLE_SERIALIZATION
