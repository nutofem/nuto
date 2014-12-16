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
    for (int theNode=0; theNode<GetNumNodesGeometry(); theNode++)
    {
        GetNodeGeometry(theNode)->GetCoordinates1D(&(rLocalCoordinates[theNode]));
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
    assert((int)rLocalDisplacements.size()==GetNumNodesField());
    for (int theNode=0; theNode<GetNumNodesField(); theNode++)
    {
        GetNodeField(theNode)->GetDisplacements1D(rTimeDerivative, &(rLocalDisplacements[theNode]));
    }
}

// interpolate geometry
void NuTo::Truss1D::InterpolateCoordinatesFrom1D(double rLocalCoordinates, double rGlobalCoordinates[3]) const
{
    // calculate shape functions
    std::vector<double> ShapeFunctions(this->GetNumNodesGeometry());
    this->CalculateShapeFunctionsGeometry(rLocalCoordinates, ShapeFunctions);

    // start interpolation
    rGlobalCoordinates[0] = 0.0;
    rGlobalCoordinates[1] = 0.0;
    rGlobalCoordinates[2] = 0.0;
    for (int theNode = 0; theNode < this->GetNumNodesGeometry(); theNode++)
    {
        // get node coordinate
        double NodeCoordinate;
        GetNode(theNode)->GetCoordinates1D(&NodeCoordinate);

        // add node contribution
        rGlobalCoordinates[0] += ShapeFunctions[theNode] *  NodeCoordinate;
    }
}

// interpolate displacements
void NuTo::Truss1D::InterpolateDisplacementsFrom1D(int rTimeDerivatrive, double rLocalCoordinates, double rGlobalDisplacements[3]) const
{
    // calculate shape functions
    std::vector<double> ShapeFunctions(this->GetNumNodesField());
    this->CalculateShapeFunctionsField(rLocalCoordinates, ShapeFunctions);

    // start interpolation
    rGlobalDisplacements[0] = 0.0;
    rGlobalDisplacements[1] = 0.0;
    rGlobalDisplacements[2] = 0.0;
    for (int theNode = 0; theNode < this->GetNumNodesField(); theNode++)
    {
        // get node displacements
        double NodeDisplacement;
       	GetNodeField(theNode)->GetDisplacements1D(rTimeDerivatrive,&NodeDisplacement);

        // add node contribution
        rGlobalDisplacements[0] += ShapeFunctions[theNode] *  NodeDisplacement;
    }
}

// build global row dofs
void NuTo::Truss1D::CalculateGlobalRowDofs(
        std::vector<int>& rGlobalRowDofs,
        int rNumDispDofs,
        int rNumTempDofs,
		int rNumNonlocalEqPlasticStrainDofs,
		int rNumNonlocalTotalStrainDofs,
        int rNumNonlocalEqStrainDofs,
        int rNumRelativeHumidityDofs,
        int rNumWaterPhaseFractionDofs) const
{
    rGlobalRowDofs.resize(rNumDispDofs+rNumTempDofs+rNumNonlocalEqPlasticStrainDofs+rNumNonlocalTotalStrainDofs+rNumNonlocalEqStrainDofs+rNumRelativeHumidityDofs+rNumWaterPhaseFractionDofs);



    // indices for the current dof type
    int iGlobalDisp                    = 0;
    int iGlobalTemp                    = 0;
    int iGlobalNonlocalEqPlasticStrain = 0;
    int iGlobalNonlocalTotalStrain     = 0;
    int iGlobalNonlocalEqStrain        = 0;
    int iGlobalWaterPhaseFraction      = 0;
    int iGlobalRelativeHumidity        = 0;

    int numNodes(this->GetNumNodes());

    for (int iNode = 0; iNode < numNodes; iNode++)
    {
        const NodeBase * node = GetNode(iNode);

        if (rNumDispDofs>0)
        {
            for (int iLocalLocal = 0; iLocalLocal < node->GetNumDisplacements(); iLocalLocal++)
            {
        	    rGlobalRowDofs[iGlobalDisp] = node->GetDofDisplacement(iLocalLocal);
        	    iGlobalDisp++;
            }
        }

        int iStartTemp = rNumDispDofs;
        if (node->GetNumTemperatures()>0 && rNumTempDofs>0)
        {
            rGlobalRowDofs[iStartTemp + iGlobalTemp] = node->GetDofTemperature();

            // is only increased if the node actually has the corresponding dofs
            iGlobalTemp++;

        }

        int iStartNonlocalEqPlasticStrain = iStartTemp + rNumTempDofs;
        if (node->GetNumNonlocalEqPlasticStrain()>0 && rNumNonlocalEqPlasticStrainDofs>0)
        {
            for (int iNonlocalEqPlasticStrain=0; iNonlocalEqPlasticStrain < node->GetNumNonlocalEqPlasticStrain(); iNonlocalEqPlasticStrain++)
            {
                // abr: NonlocalEqPlasticStrain --> e*[]
                // the final vector
                // should look like
                // [ dispDofs ... ]
                //
                // [ tempDofs ... ]
                //
                // [ e*[iNonlocalEqPlasticStrain = 0] of node 0
                // [ e*[iNonlocalEqPlasticStrain = 0] of node 1
                // [  ... ]
                // [ e*[iNonlocalEqPlasticStrain = 0] of node N
                //
                // [ e*[iNonlocalEqPlasticStrain = 1] of node 0
                // [ e*[iNonlocalEqPlasticStrain = 1] of node 1
                // [  ... ]
                // [ e*[iNonlocalEqPlasticStrain = 1] of node N

                // [ ...]
                // --> an offset of iNonlocalEqPlasticStrain * N is applied
                int iGlobalOffset =  numNodes*iNonlocalEqPlasticStrain;

				rGlobalRowDofs[iStartNonlocalEqPlasticStrain + iGlobalNonlocalEqPlasticStrain + iGlobalOffset] = node->GetDofNonlocalEqPlasticStrain(iNonlocalEqPlasticStrain);
            }
            // is only increased if the node actually has the corresponding dofs
            iGlobalNonlocalEqPlasticStrain++;
        }

        int iStartNonlocalTotalStrain = iStartNonlocalEqPlasticStrain + rNumNonlocalEqPlasticStrainDofs;
        if (node->GetNumNonlocalTotalStrain()>0 && rNumNonlocalTotalStrainDofs>0)
        {
            // arranged like the displacements
            for (int iNonlocalTotalStrain=0; iNonlocalTotalStrain < node->GetNumNonlocalTotalStrain(); iNonlocalTotalStrain++)
            {
				rGlobalRowDofs[iStartNonlocalTotalStrain + iGlobalNonlocalTotalStrain] = node->GetDofNonlocalTotalStrain(iNonlocalTotalStrain);
				iGlobalNonlocalTotalStrain++;
            }
        }

        int iStartNonlocalEqStrain = iStartNonlocalTotalStrain + rNumNonlocalTotalStrainDofs;
        if (node->GetNumNonlocalEqStrain()>0 && rNumNonlocalEqStrainDofs>0)
        {
            rGlobalRowDofs[iStartNonlocalEqStrain + iGlobalNonlocalEqStrain] = node->GetDofNonlocalEqStrain();
            iGlobalNonlocalEqStrain ++;
        }

        int iStartWaterPhaseFraction = iStartNonlocalEqStrain +rNumNonlocalEqStrainDofs;
        if (node->GetNumWaterPhaseFraction()>0 && rNumWaterPhaseFractionDofs>0)
        {
            rGlobalRowDofs[iStartWaterPhaseFraction+iGlobalWaterPhaseFraction] = node->GetDofWaterPhaseFraction();
            iGlobalWaterPhaseFraction++;
        }

        int iStartRelativeHumidity = iStartWaterPhaseFraction + rNumWaterPhaseFractionDofs;
        if (node->GetNumRelativeHumidity()>0 && rNumRelativeHumidityDofs>0)
        {
            rGlobalRowDofs[iStartRelativeHumidity+iGlobalRelativeHumidity]  = node->GetDofRelativeHumidity();
            iGlobalRelativeHumidity++;
        }

        //int iStartNextDof = iStartNonlocalEqStrain + rNumNonlocalEqStrain;


    }

//    std::cout << Eigen::Map<Eigen::VectorXi>(rGlobalRowDofs.data(), rGlobalRowDofs.size()) << std::endl;
//    std::cout << "==========================================" << std::endl;

}



// check element definition
void NuTo::Truss1D::CheckElement()
{
    // check nodes
    for (int nodeCount = 0; nodeCount < this->GetNumNodesGeometry(); nodeCount++)
    {
        int numCoordinates(GetNodeGeometry(nodeCount)->GetNumCoordinates());
    	if (numCoordinates<1 || numCoordinates>3)
        {
            throw MechanicsException("[NuTo::Truss1D::CheckElement] invalid node type (check node definition for coordinates).");
        }
    }

    // check node ordering (element length must be positive) and for changing sign in jacobian determinant
    //calculate local coordinates
    std::vector<double> localNodeCoord(this->GetNumNodesGeometry());
    this->CalculateLocalCoordinates(localNodeCoord);

    // check number of integration points
    if (this->GetNumIntegrationPoints() < 1)
    {
        throw MechanicsException("[NuTo::Truss1D::CheckElement] invalid integration type.");
    }

    // check sign of the jacobian determinant of the first integration point
    double localIPCoord;
    this->GetLocalIntegrationPointCoordinates(0, localIPCoord);

    std::vector<double> derivativeShapeFunctions(this->GetNumNodesGeometry());
    this->CalculateDerivativeShapeFunctionsGeometry(localIPCoord, derivativeShapeFunctions);

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
        this->CalculateDerivativeShapeFunctionsGeometry(localIPCoord, derivativeShapeFunctions);
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
