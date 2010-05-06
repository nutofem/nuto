// $Id: $
#ifndef Voxel8N_H
#define Voxel8N_H

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION
#include <vector>
#include "nuto/mechanics/elements/Solid.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"

namespace NuTo
{
class ConstitutiveTangentLocal6x6;
class EngineeringStress3D;

//! @author Andrea Keßler, ISM
//! @date April 2010
//! @brief ... voxel element without nodes
class Voxel8N : public Solid
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
    //! @brief constructor
    Voxel8N(NuTo::StructureBase* rStructure,unsigned int rElementID,NuTo::SparseMatrixCSRGeneral<double>& rLocalCoefficientMatrix0,
            NuTo::ElementData::eElementDataType rElementDataType , IpData::eIpDataType rIpDataType):
            	NuTo::Solid::Solid(rStructure, rElementDataType, GetStandardIntegrationType(), rIpDataType), 
            	mLocalCoefficientMatrix0 (rLocalCoefficientMatrix0)            	
    {
        mVoxelID = (int) rElementID;
    }
    //! @brief constructor
    Voxel8N(NuTo::StructureBase* rStructure,NuTo::SparseMatrixCSRGeneral<double>& rLocalCoefficientMatrix0,
    		NuTo::ElementData::eElementDataType rElementDataType, IpData::eIpDataType rIpDataType):
        NuTo::Solid::Solid(rStructure, rElementDataType, GetStandardIntegrationType(), rIpDataType),
        mLocalCoefficientMatrix0 (rLocalCoefficientMatrix0)
    {
        mVoxelID = -1;
    }

    #ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive> void serialize(Archive & ar,
                                           const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Solid);
    }
#endif  // ENABLE_SERIALIZATION

    //! @brief returns the enum (type of the element)
    //! @return enum
    NuTo::Element::eElementType GetEnumType()const
    {
        return NuTo::Element::VOXEL8N;
    }

    //! @brief returns the number of nodes in this element
    //! @return number of nodes
    int GetNumNodes() const
    {
        return 8;
    }

    //! @brief returns the number of shape functions
    //! this is required for the calculation of the derivatives of the shape functions
    //! whose size is GetLocalDimension*GetNumShapeFunctions
    //! @return local dimension
    virtual int GetNumShapeFunctions() const
    {
        return 8;
    }

    //! @brief calculates the shape functions
    //! @param rLocalCoordinates local coordinates of the integration point
    //! @param shape functions for all the nodes
    void CalculateShapeFunctions(const double rLocalCoordinates[3],
                                 std::vector<double>& rShapeFunctions) const;

    //! @brief calculates the derivative of the shape functions with respect to local coordinates
    //! @param rLocalCoordinates local coordinates of the integration point
    //! @param derivative of the shape functions for all the nodes,
    //! first all the directions for a single node, and then for the next node
    void CalculateDerivativeShapeFunctionsLocal(
        const double rLocalCoordinates[3],
        std::vector<double>& rDerivativeShapeFunctions) const;

    //! @brief returns the number of local degrees of freedom
    //! @return number of local degrees of freedom
    inline int GetNumDofs() const
    {
        return 24;
    }

    //! @TODO change routine to find the nodes

    //! @brief returns a pointer to the i-th node of the element
    //! @param local node number
    //! @return pointer to the node
    NodeBase* GetNode(int rLocalNodeNumber)
    {
        throw NuTo::MechanicsException("[NuTo::Voxel8N::GetNode] implement first.");
        /*
        assert(rLocalNodeNumber>=0 && rLocalNodeNumber<=7);
        return mNodes[rLocalNodeNumber];
         */
    }

    //! @brief returns a pointer to the i-th node of the element
    //! @param local node number
    //! @return pointer to the node
    const NodeBase* GetNode(int rLocalNodeNumber) const
    {
        throw NuTo::MechanicsException("[NuTo::Voxel8N::GetNode] implement first.");
        /*
        assert(rLocalNodeNumber>=0 && rLocalNodeNumber<=7);
        return mNodes[rLocalNodeNumber];
        */
    }

    //! @brief sets the rLocalNodeNumber-th node of the element
    //! @param local node number
    //! @param pointer to the node
    void SetNode(int rLocalNodeNumber, NodeBase* rNode)
    {
        throw NuTo::MechanicsException("[NuTo::Voxel8N::SetNode] implement first.");
        /*
       assert(rLocalNodeNumber>=0 && rLocalNodeNumber<=7);
        mNodes[rLocalNodeNumber] = rNode;
        */
    }
    //! @brief ... reorder nodes is not possible in gris structure
    void ReorderNodes()
    {
        throw NuTo::MechanicsException("[NuTo::Voxel8N::ReorderNodes] GridElements do not have nodes for reordering.");

    }
    //! @brief calculate list of global dofs related to the entries in the element stiffness matrix
    //! @param rGlobalDofsRow global dofs corresponding to the rows of the matrix
    //! @param rGlobalDofsColumn global dofs corresponding to the columns of the matrix
    void CalculateGlobalDofs(std::vector<int>& rGlobalDofsRow,
                              std::vector<int>& rGlobalDofsColumn) const
    {
        throw MechanicsException("[NuTo::Voxel8N::CalculateGlobalDofs] to be implemented.");
    }


    //! @brief returns the enum of the standard integration type for this element
    NuTo::IntegrationType::eIntegrationType GetStandardIntegrationType();

protected:
    int mVoxelID;
    NuTo::SparseMatrixCSRGeneral<double>& mLocalCoefficientMatrix0;
};
}
#endif //Voxel8N_H
