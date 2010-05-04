// $Id: $
#ifndef PLANE2D_H
#define PLANE2D_H

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/elements/Plane.h"
namespace NuTo
{
class NodeBase;
class SectionBase;
class Plane2D : public Plane
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION

public:
    Plane2D(NuTo::StructureBase* rStructure,
    		ElementData::eElementDataType rElementDataType, IntegrationType::eIntegrationType rIntegrationType,
    		IpData::eIpDataType rIpDataType);

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Plane);
    }
#endif  // ENABLE_SERIALIZATION

    //! @brief returns the global dimension of the element
    //! this is required to check, if an element can be used in a 1d, 2D or 3D Structure
    //! there is also a routine GetLocalDimension, which is e.g. 2 for plane elements and 1 for truss elements
    //! @return global dimension
    int GetGlobalDimension()const
    {
        return 2;
    }

    //! @brief calculates the local coordinates of the nodes
    //! @param localCoordinates vector with already correct size allocated
    //! this can be checked with an assertation
    void CalculateLocalCoordinates(std::vector<double>& rLocalCoordinates)const;

    //! @brief calculates the local displacements of the nodes
    //! @param localDisplacements vector with already correct size allocated
    //! this can be checked with an assertation
    void CalculateLocalDisplacements(std::vector<double>& rLocalDisplacements)const;

    //! @brief transforms the local matrix to the global system
    //! relevant only for 2D and 3D truss elements
    void BlowLocalMatrixToGlobal(NuTo::FullMatrix<double>& FullrCoefficientMatrix)const
    {}

    //! @brief transforms the local vector to the global system
    //! relevant only for 2D and 3D truss elements
    void BlowLocalVectorToGlobal(NuTo::FullMatrix<double>& rFullVector)const
    {}

protected:
    //! @brief ... extract global dofs from nodes (mapping of local row ordering of the element matrices to the global dof ordering)
    //! @param rGlobalRowDofs ... vector of global row dofs
    void CalculateGlobalRowDofs(std::vector<int>& rGlobalRowDofs) const;

    //! @brief ... extract global dofs from nodes (mapping of local column ordering of the element matrices to the global dof ordering)
    //! @param rGlobalColumnDofs ... vector of global column dofs
    void CalculateGlobalColumnDofs(std::vector<int>& rGlobalColumnDofs) const;


};
} // namespace NuTo
#endif //PLANE2D_H
