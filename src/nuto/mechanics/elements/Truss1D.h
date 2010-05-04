// $Id: $
#ifndef Truss1D_H
#define Truss1D_H

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/elements/Truss.h"
namespace NuTo
{
class NodeBase;
class SectionBase;
class Truss1D : public Truss
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION

public:
    Truss1D(NuTo::StructureBase* rStructure,
    		ElementData::eElementDataType rElementDataType, IntegrationType::eIntegrationType rIntegrationType,
    		IpData::eIpDataType rIpDataType);

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Truss);
    }
#endif  // ENABLE_SERIALIZATION

    //! @brief returns the global dimension of the element
    //! this is required to check, if an element can be used in a 1d, 2D or 3D Structure
    //! there is also a routine GetLocalDimension, which is e.g. 2 for plane elements and 1 for truss elements
    //! @return global dimension
    int GetGlobalDimension()const
    {
        return 1;
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

    //! @brief ... interpolate three-dimensional global point coordinates from one-dimensional local point coordinates (element coordinates system)
    //! @param rLocalCoordinates ... one-dimensional local point coordinates
    //! @param rGlobalCoordinates ... three-dimension global point coordinates
    void InterpolateCoordinatesFrom1D(double rLocalCoordinates, double rGlobalCoordinates[3]) const;

    //! @brief ... interpolate three-dimensional global point displacements from one-dimensional local point coordinates (element coordinates system)
    //! @param rLocalCoordinates ... one-dimensional local point coordinates
    //! @param rGlobalDisplacements ... three-dimension global point displacements
    void InterpolateDisplacementsFrom1D(double rLocalCoordinates, double rGlobalDisplacements[3]) const;

protected:
    //! @brief ... check if the element is properly defined (check node dofs, nodes are reordered if the element length is negative)
    void CheckElement();

    //! @brief ... extract global dofs from nodes (mapping of local row ordering of the element matrices to the global dof ordering)
    //! @param rGlobalRowDofs ... vector of global row dofs
    void CalculateGlobalRowDofs(std::vector<int>& rGlobalRowDofs) const;

    //! @brief ... extract global dofs from nodes (mapping of local column ordering of the element matrices to the global dof ordering)
    //! @param rGlobalColumnDofs ... vector of global column dofs
    void CalculateGlobalColumnDofs(std::vector<int>& rGlobalColumnDofs) const;


};
} // namespace NuTo
#endif //Truss1D_H
