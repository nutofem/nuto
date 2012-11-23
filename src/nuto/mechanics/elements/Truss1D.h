// $Id$
#ifndef Truss1D_H
#define Truss1D_H

#include "nuto/mechanics/elements/Truss.h"
#include <iostream>
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
    void serialize(Archive & ar, const unsigned int version);
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
    //! @param time derivative (0 temperature, 1 temperature rate, 2 second time derivative of temperature)
    //! @param localDisplacements vector with already correct size allocated
    //! this can be checked with an assertation
    void CalculateLocalDisplacements(int rTimeDerivative, std::vector<double>& rLocalDisplacements)const;

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
    //! @brief ... just for serialization
    Truss1D(){}

    //! @brief ... check if the element is properly defined (check node dofs, nodes are reordered if the element length is negative)
    void CheckElement();

    //! @brief ... extract global dofs from nodes (mapping of local row ordering of the element matrices to the global dof ordering)
    //! @param rGlobalRowDofs ... vector of global row dofs
    void CalculateGlobalRowDofs(std::vector<int>& rGlobalRowDofs,int rNumDispDofs, int rNumTempDofs) const;

    //! @brief ... extract global dofs from nodes (mapping of local column ordering of the element matrices to the global dof ordering)
    //! @param rGlobalColumnDofs ... vector of global column dofs
    void CalculateGlobalColumnDofs(std::vector<int>& rGlobalColumnDofs,int rNumDispDofs, int rNumTempDofs) const;


};
} // namespace NuTo
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::Truss1D)
#endif // ENABLE_SERIALIZATION
#endif //Truss1D_H
