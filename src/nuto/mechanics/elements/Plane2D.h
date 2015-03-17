// $Id$
#ifndef PLANE2D_H
#define PLANE2D_H

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
    Plane2D(const NuTo::StructureBase* rStructure,
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

    //! @brief calculates the nonlocal eq strains of the nodes
    //! @param rNodeNonlocalEqStrains vector with already correct size allocated
    //! this can be checked with an assertation
    void CalculateLocalNonlocalEqStrains(std::vector<double>& rNodeNonlocalEqStrains)const;

    //! @brief transforms the local matrix to the global system
    //! relevant only for 2D and 3D truss elements
    void BlowLocalMatrixToGlobal(NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& FullrCoefficientMatrix)const
    {}

    //! @brief transforms the local vector to the global system
    //! relevant only for 2D and 3D truss elements
    void BlowLocalVectorToGlobal(NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rFullVector)const
    {}

    //! @brief calculates the area of a plane element via the nodes (probably faster than sum over integration points)
    //! @return Area
    double CalculateArea()const;

    //! @brief checks if a node is inside a polygon
    //! @param rPoint (input) ... a pointer to a 2D tuple containing the coordinates
    //! @param rPoint (input) ... a pointer to a vector of 2D tuples representing the polyline
    bool CheckPointInsidePolygon( const std::tuple<double,double> *rPoint, const std::vector<std::tuple<double,double> > * rPolygon)const;

    //! @brief calculate the length of a crack passing through the center of gravity and intersecting two edges
    //! @parameter alpha... angle of the crack
    virtual double CalculateCrackLength2D(double rAlpha)const;

    //! @brief calculates the shape functions for the surfaces (required for surface loads)
    //! @param rLocalCoordinates local coordinates of the integration point (in the local surface coordinate system)
    //! @param shape functions for all the nodes, size should already be correct, but can be checked with an assert
    virtual void CalculateShapeFunctionsSurface(double rLocalCoordinates, std::vector<double>& rShapeFunctions)const=0;

    //! @brief calculates the derivative of the shape functions with respect to local coordinates for the surfaces (required for surface loads)
    //! @param rLocalCoordinates local coordinates of the integration point
    //! @param derivative of the shape functions for all the nodes, size should already be correct, but can be checked with an assert
    //! first all the directions for a single node, and then for the next node
    virtual void CalculateDerivativeShapeFunctionsLocalSurface(double rLocalCoordinates, std::vector<double>& rDerivativeShapeFunctions)const=0;

    //! @brief returns the surface nodes
    //! @param surface (numbering so that the normal (right hand /thumb rule) is pointing outwards)
    //! @param surface nodes
    virtual void GetSurfaceNodes(int rSurface, std::vector<const NodeBase*>& rSurfaceNodes)const=0;

    //! @brief returns the number of external surfaces
    //! @param surface (numbering so that the normal (right hand /thumb rule) is pointing outwards)
    //! @param surface nodes
    virtual int GetNumSurfaces()const=0;

    //! @brief cast the base pointer to an Plane, otherwise throws an exception
    const Plane2D* AsPlane2D()const override;

    //! @brief cast the base pointer to an Plane, otherwise throws an exception
    Plane2D* AsPlane2D() override;

protected:
    //! @brief ... just for serialization
    Plane2D(){};

    //! @brief ... extract global dofs from nodes (mapping of local row ordering of the element matrices to the global dof ordering)
    //! @param rGlobalRowDofs ... vector of global row dofs
    //! @param rNumDispDofs number of displacement dofs
    //! @param rNumTempDofs number of temperature dofs
    //! @param rNumNonlocalEqStrainDofs number of nonlocal eq strain dofs
    void CalculateGlobalRowDofs(std::vector<int>& rGlobalRowDofs, int rNumDispDofs, int rNumTempDofs, int rNumNonlocalEqStrainDofs) const override;

    //! @brief ... extract global dofs from nodes (mapping of local column ordering of the element matrices to the global dof ordering)
    //! @param rGlobalColumnDofs ... vector of global column dofs
    //! @param rNumDispDofs number of displacement dofs
    //! @param rNumTempDofs number of temperature dofs
    //! @param rNumNonlocalEqStrainDofs number of nonlocal eq strain dofs
    void CalculateGlobalColumnDofs(std::vector<int>& rGlobalColumnDofs, int rNumDispDofs, int rNumTempDofs, int rNumNonlocalEqStrainDofs) const override;


};
} // namespace NuTo
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::Plane2D)
#endif // ENABLE_SERIALIZATION
#endif //PLANE2D_H
