// $Id:$
#ifndef LATTICE2D_H
#define LATTICE2D_H

#ifdef ENABLE_SERIALIZATION
#include "boost/array.hpp"
#endif

#include "nuto/mechanics/elements/ElementBase.h"
#include <eigen3/Eigen/Core>

namespace NuTo
{
class StructureBase;
class DeformationGradient2D;
class ConstitutiveTangentLocal3x3;
class EngineeringStress2D;
class NodeBase;
class LatticeStrain2D;
//! @author Jörg F. Unger, ISM
//! @date March 2010
//! @brief ... lattice model in 2D (triangle edges based)
class Lattice2D : public ElementBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION

public:
    //! @brief constructor
    Lattice2D(const StructureBase* rStructure,
    		std::vector<NuTo::NodeBase* >& rNodes,
    		ElementData::eElementDataType rElementDataType,
    		IpData::eIpDataType rIpDataType
    		);

    //! @brief returns the global dimension of the element
    //! this is required to check, if an element can be used in a 1d, 2D or 3D Structure
    //! there is also a routine GetLocalDimension, which is e.g. 2 for PlaneLattice elements and 1 for truss elements
    //! @return global dimension
    int GetGlobalDimension()const
    {
        return 3;
    }

    //! @brief calculates the coefficient matrix for the 0-th derivative in the differential equation
    //! for a mechanical problem, this corresponds to the stiffness matrix
    //! @param rResult ... coefficient matrix
    //! @param rGlobalDofsRow ... row numbers in global system
    //! @param rGlobalDofsColumn ... column numbers in global system
    //! @param rSymmetry ... matrix is symmetric or not (in the symmetric case the full matrix is also stored
    Error::eError CalculateCoefficientMatrix_0(NuTo::FullMatrix<double>& rResult,
            std::vector<int>& rGlobalDofsRow, std::vector<int>& rGlobalDofsColumn, bool& rSymmetry)const;

    //! @brief calculates the coefficient matrix for the 1-th derivative in the differential equation
    //! for a mechanical problem, this corresponds to the damping matrix
    Error::eError CalculateCoefficientMatrix_1(NuTo::FullMatrix<double>& rResult,
            std::vector<int>& rGlobalDofsRow, std::vector<int>& rGlobalDofsColumn)const;

    //! @brief calculates the coefficient matrix for the 2-th derivative in the differential equation
    //! for a mechanical problem, this corresponds to the Mass matrix
    Error::eError CalculateCoefficientMatrix_2(NuTo::FullMatrix<double>& rResult,
            std::vector<int>& rGlobalDofsRow, std::vector<int>& rGlobalDofsColumn)const;

    //! @brief calculates the gradient of the internal potential
    //! for a mechanical problem, this corresponds to the internal force vector
    Error::eError CalculateGradientInternalPotential(NuTo::FullMatrix<double>& rResult,
                                            std::vector<int>& rGlobalDofs)const;

    //! @brief sets the section of an element
    //! implemented with an exception for all elements, reimplementation required for those elements
    //! which actually need a section
    //! @param rSection pointer to section
    //! @return pointer to constitutive law
    void SetSection(const SectionBase* rSection);

    //! @brief returns a pointer to the section of an element
    //! implemented with an exception for all elements, reimplementation required for those elements
    //! which actually need a section
    //! @return pointer to section
    const SectionBase* GetSection()const;

    //! @brief Update the static data of an element
    Error::eError UpdateStaticData(NuTo::Element::eUpdateType rUpdateType);

    //! @brief calculates the shape functions
    //! @param rLocalCoordinates local coordinates of the integration point
    //! @param shape functions for all the nodes
    void CalculateShapeFunctionsCoordinates(const boost::array<double,2>& rNaturalCoordinates, boost::array<double,3>& rShapeFunctions)const;

    //! @brief calculates the shape functions for rotations
    //! @param rGlobalCoordinatesFacet coordinates of the integration point
    //! @param rCoordNode coordinate of the node
    //! @return shape functions for the node
    boost::array<double,2> CalculateShapeFunctionsRotations(const boost::array<double,3>& rGlobalCoordinatesFacet,const boost::array<double,2>& rCoordNode)const;

    //! @brief calculates the integration point data with the current displacements applied
    //! @param rIpDataType data type to be stored for each integration point
    //! @param rIpData return value with dimension (dim of data type) x (numIp)
    Error::eError GetIpData(NuTo::IpData::eIpStaticDataType rIpDataType, FullMatrix<double>& rIpData)const;

    //! @brief Allocates static data for an integration point of an element
    //! @param rConstitutiveLaw constitutive law, which is called to allocate the static data object
    ConstitutiveStaticDataBase* AllocateStaticData(const ConstitutiveBase* rConstitutiveLaw)const;

    //! @brief ... interpolate three-dimensional global point coordinates from two-dimensional local point coordinates (element coordinates system)
    //! @param rNaturalCoordinates ... two-dimensional point coordinates in natural coordinate system
    //! @param rGlobalCoordinates ... three-dimension global point coordinates
    void InterpolateCoordinatesFrom2D(double rNaturalCoordinates[2], double rGlobalCoordinates[3]) const;

    // interpolate geometry
    //! @param rNaturalCoordinates ... two-dimensional point coordinates in natural coordinate system
    //! @param rGlobalCoordinates ... three-dimension global point coordinates
    void InterpolateCoordinatesFrom2D(const boost::array<double, 2> rNaturalCoordinates, boost::array<double, 3>& rGlobalCoordinates) const;

    //! @brief ... interpolate three-dimensional global point displacements from three-dimensional local point coordinates (element coordinates system)
    //! @param rNaturalCoordinates ... two-dimensional point coordinates in natural coordinate system
    //! @param rGlobalDisplacements ... three-dimension global point displacements
    void InterpolateDisplacementsFrom2D(double rNaturalCoordinates[2], double rGlobalDisplacements[3]) const;

    // interpolate displacements
    void InterpolateDisplacementsFrom2D(const boost::array<double, 2> rNaturalCoordinates,
    		boost::array<double, 3>& rGlobalDisplacements) const;

    //! @brief returns the enum (type of the element)
    //! @return enum
    NuTo::Element::eElementType GetEnumType()const
    {
    	return NuTo::Element::LATTICE2D;
    }

    //! @brief returns the number of nodes in this element
    //! @return number of nodes
    int GetNumNodes()const;

    //! @brief returns a pointer to the i-th node of the element
    //! @param local node number
    //! @return pointer to the node
    NuTo::NodeBase* GetNode(int rLocalNodeNumber);

    //! @brief returns a pointer to the i-th node of the element
    //! @param local node number
    //! @return pointer to the node
    const NuTo::NodeBase* GetNode(int rLocalNodeNumber)const;

    //! @brief sets the rLocalNodeNumber-th node of the element
    //! @param local node number
    //! @param pointer to the node
    void SetNode(int rLocalNodeNumber, NuTo::NodeBase* rNode);

    //! brief exchanges the node ptr in the full data set (elements, groups, loads, constraints etc.)
    //! this routine is used, if e.g. the data type of a node has changed, but the restraints, elements etc. are still identical
    //! @param rOldPtr old node ptr
    //! @param rNewPtr new node ptr
    void ExchangeNodePtr(NodeBase* rOldPtr, NodeBase* rNewPtr);

    // reorder nodes such that the sign of the length/area/volume of the element changes
    void ReorderNodes();

    //! @brief calculates the volume of an integration point (weight * detJac)
    //! @param rVolume  vector for storage of the ip volumes (area in 2D)
    void GetIntegrationPointVolume(std::vector<double>& rVolume)const;

    //! @brief returns the coordinates of an integration point
    //! @param rIpNum integration point
    //! @param rCoordinates coordinates to be returned
    void  GetGlobalIntegrationPointCoordinates(int rIpNum, double rCoordinates[3])const;

    //! @brief cast the base pointer to an PlaneLattice, otherwise throws an exception
    const Lattice2D* AsLattice2D()const;

    //! @brief cast the base pointer to an PlaneLattice, otherwise throws an exception
    Lattice2D* AsLattice2D();

#ifdef ENABLE_VISUALIZE
    void GetVisualizationCells(
        unsigned int& NumVisualizationPoints,
        std::vector<double>& VisualizationPointLocalCoordinates,
        unsigned int& NumVisualizationCells,
        std::vector<NuTo::CellBase::eCellTypes>& VisualizationCellType,
        std::vector<unsigned int>& VisualizationCellsIncidence,
        std::vector<unsigned int>& VisualizationCellsIP) const;
#endif// ENABLE_VISUALIZE

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif  // ENABLE_SERIALIZATION

protected:
    //! @brief ... just for serialization
    Lattice2D(){};

    //! @brief ... check if the element is properly defined (check node dofs, nodes are reordered if the element volume is negative)
    void CheckElement();

    //! @brief adds the stiffness of one lattice (one edge of the triangle to the element matrix
    void AddLatticeStiffnessToElementMatrix(Eigen::Matrix<double,6,6>& latticeStiffness,
    		int rLocalNode1, int rLocalNode2,NuTo::FullMatrix<double>& rCoefficientMatrix)const;

    //! @brief calculates the local B-matrix for an edge connecting two nodes (with one ip)
    void CalculateBMatrixAndLatticeStrain(int rIP, Eigen::Matrix<double,2,6>& rBMatrix, LatticeStrain2D& rLatticeStrain)const;

    //! @brief calculate the length of an edge (belonging to an integration point
    //! @param rIp integration point
    //! @return edge length
    double GetIpEdgeLength(int rIp)const;

    void CalculateGlobalCoordinatesEdgeFacePoints(boost::array<boost::array<double, 2>,3 >& rCoordinatesEdgePoints,
    		boost::array<double, 2>& rCoordinateFacePoint)const;

    void CalculateHelpTransformationGlobalLocal(
    		Eigen::Matrix<double,2,2>& rAinv,
    		boost::array<double, 2>& rCoordinatesNode0)const;

    void CalculateLocalCoordinateFromGlobal(
    		const Eigen::Matrix<double,2,2>& Ainv,
    		const boost::array<double, 2>& rCoordinatesNode0,
    		const boost::array<double, 2>& rGlobalCoordinates,
    		boost::array<double, 2>& rLocalCoordinates)const;

    //! @brief ... extract global dofs from nodes (mapping of local row ordering of the element matrices to the global dof ordering)
    //! @param rGlobalRowDofs ... vector of global row dofs
    void CalculateGlobalRowDofs(std::vector<int>& rGlobalRowDofs) const;

    //! @brief ... extract global dofs from nodes (mapping of local column ordering of the element matrices to the global dof ordering)
    //! @param rGlobalColumnDofs ... vector of global column dofs
    void CalculateGlobalColumnDofs(std::vector<int>& rGlobalColumnDofs) const;

    const SectionBase *mSection;

    //! @brief element nodes
    NodeBase* mNodes[3];
};

} // namespace NuTo

#endif //LATTICE2D_H
