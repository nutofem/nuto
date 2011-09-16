// $Id$
#ifndef ELEMENT_BASE_H
#define ELEMENT_BASE_H

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#include <boost/serialization/vector.hpp>
#else
#include <vector>
#endif //ENABLE_SERIALIZATION

#include "nuto/base/ErrorEnum.h"
#include "nuto/mechanics/elements/ElementDataEnum.h"
#include "nuto/mechanics/elements/ElementEnum.h"
#include "nuto/mechanics/elements/IpDataEnum.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeEnum.h"

#ifdef ENABLE_VISUALIZE
#include "nuto/visualize/VisualizeBase.h"
#include "nuto/visualize/VisualizeUnstructuredGrid.h"
#include <boost/ptr_container/ptr_list.hpp>
#endif // ENABLE_VISUALIZE

namespace NuTo
{
class CrackBase;
class ConstitutiveBase;
class ConstitutiveStaticDataBase;
class ElementDataBase;
class IntegrationTypeBase;
class NodeBase;
class Plane;
class SectionBase;
template<class T>
class SparseMatrix;
template<class T>
class FullMatrix;
class Solid;
class Structure;
class StructureBase;
class Truss;
class VisualizeComponentBase;
class IpDataBase;


//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... standard abstract class for all groups
class ElementBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
    friend class NuTo::Structure;

public:
    //! @brief constructor
    //! @param rStructure ... structure to which the element belongs
    ElementBase(const StructureBase* rStructure, ElementData::eElementDataType rElementDataType,
    		IntegrationType::eIntegrationType rIntegrationType, IpData::eIpDataType rIpDataType);

    virtual ~ElementBase();

    //! @todo copy-constructor

    //! @brief returns the enum (type of the element)
    //! @return enum
    virtual NuTo::Element::eElementType GetEnumType()const=0;

    //! @brief returns the enum of element data type
    //! @return enum of ElementDataType
    const NuTo::ElementData::eElementDataType GetElementDataType()const;

    //! @brief returns the id number of the element
    //! @return id
    int ElementGetId()const;

    //! @brief returns the global dimension of the element
    //! this is required to check, if an element can be used in a 1d, 2D or 3D Structure
    //! there is also a routine GetLocalDimension, which is e.g. 2 for plane elements and 1 for truss elements
    //! @return global dimension
    virtual int GetGlobalDimension()const=0;

    //! @brief returns the number of nodes in this element
    //! @return number of nodes
    virtual int GetNumNodes()const=0;

    //! @brief returns a pointer to the i-th node of the element
    //! @param local node number
    //! @return pointer to the node
    virtual NodeBase* GetNode(int rLocalNodeNumber)=0;

    //! @brief sets the rLocalNodeNumber-th node of the element
    //! @param local node number
    //! @param pointer to the node
    virtual void SetNode(int rLocalNodeNumber, NodeBase* rNode)=0;

    //! @brief returns a pointer to the i-th node of the element
    //! @param local node number
    //! @return pointer to the node
    virtual const NodeBase* GetNode(int rLocalNodeNumber)const=0;

    //! brief exchanges the node ptr in the full data set (elements, groups, loads, constraints etc.)
    //! this routine is used, if e.g. the data type of a node has changed, but the restraints, elements etc. are still identical
    virtual void ExchangeNodePtr(NodeBase* rOldPtr, NodeBase* rNewPtr)=0;

    //! @brief sets the constitutive law for an element
    //! @param rConstitutiveLaw Pointer to constitutive law entry
    virtual void SetConstitutiveLaw(ConstitutiveBase* rConstitutiveLaw);

    //! @brief returns a pointer to the constitutive law for an integration point
    //! @param integration point number (counting from zero)
    //! @return pointer to constitutive law
    const ConstitutiveBase* GetConstitutiveLaw(int rIp)const;

    //! @brief returns a pointer to the constitutive law for an integration point
    //! @param integration point number (counting from zero)
    //! @return pointer to constitutive law
    ConstitutiveBase* GetConstitutiveLaw(int rIp);

    //! @brief sets the fine scale model (deserialization from a binary file)
    void SetFineScaleModel(int rIp, std::string rFileName, double rLengthCoarseScale, std::string rIPName);

    //! @brief sets the fine scale parameter for all ips
    //! @parameter rName name of the parameter, e.g. YoungsModulus
    //! @parameter rParameter value of the parameter
    void SetFineScaleParameter(int rIp, const std::string& rName, double rParameter);

    //! @brief sets the fine scale parameter for all ips
    //! @parameter rName name of the parameter, e.g. YoungsModulus
    //! @parameter rParameter value of the parameter
    void SetFineScaleParameter(int rIp, const std::string& rName, std::string rParameter);

    //! @brief sets the section of an element
    //! implemented with an exception for all elements, reimplementation required for those elements
    //! which actually need a section
    //! @param rSection pointer to section
    virtual void SetSection(const SectionBase* rSection);

    //! @brief returns a pointer to the section of an element
    //! implemented with an exception for all elements, reimplementation required for those elements
    //! which actually need a section
    //! @return pointer to section
    virtual const SectionBase* GetSection()const;

    //! @brief sets the integration type of an element
    //! implemented with an exception for all elements, reimplementation required for those elements
    //! which actually need an integration type
    //! @param rIntegrationType pointer to integration type
    virtual void SetIntegrationType(const IntegrationTypeBase* rIntegrationType, NuTo::IpData::eIpDataType rIpDataType);

    //! @brief returns a pointer to the integration type of an element
    //! implemented with an exception for all elements, reimplementation required for those elements
    //! which actually need an integration type
    //! @return pointer to integration type
    virtual const IntegrationTypeBase* GetIntegrationType()const;

    //! @brief returns ip data type of the element
    //! implemented with an exception for all elements, reimplementation required for those elements
    //! which actually need an integration type
    //! @return pointer to integration type
    NuTo::IpData::eIpDataType GetIpDataType(int  rIp)const;

    //! @brief returns the number of integration points
    //! @return number of integration points
    int GetNumIntegrationPoints()const;

    //! @brief returns the weight of an integration point
    //! @param rIpNum integration point
    //! @return weight
    double GetIntegrationPointWeight(int rIpNum)const;

    //! @brief calculates the coefficient matrix for the 0-th derivative in the differential equation
    //! for a mechanical problem, this corresponds to the stiffness matrix
    //! @param rResult ... coefficient matrix
    //! @param rGlobalDofsRow ... row numbers in global system
    //! @param rGlobalDofsColumn ... column numbers in global system
    //! @param rSymmetry ... matrix is symmetric or not (in the symmetric case the full matrix is also stored
    virtual Error::eError CalculateCoefficientMatrix_0(NuTo::FullMatrix<double>& rResult,
            std::vector<int>& rGlobalDofsRow, std::vector<int>& rGlobalDofsColumn, bool& rSymmetry)const=0;

    //! @brief calculates the coefficient matrix for the 1-th derivative in the differential equation
    //! for a mechanical problem, this corresponds to the damping matrix
    virtual Error::eError CalculateCoefficientMatrix_1(NuTo::FullMatrix<double>& rResult,
            std::vector<int>& rGlobalDofsRow, std::vector<int>& rGlobalDofsColumn)const=0;

    //! @brief calculates the coefficient matrix for the 2-th derivative in the differential equation
    //! for a mechanical problem, this corresponds to the Mass matrix
    virtual Error::eError CalculateCoefficientMatrix_2(NuTo::FullMatrix<double>& rResult,
            std::vector<int>& rGlobalDofsRow, std::vector<int>& rGlobalDofsColumn)const=0;

    //! @brief calculates the gradient of the internal potential
    //! for a mechanical problem, this corresponds to the internal force vector
    virtual Error::eError CalculateGradientInternalPotential(NuTo::FullMatrix<double>& rResult,
            std::vector<int>& rGlobalDofs)const=0;

    //! @brief calculates the integration point data with the current displacements applied
    //! @param rIpDataType data type to be stored for each integration point
    //! @param rIpData return value with dimension (dim of data type) x (numIp)
    virtual Error::eError GetIpData(NuTo::IpData::eIpStaticDataType rIpDataType, FullMatrix<double>& rIpData)const=0;

    //! @brief integrates the stress over the element
    //! @param rStress integrated stress
    void GetIntegratedStress(FullMatrix<double>& rStress)const;

    //! @brief integrates the strain over the element
    //! @param rStrain integrated strain
    void GetIntegratedStrain(FullMatrix<double>& rStress)const;

    //! @brief Allocates static data for an integration point of an element
    //! @param rConstitutiveLaw constitutive law, which is called to allocate the static data object
    virtual  ConstitutiveStaticDataBase* AllocateStaticData(const ConstitutiveBase* rConstitutiveLaw)const=0;

    //! @brief Returns the static data for an integration point of an element
    //! @param rIp integration point
    //! @return static data
    ConstitutiveStaticDataBase* GetStaticData(int rIp);

    //! @brief Returns the static data for an integration point of an element
    //! @param rIp integration point
    //! @return static data
    const ConstitutiveStaticDataBase* GetStaticData(int rIp)const;

    //! @brief sets the static data for an integration point of an element
    //! @param rIp integration point
    //! @param rStaticData static data
    void SetStaticData(int rIp, ConstitutiveStaticDataBase* rStaticData);

    //! @brief Update the static data of an element
    virtual Error::eError UpdateStaticData(NuTo::Element::eUpdateType rUpdateType)=0;

    //! @brief ... interpolate three-dimensional global point coordinates from one-dimensional local point coordinates (element coordinates system)
    //! @param rLocalCoordinates ... one-dimensional local point coordinates
    //! @param rGlobalCoordinates ... three-dimension global point coordinates
    virtual void InterpolateCoordinatesFrom1D(double rLocalCoordinates, double rGlobalCoordinates[3]) const;

    //! @brief ... interpolate three-dimensional global point coordinates from two-dimensional local point coordinates (element coordinates system)
    //! @param rLocalCoordinates ... two-dimensional local point coordinates
    //! @param rGlobalCoordinates ... three-dimension global point coordinates
    virtual void InterpolateCoordinatesFrom2D(double rLocalCoordinates[2], double rGlobalCoordinates[3]) const;

    //! @brief ... interpolate three-dimensional global point coordinates from three-dimensional local point coordinates (element coordinates system)
    //! @param rLocalCoordinates ... three-dimensional local point coordinates
    //! @param rGlobalCoordinates ... three-dimension global point coordinates
    virtual void InterpolateCoordinatesFrom3D(double rLocalCoordinates[3], double rGlobalCoordinates[3]) const;

    //! @brief ... interpolate three-dimensional global point displacements from one-dimensional local point coordinates (element coordinates system)
    //! @param rLocalCoordinates ... one-dimensional local point coordinates
    //! @param rGlobalDisplacements ... three-dimension global point displacements
    virtual void InterpolateDisplacementsFrom1D(double rLocalCoordinates, double rGlobalDisplacements[3]) const;

    //! @brief ... interpolate three-dimensional global point displacements from two-dimensional local point coordinates (element coordinates system)
    //! @param rLocalCoordinates ... two-dimensional local point coordinates
    //! @param rGlobalDisplacements ... three-dimension global point displacements
    virtual void InterpolateDisplacementsFrom2D(double rLocalCoordinates[2], double rGlobalDisplacements[3]) const;

    //! @brief ... interpolate three-dimensional global point displacements from three-dimensional local point coordinates (element coordinates system)
    //! @param rLocalCoordinates ... three-dimensional local point coordinates
    //! @param rGlobalDisplacements ... three-dimension global point displacements
    virtual void InterpolateDisplacementsFrom3D(double rLocalCoordinates[3], double rGlobalDisplacements[3]) const;

    //! @brief adds the nonlocal weight to an integration point
    //! @param rLocalIpNumber local Ip
    //! @param rNonlocalElement element of the nonlocal ip
    //! @param rNonlocalIp local ip number of the nonlocal ip
    //! @param rWeight weight
    void SetNonlocalWeight(int rLocalIpNumber, const ElementBase* rNonlocalElement, int rNonlocalIp, double rWeight);

    //! @brief returns a vector of weights for an ip
    //! @param rIp local Ip
    //! @param rNonlocalElement nonlocal element (must be in the range of the nonlocal element size stored at the element data level)
    //! @retrun weights for each integration point of the nonlocal element
    const std::vector<double>& GetNonlocalWeights(int rIp, int rNonlocalElement)const;

    //! @brief returns a vector of the nonlocal elements
    //! @retrun nonlocal elements
    const std::vector<const NuTo::ElementBase*>& GetNonlocalElements()const;

    //! @brief returns the number of nonlocal elements
    //! @param rConstitutive constitutive model for the nonlocale elements
    //! @rerun number of nonlocal elements
    int GetNumNonlocalElements()const;

    //! @brief delete the nonlocal elements
    void DeleteNonlocalElements();

    //! @brief calculates the area of a plane element via the nodes (probably faster than sum over integration points)
    //! @return Area
    virtual double CalculateArea()const;

    //! @brief calculates the volume of an integration point (weight * detJac)
    //! @param rVolume  vector for storage of the ip volumes (area in 2D, length in 1D)
    virtual void GetIntegrationPointVolume(std::vector<double>& rVolume)const=0;

    //! @brief returns the coordinates of an integration point
    //! @param rIpNum integration point
    //! @param rCoordinates coordinates to be returned
    virtual void GetGlobalIntegrationPointCoordinates(int rIpNum, double rCoordinates[3])const=0;

    //! @brief computes the natural coordinates of an given point
    //! implemented with an exception for all elements, reimplementation required for those elements
    //! @param rGlobCoords (input) ... pointer to the array of coordinates
    //! @param rLocCoords (output) ... coordinates to be returned
    //! @return True if coordinates are within the element, False otherwise
    virtual bool GetLocalPointCoordinates(const double* rGlobCoords,  double* rLocCoords)const;

    //! @brief checks if a node is inside of an element
    //! implemented with an exception for all elements, reimplementation required for those elements
    //! @param rGlobCoords (input) ... pointer to the array of coordinates
    //! @return True if coordinates are within the element, False otherwise
    virtual bool CheckPointInside(const double* rGlobCoords)const;

    //! @brief Returns the vector of crack pointers of an element
    //! @return crack pointer vector
    virtual const std::vector<CrackBase*>  GetCracks() const;

    //! @brief Set the information that the element is already cracked or not
    //! @param bool (Input) cracked or not
    void IsCracked(const bool rIsCracked);

    //! @brief Give the information if the element is already cracked or not
    //! @return bool cracked or not
    const bool IsCracked() const;

    //! @brief cast the base pointer to an ElementPlane, otherwise throws an exception
    virtual const Plane* AsPlane()const;

    //! @brief cast the base pointer to an ElementPlane, otherwise throws an exception
    virtual Plane* AsPlane();

    //! @brief cast the base pointer to an ElementSolid, otherwise throws an exception
    virtual const Solid* AsSolid()const;

    //! @brief cast the base pointer to an ElementSolid, otherwise throws an exception
    virtual Solid* AsSolid();

    //! @brief cast the base pointer to an ElementTruss, otherwise throws an exception
    virtual const Truss* AsTruss()const;

    //! @brief cast the base pointer to an ElementTruss, otherwise throws an exception
    virtual Truss* AsTruss();

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif  // ENABLE_SERIALIZATION

#ifdef ENABLE_VISUALIZE
    virtual void Visualize(VisualizeUnstructuredGrid& rVisualize, const boost::ptr_list<NuTo::VisualizeComponentBase>& rWhat) const;

    //Visualize for all integration points the fine scale structure
    void VisualizeIpMultiscale(VisualizeUnstructuredGrid& rVisualize, const boost::ptr_list<NuTo::VisualizeComponentBase>& rWhat, bool rVisualizeDamage)const;

#endif // ENABLE_VISUALIZE

    //! @brief returns the structure
    const StructureBase* GetStructure()const
    {
    	return mStructure;
    }

private:
    //! @brief returns the Element Data Vector
    //! this was necessary due to recursive problems for serialization (nonlocal data)
    //! this method should only be called from the serialization routine of the structure
    NuTo::ElementDataBase* GetDataPtr()const;

    //! @brief sets the Element Data Vector
    //! this was necessary due to recursive problems for serialization (nonlocal data)
    //! this method should only be called from the serialization routine of the structure
    void SetDataPtr(NuTo::ElementDataBase* rElementData);

protected:
    //! @brief ... just for serialization
    ElementBase()
    {
        mElementData = 0;
    };

    //! @brief ... reorder nodes such that the sign of the length/area/volume of the element changes
    virtual void ReorderNodes() = 0;

    //! @brief ... check if the element is properly defined (check node dofs, nodes are reordered if the element length/area/volum is negative)
    virtual void CheckElement() = 0;

    //! @brief ... extract global dofs from nodes (mapping of local row ordering of the element matrices to the global dof ordering)
    //! @param rGlobalRowDofs ... vector of global row dofs
    virtual void CalculateGlobalRowDofs(std::vector<int>& rGlobalRowDofs) const = 0;

    //! @brief ... extract global dofs from nodes (mapping of local column ordering of the element matrices to the global dof ordering)
    //! @param rGlobalColumnDofs ... vector of global column dofs
    virtual void CalculateGlobalColumnDofs(std::vector<int>& rGlobalColumnDofs) const = 0;

    //the base class of the elements must not contain any data apart from a const pointer to the structure and a data pointer
    const StructureBase* mStructure;

    //the base class of the elements data
    ElementDataBase *mElementData;
};
}//namespace NuTo
#endif //ELEMENT_BASE_H

