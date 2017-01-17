#pragma once

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#include <boost/serialization/vector.hpp>
#else
#include <vector>
#endif //ENABLE_SERIALIZATION

#include <list>
#include <map>
#include <eigen3/Eigen/Core>
#include "nuto/mechanics/MechanicsException.h"
#include <memory>


namespace NuTo
{
class BoundaryGradientDamage1D;
class ConstitutiveBase;
class ConstitutiveStaticDataBase;
class ElementDataBase;
class IntegrationTypeBase;
class InterpolationType;
class NodeBase;
class Lattice2D;
class SectionBase;
template<class T>
class SparseMatrix;
class Structure;
class StructureBase;
class Truss;
class VisualizeComponentBase;
class VisualizeComponent;
class IpDataBase;
class ElementOutputBase;
class IpDataStaticDataBase;
enum class eError;
template<typename IOEnum> class ConstitutiveIOMap;
template <int TDim> class ContinuumElement;
template <int TDim> class ContinuumElementIGA;
template <int TDim> class ContinuumElementIGALayer;
template <int TDim> class ContinuumBoundaryElement;
template <int TDimSlave, int TDimMaster> class ContinuumContactElement;
template <class T, int rows, int cols> class FullMatrix;
template <class T, int rows> class FullVector;

#ifdef ENABLE_VISUALIZE
class VisualizeUnstructuredGrid;

class CellBase;
enum class eCellTypes;
#endif // ENABLE_VISUALIZE

namespace Constitutive
{
    enum class eInput;
    enum class eOutput;
}

using ConstitutiveInputMap = ConstitutiveIOMap<Constitutive::eInput>;
using ConstitutiveOutputMap = ConstitutiveIOMap<Constitutive::eOutput>;

namespace Element
{
    enum class eElementType;
    enum class eOutput;
}// namespace Element

namespace ElementData
{
    enum class eElementDataType;
}// namespace ElementData

namespace IpData
{
    enum class eIpDataType;
}// namespace IpData

namespace Node
{
    enum class eDof : unsigned char;
}// namespace Node

//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... standard abstract class for all elements
class ElementBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
    friend class NuTo::Structure;

public:
    //! @brief constructor
    //! @param rStructure ... structure to which the element belongs
    //! @param rElementDataType ... element data type
    //! @param rIntegrationType ... integration type (local coordinates are stored as a type, e.g. Gauss 2x2
    //! @param rIpDataType ... data type to decide what is stored at the integration point level
    //! @param rInterpolationType ... interpolation type
    ElementBase(const StructureBase* rStructure, ElementData::eElementDataType rElementDataType, IpData::eIpDataType rIpDataType, const InterpolationType* rInterpolationType);

    //! @brief constructor
    //! @param rStructure ... structure to which the element belongs
    //! @param rElementDataType ... element data type
    //! @param rNumIp ... number of integration points (here local coordinates should be stored at the ip (e.g. XFEM)
    //! @param rIpDataType ... data type to decide what is stored at the integration point level
    //! @param rInterpolationType ... interpolation type
    ElementBase(const StructureBase* rStructure, ElementData::eElementDataType rElementDataType, int rNumIp, IpData::eIpDataType rIpDataType, const InterpolationType* rInterpolationType);

    ElementBase(const ElementBase&) = default;

    virtual ~ElementBase();

    //! @brief returns the enum (type of the element)
    //! @return enum
    virtual NuTo::Element::eElementType GetEnumType() const=0;

    //! @brief returns the enum of element data type
    //! @return enum of ElementDataType
    const NuTo::ElementData::eElementDataType GetElementDataType() const;

    //! @brief returns the id number of the element
    //! @return id
    int ElementGetId() const;

    //! @brief returns the local dimension of the element
    //! this is required to check, if an element can be used in a 1d, 2D or 3D Structure
    //! @return local dimension
    virtual int GetLocalDimension() const=0;

    //! @brief returns the number of nodes in this element
    //! @return number of nodes
    virtual int GetNumNodes() const;

    //! @brief returns the number of nodes in this element that are influenced by it
    //! @remark overridden by boundary elements
    //! @return number of nodes
    virtual int GetNumInfluenceNodes() const
    {
        return GetNumNodes();
    }

    //! @brief returns a pointer to the i-th node of the element
    //! @param local node number
    //! @return pointer to the node
    virtual NodeBase* GetNode(int rLocalNodeNumber)=0;

    //! @brief returns a pointer to the i-th node of the element
    //! @param local node number
    //! @return pointer to the node
    virtual const NodeBase* GetNode(int rLocalNodeNumber) const=0;

    //! @brief returns a pointer to the i-th node of the element
    //! @remark overridden by boundary elements
    //! @param local node number
    //! @return pointer to the node
    virtual const NodeBase* GetInfluenceNode(int rLocalNodeNumber) const
    {
        return GetNode(rLocalNodeNumber);
    }

    //! @brief returns the number of nodes in this element of a specific dof
    //! @brief rDofType dof type
    //! @return number of nodes
    virtual int GetNumNodes(Node::eDof rDofType) const;

    //! @brief returns a pointer to the i-th node of the element
    //! @param local node number
    //! @brief rDofType dof type
    //! @return pointer to the node
    virtual NodeBase* GetNode(int rLocalNodeNumber, Node::eDof rDofType)=0;

    //! @brief returns a pointer to the i-th node of the element
    //! @param local node number
    //! @brief rDofType dof type
    //! @return pointer to the node
    virtual const NodeBase* GetNode(int rLocalNodeNumber, Node::eDof rDofType) const=0;

    //! @brief sets the rLocalNodeNumber-th node of the element
    //! @param local node number
    //! @param pointer to the node
    virtual void SetNode(int rLocalNodeNumber, NodeBase* rNode)=0;

    //! @brief resizes the node vector
    //! @param rNewNumNodes new number of nodes
    virtual void ResizeNodes(int rNewNumNodes) = 0;

    //! brief exchanges the node ptr in the full data set (elements, groups, loads, constraints etc.)
    //! this routine is used, if e.g. the data type of a node has changed, but the restraints, elements etc. are still identical
    virtual void ExchangeNodePtr(NodeBase* rOldPtr, NodeBase* rNewPtr)=0;

    //! @brief sets the constitutive law for an element
    //! @param rConstitutiveLaw Pointer to constitutive law entry
    virtual void SetConstitutiveLaw(ConstitutiveBase* rConstitutiveLaw);

    //! @brief sets the constitutive law for an element
    //! @param rIp id of integration point
    //! @param rConstitutiveLaw Pointer to constitutive law entry
    virtual void SetConstitutiveLaw(int rIp, ConstitutiveBase* rConstitutiveLaw);

    //! @brief returns a pointer to the constitutive law for an integration point
    //! @param integration point number (counting from zero)
    //! @return pointer to constitutive law
    const ConstitutiveBase* GetConstitutiveLaw(int rIp) const;

    //! @brief returns a pointer to the constitutive law for an integration point
    //! @param integration point number (counting from zero)
    //! @return pointer to constitutive law
    ConstitutiveBase* GetConstitutiveLaw(int rIp);

    //! @brief returns true, if the constitutive law has been assigned
    bool HasConstitutiveLawAssigned(int rIp) const;

    //! @brief sets the section of an element
    //! implemented with an exception for all elements, reimplementation required for those elements
    //! which actually need a section
    //! @param rSection pointer to section
    virtual void SetSection(const SectionBase* rSection)=0;

    //! @brief returns a pointer to the section of an element
    //! implemented with an exception for all elements, reimplementation required for those elements
    //! which actually need a section
    //! @return pointer to section
    virtual const SectionBase* GetSection() const=0;

    //! @brief sets the integration type of an element
    //! implemented with an exception for all elements, reimplementation required for those elements
    //! which actually need an integration type
    //! @param rIntegrationType pointer to integration type
    virtual void SetIntegrationType(const IntegrationTypeBase* rIntegrationType, NuTo::IpData::eIpDataType rIpDataType);

    //! @brief returns a pointer to the integration type of an element
    //! implemented with an exception for all elements, reimplementation required for those elements
    //! which actually need an integration type
    //! @return pointer to integration type
    virtual const IntegrationTypeBase* GetIntegrationType() const;

    //! @brief sets the interpolation type of an element
    //! @param rInterpolationType interpolation type
    virtual void SetInterpolationType(const InterpolationType* rInterpolationType);

    //! @brief returns a pointer to the interpolation type of an element
    //! @return pointer to interpolation type
    virtual const InterpolationType* GetInterpolationType() const;

    //! @brief returns ip data type of the element
    //! implemented with an exception for all elements, reimplementation required for those elements
    //! which actually need an integration type
    //! @return pointer to integration type
    NuTo::IpData::eIpDataType GetIpDataType(int rIp) const;

    //! @brief returns the number of integration points
    //! @return number of integration points
    int GetNumIntegrationPoints() const;

    //! @brief returns the weight of an integration point
    //! @param rIpNum integration point
    //! @return weight
    double GetIntegrationPointWeight(int rIpNum) const;

    //! @brief calculates output data for the element
    //! @param rInput ... constitutive input map for the constitutive law
    //! @param rOutput ...  coefficient matrix 0 1 or 2  (mass, damping and stiffness) and internal force (which includes inertia terms)
    virtual eError Evaluate(const ConstitutiveInputMap& rInput, std::map<Element::eOutput, std::shared_ptr<ElementOutputBase>>& rOutput)=0;

    template<int TDim>
    eError EvaluateConstitutiveLaw(const ConstitutiveInputMap &rConstitutiveInput, ConstitutiveOutputMap &rConstitutiveOutput, int rIP);

    //! @brief calculates output data for the element with a standard input (EULER_BACKWARD static data)
    //! @param rOutput ...  coefficient matrix 0 1 or 2  (mass, damping and stiffness) and internal force (which includes inertia terms)
    eError Evaluate(std::map<Element::eOutput, std::shared_ptr<ElementOutputBase>>& rOutput);

    //! @brief integrates the stress over the element
    //! @param rStress integrated stress
    void GetIntegratedStress(FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rStress);

    //! @brief integrates the strain over the element
    //! @param rStrain integrated strain
    void GetIntegratedStrain(FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rStress);

    //! @brief Allocates static data for an integration point of an element
    //! @param rConstitutiveLaw constitutive law, which is called to allocate the static data object
    virtual ConstitutiveStaticDataBase* AllocateStaticData(const ConstitutiveBase* rConstitutiveLaw) const=0;

    //! @brief Returns the static data for an integration point of an element
    //! @param rIp integration point
    //! @return static data
    ConstitutiveStaticDataBase* GetStaticData(int rIp);

    //! @brief Returns the static data for an integration point of an element
    //! @param rIp integration point
    //! @return static data
    const ConstitutiveStaticDataBase* GetStaticData(int rIp) const;

    IpDataStaticDataBase& GetStaticDataBase(int rIp);

    const IpDataStaticDataBase& GetStaticDataBase(int rIp) const;

    //! @brief sets the static data for an integration point of an element
    //! @param rIp integration point
    //! @param rStaticData static data
    void SetStaticData(int rIp, ConstitutiveStaticDataBase* rStaticData);

    //! @brief Update the static data of an element
    //virtual eError UpdateStaticData(NuTo::Element::eUpdateType rUpdateType)=0;

    Eigen::VectorXd ExtractNodeValues(Node::eDof rDofType) const
    {
        return this->ExtractNodeValues(0, rDofType);
    }

    virtual Eigen::VectorXd ExtractNodeValues(int rTimeDerivative, Node::eDof rDofType) const;

    virtual Eigen::VectorXd InterpolateDofGlobal(const Eigen::VectorXd& rNaturalCoordinates, Node::eDof rDofType) const;

    virtual Eigen::VectorXd InterpolateDofGlobal(int rTimeDerivative, const Eigen::VectorXd& rNaturalCoordinates, Node::eDof rDofType) const;

    Eigen::Vector3d InterpolateDof3D(const Eigen::VectorXd& rNaturalCoordinates, Node::eDof rDofType) const;

    //! @brief interpolates a vector dof to 3D, mainly for visualization
    //! @remark is overridden by the boundary element
    //! @param rTimeDerivative ... time derivative (0..2)
    //! @param rNaturalCoordinates ... coordinates of the point in natural element coordinates
    //! @param rDofType ... dof type
    Eigen::Vector3d InterpolateDof3D(int rTimeDerivative, const Eigen::VectorXd& rNaturalCoordinates, Node::eDof rDofType) const;

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
    const std::vector<double>& GetNonlocalWeights(int rIp, int rNonlocalElement) const;

    //! @brief returns a vector of the nonlocal elements
    //! @retrun nonlocal elements
    const std::vector<const NuTo::ElementBase*>& GetNonlocalElements() const;

    //! @brief returns the number of nonlocal elements
    //! @param rConstitutive constitutive model for the nonlocale elements
    //! @rerun number of nonlocal elements
    int GetNumNonlocalElements() const;

    //! @brief delete the nonlocal elements
    void DeleteNonlocalElements();

    //! @brief calculates the area of a plane element via the nodes (probably faster than sum over integration points)
    //! @return Area
    virtual double CalculateArea() const;

    //! @brief calculates the volume of an integration point (weight * detJac)
    //! @return rVolume  vector for storage of the ip volumes (area in 2D, length in 1D)
    virtual const Eigen::VectorXd GetIntegrationPointVolume() const=0;

    //! @brief returns the coordinates of an integration point
    //! @param rIpNum integration point
    //! @return rCoordinates coordinates to be returned
    virtual const Eigen::Vector3d GetGlobalIntegrationPointCoordinates(int rIpNum) const;

    //! @brief computes the natural coordinates of an given point
    //! implemented with an exception for all elements, reimplementation required for those elements
    //! @param rGlobCoords (input) ... pointer to the array of coordinates
    //! @param rLocCoords (output) ... coordinates to be returned
    //! @return True if coordinates are within the element, False otherwise
    virtual bool GetLocalPointCoordinates(const double* rGlobCoords, double* rLocCoords) const;

    //! @brief Gets the control node of an boundary element, if it has one
    //! @return boundary control node
    virtual NodeBase* GetBoundaryControlNode() const;

    //VHIRTHAMTODO remove???
//    //! @brief sets the water volume fraction at the boundary surface
//    //! @return water volume fraction at the boundary surface
//    virtual double GetBoundaryWaterVolumeFraction() const;

//    //! @brief sets the water volume fraction at the boundary surface
//    //! @param water volume fraction at the boundary surface
//    virtual void SetBoundaryWaterVolumeFraction(double rBoundaryWaterVolumeFraction);

//    //! @brief sets the relative humidity at the boundary surface
//    //! @param relative humidity at the boundary surface
//    virtual double GetBoundaryRelativeHumidity() const;

//    //! @brief sets the relative humidity at the boundary surface
//    //! @param relative humidity at the boundary surface
//    virtual void SetBoundaryRelativeHumidity(double rBoundaryRelativeHumidity);

    //! @brief checks if a node is inside of an element
    //! implemented with an exception for all elements, reimplementation required for those elements
    //! @param rGlobCoords (input) ... pointer to the array of coordinates
    //! @return True if coordinates are within the element, False otherwise
    virtual bool CheckPointInside(const double* rGlobCoords) const;

    void Info() const;

    //! @brief returns the knots of the element
    //! @return reference on the matrix containing the knots
    virtual const Eigen::MatrixXd& GetKnots() const
    {throw NuTo::MechanicsException(std::string("[") + __PRETTY_FUNCTION__ +"] Only implemented in ContinuumElementIGA.");}

    //! @brief returns the knotIDs of the element
    //! @return reference on the vector containing the knotIDs
    virtual const Eigen::VectorXi& GetKnotIDs() const
    {throw NuTo::MechanicsException(std::string("[") + __PRETTY_FUNCTION__ +"] Only implemented in ContinuumElementIGA.");}

    virtual Eigen::VectorXd InterpolateDofGlobalCurrentConfiguration(int rTimeDerivative, const Eigen::VectorXd& rNaturalCoordinates, Node::eDof rDofTypeInit, Node::eDof rDofTypeCurrent) const
    {throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Element is not of type ContinuumElement<TDim> or ContinuumElementIGA<TDim>.");}

    virtual Eigen::VectorXd InterpolateDofGlobalSurfaceDerivative(int rTimeDerivative, const Eigen::VectorXd& rParameter, int rDerivative, int rDirection) const
    {throw NuTo::MechanicsException(std::string("[") + __PRETTY_FUNCTION__ +"] Only implemented in ContinuumElementIGA.");}

    virtual Eigen::Matrix<Eigen::VectorXd, Eigen::Dynamic, Eigen::Dynamic> InterpolateDofGlobalSurfaceDerivativeTotal(int rTimeDerivative, const Eigen::VectorXd& rParameter, int rDerivative, int rSurface) const
    {throw NuTo::MechanicsException(std::string("[") + __PRETTY_FUNCTION__ +"] Only implemented in ContinuumElementIGA.");}

    virtual Eigen::VectorXd InterpolateDofGlobalSurfaceNormal(const Eigen::VectorXd& rParameter) const
    {throw NuTo::MechanicsException(std::string("[") + __PRETTY_FUNCTION__ +"] Only implemented in ContinuumElementIGA.");}

    virtual const ContinuumElement<1>& AsContinuumElement1D() const
    {throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Element is not of type ContinuumElement<1>.");}

    virtual const ContinuumElement<2>& AsContinuumElement2D() const
    {throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Element is not of type ContinuumElement<2>.");}

    virtual const ContinuumElement<3>& AsContinuumElement3D() const
    {throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Element is not of type ContinuumElement<3>.");}

    virtual ContinuumElement<1>& AsContinuumElement1D()
    {throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Element is not of type ContinuumElement<1>.");}

    virtual ContinuumElement<2>& AsContinuumElement2D()
    {throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Element is not of type ContinuumElement<2>.");}

    virtual ContinuumElement<3>& AsContinuumElement3D()
    {throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Element is not of type ContinuumElement<3>.");}

    virtual const ContinuumElementIGA<1>& AsContinuumElementIGA1D() const
    {throw NuTo::MechanicsException(__PRETTY_FUNCTION__ , "Element is not of type ContinuumElementIGA<1>.");}

    virtual const ContinuumElementIGA<2>& AsContinuumElementIGA2D() const
    {throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Element is not of type ContinuumElementIGA<2>.");}

    virtual const ContinuumElementIGA<3>& AsContinuumElementIGA3D() const
    {throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Element is not of type ContinuumElementIGA<3>.");}

    virtual ContinuumElementIGA<1>& AsContinuumElementIGA1D()
    {throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Element is not of type ContinuumElementIGA<1>.");}

    virtual ContinuumElementIGA<2>& AsContinuumElementIGA2D()
    {throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Element is not of type ContinuumElementIGA<2>.");}

    virtual ContinuumElementIGA<3>& AsContinuumElementIGA3D()
    {throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Element is not of type ContinuumElementIGA<3>.");}

    virtual const ContinuumElementIGALayer<1>& AsContinuumElementIGALayer1D() const
    {throw NuTo::MechanicsException(__PRETTY_FUNCTION__ , "Element is not of type ContinuumElementIGALayer<1>.");}

    virtual const ContinuumElementIGALayer<2>& AsContinuumElementIGALayer2D() const
    {throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Element is not of type ContinuumElementIGALayer<2>.");}

    virtual ContinuumElementIGALayer<1>& AsContinuumElementIGALayer1D()
    {throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Element is not of type ContinuumElementIGALayer<1>.");}

    virtual ContinuumElementIGALayer<2>& AsContinuumElementIGALayer2D()
    {throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Element is not of type ContinuumElementIGALayer<2>.");}

    virtual const ContinuumBoundaryElement<1>& AsContinuumBoundaryElement1D() const
    {throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Element is not of type ContinuumBoundaryElement<1>.");}

    virtual const ContinuumBoundaryElement<2>& AsContinuumBoundaryElement2D() const
    {throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Element is not of type ContinuumBoundaryElement<2>.");}

    virtual const ContinuumBoundaryElement<3>& AsContinuumBoundaryElement3D() const
    {throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Element is not of type ContinuumBoundaryElement<3>.");}

    virtual ContinuumBoundaryElement<1>& AsContinuumBoundaryElement1D()
    {throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Element is not of type ContinuumBoundaryElement<1>.");}

    virtual ContinuumBoundaryElement<2>& AsContinuumBoundaryElement2D()
    {throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Element is not of type ContinuumBoundaryElement<2>.");}

    virtual ContinuumBoundaryElement<3>& AsContinuumBoundaryElement3D()
    {throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Element is not of type ContinuumBoundaryElement<3>.");}


#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);

    //! @brief NodeBase-Pointer are not serialized to avoid cyclic dependencies, but are serialized as Pointer-Address (uintptr_t)
    //! Deserialization of the NodeBase-Pointer is done by searching and casting back the Address in the map
    //! @param mNodeMapCast   std::map containing the old and new Addresses
    virtual void SetNodePtrAfterSerialization(const std::map<std::uintptr_t, std::uintptr_t>& mNodeMapCast)
    {
        (void)mNodeMapCast;
        /* Do nothing until needed, see e.g. ConstraintNode-class*/
    }
#endif  // ENABLE_SERIALIZATION

#ifdef ENABLE_VISUALIZE

    //! @brief Computes all data in rVisualizationList for the visualization. Decomposes the element into small cells for the cisualization.
    //! @param rVisualize
    //! @param rVisualizationList: a list of visualization components to be visualized
    virtual void Visualize(VisualizeUnstructuredGrid& rVisualize, const std::list<std::shared_ptr<NuTo::VisualizeComponent>>& rVisualizationList);

    //! @brief Computes all data in rVisualizationList for the visualization. Extrapolates integration point data to element nodes
    //! @param rVisualize
    //! @param rVisualizationList: a list of visualization components to be visualized
    virtual void VisualizeExtrapolateToNodes(VisualizeUnstructuredGrid& rVisualize, const std::list<std::shared_ptr<NuTo::VisualizeComponent>>& rVisualizationList);

    //! @brief Computes all data in rVisualizationList for the visualization. Visualizes integration point data as vertiex elements
    //! @param rVisualize
    //! @param rVisualizationList: a list of visualization components to be visualized
    virtual void VisualizeIntegrationPointData(VisualizeUnstructuredGrid& rVisualize, const std::list<std::shared_ptr<NuTo::VisualizeComponent>>& rVisualizationList);

    virtual void GetVisualizationCells(unsigned int& NumVisualizationPoints, std::vector<double>& VisualizationPointLocalCoordinates, unsigned int& NumVisualizationCells, std::vector<NuTo::eCellTypes>& VisualizationCellType, std::vector<unsigned int>& VisualizationCellsIncidence,
            std::vector<unsigned int>& VisualizationCellsIP) const;

#endif // ENABLE_VISUALIZE

    //! @brief returns the structure
    const StructureBase* GetStructure() const
    {
        return mStructure;
    }

private:
    //! @brief returns the Element Data Vector
    //! this was necessary due to recursive problems for serialization (nonlocal data)
    //! this method should only be called from the serialization routine of the structure
    NuTo::ElementDataBase* GetDataPtr() const;

    //! @brief sets the Element Data Vector
    //! this was necessary due to recursive problems for serialization (nonlocal data)
    //! this method should only be called from the serialization routine of the structure
    void SetDataPtr(NuTo::ElementDataBase* rElementData);

protected:
    //! @brief ... just for serialization
    ElementBase()
    {
        mElementData = 0;
    }

    //! @brief ... reorder nodes such that the sign of the length/area/volume of the element changes
    virtual void ReorderNodes();

    //! @brief ... check if the element is properly defined (check node dofs, nodes are reordered if the element length/area/volum is negative)
    virtual void CheckElement() = 0;

    //! @brief ... extract global dofs from nodes (mapping of local row ordering of the element matrices to the global dof ordering)
    //! @param rGlobalRowDofs ... vector of global row dofs
    //virtual void CalculateGlobalRowDofs(std::vector<int>& rGlobalRowDofs) const = 0;

    //! @brief ... extract global dofs from nodes (mapping of local column ordering of the element matrices to the global dof ordering)
    //! @param rGlobalColumnDofs ... vector of global column dofs
    //virtual void CalculateGlobalColumnDofs(std::vector<int>& rGlobalColumnDofs) const = 0;

    //the base class of the elements must not contain any data apart from a const pointer to the structure and a data pointer
    const StructureBase* mStructure;

    //the base class of the elements data
    ElementDataBase *mElementData;

    const InterpolationType* mInterpolationType;
};
}    //namespace NuTo

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::ElementBase)
#endif // ENABLE_SERIALIZATION


