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
#include "mechanics/MechanicsException.h"
#include "mechanics/elements/IPData.h"
#include <memory>


namespace NuTo
{

class BoundaryGradientDamage1D;
class ConstitutiveBase;
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
class VisualizeComponentBase;
class VisualizeComponent;
class ElementOutputBase;
enum class eError;
template<typename IOEnum> class ConstitutiveIOMap;
template <int TDim> class ContinuumElement;
template <int TDim> class ContinuumElementIGA;
template <int TDim> class ContinuumBoundaryElement;
template <int TDim> class ContinuumContactElement;

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
    //! @param rInterpolationType ... interpolation type
    ElementBase(const StructureBase* rStructure, const InterpolationType& rInterpolationType);

    ElementBase(const ElementBase& ) = default;
    ElementBase(      ElementBase&&) = default;

    ElementBase& operator=(const ElementBase& ) = default;
    ElementBase& operator=(      ElementBase&&) = default;

    virtual ~ElementBase() = default;

    //! @brief returns the enum (type of the element)
    //! @return enum
    virtual NuTo::Element::eElementType GetEnumType() const=0;


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
    //! @param rConstitutiveLaw reference to constitutive law entry
    void SetConstitutiveLaw(ConstitutiveBase& rConstitutiveLaw);

    //! @brief returns a reference to the constitutive law for an integration point
    //! @param rIP integration point number (counting from zero)
    //! @return reference to constitutive law
    const ConstitutiveBase& GetConstitutiveLaw(unsigned int rIP) const;

    //! @brief returns a reference to the constitutive law for an integration point
    //! @param rIP integration point number (counting from zero)
    //! @return reference to constitutive law
    ConstitutiveBase& GetConstitutiveLaw(unsigned int rIP);

    //! @brief returns a reference to the IPData object
    //! @return reference to the IPData object
    IPData& GetIPData();

    //! @brief returns true, if the constitutive law has been assigned
    //! @param rIP integration point number (counting from zero)
    bool HasConstitutiveLawAssigned(unsigned int rIP) const;

    //! @brief sets the section of an element
    //! @param rSection reference to section
    virtual void SetSection(const SectionBase& rSection);

    //! @brief returns a reference to the section of an element
    //! @return pointer to section
    virtual const SectionBase& GetSection() const;

    //! @brief sets the integration type of an element
    //! @param rIntegrationType reference to integration type
    void SetIntegrationType(const IntegrationTypeBase& rIntegrationType);

    //! @brief returns a pointer to the integration type of an element
    //! @return reference to integration type
    const IntegrationTypeBase& GetIntegrationType() const;

    //! @brief sets the interpolation type of an element
    //! @param rInterpolationType interpolation type
    void SetInterpolationType(const InterpolationType& rInterpolationType);

    //! @brief returns a pointer to the interpolation type of an element
    //! @return reference to interpolation type
    const InterpolationType& GetInterpolationType() const;

    //! @brief returns the number of integration points
    //! @return number of integration points
    int GetNumIntegrationPoints() const;

    //! @brief returns the weight of an integration point
    //! @param rIP integration point
    //! @return weight
    double GetIntegrationPointWeight(unsigned int rIP) const;

    //! @brief calculates output data for the element
    //! @param rInput ... constitutive input map for the constitutive law
    //! @param rOutput ...  coefficient matrix 0 1 or 2  (mass, damping and stiffness) and internal force (which includes inertia terms)
    virtual eError Evaluate(const ConstitutiveInputMap& rInput, std::map<Element::eOutput, std::shared_ptr<ElementOutputBase>>& rOutput)=0;

    //! @brief Evaluate the constitutive law attached to an integration point.
    //! @param rConstitutiveInput Input map of the constitutive law.
    //! @param rConstitutiveOuput Output map of the constitutive law.
    //! @param IP The current integration point.
    template<int TDim>
    eError EvaluateConstitutiveLaw(
            const ConstitutiveInputMap& rConstitutiveInput,
            ConstitutiveOutputMap& rConstitutiveOutput, unsigned int IP);


    //! @brief calculates output data for the element with a standard input (EULER_BACKWARD static data)
    //! @param rOutput ...  coefficient matrix 0 1 or 2  (mass, damping and stiffness) and internal force (which includes inertia terms)
    eError Evaluate(std::map<Element::eOutput, std::shared_ptr<ElementOutputBase>>& rOutput);

    //! @brief integrates the stress over the element
    //! @param rStress integrated stress
    void GetIntegratedStress(Eigen::MatrixXd& rStress);

    //! @brief integrates the strain over the element
    //! @param rStrain integrated strain
    void GetIntegratedStrain(Eigen::MatrixXd& rStress);

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


    void Info() const;

    //! @brief returns the knots of the element
    //! @return reference on the matrix containing the knots
    virtual const Eigen::MatrixXd& GetKnots() const
    {
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__,
                "Only implemented in ContinuumElementIGA.");
    }

    //! @brief returns the knotIDs of the element
    //! @return reference on the vector containing the knotIDs
    virtual const Eigen::VectorXi& GetKnotIDs() const
    {
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__,
                "Only implemented in ContinuumElementIGA.");
    }

    virtual Eigen::VectorXd InterpolateDofGlobalSurfaceDerivative(int, const Eigen::VectorXd&, int, int) const
    {
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__,
                "Only implemented in ContinuumElementIGA.");
    }

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
    {
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__,
                "Element is not of type ContinuumElementIGA<1>.");
    }

    virtual const ContinuumElementIGA<2>& AsContinuumElementIGA2D() const
    {
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__,
                "Element is not of type ContinuumElementIGA<2>.");
    }

    virtual const ContinuumElementIGA<3>& AsContinuumElementIGA3D() const
    {
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__,
                "Element is not of type ContinuumElementIGA<3>.");
    }

    virtual ContinuumElementIGA<1>& AsContinuumElementIGA1D()
    {
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__,
                "Element is not of type ContinuumElementIGA<1>.");
    }

    virtual ContinuumElementIGA<2>& AsContinuumElementIGA2D()
    {
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__,
                "Element is not of type ContinuumElementIGA<2>.");
    }

    virtual ContinuumElementIGA<3>& AsContinuumElementIGA3D()
    {
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__,
                "Element is not of type ContinuumElementIGA<3>.");
    }


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


protected:

#ifdef ENABLE_SERIALIZATION
    ElementBase() {};
#endif // ENABLE_SERIALIZATION

    //! @brief ... reorder nodes such that the sign of the length/area/volume of the element changes
    virtual void ReorderNodes();

    //! @brief ... check if the element is properly defined (check node dofs, nodes are reordered if the element length/area/volum is negative)
    virtual void CheckElement() = 0;

    void AddPlaneStateToInput(ConstitutiveInputMap& input) const;
    //! @brief ... extract global dofs from nodes (mapping of local row ordering of the element matrices to the global dof ordering)
    //! @param rGlobalRowDofs ... vector of global row dofs
    //virtual void CalculateGlobalRowDofs(std::vector<int>& rGlobalRowDofs) const = 0;

    //! @brief ... extract global dofs from nodes (mapping of local column ordering of the element matrices to the global dof ordering)
    //! @param rGlobalColumnDofs ... vector of global column dofs
    //virtual void CalculateGlobalColumnDofs(std::vector<int>& rGlobalColumnDofs) const = 0;

    //the base class of the elements must not contain any data apart from a const pointer to the structure and a data pointer
    const StructureBase* mStructure;

    const InterpolationType* mInterpolationType;

    IPData mIPData;

};
}    //namespace NuTo

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::ElementBase)
#endif // ENABLE_SERIALIZATION


