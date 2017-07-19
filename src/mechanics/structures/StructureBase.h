#pragma once

#include <functional>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>
#include <eigen3/Eigen/Core>

#include <boost/ptr_container/ptr_map.hpp>
#include <mechanics/MechanicsEnums.h>
#include "base/Logger.h"

#include "mechanics/dofSubMatrixStorage/DofStatus.h"
#include "base/Exception.h"
#include "StructureOutputBlockVector.h"


namespace NuTo
{
class Assembler;
class ConstitutiveBase;
class ConstitutiveStaticDataMultiscale2DPlaneStrain;
class ElementBase;
class EngineeringStrain2D;
class GroupBase;
class IntegrationTypeBase;
class InterpolationBase;
class InterpolationType;
class LoadBase;
class NewtonRaphsonAuxRoutinesBase;
class NodeBase;
class Section;
class SerializeStreamOut;
class SerializeStreamIn;
class StructureOutputBase;
class StructureOutputBlockMatrix;
class StructureOutputBlockVector;
class TimeIntegrationBase;
template <typename IOEnum>
class ConstitutiveIOMap;
template <class T>
class BlockFullMatrix;
template <class T>
class BlockFullVector;
template <class T>
class Group;
template <class T>
class SparseMatrixCSRSymmetric;
template <class T>
class SparseMatrixCSRGeneral;
template <class T>
class SparseMatrixCSRVector2General;
template <class T>
class SparseMatrixCSRVector2Symmetric;

enum class eError;
enum class eGroupId;
enum class eIntegrationType;
enum class eStructureOutput;
enum class eVisualizationType;
enum class eVisualizeWhat;
enum class eDirection;

namespace Constitutive
{
class DamageLaw;
enum class eConstitutiveParameter;
enum class eConstitutiveType;
enum class eInput;
enum class eOutput;
} // namespace Constitutive

namespace Constraint
{
class Constraints;
}

namespace Element
{
enum class eOutput;
} // namespace Element

namespace IpData
{
enum class eIpStaticDataType;
} // namespace IpData

namespace Visualize
{
class UnstructuredGrid;
} // namespace Visualize

typedef ConstitutiveIOMap<Constitutive::eInput> ConstitutiveInputMap;
typedef ConstitutiveIOMap<Constitutive::eOutput> ConstitutiveOutputMap;


//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... standard abstract class for all mechanical structures
class StructureBase
{
    friend class NewmarkIndirect;
    friend class NewmarkDirect;
    friend class VelocityVerlet;

public:
    //! @brief constructor
    //! @param mDimension  Structural dimension (1,2 or 3)
    StructureBase(int mDimension);

    //! @brief deconstructor
    virtual ~StructureBase();

    //! @brief gives the dimension of the Structure
    //! @return Structural dimension (1,2 or 3)
    int GetDimension() const;

    //! @brief ... clear all visualization components
    void ClearVisualizationComponents();

    //! @brief ... export the entire structure to Vtk data file
    //! @param rResultFileName ... file name
    void ExportVtkDataFileElements(const std::string& rResultFileName, bool binary = true);

    //! @brief ... export the entire structure to Vtk data file
    //! @param rResultFileName ... file name
    void ExportVtkDataFileNodes(const std::string& rResultFileName);

    //! @brief ... export an element group to Vtk data file
    //! @param rGroupIdent ... group ident
    //! @param rResultFileName ... file name
    void ElementGroupExportVtkDataFile(int rGroupIdent, const std::string& rResultFileName);

    //! @brief Add visualizeComponent to an element group for the visualization
    //! @param rElementGroup: element group
    //! @param visualizeComponent: visualization component, i.e. displacements, stresses...
    void AddVisualizationComponent(int rElementGroup, const std::string& visualizeComponent);

#ifndef SWIG
    //! @brief Add visualizeComponent to an element group for the visualization
    //! @param rElementGroup: element group
    //! @param visualizeComponent: visualization component, i.e. displacements, stresses...
    void AddVisualizationComponent(int rElementGroup, eVisualizeWhat visualizeComponent);

    //! @brief Set the visualization type for an element group
    //! @param rElementGroup: element group
    //! @param visualization type, i.e. voronoi cell, extrapolated...
    void SetVisualizationType(const int rElementGroup, const eVisualizationType rVisualizationType);

    //! @brief ... define the data sets (scalar, vector etc for the visualize routine based on the mVisualizecomponents
    //! for an element plot
    void DefineVisualizeElementData(Visualize::UnstructuredGrid& visualizer,
                                    const std::vector<eVisualizeWhat>& visualizeComponents) const;

    //! @brief ... define the data sets (scalar, vector etc for the visualize routine based on the mVisualizecomponents
    //! for a node plot
    void DefineVisualizeNodeData(Visualize::UnstructuredGrid& visualizer,
                                 const std::vector<eVisualizeWhat>& visualizeComponents) const;

    //! @brief ... adds all the elements in the vector to the data structure that is finally visualized
    void ElementVectorAddToVisualize(Visualize::UnstructuredGrid& visualizer,
                                     const std::vector<eVisualizeWhat>& visualizeComponents,
                                     const std::vector<ElementBase*>& elements);

    //! @brief ... adds all the elements in the vector to the data structure that is finally visualized
    void ElementVectorAddToVisualize(Visualize::UnstructuredGrid& visualizer,
                                     const std::vector<eVisualizeWhat>& visualizeComponents,
                                     const std::vector<ElementBase*>& elements,
                                     const eVisualizationType rVisualizationType);

    //! @brief ... adds all the elements in the vector to the data structure that is finally visualized
    void ElementTotalAddToVisualize(Visualize::UnstructuredGrid& visualizer,
                                    const std::vector<eVisualizeWhat>& visualizeComponents);

    //! @brief ... adds all the elements in a group to the data structure that is finally visualized
    void ElementGroupAddToVisualize(int rGroupId, Visualize::UnstructuredGrid& visualizer,
                                    const std::vector<eVisualizeWhat>& visualizeComponents);

    //! @brief Get all element groups that are supposed to be visualized
    std::vector<int> GetVisualizationGroups();
#endif // SWIG


#ifndef SWIG

    //! @brief Calculates the initial value rates (velocities) of the system to meet equilibrium
    virtual void CalculateInitialValueRates(TimeIntegrationBase& rTimeIntegrationScheme);

    //! @brief ... evaluates the structure
    virtual void Evaluate(const NuTo::ConstitutiveInputMap& rInput,
                          std::map<eStructureOutput, StructureOutputBase*>& rStructureOutput) = 0;


#endif // SWIG

    //@brief determines the maximum independent sets and stores it at the structure
    // is only relevant for openmp, otherwise the routine is just empty
    void CalculateMaximumIndependentSets();

    //@brief set the number of processors for openmp parallelization
    void SetNumProcessors(int rNumProcessors);

    //@brief get the number of processors for openmp parallelization
    int GetNumProcessors() const;

    //@brief set the number of processors for openmp parallelization
    void SetOMPNested(bool rNested);

#ifndef SWIG

    NuTo::StructureOutputBlockMatrix BuildGlobalHessian(eStructureOutput rOutput);

#endif // SWIG

    NuTo::StructureOutputBlockMatrix BuildGlobalHessian0();
    NuTo::StructureOutputBlockMatrix BuildGlobalHessian1();
    NuTo::StructureOutputBlockMatrix BuildGlobalHessian2();
    NuTo::StructureOutputBlockMatrix BuildGlobalHessian2Lumped();

    NuTo::StructureOutputBlockVector BuildGlobalInternalGradient();

    //! @brief ... build global external load vector (currently for displacements only)
    //! @return  ... StructureOutputBlockVector containing the external loads
    NuTo::StructureOutputBlockVector BuildGlobalExternalLoadVector();

    NuTo::BlockFullVector<double> SolveBlockSystem(const NuTo::BlockSparseMatrix& rMatrix,
                                                   const NuTo::BlockFullVector<double>& rVector) const;


    void SolveGlobalSystemStaticElastic();

    void Contact(const std::vector<int>& rElementGroups);

    NuTo::StructureOutputBlockMatrix BuildGlobalHessian0_CDF(double rDelta);

    bool CheckHessian0(double rDelta, double rRelativeTolerance, bool rPrintWrongMatrices = true);

private: // keep the method definitions close together
    bool CheckHessian0_Submatrix(const BlockSparseMatrix& rHessian0, BlockSparseMatrix& rDiff,
                                 double rRelativeTolerance, bool rPrintWrongMatrices);

public:
    //*************************************************
    //************ Node routines        ***************
    //***  defined in structures/StructureNode.cpp  ***
    //*************************************************
    //! @brief returns the number of nodes
    //! @return number of nodes
    virtual int GetNumNodes() const = 0;

    //! @brief a reference to a node
    //! @param identifier
    //! @return reference to a node
    virtual NodeBase* NodeGetNodePtr(int rIdent) = 0;

#ifndef SWIG
    //! @brief a reference to a node
    //! @param identifier
    //! @return reference to a node
    virtual const NodeBase* NodeGetNodePtr(int rIdent) const = 0;

    //! @brief ... store all elements connected to this node in a vector
    //! @param rNode (Input) 		... node pointer
    //! @param rElements (Output) 	... vector of element pointers
    virtual void NodeGetElements(const NodeBase* rNodePtr, std::vector<ElementBase*>& rElements);

    //! @brief gives the identifier of a node
    //! @param reference to a node
    //! @return identifier
    virtual int NodeGetId(const NodeBase* rNode) const = 0;
#endif // SWIG
    //! @brief ... returns the (first) node that has the specified coordinates within the range
    //! @param rCoordinates
    //! @param rRange
    //! @return ... node id
    int NodeGetIdAtCoordinate(Eigen::VectorXd rCoordinates, double rRange = 1.e-6);

    //! @brief ... returns the (first) node that has the specified coordinates within the range
    //! @param coordinate node coordinates
    //! @param tolerance spherical search range
    //! @return reference to the node
    NodeBase& NodeGetAtCoordinate(Eigen::VectorXd coordinate, double tolerance = 1.e-6);

    //! @brief see NodeGetAtCoordinate(Vector, tolerance) for 1D
    NodeBase& NodeGetAtCoordinate(double coordinate, double tolerance = 1.e-6);

    //! @brief ... store all elements connected to this node in a vector
    //! @param rNodeId (Input) 			... node id
    //! @param rElementNumbers (Output) ... vector of element ids
    void NodeGetElements(const int rNodeId, std::vector<int>& rElementNumbers);

    //! @brief delete node
    //! @param rIdent ... node identifier
    virtual void NodeDelete(const int rIdent) = 0;

    //! @brief info about the nodes in the Structure
    virtual void NodeInfo(int mVerboseLevel) const = 0;

    //! @brief numbers the dofs in the structure
    //! @param rCallerName ... if the method throws it is nice to know by whom it was called.
    virtual void NodeBuildGlobalDofs(std::string rCallerName = "") = 0;

    //! @brief sets the displacements of a node
    //! @param rIdent node identifier
    //! @param rDisplacements matrix (one column) with the displacements
    void NodeSetDisplacements(int rId, const Eigen::VectorXd& rDisplacements);

    //! @brief sets the displacements of a node
    //! @param rIdent node identifier
    //! @param rTimeDerivative time derivative (0 disp, 1 vel, 2 acc)
    //! @param rDisplacements matrix (one column) with the displacements
    void NodeSetDisplacements(int rId, int rTimeDerivative, const Eigen::VectorXd& rDisplacements);

    //! @brief sets the displacements of a group of nodes
    //! @param rIdent node group identifier
    //! @param rTimeDerivative time derivative (0 disp, 1 vel, 2 acc)
    //! @param rDisplacements matrix (one column) with the displacements
    void NodeGroupSetDisplacements(int rIdent, int rTimeDerivative, const Eigen::VectorXd& rDisplacements);

    //! @brief sets the displacements of a group of nodes
    //! @param rIdent node group identifier
    //! @param rDisplacements matrix (one column) with the displacements
    void NodeGroupSetDisplacements(int rIdent, const Eigen::VectorXd& rDisplacements);

    //! @brief Sets the temperature of a node
    //! @param rIdent Node identifier
    //! @param rTemperature Temperature to assign to node
    void NodeSetTemperature(int rId, double rTemperature);

    //! @brief Sets the temperature of a node
    //! @param rIdent Node identifier
    //! @param rTimeDerivative Time derivative (0 temperature, 1 temperature change)
    //! @param rTemperature Temperature (change) to assign to node
    void NodeSetTemperature(int rId, int rTimeDerivative, double rTemperature);

    //! @brief returns the node ids of an node group
    //! @param rGroupId  group number
    //! @param rMembers  return vector with node ids
    void NodeGroupGetMembers(int rGroupId, std::vector<int>& rMembers);

    //! @brief gets the coordinates of a node
    //! @param rNode node identifier
    //! @param rCoordinates matrix (one column) with the coordinates
    void NodeGetCoordinates(int rNode, Eigen::VectorXd& rCoordinates) const;

    //! @brief gets the coordinates of a group of nodes (be careful, the order of the nodes in a group might change
    //! between different runs)
    //! @param rNodeGroup node group identifier
    //! @param rCoordinates matrix (rows/nodes columns/coordinates)
    void NodeGroupGetCoordinates(int rNodeGroup, Eigen::MatrixXd& rCoordinates);

    //! @brief gets the displacements of a node
    //! @param rNode node identifier
    //! @param rDisplacements matrix (one column) with the displacements
    void NodeGetDisplacements(int rNode, Eigen::VectorXd& rDisplacements) const;

    //! @brief gets the displacements of a node
    //! @param rIdent node identifier
    //! @param rTimeDerivative time derivative (0 disp, 1 velocity,2 acceleration)
    //! @param rDisplacements matrix (one column) with the displacements
    void NodeGetDisplacements(int rNode, int rTimeDerivative, Eigen::VectorXd& rDisplacements) const;


#ifndef SWIG
    //! @brief gets the dof     identifiers of a node
    //! @param rNodeId node     identifier
    //! @param rDof             degree of freedom
    //! @return A vector of ids that correspond to rDof of rNodeId
    std::vector<int> NodeGetDofIds(const int rNodeId, Node::eDof rDof) const;
#endif

    //! @brief gets the displacements of a group of nodes (be careful, the order of the nodes in a group might change
    //! between different runs)
    //! @param rNodeGroup node group identifier
    //! @param rDisplacements matrix (rows/nodes columns/rDisplacements)
    void NodeGroupGetDisplacements(int rNodeGroup, Eigen::MatrixXd& rDisplacements);

    //! @brief Get the temperature at a node
    //! @param rNode Node identifier
    //! @return Temperature at the node
    double NodeGetTemperature(int rNode) const;

    //! @brief Get the temperature (change) of a node
    //! @param rIdent Node identifier
    //! @param rTimeDerivative Time derivative (0 temperature, 1 temperature change)
    //! @return Temperature (change) at the node
    double NodeGetTemperature(int rNode, int rTimeDerivative) const;

    //! @brief extract dof values (e.g. displacements, temperatures to the nodes)
    //! @param rTimeDerivative time derivative (0 disp 1 vel 2 acc)
    //! @return ... StructureBlockVector containing the dofs (J and K)
    virtual NuTo::StructureOutputBlockVector NodeExtractDofValues(int rTimeDerivative) const = 0;

    //! @brief extract dof values (e.g. displacements, temperatures to the nodes) (time derivative = 0)
    //! @return ... StructureBlockVector containing the dofs (J and K)
    NuTo::StructureOutputBlockVector NodeExtractDofValues() const;

    //! @brief write dof values (e.g. displacements, temperatures to the nodes)
    //! @param rTimeDerivative time derivative (0 disp 1 vel 2 acc)
    //! @param rActiveDofValues ... vector of independent dof values (ordering according to global dofs, size is number
    //! of active dofs)
    //! @param rDependentDofValues ... vector of dependent  dof values (ordering according to global dofs, size is
    //! number of active dofs)
    virtual void NodeMergeDofValues(int rTimeDerivative, const NuTo::BlockFullVector<double>& rActiveDofValues,
                                    const NuTo::BlockFullVector<double>& rDependentDofValues) = 0;

    //! @brief write dof values (e.g. displacements, temperatures to the nodes)
    //! @param rTimeDerivative time derivative (0 disp 1 vel 2 acc)
    //! @param rDofValues ... StructureBlockVector containing the dofs (J and K)
    virtual void NodeMergeDofValues(int rTimeDerivative, const NuTo::StructureOutputBlockVector& rDofValues) = 0;

    //! @brief write dof values (e.g. displacements, temperatures to the nodes)
    //! @param rActiveDofValues ... vector of independent dof values (ordering according to global dofs, size is number
    //! of active dofs)
    //! @param rDependentDofValues ... vector of dependent  dof values (ordering according to global dofs, size is
    //! number of active dofs)
    virtual void NodeMergeDofValues(NuTo::StructureOutputBlockVector& rDofValues) = 0;

    //! @brief calculate dependent dof values (for the zeroth time derivative)
    //! @param rActiveDofValues ... vector of independent dof values (ordering according to global dofs, size is number
    //! of active dofs)
    //! @return  ... vector of dependent  dof values (ordering according to global dofs, size is number of active dofs)
    virtual NuTo::BlockFullVector<double>
    NodeCalculateDependentDofValues(const NuTo::BlockFullVector<double>& rActiveDofValues) const = 0;

    //! @brief calculate the internal force vector for a node
    //! @param rId ... node id
    //! @param rGradientInternalPotential ...vector for all the dofs the corresponding internal force (return value)
    void NodeInternalForce(int rId, Eigen::VectorXd& rNodeForce);

    //! @brief calculate the internal force vector for a node group of nodes
    //! @param rIdent ... group id
    //! @param rGradientInternalPotential ...vector for all the dofs the corresponding internal force (return value)
    void NodeGroupInternalForce(int rIdent, Eigen::VectorXd& rNodeForce);

#ifndef SWIG
    //! @brief calculate the internal force vector for a node
    //! @param rNodePtr  node for which this has to be calculated
    //! @param rGradientInternalPotential ...vector for all the dofs the corresponding internal force (return value)
    void NodeInternalForce(const NodeBase* rNodePtr, Eigen::VectorXd& rNodeForce);

    //! @brief ... adds all the elements in the vector to the data structure that is finally visualized
    void NodeTotalAddToVisualize(Visualize::UnstructuredGrid& visualizer,
                                 const std::vector<NuTo::eVisualizeWhat>& visualizeComponents) const;

    //! @brief ... adds all the nodes in the vector to the data structure that is finally visualized
    void NodeVectorAddToVisualize(Visualize::UnstructuredGrid& visualizer,
                                  const std::vector<NuTo::eVisualizeWhat>& virualizeComponents,
                                  const std::vector<const NodeBase*>& nodes) const;

#endif // SWIG

    //*************************************************
    //************ Element routines     ***************
    //**  defined in structures/StructureElement.cpp **
    //*************************************************
    //! @brief returns the number of elements
    //! @return number of elements
    virtual int GetNumElements() const = 0;

#ifndef SWIG
    //! @brief a reference to an element
    //! @param identifier
    //! @return reference to an element
    virtual ElementBase* ElementGetElementPtr(int rIdent) = 0;

    //! @brief a reference to an element
    //! @param identifier
    //! @return reference to an element
    virtual const ElementBase* ElementGetElementPtr(int rIdent) const = 0;

    //! @brief gives the identifier of an element
    //! @param reference to an element
    //! @return identifier
    virtual int ElementGetId(const ElementBase* rElement) const = 0;

    //! @brief info about one single element
    //! @param rElement (Input) ... pointer to the element
    //! @param rVerboseLevel (Input) ... level of verbosity
    virtual void ElementInfo(const ElementBase* rElement, int rVerboseLevel) const = 0;

#endif // SWIG

    //! @brief delete element
    //! @param rIdent ... element identifier
    virtual void ElementDelete(int rIdent) = 0;

    //! @brief info about the elements in the Structure
    virtual void ElementInfo(int rVerboseLevel) const = 0;

#ifndef SWIG

    //! @brief duilds the element hessian0 using the central differences of the internal force vector
    //! @param rElement ... element
    //! @param rDelta ... step of the central differences
    //! @return BlockFullMatrix containing the hessian0_cdf
    BlockFullMatrix<double> ElementBuildHessian0_CDF(ElementBase& rElement, double rDelta);

    //! @brief builds the element internal gradient
    //! @param rElement ... element
    //! @return BlockFullVector containing the internal gradient
    BlockFullVector<double> ElementBuildInternalGradient(ElementBase& rElement);

    //! @brief builds the element hessian
    //! @param rHessianType ... hessian0 hessian1 hessian2
    //! @param rElement ... element
    //! @return BlockFullMatrix containing the hessian
    BlockFullMatrix<double> ElementBuildHessian(Element::eOutput rHessianType, ElementBase& rElement);

    //! @brief builds the element hessian0
    //! @param rElement ... element
    //! @return BlockFullMatrix containing the hessian
    BlockFullMatrix<double> ElementBuildHessian0(ElementBase& rElement);

    //! @brief builds the element hessian1
    //! @param rElement ... element
    //! @return BlockFullMatrix containing the hessian
    BlockFullMatrix<double> ElementBuildHessian1(ElementBase& rElement);

    //! @brief builds the element hessian2
    //! @param rElement ... element
    //! @return BlockFullMatrix containing the hessian
    BlockFullMatrix<double> ElementBuildHessian2(ElementBase& rElement);

    //! @brief builds the element vector of global row dofs
    //! @param rElement ... element
    //! @return BlockFullVector containing the global row dofs
    BlockFullVector<int> ElementBuildGlobalDofsRow(ElementBase& rElement);

    //! @brief builds the element vector of global column dofs
    //! @param rElement ... element
    //! @return BlockFullVector containing the global column dofs
    BlockFullVector<int> ElementBuildGlobalDofsColumn(ElementBase& rElement);

    bool ElementCheckHessian0(ElementBase& rElement, double rDelta, double rRelativeTolerance,
                              bool rPrintWrongMatrices = true);

#endif // SWIG

    NuTo::BlockFullMatrix<double> ElementBuildHessian0(int rElementId);
    NuTo::BlockFullMatrix<double> ElementBuildHessian1(int rElementId);
    NuTo::BlockFullMatrix<double> ElementBuildHessian2(int rElementId);

    NuTo::BlockFullVector<double> ElementBuildInternalGradient(int rElementId)
    {
        return ElementBuildInternalGradient(*ElementGetElementPtr(rElementId));
    }
    NuTo::BlockFullVector<int> ElementBuildGlobalDofsRow(int rElementId)
    {
        return ElementBuildGlobalDofsRow(*ElementGetElementPtr(rElementId));
    }
    NuTo::BlockFullVector<int> ElementBuildGlobalDofsColumn(int rElementId)
    {
        return ElementBuildGlobalDofsColumn(*ElementGetElementPtr(rElementId));
    }

    bool ElementCheckHessian0(int rElementId, double rDelta, double rRelativeTolerance, bool rPrintWrongMatrices = true)
    {
        return ElementCheckHessian0(*ElementGetElementPtr(rElementId), rDelta, rRelativeTolerance, rPrintWrongMatrices);
    }

    bool ElementCheckHessian0(double rDelta, double rRelativeTolerance, bool rPrintWrongMatrices = true);

    //! @brief modifies the constitutive law of a single element
    //! @param rElementIdent identifier for the element
    //! @param rConstitutiveLawIdent identifier for the material
    void ElementSetConstitutiveLaw(int rElementId, int rConstitutiveLawIdent);

    //! @brief modifies the constitutive law of a group of elements
    //! @param rGroupIdent identifier for the group of elements
    //! @param rConstitutiveLawIdent identifier for the material
    void ElementGroupSetConstitutiveLaw(int rGroupIdent, int rConstitutiveLawIdent);

    //! @brief modifies the constitutive law of a all elements
    //! @param rConstitutiveLawIdent identifier for the material
    void ElementTotalSetConstitutiveLaw(int rConstitutiveLawIdent);

#ifndef SWIG
    //! @brief modifies the constitutive law of a single element
    //! @param rElement element pointer
    //! @param rConstitutive material pointer
    void ElementSetConstitutiveLaw(ElementBase* rElement, ConstitutiveBase* rConstitutive);

#endif // SWIG

    //! @brief Modifies the section of a group of elements
    //! @param groupId Identifier for the group of elements
    //! @param section Section to add to the elements
    void ElementGroupSetSection(int rGroupIdent, std::shared_ptr<Section> section);

    //! @brief Modifies the section of all elements
    //! @param section Section to add to the elements
    void ElementTotalSetSection(std::shared_ptr<Section> section);

    //! @brief modifies the section of a single element
    //! @param rElement element pointer
    //! @param rSection section
    void ElementSetSection(ElementBase* rElement, std::shared_ptr<Section> rSection);

    //! @brief modifies the interpolation type of a group of elements
    //! @param rGroupId ... identifier for the group of elements
    //! @param rInterpolationTypeId ... interpolation type id
    void ElementGroupSetInterpolationType(int rGroupId, int rInterpolationTypeId);

#ifndef SWIG
    //! @brief modifies the interpolation type of a single element
    //! @param rElement element pointer
    //! @param rInterpolationType interpolation type
    void ElementSetInterpolationType(ElementBase* rElement, InterpolationType* rInterpolationType);


    //! @brief calculates static ip data
    //! @param rElemIdent  element number
    //! @param rType static ip data type
    //! @param rIPData matrix with (... x numIP), x varies depending on IPData type
    Eigen::MatrixXd ElementGetStaticIPData(int rElementId, IpData::eIpStaticDataType rType);

#endif // SWIG

    //! @brief calculates static ip data
    //! @param rElemIdent  element number
    //! @param rType static ip data type
    //! @param rIPData matrix with (... x numIP), x varies depending on IPData type
    Eigen::MatrixXd ElementGetStaticIPData(int rElementId, std::string rType);

    //! @brief calculates the engineering strain
    //! @param rElemIdent  element number
    //! @param rEngineerungStrain engineering strain (return value, always 6xnumIp matrix)
    Eigen::MatrixXd ElementGetEngineeringStrain(int rElementId);

    //! @brief calculates the engineering plastic strain
    //! @param rElemIdent  element number
    //! @param rEngineerungStrain engineering plastic strain (return value, always 6xnumIp matrix)
    Eigen::MatrixXd ElementGetEngineeringPlasticStrain(int rElementId);

    //! @brief calculates the engineering stress
    //! @param rElemIdent  element number
    //! @param rEingineeringStress Engineering Stress (return value, always 6xnumIp matrix)
    Eigen::MatrixXd ElementGetEngineeringStress(int rElementId);

    //! @brief calculates the damage
    //! @param rElemIdent  identifier for the element
    //! @param rDamage (return value, always 1xnumIp matrix)
    Eigen::MatrixXd ElementGetDamage(int rElementId);

    //! @brief calculates the global integration point coordinates
    //! @param rElemIdent  identifier for the element
    //! @param rCoordinates integration point coordinates (return value, always 3xnumIp matrix)
    Eigen::MatrixXd ElementGetIntegrationPointCoordinates(int rElementId);


    //! @brief calculates the maximum damage in all elements
    //! @param rElemIdent  identifier for the element
    //! @return max damage value
    double ElementTotalGetMaxDamage();

    //! @brief calculates some error in the static data, e.g. for IMPL-EX
    double ElementTotalGetStaticDataExtrapolationError();

    //! @brief allocates additional static data for an element group
    //! @param rElementGroupId ... element group id
    //! @param rNumAdditionalStaticData ... number of addidional static data objects
    void ElementGroupAllocateAdditionalStaticData(int rElementGroupId, int rNumAdditionalStaticData);

    //! @brief updates the history data of a all elements
    void ElementTotalUpdateStaticData();

    //! @brief updates the temprory static data of a all elements
    //! its is a const function, since only mutuable data (instead of const) is updated (kind of temporary data)
    void ElementTotalUpdateTmpStaticData();

    //! @brief For each integration point, shift the constitutive static data one step into the past.
    void ElementTotalShiftStaticDataToPast();

    //! @brief For each integration point, shift the constitutive static data one step into the future.
    void ElementTotalShiftStaticDataToFuture();

    //! @brief extrapolates static data of a all elements
    void ElementTotalExtrapolateStaticData();

    //! @brief calculates the average stress
    //! @param rVolume  volume of the structure in 3D /area in 2D/ length in 1D
    //! this is a parameter of the model, since holes have to be considered (zero stress, but still nonzero area)
    //! @param rEngineeringStress  average stress (return value)
    void ElementTotalGetAverageStress(double rVolume, Eigen::MatrixXd& rEngineeringStress);

    //! @brief calculates the average stress
    //! @param rGroupId  group number
    //! @param rVolume  volume of the structure in 3D /area in 2D/ length in 1D
    //! this is a parameter of the model, since holes have to be considered (zero stress, but still nonzero area)
    //! @param rEngineeringStress  average stress (return value)
    void ElementGroupGetAverageStress(int rGroupId, double rVolume, Eigen::MatrixXd& rEngineeringStress);

    //! @brief calculates the average strain
    //! @param rVolume  volume of the structure in 3D /area in 2D/ length in 1D
    //! this is a parameter of the model, since holes have to be considered (zero stress, but still nonzero area)
    //! @param rEngineeringStrain  average strain (return value)
    void ElementTotalGetAverageStrain(double rVolume, Eigen::MatrixXd& rEngineeringStrain);

    //! @brief calculates the average strain
    //! @param rGroupId  group number
    //! @param rVolume  volume of the structure in 3D /area in 2D/ length in 1D
    //! this is a parameter of the model, since holes have to be considered (zero stress, but still nonzero area)
    //! @param rEngineeringStrain  average strain (return value)
    void ElementGroupGetAverageStrain(int rGroupId, double rVolume, Eigen::MatrixXd& rEngineeringStrain);

    //! @brief returns the element ids of an element group
    //! @param rGroupId  group number
    //! @param rMembers  return vector with element ids
    void ElementGroupGetMembers(int rGroupId, std::vector<int>& rMembers);

    //! @brief calculates the volume of the elements
    //! @param rGroupId  group number
    //! @return volume of the structure in 3D /area in 2D/ length in 1D
    double ElementGroupGetVolume(int rGroupId);

    //! @brief calculate the largest element eigenvalue for all elements solving the generalized eigenvalue problem
    //! Ku=lambda Mu
    //! this is used for the estimation of the critical time step
    double ElementTotalCalculateLargestElementEigenvalue();

    //! @brief calculate the largest element eigenvalue for a group of elements solving the generalized eigenvalue
    //! problem Ku=lambda Mu
    //! this is used for the estimation of the critical time step
    double ElementGroupCalculateLargestElementEigenvalue(int rGroupId);

    //! @brief ... create a constraint equation
    //! @param rConstraint ... constraint id
    //! @param rNode ... node id in the first constraint equation term
    //! @param rDofType ... type of dof in the first constraint equation term (e.g DISPLACEMENTS, ROTATIONS,
    //! TEMPERATURES)
    //! @param rDofComponent ... dof component (0, 1, 2)
    //! @param rCoefficient ... weight factor of this term
    //! @param rRHS ... prescribed right hand side value

    //! @brief Creates a constraint equation that couples the degrees of freedom of an arbitrary node to a point in an
    //! element
    //! @param rNode ... node id in the first constraint equation term
    //! @param rElement ... element group id
    //! @param rDofType ... type of dof in the first constraint equation term (e.g DISPLACEMENTS, ROTATIONS,
    //! TEMPERATURES)
    //! @param numNearestNeighbours ... number of nearest neighbours to be found by the approximate nearest neighbour
    //! algorithm
    void ConstraintLinearEquationNodeToElementCreate(int rNode, int rElementGroup, NuTo::Node::eDof rDofType,
                                                     const double rTolerance = 1.e-6,
                                                     Eigen::Vector3d rNodeCoordOffset = Eigen::Vector3d::Zero());

public:
    //*************************************************
    //************ Load routines        ***************
    //***  defined in structures/StructureBaseLoad.cpp  ***
    //*************************************************

    //! @brief adds a scalar source for a node
    //! @param rNodeIdent ... identifier for node
    //! @param rValue ... source value
    //! @return integer id to delete or modify the load
    int LoadCreateScalarSource(int rNodeIdent, double rValue);

    //! @brief adds a force for a node
    //! @param rNodeIdent ... identifier for node
    //! @param rDirection ... direction of the force
    //! @param rValue ... force
    //! @return integer id to delete or modify the load
    int LoadCreateNodeForce(int rNodeIdent, const Eigen::MatrixXd& rDirection, double rValue);

    //! @brief adds a force for a node group
    //! @param rGroupIdent ... identifier for node group
    //! @param rDirection ... direction of the force
    //! @param rValue ... force
    //! @return integer id to delete or modify the load
    int LoadCreateNodeGroupForce(int rGroupIdent, const Eigen::MatrixXd& rDirection, double rValue);

    //! @brief adds a surface load to 2D plane elements (2D)
    //! @param rElementGroupId ... specifies the elements with surface loads
    //! @param rNodeGroupId ... specifies the surfaces (if all nodes of an elemental surface is included in this group,
    //! the
    // surface is considered to be loaded
    //! @param rLoadVector ... constant load vector (independent of position and orientation of the loading surface
    //! @return integer id to delete or modify the load
    int LoadSurfaceConstDirectionCreate2D(int rElementGroupId, int rNodeGroupId, const Eigen::VectorXd& rLoadVector);

    //! @brief adds a surface load to 3D solid elements
    //! @param rElementGroupId ... specifies the elements with surface loads
    //! @param rNodeGroupId ... specifies the surfaces (if all nodes of an elemental surface is included in this group,
    //! the
    // surface is considered to be loaded
    //! @param rLoadVector ... constant load vector (independent of position and orientation of the loading surface
    //! @return integer id to delete or modify the load
    int LoadSurfaceConstDirectionCreate3D(int rElementGroupId, int rNodeGroupId, const Eigen::VectorXd& rLoadVector);

    //! @brief adds a surface load (pressure) to 3D solid elements
    //! @param rElementGroupId ... specifies the elements with surface loads
    //! @param rNodeGroupId ... specifies the surfaces (if all nodes of an elemental surface is included in this group,
    //! the
    // surface is considered to be loaded
    //! @param rPressure value ... normal to the surface, positive for compression
    //! @return integer id to delete or modify the load
    int LoadSurfacePressureCreate2D(int rElementGroupId, int rNodeGroupId, double rPressure);

    //! @brief adds a surface load pressure-function to 2D elements
    //! @param rElementGroupId ... specifies the elements with surface loads
    //! @param rNodeGroupId ... specifies the surfaces (if all nodes of an elemental surface is included in this group,
    //! the surface is considered to be loaded
    //! @param rLoadFunction ... pressure function on the boundary
    //! @return integer id to delete or modify the load
    int LoadSurfacePressureFunctionCreate2D(int rElementGroupId, int rNodeGroupId,
                                            const std::function<Eigen::Vector2d(Eigen::Vector2d)>& rLoadFunction);

    //! @brief adds a surface load (pressure) to 3D solid elements
    //! @param rElementGroupId ... specifies the elements with surface loads
    //! @param rNodeGroupId ... specifies the surfaces (if all nodes of an elemental surface is included in this group,
    //! the
    // surface is considered to be loaded
    //! @param rPressure value ... normal to the surface, positive for compression
    //! @return integer id to delete or modify the load
    int LoadSurfacePressureCreate3D(int rElementGroupId, int rNodeGroupId, double rPressure);

#ifndef SWIG

    //! @brief adds a scalar source for a node
    //! @param rNode ... pointer to node
    //! @param rValue ... source value
    //! @return integer id to delete or modify the load
    int LoadCreateScalarSource(const NodeBase* rNode, double rValue);

    //! @brief adds a force for a node
    //! @param rNode ... pointer to node
    //! @param rDirection ... direction of the force
    //! @param rValue ... force
    //! @return integer id to delete or modify the load
    int LoadCreateNodeForce(const NodeBase* rNode, const Eigen::MatrixXd& rDirection, double rValue);

    //! @brief Adds a heat flux to a node.
    //! @param rNode Pointer to node
    //! @param rDirection Direction of the heat flux
    //! @param rValue Value of the flux
    //! @return Integer id to delete or modify the load
    int LoadCreateNodeHeatFlux(const NodeBase* rNode, const Eigen::MatrixXd& rDirection, double rValue);

    //! @brief adds a force for a node grpup
    //! @param rNodeGroup ... pointer to node group
    //! @param rDirection ... direction of the force
    //! @param rValue ... force
    //! @return integer id to delete or modify the load
    int LoadCreateNodeGroupForce(const Group<NodeBase>* rNodeGroup, const Eigen::MatrixXd& rDirection, double rValue);

    //! @brief ... get the pointer to a load from the load identifier
    //! @param rIdent ... load identifier
    //! @return ... pointer to the load
    LoadBase* LoadGetLoadPtr(int rIdent);

    //! @brief ... get the pointer to a load from the load identifier
    //! @param rIdent ... load identifier
    //! @return ... pointer to the load
    const LoadBase* LoadGetLoadPtr(int rIdent) const;
#endif

    //*************************************************
    //************ Constitutive Law routines *********
    //**  defined in structures/StructureBaseConstitutiveLaws.cpp **
    //*************************************************
    //! @brief ... returns the number of constitutive laws
    //! @return number of constitutive laws
    inline int GetNumConstitutiveLaws() const
    {
        return this->mConstitutiveLawMap.size();
    }

    //! @brief ... create a new constitutive law
    //! @param rIdent ... constitutive law identifier
    //! @param rType ... constitutive law type
    int ConstitutiveLawCreate(const std::string& rType);

    //! @brief ... delete an existing constitutive law
    //! @param rIdent ... constitutive law identifier
    void ConstitutiveLawDelete(int rIdent);

    //! @brief ... print information about all constitutive laws
    //! @param rVerboseLevel ... controls the verbosity of the information
    void ConstitutiveLawInfo(unsigned short rVerboseLevel) const;

    //! @brief ... print information of a single constitutive law
    //! @param rIdent ... constitutive law identifier
    //! @param rVerboseLevel ... controls the verbosity of the information
    void ConstitutiveLawInfo(int rIdent, unsigned short rVerboseLevel) const;


    //! @brief ... gets a variable of the constitutive law which is selected by a string
    //! @param rIdent ... constitutive law identifier
    //! @param rIdentifier ... String to identify the requested variable
    //! @return ... value of the requested variable
    bool ConstitutiveLawGetParameterBool(int rIdent, const std::string& rIdentifier) const;

    //! @brief ... sets a variable of the constitutive law which is selected by a string
    //! @param rIdent ... constitutive law identifier
    //! @param rIdentifier ... String to identify the requested variable
    //! @param rValue ... new value for requested variable
    void ConstitutiveLawSetParameterBool(int rIdent, const std::string& rIdentifier, bool rValue);

    //! @brief ... gets a variable of the constitutive law which is selected by a string
    //! @param rIdent ... constitutive law identifier
    //! @param rIdentifier ... String to identify the requested variable
    //! @return ... value of the requested variable
    double ConstitutiveLawGetParameterDouble(int rIdent, const std::string& rIdentifier) const;

    //! @brief ... sets a variable of the constitutive law which is selected by a string
    //! @param rIdent ... constitutive law identifier
    //! @param rIdentifier ... String to identify the requested variable
    //! @param rValue ... new value for requested variable
    void ConstitutiveLawSetParameterDouble(int rIdent, const std::string& rIdentifier, double rValue);

    //! @brief ... gets a variable of the constitutive law which is selected by a string
    //! @param rIdent ... constitutive law identifier
    //! @param rIdentifier ... String to identify the requested variable
    //! @return ... value of the requested variable
    Eigen::VectorXd ConstitutiveLawGetParameterFullVectorDouble(int rIdent, const std::string& rIdentifier) const;

    //! @brief ... sets a variable of the constitutive law which is selected by a string
    //! @param rIdent ... constitutive law identifier
    //! @param rIdentifier ... String to identify the requested variable
    //! @param rValue ... new value for requested variable
    void ConstitutiveLawSetParameterFullVectorDouble(int rIdent, const std::string& rIdentifier,
                                                     Eigen::VectorXd rValue);

    //! @brief ... gets a variable of the constitutive law which is selected by a string
    //! @param rIdent ... constitutive law identifier
    //! @param rIdentifier ... String to identify the requested variable
    //! @return ... value of the requested variable
    Eigen::MatrixXd ConstitutiveLawGetParameterMatrixDouble(int rIdent, const std::string& rIdentifier) const;

    //! @brief ... sets a variable of the constitutive law which is selected by a string
    //! @param rIdent ... constitutive law identifier
    //! @param rIdentifier ... String to identify the requested variable
    //! @param rValue ... new value for requested variable
    void ConstitutiveLawSetParameterMatrixDouble(int rIdent, const std::string& rIdentifier, Eigen::MatrixXd rValue);


#ifndef SWIG

    //! @brief ... gets a parameter of the constitutive law which is selected by an enum
    //! @param rIdent ... constitutive law identifier
    //! @param rIdentifier ... Enum to identify the requested parameter
    //! @return ... value of the requested variable
    bool ConstitutiveLawGetParameterBool(int rIdent, Constitutive::eConstitutiveParameter rIdentifier) const;

    //! @brief ... sets a parameter of the constitutive law which is selected by an enum
    //! @param rIdent ... constitutive law identifier
    //! @param rIdentifier ... Enum to identify the requested parameter
    //! @param rValue ... new value for requested variable
    void ConstitutiveLawSetParameterBool(int rIdent, Constitutive::eConstitutiveParameter rIdentifier, bool rValue);

    //! @brief ... gets a parameter of the constitutive law which is selected by an enum
    //! @param rIdent ... constitutive law identifier
    //! @param rIdentifier ... Enum to identify the requested parameter
    //! @return ... value of the requested variable
    double ConstitutiveLawGetParameterDouble(int rIdent, Constitutive::eConstitutiveParameter rIdentifier) const;

    //! @brief ... sets a parameter of the constitutive law which is selected by an enum
    //! @param rIdent ... constitutive law identifier
    //! @param rIdentifier ... Enum to identify the requested parameter
    //! @param rValue ... new value for requested variable
    void ConstitutiveLawSetParameterDouble(int rIdent, Constitutive::eConstitutiveParameter rIdentifier, double rValue);

    //! @brief ... gets a parameter of the constitutive law which is selected by an enum
    //! @param rIdent ... constitutive law identifier
    //! @param rIdentifier ... Enum to identify the requested parameter
    //! @return ... value of the requested variable
    Eigen::VectorXd ConstitutiveLawGetParameterFullVectorDouble(int rIdent,
                                                                Constitutive::eConstitutiveParameter rIdentifier) const;

    //! @brief ... sets a parameter of the constitutive law which is selected by an enum
    //! @param rIdent ... constitutive law identifier
    //! @param rIdentifier ... Enum to identify the requested parameter
    //! @param rValue ... new value for requested variable
    void ConstitutiveLawSetParameterFullVectorDouble(int rIdent, Constitutive::eConstitutiveParameter rIdentifier,
                                                     Eigen::VectorXd rValue);

    //! @brief ... gets a parameter of the constitutive law which is selected by an enum
    //! @param rIdent ... constitutive law identifier
    //! @param rIdentifier ... Enum to identify the requested parameter
    //! @return ... value of the requested variable
    Eigen::MatrixXd ConstitutiveLawGetParameterMatrixDouble(int rIdent,
                                                            Constitutive::eConstitutiveParameter rIdentifier) const;

    //! @brief ... sets a parameter of the constitutive law which is selected by an enum
    //! @param rIdent ... constitutive law identifier
    //! @param rIdentifier ... Enum to identify the requested parameter
    //! @param rValue ... new value for requested variable
    void ConstitutiveLawSetParameterMatrixDouble(int rIdent, Constitutive::eConstitutiveParameter rIdentifier,
                                                 Eigen::MatrixXd rValue);
#endif

    //! @brief ... set damage law
    //! @param lawId ... lawId
    //! @param damageLaw ... damage law
    void ConstitutiveLawSetDamageLaw(int lawId, std::shared_ptr<Constitutive::DamageLaw> damageLaw);

    //! @brief ... gets the equilibrium water volume fraction depend on the relative humidity
    //! @param rIdent ... constitutive law identifier
    //! @param rRelativeHumidity ... relative humidity
    //! @param rCoeffs ... polynomial coefficients of the sorption curve
    //! @return ... equilibrium water volume fraction
    double ConstitutiveLawGetEquilibriumWaterVolumeFraction(int rIdent, double rRelativeHumidity,
                                                            Eigen::VectorXd rCoeffs) const;

//    //VHIRTHAMTODO Delete???
//    //! @brief ... adds a constitutive law to a multi physics model
//    //! @param rIdentMultiPhysics ... multi physics constitutive law to which the constitutive law should be added
//    //! @param rIdentConstitutiveLaw ... constitutive law which should be added to the multi physics model
//    void ConstitutiveLawMultiPhysicsAddConstitutiveLaw(int rIdentMultiPhysics, int rIdentConstitutiveLaw);

#ifndef SWIG

    //! @brief ... create a new constitutive law
    //! @param rIdent ... constitutive law identifier
    //! @param rType ... constitutive law type
    int ConstitutiveLawCreate(Constitutive::eConstitutiveType rType);

    //! @brief ... create a new constitutive law
    //! @param rIdent ... constitutive law identifier
    //! @param rType ... constitutive law type
    void ConstitutiveLawCreate(int rIdent, Constitutive::eConstitutiveType rType);

    //! @brief Adds a constitutive law
    //! @param rConstitutiveLawPtr: Pointer to constitutive law
    int AddConstitutiveLaw(ConstitutiveBase* rConstitutiveLawPtr);

    //! @brief ... get the pointer to a constitutive law from the constitutive law identifier
    //! @param rIdent ... constitutive law identifier
    //! @return ... pointer to the constitutive law
    ConstitutiveBase* ConstitutiveLawGetConstitutiveLawPtr(int rIdent);

    //! @brief ... get the pointer to a constitutive law from the constitutive law identifier
    //! @param rIdent ... constitutive law identifier
    //! @return ... pointer to the constitutive law
    const ConstitutiveBase* ConstitutiveLawGetConstitutiveLawPtr(int rIdent) const;

    //! @brief ... get the identifier of the constitutive law from the pointer to constitutive law object
    //! @param rConstitutiveLawPtr ... pointer to constitutive law
    //! @return ... constitutive law identifier
    int ConstitutiveLawGetId(const ConstitutiveBase* rConstitutiveLawPtr) const;
#endif

    //*************************************************
    //************ Group routines     ***************
    //**  defined in structures/StructureBaseGroup.cpp **
    //*************************************************
    //! @brief ... Info routine that prints general information about the groups
    //! @param rVerboseLevel describes how detailed the information is
    void GroupInfo(int rVerboseLevel) const;

    //! @brief gives the identifier of an element
    //! @param reference to an element
    //! @return identifier
    int GroupGetId(GroupBase* rGroup) const;

    //! @brief ... get the pointer to a group from the group identifier
    //! @param rIdent ... group identifier
    //! @return ... pointer to the group
    GroupBase* GroupGetGroupPtr(int rIdent);

    //! @brief ... get the pointer to a group from the group identifier
    //! @param rIdent ... group identifier
    //! @return ... pointer to the group
    const GroupBase* GroupGetGroupPtr(int rIdent) const;

#ifndef SWIG

    //! @brief Creates a group of node ids for the structure
    //! @return rIdent identifier for the node group
    int GroupCreateNodeGroup()
    {
        return GroupCreate(NuTo::eGroupId::Nodes);
    }

    //! @brief Creates a group element ids for the structure
    //! @return rIdent identifier for the element group
    int GroupCreateElementGroup()
    {
        return GroupCreate(NuTo::eGroupId::Elements);
    }

    //! @brief ... Creates a group for the structure
    //! @param rType  type of the group, e.g. "NODES" or "ELEMENTS"
    //! @return ... rIdent identifier for the group
    int GroupCreate(NuTo::eGroupId rEnumType);

    //! @brief ... Creates a group for the structure
    //! @param rIdent identifier for the group
    //! @param rType  type of the group
    void GroupCreate(int id, NuTo::eGroupId rEnumType);
#endif

    //! @brief ... Creates a group for the structure
    //! @param rIdent identifier for the group
    //! @param rType  type of the group, e.g. "NODES" or "ELEMENTS"
    int GroupCreate(const std::string& rType);


    //! @brief ... Deletes a group from the structure
    //! @param rIdent identifier for the group
    void GroupDelete(int rIdent);

    //! @brief ... Unites two groups and stores the result in a new group
    //! @param rIdentGroup1 identifier for the first group
    //! @param rIdentGroup2 identifier for the second group
    //! @result ... rIdentGroupResult identifier for the created result group
    int GroupUnion(int rIdentGroup1, int rIdentGroup2);
    //! @brief ... Difference between two groups and stores the result in a new group
    //! @param rIdentGroup1 identifier for the first group
    //! @param rIdentGroup2 identifier for the second group
    //! @result ... rIdentGroupResult identifier for the created result group
    int GroupDifference(int rIdentGroup1, int rIdentGroup2);

    //! @brief ... Calculates the intersection between two groups and stores the result in a new group
    //! @param rIdentGroup1 identifier for the first group
    //! @param rIdentGroup2 identifier for the second group
    //! @result ... rIdentGroupResult identifier for the created result group
    int GroupIntersection(int rIdentGroup1, int rIdentGroup2);

    //! @brief ... Calculates the symmetric difference between two groups and stores the result in a new group
    //! @param rIdentGroup1 identifier for the first group
    //! @param rIdentGroup2 identifier for the second group
    //! @result ... rIdentGroupResult identifier for the created result group
    int GroupSymmetricDifference(int rIdentGroup1, int rIdentGroup2);

    //! @brief ... Adds a node to a node group
    //! @param rIdentGroup identifier for the group
    //! @param rIdentNode  identifier for the node
    void GroupAddNode(int rIdentGroup, int rIdNode);

    //! @brief ... Adds all nodes to a group whose coordinates are in the specified range
    //! @param rIdentGroup identifier for the group
    //! @param rDirection either 0,1,2 for x,y, or z
    //! @param rMin ... minimum value
    //! @param rMax ... maximum value
    virtual void GroupAddNodeCoordinateRange(int rIdentGroup, int rDirection, double rMin, double rMax);

    //! @brief Adds all nodes to a group whose coordinates are in the specified range
    //! @param rIdentGroup identifier for the group
    //! @param rDirection either X, Y, Z
    //! @param rMin ... minimum value
    //! @param rMax ... maximum value
    virtual void GroupAddNodeCoordinateRange(int rIdentGroup, eDirection rDirection, double rMin, double rMax);

    //! @brief creates a new node group and selects all nodes whose coordinates are in the specified range [min .. max]
    //! @param direction either X, Y, Z
    //! @param min minimum value
    //! @param max maximum value
    //! return reference to the node group
    //! @remark This should be rewritten in a method that returns a Group by value without storing it
    //! at all. But this collides with some other features (mainly postprocessing) that rely on a
    //! "GroupId". Here, you could obtain the GroupId by calling GroupGetId() if you need it.
    Group<NodeBase>& GroupGetNodeCoordinateRange(eDirection direction, double min, double max);


    //! @brief creates a new node group and selects all nodes whose coordinates at value (+- tolerance)
    //! Equal to GroupGetNodeCoordinateRange(direction, value - tolerance, value + tolerance)
    //! Example: direction = Y, value = 42. --> selects all nodes with Y = 42
    //! @param direction either X, Y, Z
    //! @param value specific coordinate value
    //! @param tolerance tolerance -+
    //! return reference to the node group
    //! @remark This should be rewritten in a method that returns a Group by value without storing it
    //! at all. But this collides with some other features (mainly postprocessing) that rely on a
    //! "GroupId". Here, you could obtain the GroupId by calling GroupGetId() if you need it.
    Group<NodeBase>& GroupGetNodesAtCoordinate(eDirection direction, double value, double tolerance = 1.e-6);


#ifndef SWIG
    //! @brief ... Adds all nodes which fulfill the conditions specified in a std::function
    //! @param rIdentGroup identifier for the group
    //! @param rFunction std::function
    void GroupAddNodeFunction(int rIdentGroup, std::function<bool(NodeBase*)> rFunction);

    //! @brief ... Adds all nodes which fulfill the conditions specified in a std::function
    //! @param rIdentNewGroup identifier for the group where to add the nodes
    //! @param rIdentOldGroup identifier for the group where the ids are searched
    //! @param rFunction std::function
    void GroupAddNodeFunction(int rIdentNewGroup, int rIdentOldGroup, std::function<bool(NuTo::NodeBase*)> rFunction);
#endif

    //! @brief ... Adds an element to an element group
    //! @param rIdentGroup identifier for the group
    //! @param rIdentNode  identifier for the element
    void GroupAddElement(int rIdentGroup, int rIdElement);

    //! @brief ... Adds all elements to a group whose nodes are in the given node group
    //! @param rElementGroupId identifier for the element group
    //! @param rNodeGroupId identifier for the node group
    //! @param rHaveAllNodes if set to true, the element is only selected when all element nodes are in the node group,
    //! if set
    //! to false, the element is select if at least one node is in the node group
    void GroupAddElementsFromNodes(int rElementGroupId, int rNodeGroupId, bool rHaveAllNodes);

    //! @brief ... Adds all the nodes from the group-rElementGroupId to the group rNodeGroupId
    //! @param rNodeGroupId id for the node group
    //! @param rElementGroupId id for the element group
    void GroupAddNodesFromElements(int rNodeGroupId, int rElementGroupId);


    //! @brief ... Adds all the nodes from the group-rElementGroupId to a new group
    //! @param rElementGroupId id for the element group
    //! @return rNodeGroupId id for the node group
    int GroupCreateNodeGroupFromElements(int rElementGroupId);

    //! @brief ... Adds all nodes to a group whose coordinates are in the specified range
    //! @param rIdentGroup identifier for the group
    //! @param rCenter center of the selection circle
    //! @param rMin ... minimum radius
    //! @param rMax ... maximum radius
    void GroupAddNodeRadiusRange(int rIdentGroup, Eigen::VectorXd rCenter, double rMin, double rMax);

    Group<NodeBase>& GroupGetNodeRadiusRange(Eigen::VectorXd center, double min = 0., double max = 1.e-6);

    //! @brief ... Adds all nodes to a group whose coordinates are on a cylinder with the radius in the in the specified
    //! range
    //! @param rIdentGroup identifier for the group
    //! @param rCenter center of the cylinder
    //! @param rAxis axis of the cylinder
    //! @param rMin ... minimum radius
    //! @param rMax ... maximum radius
    void GroupAddNodeCylinderRadiusRange(int rIdentGroup, Eigen::VectorXd rCenter, Eigen::VectorXd rDirection,
                                         double rMin, double rMax);

    //! @brief ... Returns the number of members in a group
    //! @param rIdentGroup identifier for the group
    //! @return ... number of members
    int GroupGetNumMembers(int rIdentGroup) const;

    //! @brief ... Returns a vector with the members of a group
    //! @param rIdentGroup identifier for the group
    //! @return ... vector of members
    std::vector<int> GroupGetMemberIds(int rIdentGroup) const;

    //! @brief ... checks for a member in a group
    //! @param rIdentGroup identifier for the group
    //! @return ... rMember id (element id, node id etc.)
    bool GroupContainsMember(int rIdentGroup, int rMember) const;

    //*************************************************************
    //************ Integration type routines     ******************
    //**  defined in structures/StructureBaseIntegrationType.cpp **
    //*************************************************************
    //! @brief ... Info routine that prints general information about the allocated integration types
    //! an integration type is only allocated if required (from created elements)
    //! @param rVerboseLevel describes how detailed the information is
    void IntegrationTypeInfo(int rVerboseLevel) const;

#ifndef SWIG
    //! @brief ... Returns a pointer to an integration type
    //! if the integration type does not exist in the map, an exception is thrown
    //! @param identIntegrationType Identifier for an integration type
    NuTo::IntegrationTypeBase* GetPtrIntegrationType(std::string rIdentIntegrationType);

    //! @brief ... Returns a pointer to an integration type
    //! if the integration type does not exist (in the map), the integration type is created
    //! @param identIntegrationType Identifier for an integration type
    NuTo::IntegrationTypeBase* GetPtrIntegrationType(NuTo::eIntegrationType rIdentIntegrationType);
#endif // SWIG

    //*************************************************
    //************    Logger routines    ***************
    //*************************************************
    //! @brief opens a logger file for the output to a log file
    //! @param rFileName file name
    void LoggerOpenFile(std::string rFileName);

    //! @brief set the logger to be quiet (output only to file, if set)
    //! @param rQuiet (true for quiet logger, false for output to standard output)
    void LoggerSetQuiet(bool rQuiet);

#ifndef SWIG
    //! @brief returns the logger
    inline Logger& GetLogger() const
    {
        return mLogger;
    }
#endif // SWIG


    //*************************************************
    //************ Basic routines     ***************
    //**  defined in structures/StructureBase.cpp **
    //*************************************************
    //! @brief ... number of time derivatives (0 : static, 1: velocities, 2: accelerations)
    void SetNumTimeDerivatives(int rNumTimeDerivatives);

    //! @brief ... return number of time derivatives (0 : static, 1: velocities, 2: accelerations)
    int GetNumTimeDerivatives() const;

    //! @brief ... Info routine that prints general information about the object (detail according to verbose level)
    virtual void Info() const;

    //! @brief absolute tolerance for entries of the global stiffness matrix (coefficientMatrix0)
    //! values smaller than that one will not be added to the global matrix
    void SetToleranceStiffnessEntries(double rToleranceStiffnessEntries);

    //! @brief absolute tolerance for entries of the global stiffness matrix (coefficientMatrix0)
    //! values smaller than that one will not be added to the global matrix
    double GetToleranceStiffnessEntries() const;

    //! @brief returns the number of degrees of freedom
    //! @return ... number of degrees of freedom
    int GetNumTotalDofs() const;

    //! @brief returns the number of active degrees of freedom
    //! @return ... number of active degrees of freedom
    int GetNumTotalActiveDofs() const;

    //! @brief returns the number of dependent degrees of freedom
    //! @return ... number of dependent degrees of freedom
    int GetNumTotalDependentDofs() const;

#ifndef SWIG
    //! @brief returns the number of degrees of freedom for dof type rDof
    //! @param rDofType ... dof type
    //! @return ... number of degrees of freedom
    int GetNumDofs(Node::eDof rDofType) const;

    //! @brief returns the number of active degrees of freedom for dof type rDof
    //! @param rDofType ... dof type
    //! @return ... number of active degrees of freedom
    int GetNumActiveDofs(Node::eDof rDofType) const;

    //! @brief returns the number of dependent degrees of freedom for dof type rDof
    //! @param rDofType ... dof type
    //! @return ... number of dependent degrees of freedom
    int GetNumDependentDofs(Node::eDof rDofType) const;

    //! @brief returns all dof types of the structure
    //! @return ... set of active dof types
    std::set<Node::eDof> DofTypesGet() const;

    //! @brief Sets all dofs inactive
    void DofTypeDeactivateAll();

    //! @brief Sets all dofs active
    void DofTypeActivateAll();

    //! @brief forwards the property to all interpolation types. Inactive dofs are not solved for
    //! @param rDofType ... dof type
    //! @param rIsActive ... active/inactive
    void DofTypeSetIsActive(Node::eDof rDofType, bool rIsActive);

    //! @brief Sets a set of active dofs, deactivates others
    //! @param rActiveDofTypes ... dof type
    void DofTypeSetIsActive(const std::set<Node::eDof>& rActiveDofTypes);

    //! @brief returns true if rDofType is an active dof type
    //! @param rDofType ... dof type
    bool DofTypeIsActive(Node::eDof rDofType) const;

    //! @brief forwards the property to all interpolation types
    //! @param rDofType ... dof type
    //! @param rIsConstitutiveInput ... is/is not constitutive input
    void DofTypeSetIsConstitutiveInput(Node::eDof rDofType, bool rIsConstitutiveInput);

    //! @brief adds/removes the dof type from the mDofTypesSymmetric-set
    //! @param rDofType ... dof type
    //! @param rIsSymmetric ... is/is not symmetric input
    //! TODO: move to constitutive law somehow
    void DofTypeSetIsSymmetric(Node::eDof rDofType, bool rIsSymmetric);

    //! @brief returns if rDof is symmetric
    //! @param rDofType ... dof type
    //! @return mSymmetricDof contains rDof
    bool DofTypeIsSymmetric(Node::eDof rDofType) const;

    //! @brief returns all active dof types
    //! @return ... set of active dof types
    std::set<Node::eDof> DofTypesGetActive() const;

    //! @brief updates the mDofStatus with current information from the interpolation types
    void UpdateDofStatus();

#endif // SWIG

    /// !!!!!! IMPORTANT ------------------------------------------
    /// DofStatus-Geter should never get a non-const version
    /// Use Structure to change dofStatus
    const NuTo::DofStatus& GetDofStatus() const;


    //! @brief returns the number of degrees of freedom for dof type rDof
    //! @param rDofType ... dof type
    //! @return ... number of degrees of freedom
    int GetNumDofs(std::string rDofType) const;

    //! @brief returns the number of active degrees of freedom for dof type rDof
    //! @param rDofType ... dof type
    //! @return ... number of active degrees of freedom
    int GetNumActiveDofs(std::string rDofType) const;

    //! @brief returns the number of dependent degrees of freedom for dof type rDof
    //! @param rDofType ... dof type
    //! @return ... number of dependent degrees of freedom
    int GetNumDependentDofs(std::string rDofType) const;

    //! @brief adds/removes the dof type from the mDofTypesActive-set. Inactive dofs are not solved for
    //! @param rDofType ... dof type
    //! @param rIsActive ... active/inactive
    void DofTypeSetIsActive(std::string rDofType, bool rIsActive);

    //! @brief adds/removes the dof type from the mDofTypesConstitutiveInput-set
    //! @param rDofType ... dof type
    //! @param rIsConstitutiveInput ... is/is not constitutive input
    void DofTypeSetIsConstitutiveInput(std::string rDofType, bool rIsConstitutiveInput);

    //! @brief writes the current state of the structure (nodal values, history variables and global time) as a binary
    //! file
    //! @param filename ... file name
    //! @param globalTime ... global time
    void WriteRestartFile(std::string filename, double globalTime);

    //! @brief reads a state of the structure (nodal values and history variables) from a restart file
    //! @param filename ... file name
    //! @return globalTime
    double ReadRestartFile(std::string filename);

    //! @brief defines the serialization of this class
    //! @param rStream serialize output stream
    virtual void NuToSerializeSave(SerializeStreamOut&){/* currently no members to serialize */};

    //! @brief defines the serialization of this class
    //! @param rStream serialize input stream
    virtual void NuToSerializeLoad(SerializeStreamIn&){/* currently no members to serialize */};

    //! @brief this routine is only relevant for the multiscale model, since an update on the fine scale should only be
    //! performed
    // for an update on the coarse scale
    // as a consequence, in an iterative solution with updates in between the initial state has to be restored after
    // leaving the routine
    // this routine saves the current state before an update in the Newton Raphson iteration is performed
    // this only happens for more than one load step (either prescibed or with automatic load control)
    virtual void SaveStructure(std::stringstream&) const
    {
        throw Exception(
                "[StructureBase::SaveStructure] Saving of the structure not implemented in derived class.");
    }

    virtual void RestoreStructure(std::stringstream&)
    {
        throw Exception(
                "[StructureBase::RestoreStructure] Saving of the structure not implemented in derived class.");
    }

    void SetUpdateTmpStaticDataRequired()
    {
        mUpdateTmpStaticDataRequired = true;
    }

    //! @brief calculate the critical time step for a vector of elements solving the generalized eigenvalue problem
    //! Ku=lambda Mu
    double ElementCalculateLargestElementEigenvalue(const std::vector<ElementBase*>& rElementVector);

    //! @brief returns whether or not the dof is constitutive input at least in one InrepolationType
    //! @param rInterpolationTypeId ... interpolation type id
    //! @param rDofType ... dof type
    bool InterpolationTypeIsConstitutiveInput(NuTo::Node::eDof rDofType);


    ///
    /// \brief DofStatusSetHasInteractingConstraints
    /// \param rHasInteractingConstraints
    ///
    /// Sets the member variable mHasInteractingConstraints of mDofStatus to either true or false.
    /// The mHasInteractingConstraints determines whether K_KJ and K_KK are assembled during the Evaluate function
    ///
    void DofStatusSetHasInteractingConstraints(bool rHasInteractingConstraints);

    bool GetShowTime() const;

    void SetShowTime(bool showTime);

    unsigned short GetVerboseLevel() const;

    void SetVerboseLevel(unsigned short verboseLevel);

    const Constraint::Constraints& Constraints() const;
    Constraint::Constraints& Constraints();

    const Assembler& GetAssembler() const
    {
        return *mAssembler;
    }
    Assembler& GetAssembler()
    {
        return *mAssembler;
    }

protected:
    //! @brief finds an unused ID in rMap
    //! @param rMap hopefully any kind of map, boost, std, ...
    //! @return unused ID
    template <typename T>
    int GetUnusedId(const T& rMap)
    {
        int unused = rMap.size();
        auto it = rMap.find(unused);
        while (it != rMap.end())
        {
            unused++;
            it = rMap.find(unused);
        }
        return unused;
    }


    //! @brief ... number of time derivatives (0 : static, 1: velocities, 2: accelerations)
    int mNumTimeDerivatives;

    int mDimension;

    std::unique_ptr<NuTo::Assembler> mAssembler;

    //! @brief ... map storing the name and the pointer to the constitutive law
    //! @sa ConstitutiveBase
    boost::ptr_map<int, ConstitutiveBase> mConstitutiveLawMap;

    //! @brief ... map storing node loads
    //! @sa LoadBase
    boost::ptr_map<int, LoadBase> mLoadMap;

    //! @brief ... map storing the groups and a pointer to the objects
    //! @sa GroupBase
    boost::ptr_map<int, GroupBase> mGroupMap;

    //! @brief ... map storing the name and the pointer to the integration types
    //! @sa IntegrationTypeBase
    boost::ptr_map<std::string, IntegrationTypeBase> mIntegrationTypeMap;

    //! @brief ... map storing the interpolation types
    //! @sa InterpolationType
    boost::ptr_map<int, InterpolationType> mInterpolationTypeMap;

    //! @brief ... map storing the components (displacements, strains, nonlocal weights etc) to be included in the
    //! output (VTK) file
    std::map<int, std::vector<eVisualizeWhat>> mGroupVisualizeComponentsMap;

    //! @brief ... map storing the type of visualization for the output (VTK) file
    std::map<int, eVisualizationType> mGroupVisualizationType;

    //! @brief is set to true, if at least one constitutive model requires an update of tmpStaticData before stress and
    //! stiffness routines are called
    bool mHaveTmpStaticData;

    //! @brief is set to false, if the structure is changed (nodes, elements) or (DOFs at the nodes)
    bool mUpdateTmpStaticDataRequired;

    //! @brief absolute tolerance for entries of the global stiffness matrix (coefficientMatrix0)
    //! values smaller than that one will not be added to the global matrix
    double mToleranceStiffnessEntries;

#ifdef _OPENMP
    //@brief maximum independent sets used for parallel assembly of the stiffness resforce etc.
    mutable std::vector<std::vector<ElementBase*>> mMIS;
    //@brief number of processors used in an openmp simulation
    int mNumProcessors;
#endif

    //! @brief logger class to redirect the output to some file or the console (or both), can be changed even for const
    //! routines
    mutable Logger mLogger;

    bool mShowTime;

    unsigned short mVerboseLevel;


#ifndef SWIG
    //! @brief ... store all elements of a structure in a vector
    //! @param rElements ... vector of element pointer
    virtual void GetElementsTotal(std::vector<const ElementBase*>& rElements) const = 0;

    //! @brief ... store all elements of a structure in a vector
    //! @param rElements ... vector of element pointer
    virtual void GetElementsTotal(std::vector<std::pair<int, const ElementBase*>>& rElements) const = 0;

    //! @brief ... store all elements of a structure in a vector
    //! @param rElements ... vector of element pointer
    virtual void GetElementsTotal(std::vector<ElementBase*>& rElements) = 0;

    //! @brief ... store all elements of a structure in a vector
    //! @param rElements ... vector of element pointer
    virtual void GetElementsTotal(std::vector<std::pair<int, ElementBase*>>& rElements) = 0;

    //! @brief ... store all nodes of a structure in a vector
    //! @param rNodes ... vector of element pointer
    virtual void GetNodesTotal(std::vector<const NodeBase*>& rNodess) const = 0;

    //! @brief ... store all nodes of a structure in a vector
    //! @param rNodes ... vector of element pointer
    virtual void GetNodesTotal(std::vector<std::pair<int, const NodeBase*>>& rNodes) const = 0;

    //! @brief ... store all nodes of a structure in a vector
    //! @param rNodes ... vector of element pointer
    virtual void GetNodesTotal(std::vector<NodeBase*>& rNodes) = 0;

    //! @brief ... store all nodes of a structure in a vector
    //! @param rNodes ... vector of element pointer
    virtual void GetNodesTotal(std::vector<std::pair<int, NodeBase*>>& rNodes) = 0;
#endif

    //! @brief ... get all elements of a group in a vector
    //! @param rElementGroup ... element group
    //! @param rElements ... vector of element pointer
    void GetElementsByGroup(Group<ElementBase>* rElementGroup, std::vector<ElementBase*>& rElements);
};
} // namespace NuTo
