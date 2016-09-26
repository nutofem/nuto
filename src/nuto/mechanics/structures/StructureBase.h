// $Id$

#pragma once

#include <functional>
#include <list>
#include <map>
#include <memory>
#include <set>
#include <string>

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION


// parent class
#include "nuto/base/NuToObject.h"

// class member
#include <boost/ptr_container/ptr_map.hpp>
#include "nuto/base/Logger.h"
#include "nuto/math/FullVector_Def.h"
#include "nuto/mechanics/dofSubMatrixStorage/BlockFullVector.h"
#include "nuto/mechanics/dofSubMatrixStorage/BlockSparseMatrix.h"
#include "nuto/mechanics/dofSubMatrixStorage/DofStatus.h"


#include "nuto/mechanics/MechanicsException.h"



namespace NuTo
{
class ConstitutiveBase;
class ConstitutiveStaticDataMultiscale2DPlaneStrain;
class ConstraintBase;
class CrackBase;
class ElementBase;
class EngineeringStrain2D;
class GroupBase;
class IntegrationTypeBase;
class InterpolationBase;
class InterpolationType;
class LoadBase;
class NewtonRaphsonAuxRoutinesBase;
class NodeBase;
class SectionBase;
class StructureOutputBase;
class StructureOutputBlockMatrix;
class StructureOutputBlockVector;
class VisualizeComponent;
class VisualizeUnstructuredGrid;
template<typename T> class BlockFullMatrix;
template<typename T> class BlockFullVector;
template<typename IOEnum> class ConstitutiveIOMap;
template<class T> class Group;
template<class T, int rows, int cols> class FullMatrix;
template<class T, int rows> class FullVector;
template<class T> class SparseMatrixCSRSymmetric;
template<class T> class SparseMatrixCSRGeneral;
template<class T> class SparseMatrixCSRVector2General;
template<class T> class SparseMatrixCSRVector2Symmetric;

enum class eError;
enum class eGroupId;
enum class eIntegrationType;
enum class eSectionType;
enum class eStructureOutput;
enum class eVisualizationType;
enum class eVisualizeWhat;

namespace Constitutive
{
    enum class eConstitutiveParameter;
    enum class eConstitutiveType;
    enum class eDamageLawType;
    enum class eInput;
    enum class eOutput;
}// namespace Constitutive


namespace Element
{
    enum class eOutput;
}// namespace Element

namespace IpData
{
    enum class eIpDataType;
    enum class eIpStaticDataType;
}// namespace IpData


using ConstitutiveInputMap = ConstitutiveIOMap<Constitutive::eInput>;
using ConstitutiveOutputMap = ConstitutiveIOMap<Constitutive::eOutput>;


//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... standard abstract class for all mechanical structures
class StructureBase : public NuToObject
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
    friend class NewmarkIndirect;
    friend class NewmarkDirect;
    friend class VelocityVerlet;
public:
    //! @brief constructor
    //! @param mDimension  Structural dimension (1,2 or 3)
    StructureBase(int mDimension);

    //! @brief deconstructor
    virtual ~StructureBase();

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);


    //! @brief ... save the object to a file
    //! @param filename ... filename
    //! @param aType ... type of file, either BINARY, XML or TEXT
    virtual void SaveUpdate (const std::string &filename, std::string rType ) const
    {
        (void)filename;
        (void)rType;
    }

    //! @brief ... restore the object from a file
    //! @param filename ... filename
    //! @param aType ... type of file, either BINARY, XML or TEXT
    virtual void RestoreUpdate (const std::string &filename, std::string rType )
    {
        (void)filename;
        (void)rType;
    }

#endif  // ENABLE_SERIALIZATION

    //! @brief gives the dimension of the Structure
    //! @return Structural dimension (1,2 or 3)
    int GetDimension()const;

    //! @brief ... clear all visualization components
    void ClearVisualizationComponents();

    //! @brief ... export the entire structure to Vtk data file
    //! @param rResultFileName ... file name
    //! @param rTimeStep ... time step for the output files
    //! @param rXML ... if true, a vtu file is exported, otherwise a legacy vtk file is produced
    void ExportVtkDataFileElements(const std::string& rResultFileName, bool rXML = false);

    //! @brief ... export the entire structure to Vtk data file
    //! @param rResultFileName ... file name
    //! @param rTimeStep ... time step for the output files
    //! @param rXML ... if true, a vtu file is exported, otherwise a legacy vtk file is produced
    void ExportVtkDataFileNodes(const std::string& rResultFileName, bool rXML = false);

    //! @brief ... export an element group to Vtk/xml data file
    //! @param rGroupIdent ... group ident
    //! @param rResultFileName ... file name
    //! @param rTimeStep ... time step for the output files
    //! @param rXML ... use xml or vtk format (true for xml)
    void ElementGroupExportVtkDataFile(int rGroupIdent, const std::string& rResultFileName, bool rXML);

    //! @brief Add rVisualizeComponent to an element group for the visualization
    //! @param rElementGroup: element group
    //! @param rVisualizeComponent: visualization component, i.e. displacements, stresses...
    void AddVisualizationComponent(int rElementGroup, const std::string& rVisualizeComponent);

    //! @brief Add nonlocal weights to an element group for the visualization
    //! @param rElementGroup: element group
    //! @param rElementId: element id
    //! @param rIp: integration point
    void AddVisualizationComponentNonlocalWeights(int rElementGroup, int rElementId, int rIp);

#ifndef SWIG
    //! @brief Add rVisualizeComponent to an element group for the visualization
    //! @param rElementGroup: element group
    //! @param rVisualizeComponent: visualization component, i.e. displacements, stresses...
    void AddVisualizationComponent(int rElementGroup, eVisualizeWhat rVisualizeComponent);

    //! @brief Set tje visualization type for an element group
    //! @param rElementGroup: element group
    //! @param rVisualizeComponent: visualization type, i.e. voronoi cell, extrapolated...
    void SetVisualizationType(const int rElementGroup, const eVisualizationType rVisualizationType);

    //! @brief ... define the data sets (scalar, vector etc for the visualize routine based on the mVisualizecomponents for an element plot
    void DefineVisualizeElementData(VisualizeUnstructuredGrid& rVisualize, const std::list<std::shared_ptr<NuTo::VisualizeComponent>>& rVisualizationList)const;

    //! @brief ... define the data sets (scalar, vector etc for the visualize routine based on the mVisualizecomponents for a node plot
    void DefineVisualizeNodeData(VisualizeUnstructuredGrid& rVisualize, const std::list<std::shared_ptr<NuTo::VisualizeComponent>>& rVisualizationList)const;

    //! @brief ... adds all the elements in the vector to the data structure that is finally visualized
    void ElementVectorAddToVisualize(VisualizeUnstructuredGrid& rVisualize, const std::list<std::shared_ptr<NuTo::VisualizeComponent>>& rVisualizationList, const std::vector<ElementBase*>& rElements);

    //! @brief ... adds all the elements in the vector to the data structure that is finally visualized
    void ElementVectorAddToVisualize(VisualizeUnstructuredGrid& rVisualize, const std::list<std::shared_ptr<NuTo::VisualizeComponent>>& rVisualizationList, const std::vector<ElementBase*>& rElements, const eVisualizationType rVisualizationType);

    //! @brief ... adds all the elements in the vector to the data structure that is finally visualized
    void ElementTotalAddToVisualize(VisualizeUnstructuredGrid& rVisualize, const std::list<std::shared_ptr<NuTo::VisualizeComponent>>& rVisualizationList);

    //! @brief ... adds all the elements in a group to the data structure that is finally visualized
    void ElementGroupAddToVisualize(int rGroupId, VisualizeUnstructuredGrid& rVisualize, const std::list<std::shared_ptr<NuTo::VisualizeComponent>>& rVisualizationList);

    //! @brief ... returns the map that contains the visualization components to be exported for each element group
    std::map<int, std::list<std::shared_ptr<VisualizeComponent>>>& GetGroupVisualizeComponentsMap(void);
#endif //SWIG


#ifndef SWIG

    //! @brief ... evaluates the structure
    virtual void Evaluate(const NuTo::ConstitutiveInputMap& rInput, std::map<eStructureOutput, StructureOutputBase*> &rStructureOutput);


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

#endif //SWIG

    NuTo::StructureOutputBlockMatrix BuildGlobalHessian0();
    NuTo::StructureOutputBlockMatrix BuildGlobalHessian1();
    NuTo::StructureOutputBlockMatrix BuildGlobalHessian2();
    NuTo::StructureOutputBlockMatrix BuildGlobalHessian2Lumped();

    NuTo::StructureOutputBlockVector BuildGlobalInternalGradient();

    //! @brief ... build global external load vector (currently for displacements only)
    //! @param rLoadCase ... load case
    //! @return  ... StructureOutputBlockVector containing the external loads
    NuTo::StructureOutputBlockVector BuildGlobalExternalLoadVector(int rLoadCase);

    NuTo::BlockFullVector<double> SolveBlockSystem(const NuTo::BlockSparseMatrix& rMatrix, const NuTo::BlockFullVector<double>& rVector) const;


    void SolveGlobalSystemStaticElastic(int rLoadCase = 0);

    NuTo::StructureOutputBlockMatrix BuildGlobalHessian0_CDF(double rDelta);

    bool CheckHessian0(double rDelta, double rRelativeTolerance, bool rPrintWrongMatrices = true);

private: // keep the method definitions close together
    bool CheckHessian0_Submatrix(const BlockSparseMatrix& rHessian0, BlockSparseMatrix& rDiff, double rRelativeTolerance, bool rPrintWrongMatrices);
public:

//*************************************************
//************ Node routines        ***************
//***  defined in structures/StructureNode.cpp  ***
//*************************************************
    //! @brief returns the number of nodes
    //! @return number of nodes
    virtual int GetNumNodes() const =0;

#ifndef SWIG
    //! @brief a reference to a node
    //! @param identifier
    //! @return reference to a node
    virtual NodeBase* NodeGetNodePtr(int rIdent)=0;

    //! @brief a reference to a node
    //! @param identifier
    //! @return reference to a node
    virtual const NodeBase* NodeGetNodePtr(int rIdent)const=0;

    //! @brief ... store all elements connected to this node in a vector
    //! @param rNode (Input) 		... node pointer
    //! @param rElements (Output) 	... vector of element pointers
    virtual void NodeGetElements(const NodeBase* rNodePtr, std::vector<ElementBase*>& rElements);

    //! @brief gives the identifier of a node
    //! @param reference to a node
    //! @return identifier
    virtual int NodeGetId(const NodeBase* rNode)const=0;
#endif //SWIG
    //! @brief ... returns the (first) node that has the specified coordinates within the range
    //! @param ... rCoordinates
    //! @param ... rRange
    //! @return ... node id
    int NodeGetIdAtCoordinate(FullVector<double, Eigen::Dynamic> rCoordinates, double rRange);

    //! @brief ... store all elements connected to this node in a vector
    //! @param rNodeId (Input) 			... node id
    //! @param rElementNumbers (Output) ... vector of element ids
    void NodeGetElements(const int rNodeId, NuTo::FullVector<int,Eigen::Dynamic>& rElementNumbers);

    //! @brief delete node
    //! @param rIdent ... node identifier
    virtual void NodeDelete(const int rIdent)=0;

    //! @brief info about the nodes in the Structure
    virtual void NodeInfo(int mVerboseLevel)const=0;

    //! @brief numbers the dofs in the structure
    //! @param rCallerName ... if the method throws it is nice to know by whom it was called.
    virtual void NodeBuildGlobalDofs(std::string rCallerName = "")=0;

    //! @brief sets the displacements of a node
    //! @param rIdent node identifier
    //! @param rDisplacements matrix (one column) with the displacements
    void NodeSetDisplacements(int rId,const NuTo::FullVector<double,Eigen::Dynamic>& rDisplacements);

    //! @brief sets the displacements of a node
    //! @param rIdent node identifier
    //! @param rTimeDerivative time derivative (0 disp, 1 vel, 2 acc)
    //! @param rDisplacements matrix (one column) with the displacements
    void NodeSetDisplacements(int rId, int rTimeDerivative, const NuTo::FullVector<double,Eigen::Dynamic>& rDisplacements);

    //! @brief sets the displacements of a node
    //! @param rIdent node identifier
    //! @param rRotations matrix (one column) with the rotations
    void NodeSetRotations(int rId,const NuTo::FullVector<double,Eigen::Dynamic>& rRotations);

    //! @brief sets the displacements of a group of nodes
    //! @param rIdent node group identifier
    //! @param rTimeDerivative time derivative (0 disp, 1 vel, 2 acc)
    //! @param rDisplacements matrix (one column) with the displacements
    void NodeGroupSetDisplacements(int rIdent, int rTimeDerivative, const NuTo::FullVector<double,Eigen::Dynamic>& rDisplacements);

    //! @brief sets the displacements of a group of nodes
    //! @param rIdent node group identifier
    //! @param rDisplacements matrix (one column) with the displacements
    void NodeGroupSetDisplacements(int rIdent, const NuTo::FullVector<double,Eigen::Dynamic>& rDisplacements);

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
    void NodeGroupGetMembers(int rGroupId, NuTo::FullVector<int,Eigen::Dynamic>& rMembers);

    //! @brief gets the coordinates of a node
    //! @param rNode node identifier
    //! @param rCoordinates matrix (one column) with the coordinates
    void NodeGetCoordinates(int rNode, NuTo::FullVector<double,Eigen::Dynamic>& rCoordinates)const;

    //! @brief gets the coordinates of a group of nodes (be careful, the order of the nodes in a group might change between different runs)
    //! @param rNodeGroup node group identifier
    //! @param rCoordinates matrix (rows/nodes columns/coordinates)
    void NodeGroupGetCoordinates(int rNodeGroup, NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rCoordinates);

    //! @brief gets the displacements of a node
    //! @param rNode node identifier
    //! @param rDisplacements matrix (one column) with the displacements
    void NodeGetDisplacements(int rNode, NuTo::FullVector<double,Eigen::Dynamic>& rDisplacements)const;

    //! @brief gets the displacements of a node
    //! @param rIdent node identifier
    //! @param rTimeDerivative time derivative (0 disp, 1 velocity,2 acceleration)
    //! @param rDisplacements matrix (one column) with the displacements
    void NodeGetDisplacements(int rNode, int rTimeDerivative, FullVector<double,Eigen::Dynamic>& rDisplacements)const;

    //! @brief gets the displacement dofs of a node
    //! @param rIdent node identifier
    //! @param rDisplacements matrix (one column) with the displacements
    void NodeGetDisplacementDofs(int rNode, FullVector<int,Eigen::Dynamic>& rDisplacementDofs)const;

    //! @brief gets the rotations of a node
    //! @param rNode node identifier
    //! @param rRotation matrix (one column) with the rotations
    void NodeGetRotations(int rNode, NuTo::FullVector<double,Eigen::Dynamic>& rRotations)const;

    //! @brief gets the global nonlocal eq plastic strain variables of a node
    //! @param rNode node identifier
    //! @return global (nodal) nonlocal eq plastic strain
    void NodeGetNonlocalEqPlasticStrain(int rNode, NuTo::FullVector<double,Eigen::Dynamic>& rNonlocalEqPlasticStrain)const;

    //! @brief gets the global nonlocal total strain variables of a node
    //! @param rNode node identifier
    //! @return global (nodal) nonlocal total strain
    void NodeGetNonlocalTotalStrain(int rNode, NuTo::FullVector<double,Eigen::Dynamic>& rNonlocalTotalStrain)const;

    //! @brief gets the displacements of a group of nodes (be careful, the order of the nodes in a group might change between different runs)
    //! @param rNodeGroup node group identifier
    //! @param rDisplacements matrix (rows/nodes columns/rDisplacements)
    void NodeGroupGetDisplacements(int rNodeGroup, NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rDisplacements);

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
    //! @param rActiveDofValues ... vector of independent dof values (ordering according to global dofs, size is number of active dofs)
    //! @param rDependentDofValues ... vector of dependent  dof values (ordering according to global dofs, size is number of active dofs)
    virtual void NodeMergeDofValues(int rTimeDerivative, const NuTo::BlockFullVector<double>& rActiveDofValues, const NuTo::BlockFullVector<double>& rDependentDofValues) = 0;

    //! @brief write dof values (e.g. displacements, temperatures to the nodes)
    //! @param rTimeDerivative time derivative (0 disp 1 vel 2 acc)
    //! @param rDofValues ... StructureBlockVector containing the dofs (J and K)
    void NodeMergeDofValues(int rTimeDerivative, const NuTo::StructureOutputBlockVector& rDofValues);

    //! @brief write dof values (e.g. displacements, temperatures to the nodes)
    //! @param rActiveDofValues ... vector of independent dof values (ordering according to global dofs, size is number of active dofs)
    //! @param rDependentDofValues ... vector of dependent  dof values (ordering according to global dofs, size is number of active dofs)
    inline void NodeMergeDofValues(NuTo::StructureOutputBlockVector& rDofValues)
    {
        NodeMergeDofValues(0,rDofValues);
    }

    //! @brief calculate dependent dof values (for the zeroth time derivative)
    //! @param rActiveDofValues ... vector of independent dof values (ordering according to global dofs, size is number of active dofs)
    //! @return  ... vector of dependent  dof values (ordering according to global dofs, size is number of active dofs)
    virtual NuTo::BlockFullVector<double> NodeCalculateDependentDofValues(const NuTo::BlockFullVector<double>& rActiveDofValues) const = 0;

    //! @brief calculate the internal force vector for a node
    //! @param rId ... node id
    //! @param rGradientInternalPotential ...vector for all the dofs the corresponding internal force (return value)
    void NodeInternalForce(int rId, NuTo::FullVector<double,Eigen::Dynamic>& rNodeForce);

    //! @brief calculate the internal force vector for a node group of nodes
    //! @param rIdent ... group id
    //! @param rGradientInternalPotential ...vector for all the dofs the corresponding internal force (return value)
    void NodeGroupInternalForce(int rIdent, NuTo::FullVector<double,Eigen::Dynamic>& rNodeForce);

#ifndef SWIG
    //! @brief calculate the internal force vector for a node
    //! @param rNodePtr  node for which this has to be calculated
    //! @param rGradientInternalPotential ...vector for all the dofs the corresponding internal force (return value)
    void NodeInternalForce(const NodeBase* rNodePtr, NuTo::FullVector<double,Eigen::Dynamic>& rNodeForce);

    //! @brief ... adds all the elements in the vector to the data structure that is finally visualized
    void NodeTotalAddToVisualize(VisualizeUnstructuredGrid& rVisualize, const std::list<std::shared_ptr<NuTo::VisualizeComponent>>& rVisualizationList) const;

    //! @brief ... adds all the nodes in the vector to the data structure that is finally visualized
    void NodeVectorAddToVisualize(VisualizeUnstructuredGrid& rVisualize, const std::list<std::shared_ptr<NuTo::VisualizeComponent>>& rVisualizationList, const std::vector<const NodeBase*>& rNodes) const;

#endif //SWIG

//*************************************************
//************ Element routines     ***************
//**  defined in structures/StructureElement.cpp **
//*************************************************
    //! @brief returns the number of elements
    //! @return number of elements
    virtual int GetNumElements() const =0;

#ifndef SWIG
    //! @brief a reference to an element
    //! @param identifier
    //! @return reference to an element
    virtual ElementBase* ElementGetElementPtr(int rIdent)=0;

    //! @brief a reference to an element
    //! @param identifier
    //! @return reference to an element
    virtual const ElementBase* ElementGetElementPtr(int rIdent)const=0;

    //! @brief gives the identifier of an element
    //! @param reference to an element
    //! @return identifier
    virtual int ElementGetId(const ElementBase* rElement)const=0;

    //! @brief info about one single element
    //! @param rElement (Input) ... pointer to the element
    //! @param rVerboseLevel (Input) ... level of verbosity
    virtual void ElementInfo(const ElementBase* rElement, int rVerboseLevel)const=0;

#endif //SWIG

    //! @brief delete element
    //! @param rIdent ... element identifier
    virtual void ElementDelete(int rIdent)=0;

    //! @brief info about the elements in the Structure
    virtual void ElementInfo(int rVerboseLevel)const=0;

#ifndef SWIG

    //! @brief duilds the element hessian0 using the central differences of the internal force vector
    //! @param rElement ... element
    //! @param rDelta ... step of the central differences
    //! @return BlockFullMatrix containing the hessian0_cdf
    BlockFullMatrix<double> ElementBuildHessian0_CDF(ElementBase* rElement, double rDelta);

    //! @brief builds the element internal gradient
    //! @param rElement ... element
    //! @return BlockFullVector containing the internal gradient
    BlockFullVector<double> ElementBuildInternalGradient(ElementBase* rElement);

    //! @brief builds the element hessian
    //! @param rHessianType ... hessian0 hessian1 hessian2
    //! @param rElement ... element
    //! @return BlockFullMatrix containing the hessian
    BlockFullMatrix<double> ElementBuildHessian(Element::eOutput rHessianType, ElementBase* rElement);

    //! @brief builds the element hessian0
    //! @param rElement ... element
    //! @return BlockFullMatrix containing the hessian
    BlockFullMatrix<double> ElementBuildHessian0(ElementBase* rElement);

    //! @brief builds the element hessian1
    //! @param rElement ... element
    //! @return BlockFullMatrix containing the hessian
    BlockFullMatrix<double> ElementBuildHessian1(ElementBase* rElement);

    //! @brief builds the element hessian2
    //! @param rElement ... element
    //! @return BlockFullMatrix containing the hessian
    BlockFullMatrix<double> ElementBuildHessian2(ElementBase* rElement);

    //! @brief builds the element vector of global row dofs
    //! @param rElement ... element
    //! @return BlockFullVector containing the global row dofs
    BlockFullVector<int> ElementBuildGlobalDofsRow(ElementBase* rElement);

    //! @brief builds the element vector of global column dofs
    //! @param rElement ... element
    //! @return BlockFullVector containing the global column dofs
    BlockFullVector<int> ElementBuildGlobalDofsColumn(ElementBase* rElement);

    bool ElementCheckHessian0(ElementBase* rElement, double rDelta, double rRelativeTolerance, bool rPrintWrongMatrices = true);

#endif // SWIG

    NuTo::BlockFullMatrix<double> ElementBuildHessian0        (int rElementId);
    NuTo::BlockFullMatrix<double> ElementBuildHessian1        (int rElementId);
    NuTo::BlockFullMatrix<double> ElementBuildHessian2        (int rElementId);

    NuTo::BlockFullVector<double> ElementBuildInternalGradient(int rElementId) { return ElementBuildInternalGradient(ElementGetElementPtr(rElementId)); }
    NuTo::BlockFullVector<int>    ElementBuildGlobalDofsRow   (int rElementId) { return ElementBuildGlobalDofsRow   (ElementGetElementPtr(rElementId)); }
    NuTo::BlockFullVector<int>    ElementBuildGlobalDofsColumn(int rElementId) { return ElementBuildGlobalDofsColumn(ElementGetElementPtr(rElementId)); }

    bool ElementCheckHessian0(int rElementId, double rDelta, double rRelativeTolerance, bool rPrintWrongMatrices = true)
    {
        return ElementCheckHessian0(ElementGetElementPtr(rElementId), rDelta, rRelativeTolerance, rPrintWrongMatrices);
    }

    bool ElementCheckHessian0(double rDelta, double rRelativeTolerance, bool rPrintWrongMatrices = true);

    //! @brief modifies the constitutive law of a single element
    //! @param rElementIdent identifier for the element
    //! @param rConstitutiveLawIdent identifier for the material
    void ElementSetConstitutiveLaw(int rElementId, int rConstitutiveLawIdent);

    //! @brief sets the constitutive law of a single element
    //! @param rElementIdent identifier for the element
    //! @param rIp  id of integration point
    //! @param rConstitutiveLawIdent identifier for the material
    void ElementSetConstitutiveLaw(int rElementId,int rIp, int rConstitutiveLawIdent);

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

    //! @brief sets the constitutive law of a single ip at an element
    //! @param rElement element pointer
    //! @param rIp number of integration point
    //! @param rConstitutive material pointer
    void ElementSetConstitutiveLaw(ElementBase* rElement,int rIp, ConstitutiveBase* rConstitutive);
#endif //SWIG

     //! @brief modifies the section of a single element
    //! @param rElementIdent element number
    //! @param rSectionIdent identifier for the section
    void ElementSetSection(int rElementId, int rSectionId);

    //! @brief modifies the section of a group of elements
    //! @param rGroupIdent identifier for the group of elements
    //! @param rSectionId identifier for the section
    void ElementGroupSetSection(int rGroupIdent, int rSectionId);

    //! @brief modifies the section of a all elements
    //! @param rSectionIdent identifier for the section
    void ElementTotalSetSection(int rSectionId);

#ifndef SWIG
    //! @brief returns the enum of string identifier for an integration type
    //! @param rIpDataTypeStr string
    //! @return enum
    NuTo::IpData::eIpDataType ElementGetEnumIntegrationType(const std::string& rIpDataTypeStr);

    //! @brief modifies the section of a single element
    //! @param rElement element pointer
    //! @param rSection section
    void ElementSetSection(ElementBase* rElement, SectionBase* rSection);
#endif //SWIG

    //! @brief modifies the interpolation type of a single element
    //! @param rElementId ... element number
    //! @param rInterpolationTypeId ... interpolation type id
    void ElementSetInterpolationType(int rElementId, int rInterpolationTypeId);

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
    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> ElementGetStaticIPData(int rElementId, IpData::eIpStaticDataType rType);

#endif //SWIG

    //! @brief calculates static ip data
    //! @param rElemIdent  element number
    //! @param rType static ip data type
    //! @param rIPData matrix with (... x numIP), x varies depending on IPData type
    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> ElementGetStaticIPData(int rElementId, std::string rType);

    //! @brief calculates the engineering strain
    //! @param rElemIdent  element number
    //! @param rEngineerungStrain engineering strain (return value, always 6xnumIp matrix)
    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> ElementGetEngineeringStrain(int rElementId);

    //! @brief calculates the engineering plastic strain
    //! @param rElemIdent  element number
    //! @param rEngineerungStrain engineering plastic strain (return value, always 6xnumIp matrix)
    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> ElementGetEngineeringPlasticStrain(int rElementId);

    //! @brief calculates the engineering stress
    //! @param rElemIdent  element number
    //! @param rEingineeringStress Engineering Stress (return value, always 6xnumIp matrix)
    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> ElementGetEngineeringStress(int rElementId);

    //! @brief calculates the damage
    //! @param rElemIdent  identifier for the element
    //! @param rDamage (return value, always 1xnumIp matrix)
    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> ElementGetDamage(int rElementId);

    //! @brief calculates the global integration point coordinates
    //! @param rElemIdent  identifier for the element
    //! @param rCoordinates integration point coordinates (return value, always 3xnumIp matrix)
    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> ElementGetIntegrationPointCoordinates(int rElementId);


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
    NuTo::eError ElementTotalUpdateStaticData();

    //! @brief updates the temprory static data of a all elements
    //! its is a const function, since only mutuable data (instead of const) is updated (kind of temporary data)
    NuTo::eError ElementTotalUpdateTmpStaticData();

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
    void ElementTotalGetAverageStress(double rVolume, NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rEngineeringStress);

    //! @brief calculates the average stress
    //! @param rGroupId  group number
    //! @param rVolume  volume of the structure in 3D /area in 2D/ length in 1D
    //! this is a parameter of the model, since holes have to be considered (zero stress, but still nonzero area)
    //! @param rEngineeringStress  average stress (return value)
    void ElementGroupGetAverageStress(int rGroupId, double rVolume, NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rEngineeringStress);

    //! @brief calculates the average strain
    //! @param rVolume  volume of the structure in 3D /area in 2D/ length in 1D
    //! this is a parameter of the model, since holes have to be considered (zero stress, but still nonzero area)
    //! @param rEngineeringStrain  average strain (return value)
    void ElementTotalGetAverageStrain(double rVolume, NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rEngineeringStrain);

    //! @brief calculates the average strain
    //! @param rGroupId  group number
    //! @param rVolume  volume of the structure in 3D /area in 2D/ length in 1D
    //! this is a parameter of the model, since holes have to be considered (zero stress, but still nonzero area)
    //! @param rEngineeringStrain  average strain (return value)
    void ElementGroupGetAverageStrain(int rGroupId, double rVolume, NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rEngineeringStrain);

    //! @brief returns the element ids of an element group
    //! @param rGroupId  group number
    //! @param rMembers  return vector with element ids
    void ElementGroupGetMembers(int rGroupId, NuTo::FullVector<int,Eigen::Dynamic>& rMembers);

    //! @brief calculates the volume of the elements
    //! @param rGroupId  group number
    //! @return volume of the structure in 3D /area in 2D/ length in 1D
    double ElementGroupGetVolume(int rGroupId);

    //! @brief calculates the internal energy of the system
    //! @return total energy
    virtual double ElementTotalGetInternalEnergy();

    //! @brief calculates the total energy of the system
    //! @return total energy
    virtual double ElementGroupGetTotalEnergy(int rGroupId);

    //! @brief calculates the elastic energy of the system
    //! @return elastic energy
    double ElementTotalGetElasticEnergy();

    //! @brief calculate the largest element eigenvalue for all elements solving the generalized eigenvalue problem Ku=lambda Mu
    //! this is used for the estimation of the critical time step
    double ElementTotalCalculateLargestElementEigenvalue();

    //! @brief calculate the largest element eigenvalue for a group of elements solving the generalized eigenvalue problem Ku=lambda Mu
    //! this is used for the estimation of the critical time step
    double ElementGroupCalculateLargestElementEigenvalue(int rGroupId);

    //*************************************************
    //************ Constraint routines     ***************
    //**  defined in StructureBaseConstraints.cpp **
    //*************************************************
    //! @brief deletes a constraint equation
    //! @param ConstraintId constraint id
    void ConstraintDelete(int ConstraintId);

#ifndef SWIG
    //! @brief adds a constraint to the map
    //! @param ConstraintId constraint id
    //! @param
    void ConstraintAdd(int rConstraintId, NuTo::ConstraintBase* rConstraint);
#endif

    //! @brief releases a constraint, (remove from the list but don't delete it)
    //!@param rConstraintEquation id of the constraint equation
    //! @return ptr to constraint
    //! @brief releases a constraint, (remove from the list but don't delete it)
    NuTo::ConstraintBase* ConstraintRelease(int rConstraintId);

#ifndef SWIG

    //! @brief adds a displacement constraint equation for a node
    //! @param rDOFType Type of the DOF that should be constrained (displacements, relativehumidity etc.)
    //! @param rNode pointer to node
    //! @param rValue prescribed value (e.g. zero to fix a displacement to zero)
    //! @return integer id to delete or modify the constraint
    int ConstraintLinearSetNode(Node::eDof rDOFType,NodeBase* rNode, double rValue);

    //! @brief adds a displacement constraint equation for a node
    //! @param rDOFType Type of the DOF that should be constrained (displacements, relativehumidity etc.)
    //! @param rNode pointer to node
    //! @param rDirection direction of the constraint (in 2D a point with 2 entries, in 3D 3 entries, in 1D not used)
    //! @param rValue prescribed value (e.g. zero to fix a displacement to zero)
    //! @return integer id to delete or modify the constraint
    int ConstraintLinearSetNode(Node::eDof rDOFType,NodeBase* rNode, const NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rDirection, double rValue);

    //! @brief adds a displacement constraint equation for a node
    //! @param rNode pointer to node
    //! @param rDirection direction of the constraint (in 2D a point with 2 entries, in 3D 3 entries, in 1D not used)
    //! @param rValue prescribed value (e.g. zero to fix a displacement to zero)
    //! @return integer id to delete or modify the constraint
    int ConstraintLinearSetDisplacementNode(NodeBase* rNode, const NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rDirection, double rValue);
#endif

    //! @brief adds a displacement constraint equation for a node
    //! @param rNode identifier for node
    //! @param rComponent e.g. the first (count from zero) displacement component
    //! @param rValue prescribed value (e.g. zero to fix a displacement to zero)
    int ConstraintLinearSetDisplacementNode(int rIdent, const NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rDirection, double rValue);

    //! @brief adds a relative humidity constraint equation for node
    //! @param rNode pointer to node
    //! @param rValue prescribed value (e.g. zero to fix a displacement to zero)
    //! @return integer id to delete or modify the constraint
    int ConstraintLinearSetRelativeHumidityNode(NodeBase* rNode, double rValue);

    //! @brief adds a relative humidity constraint for a node
    //! @param rIdent identifier for node
    //! @param rValue prescribed value (e.g. zero to fix a displacement to zero)
    //! @return integer id to delete or modify the constraint
    int ConstraintLinearSetRelativeHumidityNode(int rIdent, double rValue);

#ifndef SWIG
    //! @brief adds a rotation constraint equation for a node
    //! @param rNode pointer to node
    //! @param rValue prescribed value (e.g. zero to fix a displacement to zero)
    //! @return integer id to delete or modify the constraint
    int ConstraintLinearSetRotationNode(NodeBase* rNode, double rValue);
#endif

    //! @brief adds a rotation constraint equation for a node
    //! @param rNode identifier for node
    //! @param rValue prescribed value (e.g. zero to fix a displacement to zero)
    int  ConstraintLinearSetRotationNode(int rIdent, double rValue);

    //! @brief Adds a temperature constraint equation for node.
    //! @param rNode Pointer to node
    //! @param rValue Prescribed value (e.g. zero to fix a displacement to zero)
    //! @return Integer id to delete or modify the constraint
    int ConstraintLinearSetTemperatureNode(NodeBase* rNode, double rValue);

    //! @brief Adds a relative humidity constraint for a node.
    //! @param rIdent Identifier for node
    //! @param rValue Prescribed value (e.g. zero to fix a displacement to zero)
    //! @return Integer id to delete or modify the constraint
    int ConstraintLinearSetTemperatureNode(int rIdent, double rValue);

    //! @brief adds a water volume fraction constraint for a node
    //! @param rNode pointer to node
    //! @param rValue prescribed value (e.g. zero to fix a displacement to zero)
    //! @return integer id to delete or modify the constraint
    int ConstraintLinearSetWaterVolumeFractionNode(NodeBase* rNode, double rValue);

    //! @brief adds a water volume fraction constraint for a node
    //! @param rIdent identifier for node
    //! @param rValue prescribed value (e.g. zero to fix a displacement to zero)
    //! @return integer id to delete or modify the constraint
    int ConstraintLinearSetWaterVolumeFractionNode(int rIdent, double rValue);

#ifndef SWIG
    //! @brief adds a displacement constraint equation for a group of node
    //! @param rNode pointer to group of nodes
    //! @param rDirection direction of the constraint (in 2D a point with 2 entries, in 3D 3 entries, in 1D not used)
    //! @param rValue prescribed value (e.g. zero to fix a displacement to zero)
    //! @return integer id to delete or modify the constraint
    int ConstraintLinearSetDisplacementNodeGroup(Group<NodeBase>* rGroup, const NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rDirection, double rValue);
#endif

    //! @brief adds a constraint equation for a group of nodes
    //! @param rGroupIdent identifier for group of nodes
    //! @param rDof displacements, rotations, temperatures
    //! @param rComponent e.g. the first (count from zero) displacement component
    //! @param rValue prescribed value (e.g. zero to fix a displacement to zero)
    int ConstraintLinearSetDisplacementNodeGroup(int rGroupIdent, const NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rDirection, double rValue);

#ifndef SWIG
    //! @brief adds a rotation constraint equation for a group of node
    //! @param rNode pointer to group of nodes
    //! @param rValue prescribed value (e.g. zero to fix a rotation to zero)
    //! @return integer id to delete or modify the constraint
    int ConstraintLinearSetRotationNodeGroup(Group<NodeBase>* rGroup, double rValue);
#endif

    //! @brief adds a constraint equation for a group of nodes
    //! @param rGroupIdent identifier for group of nodes
    //! @param rValue prescribed value (e.g. zero to fix a displacement to zero)
    int ConstraintLinearSetRotationNodeGroup(int rGroupIdent, double rValue);

#ifndef SWIG
    //! @brief adds a temperature constraint equation for a group of nodes
    //! @param rNode pointer to group of nodes
    //! @param rValue prescribed value (e.g. zero to fix a temperature to zero)
    //! @return integer id to delete or modify the constraint
    int ConstraintLinearSetTemperatureNodeGroup(Group<NodeBase>* rGroup, double rValue);
#endif

    //! @brief adds a constraint equation for a group of nodes
    //! @param rGroupIdent identifier for group of nodes
    //! @param rValue prescribed value (e.g. zero to fix a temperature to zero)
    int ConstraintLinearSetTemperatureNodeGroup(int rGroupIdent, double rValue);


    // ######################################
    // ##                                  ##
    // ##    CONSTRAINT CALCULATION        ##
    // ##                                  ##
    // ######################################

#ifndef SWIG
    //! @brief returns the number of constraint equations for a specific dof type
    //! @return number of constraints
    //! @param rDofType  dof type
    int ConstraintGetNumLinearConstraints(Node::eDof rDof) const;

#endif

    //! @brief returns the number of constraint equations for a specific dof type
    //! @return number of constraints
    //! @param rDofType  dof type
    int ConstraintGetNumLinearConstraints(std::string rDof) const;


    //! @brief calculates the constraint matrix that builds relations between the nodal degrees of freedom (before gauss elimination)
    //! @param rConstraintMatrix constraint matrix
    NuTo::BlockSparseMatrix ConstraintGetConstraintMatrixBeforeGaussElimination() const;

    //! @brief returns the constraint vector after gauss elimination
    //! rConstraintMatrix*DOFS = RHS
    //! @return rhs
    const NuTo::BlockFullVector<double>& ConstraintGetRHSAfterGaussElimination() const;

    //! @brief returns the constraint vector after gauss elimination
    //! rConstraintMatrix*DOFS = RHS
    NuTo::BlockFullVector<double> ConstraintGetRHSBeforeGaussElimination();


    //! @brief calculates the right hand side of the constraint equations based on the mapping matrix and the rhs before the gauss elimination
    //! the result is stored internally in mConstraintRHS
    void ConstraintUpdateRHSAfterGaussElimination();

    //!@brief sets/modifies the right hand side of the constraint equations
    //!@param rRHS new right hand side
    //!@param rRHS new right hand side
    void ConstraintSetRHS(int rConstraintEquation, double rRHS);

    //!@brief gets the right hand side of the constraint equations
    //!@param rConstraintEquation constraint equation
    //!@return rRHS
    double ConstraintGetRHS(int rConstraintEquation)const;


    //! @brief ... create a constraint equation
    //! @param rNode ... node id in the first constraint equation term
    //! @param rDof ... dof in the first constraint equation term (e.g "X_DISPLACEMENT", "Z_Rotation", "Temperature")
    //! @param rCoefficient ... weight factor of this term
    //! @param rRHS ... prescribed right hand side value
    //! @return integer id of the constraint
    int ConstraintLinearEquationCreate(int rNode, const std::string& rDof, double rCoefficient, double rRHS = 0);

    //! @brief ... create a constraint equation
    //! @param rConstraint ... constraint id
    //! @param rNode ... node id in the first constraint equation term
    //! @param rDof ... dof in the first constraint equation term (e.g "X_DISPLACEMENT", "Z_Rotation", "Temperature")
    //! @param rCoefficient ... weight factor of this term
    //! @param rRHS ... prescribed right hand side value
    void ConstraintLinearEquationCreate(int rConstraint, int rNode, const std::string& rDof, double rCoefficient, double rRHS = 0);

#ifndef SWIG
    //! @brief ... create a constraint equation
    //! @param rConstraint ... constraint id
    //! @param rNode ... node id in the first constraint equation term
    //! @param rDofType ... type of dof in the first constraint equation term (e.g DISPLACEMENTS, ROTATIONS, TEMPERATURES)
    //! @param rDofComponent ... dof component (0, 1, 2)
    //! @param rCoefficient ... weight factor of this term
    //! @param rRHS ... prescribed right hand side value
    void ConstraintLinearEquationCreate(int rConstraint, int rNode, NuTo::Node::eDof rDofType, int rDofComponent, double rCoefficient, double rRHS = 0);

    //! @brief Creates a constraint equation that couples the degrees of freedom of an arbitrary node to a point in an element
    //! @param rNode ... node id in the first constraint equation term
    //! @param rElement ... element group id
    //! @param rDofType ... type of dof in the first constraint equation term (e.g DISPLACEMENTS, ROTATIONS, TEMPERATURES)
    //! @param numNearestNeighbours ... number of nearest neighbours to be found by the approximate nearest neighbour algorithm
    void ConstraintLinearEquationNodeToElementCreate(int rNode, int rElementGroup, NuTo::Node::eDof rDofType, int rNumNearestNeighbours, double rTolerance = 1.0e-6);


#endif

    //! @brief ... add a term to a constraint equation
    //! @param rConstraint ... constraint id
    //! @param rNode ... node id
    //! @param rDof ... dof (e.g "X_DISPLACEMENT", "Z_Rotation", "Temperature")
    //! @param rCoefficient ... weight factor of this term
    void ConstraintLinearEquationAddTerm(int rConstraint, int rNode, const std::string& rDof, double rCoefficient);

    //! @brief ... set periodic boundary conditions according to a prescibed angle of a localization zone
    //! @param  rAngle... angle in deg
    //! @param  rStrain... average strain to be applied (epsilon_xx, epsilon_yy, gamma_xy)
    //! @param  rNodeGroupUpper... all nodes on the upper boundary
    //! @param  rNodeGrouplower... all nodes on the lower boundary
    //! @param  rNodeGroupLeft... all nodes on the left boundary
    //! @param  rNodeGroupRight...  all nodes on the right boundary
    int ConstraintLinearDisplacementsSetPeriodic2D(double angle, NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> rStrain,
            NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> rCrackOpening, double rRadiusToCrackWithoutConstraints,
            int rNodeGroupUpper, int rNodeGrouplower, int rNodeGroupLeft, int rNodeGroupRight);

#ifndef SWIG
    //! @brief ... add a term to a constraint equation
    //! @param rConstraint ... constraint id
    //! @param rNode ... node id
    //! @param rDofType ... type of dof (e.g DISPLACEMENTS, ROTATIONS, TEMPERATURES)
    //! @param rDofComponent ... dof component (0, 1, 2)
    //! @param rCoefficient ... weight factor of this term
    void ConstraintLinearEquationAddTerm(int rConstraint, int rNode, NuTo::Node::eDof rDofType, int rDofComponent, double rCoefficient);

    //! @brief info about the nodes in the Structure
    void ConstraintInfo(int mVerboseLevel)const;

#endif

private:
    //! @brief ... convert input string to dof type and dof component
    //! @brief rDof ... input string
    //! @param rDofType ... type of dof (e.g DISPLACEMENTS, ROTATIONS, TEMPERATURES)
    //! @param rDofComponent ... dof component (0, 1, 2)
    void ConstraintEquationGetDofInformationFromString(const std::string& rDof, NuTo::Node::eDof& rDofType, int& rDofComponent);

public:
    //*************************************************
    //************ Load routines        ***************
    //***  defined in structures/StructureBaseLoad.cpp  ***
    //*************************************************

    //! @brief adds a force for a node
    //! @param rNodeIdent ... identifier for node
    //! @param rDirection ... direction of the force
    //! @param rValue ... force
    //! @return integer id to delete or modify the load
    int LoadCreateNodeForce(int rLoadCase, int rNodeIdent, const NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rDirection, double rValue);

    //! @brief adds a force for a node group
    //! @param rGroupIdent ... identifier for node group
    //! @param rDirection ... direction of the force
    //! @param rValue ... force
    //! @return integer id to delete or modify the load
    int LoadCreateNodeGroupForce(int rLoadCase, int rGroupIdent, const NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rDirection, double rValue);

    //! @brief Adds a heat flux to a node.
    //! @param rNodeIdent Identifier for node
    //! @param rDirection Direction of the flux
    //! @param rValue Value of the flux
    //! @return Integer id to delete or modify the load
    int LoadCreateNodeHeatFlux(int rLoadCase, int rNodeIdent, const NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rDirection, double rValue);

    //! @brief adds a surface load to 2D plane elements (2D)
    //! @param rElementGroupId ... specifies the elements with surface loads
    //! @param rNodeGroupId ... specifies the surfaces (if all nodes of an elemental surface is included in this group, the
    // surface is considered to be loaded
    //! @param rLoadVector ... constant load vector (independent of position and orientation of the loading surface
    //! @return integer id to delete or modify the load
    int LoadSurfaceConstDirectionCreate2D(int rLoadCase, int rElementGroupId, int rNodeGroupId,
    		const NuTo::FullVector<double,Eigen::Dynamic>& rLoadVector);

    //! @brief adds a surface load to 3D solid elements
    //! @param rElementGroupId ... specifies the elements with surface loads
    //! @param rNodeGroupId ... specifies the surfaces (if all nodes of an elemental surface is included in this group, the
    // surface is considered to be loaded
    //! @param rLoadVector ... constant load vector (independent of position and orientation of the loading surface
    //! @return integer id to delete or modify the load
    int LoadSurfaceConstDirectionCreate3D(int rLoadCase, int rElementGroupId, int rNodeGroupId,
    		const NuTo::FullVector<double,Eigen::Dynamic>& rLoadVector);

    //! @brief adds a surface load (pressure) to 3D solid elements
    //! @param rElementGroupId ... specifies the elements with surface loads
    //! @param rNodeGroupId ... specifies the surfaces (if all nodes of an elemental surface is included in this group, the
    // surface is considered to be loaded
    //! @param rPressure value ... normal to the surface, positive for compression
    //! @return integer id to delete or modify the load
    int LoadSurfacePressureCreate2D(int rLoadCase, int rElementGroupId, int rNodeGroupId, double rPressure);

    //! @brief adds a surface load (pressure) to 3D solid elements
    //! @param rElementGroupId ... specifies the elements with surface loads
    //! @param rNodeGroupId ... specifies the surfaces (if all nodes of an elemental surface is included in this group, the
    // surface is considered to be loaded
    //! @param rPressure value ... normal to the surface, positive for compression
    //! @return integer id to delete or modify the load
    int LoadSurfacePressureCreate3D(int rLoadCase, int rElementGroupId, int rNodeGroupId, double rPressure);

    //! @brief delete load
    //! @param rIdent ... load identifier
    void LoadDelete(int rIdent);

    //! @brief returns the number of load cases
    void SetNumLoadCases(int rNumLoadCases)
    {
    	mNumLoadCases = rNumLoadCases;
    }

    //! @brief returns the number of load cases
    int GetNumLoadCases()const
    {
    	return mNumLoadCases;
    }

#ifndef SWIG
    //! @brief adds a force for a node
    //! @param rNode ... pointer to node
    //! @param rDirection ... direction of the force
    //! @param rValue ... force
    //! @return integer id to delete or modify the load
    int LoadCreateNodeForce(int rLoadCase, const NodeBase* rNode, const NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rDirection, double rValue);

    //! @brief Adds a heat flux to a node.
    //! @param rNode Pointer to node
    //! @param rDirection Direction of the heat flux
    //! @param rValue Value of the flux
    //! @return Integer id to delete or modify the load
    int LoadCreateNodeHeatFlux(int rLoadCase, const NodeBase* rNode, const NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rDirection, double rValue);

    //! @brief adds a force for a node grpup
    //! @param rNodeGroup ... pointer to node group
    //! @param rDirection ... direction of the force
    //! @param rValue ... force
    //! @return integer id to delete or modify the load
    int LoadCreateNodeGroupForce(int rLoadCase, const Group<NodeBase>* rNodeGroup, const NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rDirection, double rValue);

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
    NuTo::FullVector<double, Eigen::Dynamic> ConstitutiveLawGetParameterFullVectorDouble(int rIdent, const std::string& rIdentifier) const;

    //! @brief ... sets a variable of the constitutive law which is selected by a string
    //! @param rIdent ... constitutive law identifier
    //! @param rIdentifier ... String to identify the requested variable
    //! @param rValue ... new value for requested variable
    void ConstitutiveLawSetParameterFullVectorDouble(int rIdent, const std::string& rIdentifier, NuTo::FullVector<double, Eigen::Dynamic>  rValue);


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
    NuTo::FullVector<double, Eigen::Dynamic> ConstitutiveLawGetParameterFullVectorDouble(int rIdent, Constitutive::eConstitutiveParameter rIdentifier) const;

    //! @brief ... sets a parameter of the constitutive law which is selected by an enum
    //! @param rIdent ... constitutive law identifier
    //! @param rIdentifier ... Enum to identify the requested parameter
    //! @param rValue ... new value for requested variable
    void ConstitutiveLawSetParameterFullVectorDouble(int rIdent, Constitutive::eConstitutiveParameter rIdentifier, NuTo::FullVector<double, Eigen::Dynamic>  rValue);
#endif // SWIG


    //! @brief ... get yield strength for multilinear response
    //! @return ... first column: equivalent plastic strain
    //! @return ... second column: corresponding yield strength
    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> ConstitutiveLawGetYieldStrength(int rIdent) const;

    //! @brief ... add yield strength
    //! @param rEpsilon ...  equivalent plastic strain
    //! @param rSigma ...  yield strength
    void ConstitutiveLawAddYieldStrength(int rIdent, double rEpsilon, double rSigma);

    //! @brief ... get hardening modulus for multilinear response
    //! @return ... first column: equivalent plastic strain
    //! @return ... second column: corresponding hardening modulus
    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> ConstitutiveLawGetHardeningModulus(int rIdent) const;

    //! @brief ... add hardening modulus
    //! @param rEpsilon ...  equivalent plastic strain
    //! @param rSigma ...  hardening modulus
    void ConstitutiveLawAddHardeningModulus(int rIdent, double rEpsilon, double rH);


#ifndef SWIG
    //! @brief ... set damage law
    //! @param rDamageLaw ... damage law
    void ConstitutiveLawSetDamageLaw(int rIdent, Constitutive::eDamageLawType rDamageLaw);
#endif

    //! @brief ... set damage law
    //! @param rDamageLaw ... damage law
    void ConstitutiveLawSetDamageLaw(int rIdent, std::string rDamageLaw);




    //! @brief ... gets the equilibrium water volume fraction depend on the relative humidity
    //! @param rIdent ... constitutive law identifier
    //! @param rRelativeHumidity ... relative humidity
    //! @param rCoeffs ... polynomial coefficients of the sorption curve
    //! @return ... equilibrium water volume fraction
    double ConstitutiveLawGetEquilibriumWaterVolumeFraction(int rIdent, double rRelativeHumidity, NuTo::FullVector<double,Eigen::Dynamic> rCoeffs) const;

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
    //************ Section routines     ***************
    //**  defined in structures/StructureBaseSection.cpp **
    //*************************************************
    //! @brief ... returns the number of sections
    //! @return number of sections
    inline int GetNumSections() const
    {
        return this->mSectionMap.size();
    }

    //! @brief ... create a new section
    //! @param rIdent ... section identifier
    //! @param rType ... section type
    int SectionCreate(const std::string& rType);

    //! @brief ... delete an existing section
    //! @param rIdent ... section identifier
    void SectionDelete(int rId);

    //! @brief ... set section cross-section area
    //! @param rIdent ... section identifier
    //! @param rArea ... cross-section area
    void SectionSetArea(int rId, double rArea);

    //! @brief ... get section cross-section area
    //! @param rIdent ... section identifier
    //! @return section cross-section area
    double SectionGetArea(int rId) const;

    //! @brief ... set section thickness
    //! @param rIdent ... section identifier
    //! @param rThickness ... thickness
    void SectionSetThickness(int rId, double rThickness);

    //! @brief ... get section thickness
    //! @param rIdent ... section identifier
    //! @return section thickness
    double SectionGetThickness(int rId) const;

    //! @brief ... set section circumference for interface bewteen matrix and fibre
    //! @param rIdent ... section identifier
    //! @param rThickness ... circumference of the fibre
    void SectionSetCircumference(int rId, double rCircumference);

    //! @brief ... get section thickness
    //! @param rIdent ... section identifier
    //! @return section circumference of the fibre
    double SectionGetCircumference(int rId) const;

    //! @brief ... print information about all sections
    //! @param rVerboseLevel ... controls the verbosity of the information
    void SectionInfo(unsigned short rVerboseLevel) const;

    //! @brief ... print information of a single section
    //! @param rIdent ... section identifier
    //! @param rVerboseLevel ... controls the verbosity of the information
    void SectionInfo(int rId, unsigned short rVerboseLevel) const;

#ifndef SWIG
    //! @brief ... create a new section
    //! @param rIdent ... section identifier
    //! @param rType ... section type
    int SectionCreate(eSectionType rType);

    //! @brief ... get the pointer to a section from the section identifier
    //! @param rIdent ... section identifier
    //! @return ... pointer to the section
    SectionBase* SectionGetSectionPtr(int rId);

    //! @brief ... get the pointer to a section from the section identifier
    //! @param rIdent ... section identifier
    //! @return ... pointer to the section
    const SectionBase* SectionGetSectionPtr(int rId) const;

    //! @brief ... get the identifier of an section from the pointer to section object
    //! @param rSectionPtr ... pointer to section object
    //! @return ... section identifier
    int SectionGetId(const SectionBase* rSectionPtr) const;
#endif

    //*************************************************
    //************ Group routines     ***************
    //**  defined in structures/StructureBaseGroup.cpp **
    //*************************************************
    //! @brief ... Info routine that prints general information about the groups
    //! @param ... rVerboseLevel describes how detailed the information is
    void GroupInfo(int rVerboseLevel)const;

#ifndef SWIG
    //! @brief gives the identifier of an element
    //! @param reference to an element
    //! @return identifier
    int GroupGetId(GroupBase* rGroup)const;

    //! @brief ... get the pointer to a group from the group identifier
    //! @param rIdent ... group identifier
    //! @return ... pointer to the group
    GroupBase* GroupGetGroupPtr(int rIdent);

    //! @brief ... get the pointer to a group from the group identifier
    //! @param rIdent ... group identifier
    //! @return ... pointer to the group
    const GroupBase* GroupGetGroupPtr(int rIdent) const;

    //! @brief ... Creates a group for the structure
    //! @param ... rType  type of the group, e.g. "NODES" or "ELEMENTS"
    //! @return ... rIdent identifier for the group
    int GroupCreate(NuTo::eGroupId rEnumType);

    //! @brief ... Creates a group for the structure
    //! @param ... rIdent identifier for the group
    //! @param ... rType  type of the group
    void GroupCreate(int id, NuTo::eGroupId rEnumType);
#endif

    //! @brief ... Creates a group for the structure
    //! @param ... rIdent identifier for the group
    //! @param ... rType  type of the group, e.g. "NODES" or "ELEMENTS"
    int GroupCreate(const std::string& rType);


    //! @brief ... Deletes a group from the structure
    //! @param ... rIdent identifier for the group
    void GroupDelete(int rIdent);

    //! @brief ... Unites two groups and stores the result in a new group
    //! @param ... rIdentGroup1 identifier for the first group
    //! @param ... rIdentGroup2 identifier for the second group
    //! @result ... rIdentGroupResult identifier for the created result group
    int GroupUnion(int rIdentGroup1, int rIdentGroup2);
    //! @brief ... Difference between two groups and stores the result in a new group
    //! @param ... rIdentGroup1 identifier for the first group
    //! @param ... rIdentGroup2 identifier for the second group
    //! @result ... rIdentGroupResult identifier for the created result group
    int GroupDifference(int rIdentGroup1, int rIdentGroup2);

    //! @brief ... Calculates the intersection between two groups and stores the result in a new group
    //! @param ... rIdentGroup1 identifier for the first group
    //! @param ... rIdentGroup2 identifier for the second group
    //! @result ... rIdentGroupResult identifier for the created result group
    int GroupIntersection(int rIdentGroup1, int rIdentGroup2);

    //! @brief ... Calculates the symmetric difference between two groups and stores the result in a new group
    //! @param ... rIdentGroup1 identifier for the first group
    //! @param ... rIdentGroup2 identifier for the second group
    //! @result ... rIdentGroupResult identifier for the created result group
    int GroupSymmetricDifference(int rIdentGroup1, int rIdentGroup2);

    //! @brief ... Adds a node to a node group
    //! @param ... rIdentGroup identifier for the group
    //! @param ... rIdentNode  identifier for the node
    void GroupAddNode(int rIdentGroup, int rIdNode);

    //! @brief ... Adds all nodes to a group whose coordinates are in the specified range
    //! @param ... rIdentGroup identifier for the group
    //! @param ... rDirection either 0,1,2 for x,y, or z
    //! @param ... rMin ... minimum value
    //! @param ... rMax ... maximum value
    virtual void GroupAddNodeCoordinateRange(int rIdentGroup, int rDirection, double rMin, double rMax);

#ifndef SWIG
    //! @brief ... Adds all nodes which fulfill the conditions specified in a std::function
    //! @param ... rIdentGroup identifier for the group
    //! @param ... rFunction std::function
    void GroupAddNodeFunction(int rIdentGroup,
                              std::function<bool (NodeBase*)> rFunction);
#endif

    //! @brief ... Adds an element to an element group
    //! @param ... rIdentGroup identifier for the group
    //! @param ... rIdentNode  identifier for the element
    void GroupAddElement(int rIdentGroup, int rIdElement);

    //! @brief ... Adds all elements to a group whose nodes are in the given node group
    //! @param ... rElementGroupId identifier for the element group
    //! @param ... rNodeGroupId identifier for the node group
    //! @param ... rHaveAllNodes if set to true, the element is only selected when all element nodes are in the node group, if set
    //! to false, the element is select if at least one node is in the node group
    void GroupAddElementsFromNodes(int rElementGroupId, int rNodeGroupId, bool rHaveAllNodes);

    //! @brief ... Adds all the nodes from the group-rElementGroupId to the group rNodeGroupId
    //! @param ... rNodeGroupId id for the node group
    //! @param ... rElementGroupId id for the element group
    void GroupAddNodesFromElements(int rNodeGroupId, int rElementGroupId);

    //! @brief ... Adds all nodes to a group whose coordinates are in the specified range
    //! @param ... rIdentGroup identifier for the group
    //! @param ... rCenter center of the selection circle
    //! @param ... rMin ... minimum radius
    //! @param ... rMax ... maximum radius
    void GroupAddNodeRadiusRange(int rIdentGroup, NuTo::FullVector<double,Eigen::Dynamic> rCenter, double rMin, double rMax);

    //! @brief ... Adds all nodes to a group whose coordinates are on a cylinder with the radius in the in the specified range
    //! @param ... rIdentGroup identifier for the group
    //! @param ... rCenter center of the cylinder
    //! @param ... rAxis axis of the cylinder
    //! @param ... rMin ... minimum radius
    //! @param ... rMax ... maximum radius
    void GroupAddNodeCylinderRadiusRange(int rIdentGroup, NuTo::FullVector<double,Eigen::Dynamic> rCenter,
    		 NuTo::FullVector<double,Eigen::Dynamic> rDirection, double rMin, double rMax);

    //! @brief ... Returns the number of members in a group
    //! @param ... rIdentGroup identifier for the group
    //! @return ... number of members
    int GroupGetNumMembers(int rIdentGroup)const;

    //! @brief ... Returns a vector with the members of a group
    //! @param ... rIdentGroup identifier for the group
    //! @return ... vector of members
    NuTo::FullVector<int, Eigen::Dynamic> GroupGetMemberIds(int rIdentGroup)const;

    //! @brief ... checks for a member in a group
    //! @param ... rIdentGroup identifier for the group
    //! @return ... rMember id (element id, node id etc.)
    bool GroupContainsMember(int rIdentGroup, int rMember)const;

    //*************************************************************
    //************ Integration type routines     ******************
    //**  defined in structures/StructureBaseIntegrationType.cpp **
    //*************************************************************
    //! @brief ... Info routine that prints general information about the allocated integration types
    //! an integration type is only allocated if required (from created elements)
    //! @param ... rVerboseLevel describes how detailed the information is
    void IntegrationTypeInfo(int rVerboseLevel)const;

#ifndef SWIG
    //! @brief ... Returns a pointer to an integration type
    //! if the integration type does not exist in the map, an exception is thrown
    //! @param identIntegrationType Identifier for an integration type
    NuTo::IntegrationTypeBase* GetPtrIntegrationType(const std::string& rIdentIntegrationType);

    //! @brief ... Returns a pointer to an integration type
    //! if the integration type does not exist (in the map), the integration type is created
    //! @param identIntegrationType Identifier for an integration type
    NuTo::IntegrationTypeBase* GetPtrIntegrationType(NuTo::eIntegrationType rIdentIntegrationType);
#endif //SWIG

    //*************************************************
    //************    Crack routines    ***************
    //**  pure virtual                               **
    //*************************************************
	//! @brief returns the number of cracks
	//! @return number of cracks
	virtual unsigned int GetNumCracks()const=0;

#ifndef SWIG
	//! @brief a reference to a crack
	//! @param identifier
	//! @return reference to a cracks
	virtual CrackBase* CrackGetCrackPtr(int rIdent)=0;

	//! @brief a reference to a crack
	//! @param identifier
	//! @return reference to a crack
	virtual const CrackBase* CrackGetCrackPtr(int rIdent)const=0;

	//! @brief gives the identifier of a crack
	//! @param reference to a crack
	//! @return identifier
	virtual int CrackGetId(const CrackBase* rCrack)const=0;
#endif //SWIG

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
	inline Logger& GetLogger()const
	{
		return mLogger;
	}
#endif //SWIG


    //*************************************************
    //************ Basic routines     ***************
    //**  defined in structures/StructureBase.cpp **
    //*************************************************
    //! @brief ... number of time derivatives (0 : static, 1: velocities, 2: accelerations)
	void SetNumTimeDerivatives(int rNumTimeDerivatives);

	//! @brief ... return number of time derivatives (0 : static, 1: velocities, 2: accelerations)
	int GetNumTimeDerivatives()const;

	//! @brief ... Info routine that prints general information about the object (detail according to verbose level)
    virtual void Info()const;

    //! @brief set the beginning of the time increment to the structure
    void SetPrevTime(double rPrevTime);

    //! @brief get the beginning of the time increment of the structure
    double GetPrevTime() const;

    //! @brief set the end of the time increment to the structure (current time)
    void SetTime(double rTime);

    //! @brief get the end of the time increment of the structure (current time)
    double GetTime() const;

    //! @brief set number of cycles to be extrapolated in the cycle jump routine
    //! @brief ... rNumber[0] is the number of extrapolated cycles itself Njump
    //! @brief ... rNumber[1] is the weighting coefficient of the implicit term
    //! @brief ... rNumber[2] is the weighting coefficient of the explicit term
    //! @brief ... rNumber[3] and higher are the weighting coefficients of the terms for a higher-order extrapolation
    void SetNumExtrapolatedCycles(NuTo::FullVector<double,Eigen::Dynamic> rNumber);

    //! @brief get the number of cycles to be extrapolated in the cycle jump routine. Returns:
    //! @brief ... [0] is the number of extrapolated cycles itself Njump
    //! @brief ... [1] is the weighting coefficient of the implicit term
    //! @brief ... [2] is the weighting coefficient of the explicit term
    //! @brief ... [3] and higher are the weighting coefficients of the terms for a higher-order extrapolation
    NuTo::FullVector<double,Eigen::Dynamic> GetNumExtrapolatedCycles() const;

    //! @brief absolute tolerance for entries of the global stiffness matrix (coefficientMatrix0)
    //! values smaller than that one will not be added to the global matrix
    void SetToleranceStiffnessEntries(double rToleranceStiffnessEntries);

    //! @brief absolute tolerance for entries of the global stiffness matrix (coefficientMatrix0)
    //! values smaller than that one will not be added to the global matrix
    double GetToleranceStiffnessEntries()const;

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

    //! @brief returns the a reference to the constraint matrix
    const NuTo::BlockSparseMatrix& GetConstraintMatrix() const;

    //! @brief this routine is only relevant for the multiscale model, since an update on the fine scale should only be performed
    //for an update on the coarse scale
    //as a consequence, in an iterative solution with updates in between the initial state has to be restored after leaving the routine
    //this routine saves the current state before an update in the Newton Raphson iteration is performed
    //this only happens for more than one load step (either prescibed or with automatic load control)
    virtual void SaveStructure(std::stringstream& rSaveStringStream)const
    {
    	throw MechanicsException("[StructureBase::SaveStructure] Saving of the structure not implemented in derived class.");
    }

    virtual void RestoreStructure(std::stringstream& rSaveStringStream)
    {
    	throw MechanicsException("[StructureBase::RestoreStructure] Saving of the structure not implemented in derived class.");
    }

    void SetUpdateTmpStaticDataRequired()
    {
    	mUpdateTmpStaticDataRequired = true;
    }

    //! @brief calculate the critical time step for a vector of elements solving the generalized eigenvalue problem Ku=lambda Mu
    double ElementCalculateLargestElementEigenvalue(const std::vector< ElementBase*>& rElementVector);

	//! @brief returns whether or not the dof is constitutive input at least in one InrepolationType
	//! @param rInterpolationTypeId ... interpolation type id
	//! @param rDofType ... dof type
	bool InterpolationTypeIsConstitutiveInput(NuTo::Node::eDof rDofType);

protected:
    //! @brief ... number of time derivatives (0 : static, 1: velocities, 2: accelerations)
	int mNumTimeDerivatives;

	//! @brief ... storing the beginning of the time increment
	double mPrevTime;

    //! @brief ... storing the end of the time increment (current time)
	double mTime;

    int mDimension;

    //! @brief ... number of cycles applied for extrapolation in the cycle jump.
    //! @brief ... extrapolation is Njump * (ImplicitTerm*A1 + ExplicitTerm*A2 + HigherOrder*A3 + ... )
    //! @brief ... [0] is the number of extrapolated cycles itself Njump
    //! @brief ... [1] is the weighting coefficient of the implicit term A1
    //! @brief ... [2] is the weighting coefficient of the explicit term A2
    //! @brief ... [3] and higher are the weighting coefficients of the terms for a higher-order extrapolation
    NuTo::FullVector<double,Eigen::Dynamic> mNumExtrapolatedCycles;

    //! @brief ... map storing the name and the pointer to the constitutive law
    //! @sa ConstitutiveBase
    boost::ptr_map<int,ConstitutiveBase> mConstitutiveLawMap;

    //! @brief ... map storing the constraints
    //! @sa ConstraintBase
    boost::ptr_map<int,ConstraintBase> mConstraintMap;

    //! @brief ... map storing node loads
    //! @sa LoadBase
    int mNumLoadCases;        //number of load cases to be considered
    boost::ptr_map<int,LoadBase> mLoadMap;

    //! @brief ... map storing the groups and a pointer to the objects
    //! @sa GroupBase
    boost::ptr_map<int,GroupBase> mGroupMap;

    //! @brief ... map storing the name and the pointer to the integration types
    //! @sa IntegrationTypeBase
    boost::ptr_map<std::string,IntegrationTypeBase> mIntegrationTypeMap;

    //! @brief ... map storing the section name and the pointer to the section object
    //! @sa SectionBase
    boost::ptr_map<int,SectionBase> mSectionMap;

    //! @brief ... map storing the interpolation types
    //! @sa InterpolationType
    boost::ptr_map<int,InterpolationType> mInterpolationTypeMap;

    //! @brief ... a mapping from the enums of the predefined integration types to their corresponding string name
    std::vector<std::string> mMappingIntEnum2String;

    //! @brief ... map storing the components (displacements, strains, nonlocal weights etc) to be included in the output (VTK) file
    std::map<int, std::list<std::shared_ptr<VisualizeComponent>>> mGroupVisualizeComponentsMap;

    //! @brief ... map storing the type of visualization for the output (VTK) file
    std::map<int, eVisualizationType> mGroupVisualizationType;

    //! @brief summarizes information to dof numbering, active dof types, symmetric dof types, constant dof types
    DofStatus mDofStatus;

    //!brief ... renumbering of nodal DOFs required or not
    bool mNodeNumberingRequired;

    //! @brief constraint matrix relating the prescibed nodal unknowns to the free parameters
    BlockSparseMatrix mConstraintMatrix;

    //! @brief mapping matrix of the rhs to relate the rhs before the gauss elimination to the constraint matrix after
    // (mConstraintRHS (after elimination) = mConstraintMappingRHS *  mConstraintRHS (before elimination)
    // (the values of the RHS before elimination are stored at the individual constraints
    //the initial system is e.g.
    //[1 1 0]* [d1 d2 d3]^T = [rhs1]
    //[0 0 2]                 [rhs2]
    //this is replaced by
    //[1 1 0]* [d1 d2 d3]^T = rhs1 *[1] + rhs2 *[0]
    //[0 0 2]                       [0]         [1]
    //after gauss elimination and reordering this gives
    //[1 0 1]* [d1 d3 d2]^T = rhs1 *[1] + rhs2 *[0]
    //[0 1 0]                       [0]         [0.5]
    //as a consequence, the gauss elimination has only to be performed for a change of the constraint matrix
    //for a change of the rhs it is sufficient to recalculate the rhs from the above vectors
    //the mapping matrix [1,0; 0,0.5] is stored and the rhs is calculated from mConstraintMappingRHS*mConstraintRHSBeforGaussElimination
    BlockSparseMatrix mConstraintMappingRHS;

    //! @brief right hand side of the constraint equations
    BlockFullVector<double> mConstraintRHS;

    //! @brief is set to true, if at least one constitutive model requires an update of tmpStaticData before stress and stiffness routines are called
    bool mHaveTmpStaticData;

    //! @brief is set to false, if the structure is changed (nodes, elements) or (DOFs at the nodes)
    bool mUpdateTmpStaticDataRequired;

    //! @brief absolute tolerance for entries of the global stiffness matrix (coefficientMatrix0)
    //! values smaller than that one will not be added to the global matrix
    double mToleranceStiffnessEntries;

#ifdef _OPENMP
    //@brief maximum independent sets used for parallel assembly of the stiffness resforce etc.
    mutable std::vector<std::vector<ElementBase*> > mMIS;
    //@brief number of processors used in an openmp simulation
    int mNumProcessors;
#endif

    //! @brief logger class to redirect the output to some file or the console (or both), can be changed even for const routines
    mutable Logger mLogger;


#ifdef ENABLE_SERIALIZATION
    //! @brief ... standard constructor just for the serialization routine
    StructureBase() :
            mConstraintMatrix(mDofStatus, false),
            mConstraintMappingRHS(mDofStatus, false),
            mConstraintRHS(mDofStatus)
    {}
#endif  // ENABLE_SERIALIZATION

#ifndef SWIG
    //! @brief ... store all elements of a structure in a vector
    //! @param rElements ... vector of element pointer
    virtual void GetElementsTotal(std::vector<const ElementBase*>& rElements) const = 0;

    //! @brief ... store all elements of a structure in a vector
    //! @param rElements ... vector of element pointer
    virtual void GetElementsTotal(std::vector<std::pair<int,const ElementBase*> >& rElements) const = 0;

    //! @brief ... store all elements of a structure in a vector
    //! @param rElements ... vector of element pointer
    virtual void GetElementsTotal(std::vector<ElementBase*>& rElements) = 0;

    //! @brief ... store all elements of a structure in a vector
    //! @param rElements ... vector of element pointer
    virtual void GetElementsTotal(std::vector<std::pair<int,ElementBase*> >&  rElements) = 0;

    //! @brief ... store all nodes of a structure in a vector
    //! @param rNodes ... vector of element pointer
    virtual void GetNodesTotal(std::vector<const NodeBase*>& rNodess) const = 0;

    //! @brief ... store all nodes of a structure in a vector
    //! @param rNodes ... vector of element pointer
    virtual void GetNodesTotal(std::vector<std::pair<int,const NodeBase*> >& rNodes) const = 0;

    //! @brief ... store all nodes of a structure in a vector
    //! @param rNodes ... vector of element pointer
    virtual void GetNodesTotal(std::vector<NodeBase*>& rNodes) = 0;

    //! @brief ... store all nodes of a structure in a vector
    //! @param rNodes ... vector of element pointer
    virtual void GetNodesTotal(std::vector<std::pair<int,NodeBase*> >& rNodes) = 0;
#endif

    //! @brief ... get all elements of a group in a vector
    //! @param rElementGroup ... element group
    //! @param rElements ... vector of element pointer
    void GetElementsByGroup(const Group<ElementBase>* rElementGroup, std::vector<const ElementBase*>& rElements) const;

    //! @brief ... get all elements of a group in a vector
    //! @param rElementGroup ... element group
    //! @param rElements ... vector of element pointer
    void GetElementsByGroup(Group<ElementBase>* rElementGroup, std::vector< ElementBase*>& rElements);

};
} //namespace NuTo
#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
#include <boost/serialization/assume_abstract.hpp>
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::StructureBase)
#endif // SWIG
#endif  // ENABLE_SERIALIZATION
