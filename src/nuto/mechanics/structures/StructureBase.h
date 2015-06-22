// $Id$

#ifndef STRUCTUREBASE_H
#define STRUCTUREBASE_H

#include <ctime>
#include <array>

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

#include <boost/ptr_container/ptr_map.hpp>
#include <boost/ptr_container/ptr_list.hpp>

#include <string>
#include "nuto/base/Logger.h"
#include "nuto/base/NuToObject.h"
#include "nuto/math/FullMatrix_Def.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"
#include "nuto/mechanics/constitutive/ConstitutiveBase.h"
#include "nuto/mechanics/constraints/ConstraintBase.h"
#include "nuto/mechanics/groups/GroupBase.h"
#include "nuto/mechanics/elements/IpDataEnum.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeBase.h"
#include "nuto/mechanics/interpolationtypes/InterpolationType.h"
#include "nuto/mechanics/loads/LoadBase.h"
#include "nuto/mechanics/nodes/NodeEnum.h"
#include "nuto/mechanics/sections/SectionBase.h"
#include "nuto/mechanics/sections/SectionEnum.h"
#include "nuto/mechanics/structures/StructureBaseEnum.h"

#ifdef ENABLE_VISUALIZE
#include "nuto/visualize/VisualizeBase.h"
#include "nuto/visualize/VisualizeComponentBase.h"
#endif // ENABLE_VISUALIZE

namespace NuTo
{
class ElementBase;
class NodeBase;
template<class T, int rows> class FullVector;
template<class T, int rows, int cols> class FullMatrix;
template<class T> class SparseMatrixCSRSymmetric;
template<class T> class SparseMatrixCSRGeneral;
template<class T> class SparseMatrixCSRVector2General;
template<class T> class SparseMatrixCSRVector2Symmetric;
class EngineeringStrain2D;
class NewtonRaphsonAuxRoutinesBase;
class CrackBase;
class ConstitutiveStaticDataMultiscale2DPlaneStrain;

#ifdef ENABLE_VISUALIZE
class VisualizeUnstructuredGrid;
class VisualizeComponentBase;
#endif // ENABLE_VISUALIZE

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
    virtual ~StructureBase()
    {}

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif  // ENABLE_SERIALIZATION

    //! @brief gives the dimension of the Structure
    //! @return Structural dimension (1,2 or 3)
    int GetDimension()const;

#ifdef ENABLE_VISUALIZE
    //! @brief ... Add the damage variable to the internal list, which is finally exported via the ExportVtkDataFile command
    void AddVisualizationComponentDamage();

    //! @brief ... Add visualization displacements to the internal list, which is finally exported via the ExportVtkDataFile command
    void AddVisualizationComponentDisplacements();

    //! @brief ... Add element id to the internal list, which is finally exported via the ExportVtkDataFile command
    void AddVisualizationComponentElement();

    //! @brief ... Add engineering strains to the internal list, which is finally exported via the ExportVtkDataFile command
    void AddVisualizationComponentEngineeringStrain();

    //! @brief ... Add engineering plastic strains to the internal list, which is finally exported via the ExportVtkDataFile command
    void AddVisualizationComponentEngineeringPlasticStrain();

    //! @brief ... Add engineering stress to the internal list, which is finally exported via the ExportVtkDataFile command
    void AddVisualizationComponentEngineeringStress();

    //! @brief ... Add section id to the internal list, which is finally exported via the ExportVtkDataFile command
    void AddVisualizationComponentSection();

    //! @brief ... Add constitutive id to the internal list, which is finally exported via the ExportVtkDataFile command
    void AddVisualizationComponentConstitutive();

    //! @brief ... Add crack id vector to the internal list, which is finally exported via the ExportVtkDataFile command
    void AddVisualizationComponentCracks();

    //! @brief ... Add visualization particle radius to the internal list, which is finally exported via the ExportVtkDataFile command
    void AddVisualizationComponentParticleRadius();

    //! @brief ... Add visualization of principal stresses to the internal list, which is finally exported via the ExportVtkDataFile command
    void AddVisualizationComponentPrincipalEngineeringStress();

    //! @brief ... Add nonlocal weights to the internal list, which is finally exported via the ExportVtkDataFile command
    //! @param rElementId ... Element id
    //! @param rIp ... local ip number
    void AddVisualizationComponentNonlocalWeights(int rElementId, int rIp);

    //! @brief ... Add visualization of rotations to the internal list, which is finally exported via the ExportVtkDataFile command
    void AddVisualizationComponentRotation();

    //! @brief ... Add visualization of velocity to the internal list, which is finally exported via the ExportVtkDataFile command
    void AddVisualizationComponentVelocity();

    //! @brief ... Add visualization of acceration to the internal list, which is finally exported via the ExportVtkDataFile command
    void AddVisualizationComponentAcceleration();

    //! @brief ... Add visualization of angular velocity to the internal list, which is finally exported via the ExportVtkDataFile command
    void AddVisualizationComponentAngularVelocity();

    //! @brief ... Add visualization of angular velocity to the internal list, which is finally exported via the ExportVtkDataFile command
    void AddVisualizationComponentAngularAcceleration();

    //! @brief ... Add visualization of temperature to the internal list, which is finally exported via the ExportVtkDataFile command
    void AddVisualizationComponentTemperature();

    //! @brief ... Add visualization of heat flux to the internal list, which is finally exported via the ExportVtkDataFile command
    void AddVisualizationComponentHeatFlux();

    //! @brief ... Add visualization of nonlocal equivalent strain to the internal list, which is finally exported via the ExportVtkDataFile command
    void AddVisualizationComponentNonlocalEqStrain();

    //! @brief ... Add visualization of relative humidity to the internal list, which is finally exported via the ExportVtkDataFile command
    void AddVisualizationComponentRelativeHumidity();

    //! @brief ... Add visualization of water volume fraction to the internal list, which is finally exported via the ExportVtkDataFile command
    void AddVisualizationComponentWaterVolumeFraction();

    //! @brief ... clear all visualization components
    void ClearVisualizationComponents();

    //! @brief ... export the entire structure to Vtk data file
    //! @param rFileName ... file name
    //! @param rWhat ... string which describes what to plot
    void ExportVtkDataFile(const std::string& rFileName);

    //! @brief ... export the entire structure to Vtk data file
    //! @param rFileName ... file name
    //! @param rWhat ... string which describes what to plot
    void ExportVtkDataFileElements(const std::string& rFileName);

    //! @brief ... export the entire structure to Vtk data file
    //! @param rFileName ... file name
    //! @param rWhat ... string which describes what to plot
    void ExportVtkDataFileNodes(const std::string& rFileName);

    //! @brief ... export the entire structure to Vtk data file
    //! @param rFileName ... file name
    //! @param rWhat ... string which describes what to plot
    //! @param rXML ... if true, a vtu file is exported, otherwise a legacy vtk file is produced
    void ExportVtkDataFileElements(const std::string& rFileName, bool rXML);

    //! @brief ... export the entire structure to Vtk data file
    //! @param rFileName ... file name
    //! @param rWhat ... string which describes what to plot
    //! @param rXML ... if true, a vtu file is exported, otherwise a legacy vtk file is produced
    void ExportVtkDataFileNodes(const std::string& rFileName, bool rXML);

    //! @brief ... export an element group to Vtk/xml data file
    //! @param rGroupIdent ... group ident
    //! @param rFileName ... file name
    //! @param rXML ... use xml or vtk format (true for xml)
    void ElementGroupExportVtkDataFile(int rGroupIdent, const std::string& rFileName, bool rXML);

#ifndef SWIG
    //! @brief ... define the data sets (scalar, vector etc for the visualize routine based on the mVisualizecomponents for an element plot
    void DefineVisualizeElementData(VisualizeUnstructuredGrid& rVisualize, const boost::ptr_list<NuTo::VisualizeComponentBase>& rWhat)const;

    //! @brief ... define the data sets (scalar, vector etc for the visualize routine based on the mVisualizecomponents for a node plot
    void DefineVisualizeNodeData(VisualizeUnstructuredGrid& rVisualize, const boost::ptr_list<NuTo::VisualizeComponentBase>& rWhat)const;

    //! @brief ... adds all the elements in the vector to the data structure that is finally visualized
    void ElementVectorAddToVisualize(VisualizeUnstructuredGrid& rVisualize, const boost::ptr_list<NuTo::VisualizeComponentBase>& rWhat, const std::vector<ElementBase*>& rElements);

    //! @brief ... adds all the elements in the vector to the data structure that is finally visualized
    void ElementTotalAddToVisualize(VisualizeUnstructuredGrid& rVisualize, const boost::ptr_list<NuTo::VisualizeComponentBase>& rWhat);

    //! @brief ... adds all the elements in a group to the data structure that is finally visualized
    void ElementGroupAddToVisualize(int rGroupId, VisualizeUnstructuredGrid& rVisualize, const boost::ptr_list<NuTo::VisualizeComponentBase>& rWhat);

#endif //SWIG
#endif // ENABLE_VISUALIZE

    //@brief determines the maximum independent sets and stores it at the structure
    // is only relevant for openmp, otherwise the routine is just empty
    void CalculateMaximumIndependentSets();

    //@brief determines if in the omp parallelization the maximum independent sets are used (parallel assembly of the stiffness, generally faster)
    // or sequential insertion of the element stiffness using a barrier (faster for different load balancing of the elements)
    // is only relevant for openmp, otherwise the routine is just empty
    void UseMaximumIndependentSets(bool rUseMIS);

    //@brief set the number of processors for openmp parallelization
    void SetNumProcessors(int rNumProcessors);

    //@brief set the number of processors for openmp parallelization
    void SetOMPNested(bool rNested);

#ifndef SWIG
    //! @brief ... build global coefficient matrix (e.g stiffness) for primary dofs (e.g displacements, rotations, temperature)
    //! @param rType ... matrix type (stiffness or mass)
    //! @param rMatrix ... global coefficient matrix (nonsymmetric)
    //! @param rVector ... global equivalent load vector (e.g. due to prescribed displacements)
    NuTo::Error::eError BuildGlobalCoefficientMatrix(NuTo::StructureBaseEnum::eMatrixType rType, NuTo::SparseMatrixCSRGeneral<double>& rMatrix, NuTo::FullVector<double,Eigen::Dynamic>& rVector);

    //! @brief ... build global coefficient matrix (e.g stiffness) for primary dofs (e.g displacements, rotations, temperature)
    //! @param rMatrix ... global coefficient matrix (symmetric)
    //! @param rVector ... global equivalent load vector (e.g. due to prescribed displacements)
    NuTo::Error::eError BuildGlobalCoefficientMatrix(NuTo::StructureBaseEnum::eMatrixType rType, NuTo::SparseMatrixCSRSymmetric<double>& rMatrix, NuTo::FullVector<double,Eigen::Dynamic>& rVector);

    //! @brief ... build global coefficient matrix (e.g stiffness) for primary dofs (e.g displacements, rotations, temperature)
    //! @param rMatrix ... global coefficient matrix (nonsymmetric)
    //! @param rVector ... global equivalent load vector (e.g. due to prescribed displacements)
    NuTo::Error::eError BuildGlobalCoefficientMatrix(NuTo::StructureBaseEnum::eMatrixType rType, NuTo::SparseMatrixCSRVector2General<double>& rMatrix, NuTo::FullVector<double,Eigen::Dynamic>& rVector);

    //! @brief ... build global coefficient matrix (e.g stiffness) for primary dofs (e.g displacements, rotations, temperature)
    //! @param rMatrix ... global coefficient matrix
    //! @param rVector ... global equivalent load vector (e.g. due to prescribed displacements)
    NuTo::Error::eError BuildGlobalCoefficientMatrix(NuTo::StructureBaseEnum::eMatrixType rType, NuTo::SparseMatrixCSRVector2Symmetric<double>& rMatrix, FullVector<double,Eigen::Dynamic>& rVector);
#endif //SWIG

    //! @brief ... build global coefficient matrix (stiffness) for primary dofs (e.g displacements, rotations, temperature)
    //! @param rMatrix ... global coefficient matrix (nonsymmetric)
    //! @param rVector ... global equivalent load vector (e.g. due to prescribed displacements)
    NuTo::Error::eError BuildGlobalCoefficientMatrix0(NuTo::SparseMatrixCSRGeneral<double>& rMatrix, NuTo::FullVector<double,Eigen::Dynamic>& rVector);

    //! @brief ... build global coefficient matrix (stiffness) for primary dofs (e.g displacements, rotations, temperature)
    //! @param rMatrix ... global coefficient matrix (symmetric)
    //! @param rVector ... global equivalent load vector (e.g. due to prescribed displacements)
    NuTo::Error::eError BuildGlobalCoefficientMatrix0(NuTo::SparseMatrixCSRSymmetric<double>& rMatrix, NuTo::FullVector<double,Eigen::Dynamic>& rVector);

    //! @brief ... build global coefficient matrix (stiffness) for primary dofs (e.g displacements, rotations, temperature)
    //! @param rMatrix ... global coefficient matrix (nonsymmetric)
    //! @param rVector ... global equivalent load vector (e.g. due to prescribed displacements)
    NuTo::Error::eError BuildGlobalCoefficientMatrix0(NuTo::SparseMatrixCSRVector2General<double>& rMatrix, NuTo::FullVector<double,Eigen::Dynamic>& rVector);

    //! @brief ... build global coefficient matrix (stiffness) for primary dofs (e.g displacements, rotations, temperature)
    //! @param rMatrix ... global coefficient matrix
    //! @param rVector ... global equivalent load vector (e.g. due to prescribed displacements)
    NuTo::Error::eError BuildGlobalCoefficientMatrix0(NuTo::SparseMatrixCSRVector2Symmetric<double>& rMatrix, NuTo::FullVector<double,Eigen::Dynamic>& rVector);

    //! @brief ... build global coefficient matrix (damping) for primary dofs (e.g displacements, rotations, temperature)
    //! @param rMatrix ... global coefficient matrix (symmetric)
    //! @param rVector ... global equivalent load vector (e.g. due to prescribed displacements)
    NuTo::Error::eError BuildGlobalCoefficientMatrix1(NuTo::SparseMatrixCSRVector2General<double>& rMatrix, NuTo::FullVector<double,Eigen::Dynamic>& rVector);

    //! @brief ... build global coefficient matrix (mass) for primary dofs (e.g displacements, rotations, temperature)
    //! @param rMatrix ... global coefficient matrix (nonsymmetric)
    //! @param rVector ... global equivalent load vector (e.g. due to prescribed displacements)
    NuTo::Error::eError BuildGlobalCoefficientMatrix2(NuTo::SparseMatrixCSRGeneral<double>& rMatrix, NuTo::FullVector<double,Eigen::Dynamic>& rVector);

    //! @brief ... build global coefficient matrix (mass) for primary dofs (e.g displacements, rotations, temperature)
    //! @param rMatrix ... global coefficient matrix (symmetric)
    //! @param rVector ... global equivalent load vector (e.g. due to prescribed displacements)
    NuTo::Error::eError BuildGlobalCoefficientMatrix2(NuTo::SparseMatrixCSRSymmetric<double>& rMatrix, NuTo::FullVector<double,Eigen::Dynamic>& rVector);

    //! @brief ... build global coefficient matrix (mass) for primary dofs (e.g displacements, rotations, temperature)
    //! @param rMatrix ... global coefficient matrix (nonsymmetric)
    //! @param rVector ... global equivalent load vector (e.g. due to prescribed displacements)
    NuTo::Error::eError BuildGlobalCoefficientMatrix2(NuTo::SparseMatrixCSRVector2General<double>& rMatrix, NuTo::FullVector<double,Eigen::Dynamic>& rVector);

    //! @brief ... build global coefficient matrix (mass) for primary dofs (e.g displacements, rotations, temperature)
    //! @param rMatrix ... global coefficient matrix
    //! @param rVector ... global equivalent load vector (e.g. due to prescribed displacements)
    NuTo::Error::eError BuildGlobalCoefficientMatrix2(NuTo::SparseMatrixCSRVector2Symmetric<double>& rMatrix, NuTo::FullVector<double,Eigen::Dynamic>& rVector);

    //! @brief ... build global external load vector
    //! @param rVector ... external load vector on independent dofs - C external load vector on dependent dofs
    void BuildGlobalExternalLoadVector(int rLoadCase, NuTo::FullVector<double,Eigen::Dynamic>& rVector_j);

    //! @brief ... build global external load vector
    //! @param rVector ... external load vector on independent dofs
    //! @param rVector ... external load vector on dependent dofs
    void BuildGlobalExternalLoadVector(int rLoadCase, NuTo::FullVector<double,Eigen::Dynamic>& rVector_j, NuTo::FullVector<double,Eigen::Dynamic>& rVector_k);

    //! @brief ... build global gradient of the internal potential (e.g. the internal forces)
    //! @param rVector ... global gradient of the internal potential (e.g. internal force vector)
    NuTo::Error::eError BuildGlobalGradientInternalPotentialVector(NuTo::FullVector<double,Eigen::Dynamic>& rVector);

    //! @brief ... based on the global dofs build submatrices of the global coefficent matrix
    //! @param rType   ... matrix type (stiffness damping mass)
    //! @param rMatrixJJ ... submatrix jj (number of active dof x number of active dof)
    //! @param rMatrixJK ... submatrix jk (number of active dof x number of dependent dof)
    virtual Error::eError BuildGlobalCoefficientSubMatricesGeneral(NuTo::StructureBaseEnum::eMatrixType rType, NuTo::SparseMatrix<double>& rMatrixJJ, NuTo::SparseMatrix<double>& rMatrixJK) = 0;

    //! @brief ... based on the global dofs build submatrices of the global coefficent matrix0
    //! @param rType   ... matrix type (stiffness damping mass)
    //! @param rMatrixJJ ... submatrix jj (number of active dof x number of active dof)
    //! @param rMatrixJK ... submatrix jk (number of active dof x number of dependent dof)
    //! @param rMatrixKJ ... submatrix kj (number of dependent dof x number of active dof)
    //! @param rMatrixKK ... submatrix kk (number of dependent dof x number of dependent dof)
    virtual Error::eError BuildGlobalCoefficientSubMatricesGeneral(NuTo::StructureBaseEnum::eMatrixType rType, NuTo::SparseMatrix<double>& rMatrixJJ, NuTo::SparseMatrix<double>& rMatrixJK, NuTo::SparseMatrix<double>& rMatrixKJ, NuTo::SparseMatrix<double>& rMatrixKK) = 0;

    //! @brief ... based on the global dofs build submatrices of the global stiffness matrix0
    //! @brief ... presumes elastic deformation, that is the state variables remain constant
    //! @param rMatrixJJ ... submatrix jj (number of active dof x number of active dof)
    //! @param rMatrixJK ... submatrix jk (number of active dof x number of dependent dof)
    //! @param rMatrixKJ ... submatrix kj (number of dependent dof x number of active dof)
    //! @param rMatrixKK ... submatrix kk (number of dependent dof x number of dependent dof)
    virtual Error::eError BuildGlobalElasticStiffnessSubMatricesGeneral(NuTo::SparseMatrix<double>& rMatrixJJ, NuTo::SparseMatrix<double>& rMatrixJK, NuTo::SparseMatrix<double>& rMatrixKJ, NuTo::SparseMatrix<double>& rMatrixKK) = 0;

    //! @brief ... based on the global dofs build submatrices of the global coefficent matrix
    //! @param rType   ... matrix type (stiffness damping mass)
    //! @param rMatrixJJ ... submatrix jj (number of active dof x number of active dof)
    //! @param rMatrixJK ... submatrix jk (number of active dof x number of dependent dof)
    virtual Error::eError BuildGlobalCoefficientSubMatricesSymmetric(NuTo::StructureBaseEnum::eMatrixType rType, NuTo::SparseMatrix<double>& rMatrixJJ, NuTo::SparseMatrix<double>& rMatrixJK) = 0;

    //! @brief ... based on the global dofs build submatrices of the global coefficent matrix
    //! @param rType   ... matrix type (stiffness damping mass)
    //! @param rMatrixJJ ... submatrix jj (number of active dof x number of active dof)
    //! @param rMatrixJK ... submatrix jk (number of active dof x number of dependent dof)
    //! @param rMatrixKK ... submatrix kk (number of dependent dof x number of dependent dof)
    virtual Error::eError BuildGlobalCoefficientSubMatricesSymmetric(NuTo::StructureBaseEnum::eMatrixType rType, NuTo::SparseMatrix<double>& rMatrixJJ, NuTo::SparseMatrix<double>& rMatrixJK, NuTo::SparseMatrix<double>& rMatrixKK) = 0;

    //! @brief ... based on the global dofs build sub-vectors of the global lumped mass
    //! @param rActiveDofVector ... global lumped mass which corresponds to the active dofs
    //! @param rDependentDofVector ... global lumped mass which corresponds to the dependent dofs
    virtual Error::eError BuildGlobalLumpedHession2(NuTo::FullVector<double,Eigen::Dynamic>& rActiveDofVector,
    		NuTo::FullVector<double,Eigen::Dynamic>& rDependentDofVector) = 0;

    //! @brief ... based on the global dofs build sub-vectors of the global internal potential gradient
    //! @param rActiveDofGradientVector ... global internal potential gradient which corresponds to the active dofs
    //! @param rDependentDofGradientVector ... global internal potential gradient which corresponds to the dependent dofs
    //! @param rUpdateHistoryVariables (update history variables after having calculated the response)
    virtual Error::eError BuildGlobalGradientInternalPotentialSubVectors(NuTo::FullVector<double,Eigen::Dynamic>& rActiveDofGradientVector,
    		NuTo::FullVector<double,Eigen::Dynamic>& rDependentDofGradientVector, bool rUpdateHistoryVariables) = 0;

    //! @brief ... based on the global dofs build sub-vectors of the global internal potential gradient
    //! @brief ... presumes elastic deformation, that is the state variables remain constant
    //! @param rActiveDofGradientVector ... global internal potential gradient which corresponds to the active dofs
    //! @param rDependentDofGradientVector ... global internal potential gradient which corresponds to the dependent dofs
    //! @param update of history variables is not performed
    virtual Error::eError BuildGlobalElasticGradientInternalPotentialSubVectors(NuTo::FullVector<double,Eigen::Dynamic>& rActiveDofGradientVector,
    		NuTo::FullVector<double,Eigen::Dynamic>& rDependentDofGradientVector) = 0;

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
    void NodeGetElements(const NodeBase* rNodePtr, std::vector<ElementBase*>& rElements);

    //! @brief gives the identifier of a node
    //! @param reference to a node
    //! @return identifier
    virtual int NodeGetId(const NodeBase* rNode)const=0;
#endif //SWIG
    //! @brief ... returns the (first) node that has the specified coordinates within the range
    //! @param ... rCoordinates
    //! @param ... rRange
    //! @return ... node id
    int NodeGetIdAtCoordinate(FullVector<double, Eigen::Dynamic>& rCoordinates, double rRange);

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
    virtual void NodeBuildGlobalDofs()=0;

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

    //! @brief extract dof values (e.g. displacements, temperatures to the nodes)
    //! @param rTimeDerivative time derivative (0 disp 1 vel 2 acc)
    //! @param rActiveDofValues ... vector of global active dof values (ordering according to global dofs, size is number of active dofs)
    //! @param rDependentDofValues ... vector of global dependent dof values (ordering according to (global dofs) - (number of active dofs), size is (total number of dofs) - (number of active dofs))
    virtual void NodeExtractDofValues(int rTimeDerivative, NuTo::FullVector<double,Eigen::Dynamic>& rActiveDofValues, NuTo::FullVector<double,Eigen::Dynamic>& rDependentDofValues) const = 0;

    //! @brief extract dof values (e.g. displacements, temperatures to the nodes)
    //! @param rActiveDofValues ... vector of global active dof values (ordering according to global dofs, size is number of active dofs)
    //! @param rDependentDofValues ... vector of global dependent dof values (ordering according to (global dofs) - (number of active dofs), size is (total number of dofs) - (number of active dofs))
    inline void NodeExtractDofValues(NuTo::FullVector<double,Eigen::Dynamic>& rActiveDofValues, NuTo::FullVector<double,Eigen::Dynamic>& rDependentDofValues) const
    {
    	NodeExtractDofValues(0,rActiveDofValues,rDependentDofValues);
    }

    //! @brief write dof values (e.g. displacements, temperatures to the nodes)
    //! @param rTimeDerivative time derivative (0 disp 1 vel 2 acc)
    //! @param rActiveDofValues ... vector of independent dof values (ordering according to global dofs, size is number of active dofs)
    //! @param rDependentDofValues ... vector of dependent  dof values (ordering according to global dofs, size is number of active dofs)
    virtual void NodeMergeDofValues(int rTimeDerivative, const NuTo::FullVector<double,Eigen::Dynamic>& rActiveDofValues, const NuTo::FullVector<double,Eigen::Dynamic>& rDependentDofValues)=0;

    //! @brief write dof values (e.g. displacements, temperatures to the nodes)
    //! @param rActiveDofValues ... vector of independent dof values (ordering according to global dofs, size is number of active dofs)
    //! @param rDependentDofValues ... vector of dependent  dof values (ordering according to global dofs, size is number of active dofs)
    inline void NodeMergeDofValues(const NuTo::FullVector<double,Eigen::Dynamic>& rActiveDofValues, const NuTo::FullVector<double,Eigen::Dynamic>& rDependentDofValues)
    {
    	NodeMergeDofValues(0,rActiveDofValues,rDependentDofValues);
    }

    //! @brief write dof values (e.g. displacements, temperatures to the nodes)
    //! @param rTimeDerivative time derivative (0 disp 1 vel 2 acc)
    //! @param rActiveDofValues ... vector of global dof values (ordering according to global dofs, size is number of active dofs)
    virtual void NodeMergeActiveDofValues(int rTimeDerivative, const NuTo::FullVector<double,Eigen::Dynamic>& rActiveDofValues)=0;

    //! @brief write dof values (e.g. displacements, temperatures to the nodes)
    //! @param rActiveDofValues ... vector of global dof values (ordering according to global dofs, size is number of active dofs)
    inline void NodeMergeActiveDofValues(const NuTo::FullVector<double,Eigen::Dynamic>& rActiveDofValues)
    {
    	NodeMergeActiveDofValues(0,rActiveDofValues);
    }

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

#ifdef ENABLE_VISUALIZE
    //! @brief ... adds all the elements in the vector to the data structure that is finally visualized
    void NodeTotalAddToVisualize(VisualizeUnstructuredGrid& rVisualize, const boost::ptr_list<NuTo::VisualizeComponentBase>& rWhat) const;

    //! @brief ... adds all the nodes in the vector to the data structure that is finally visualized
    void NodeVectorAddToVisualize(VisualizeUnstructuredGrid& rVisualize, const boost::ptr_list<NuTo::VisualizeComponentBase>& rWhat, const std::vector<const NodeBase*>& rNodes) const;
#endif // ENABLE_VISUALIZE
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

    //! @brief calls ElementCoefficientMatrix_0,
    //! renaming only for clarification in mechanical problems for the end user
    void ElementStiffness(int rElementId, NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rResult ,
                          NuTo::FullVector<int,Eigen::Dynamic>& rGlobalDofsRow,
                          NuTo::FullVector<int,Eigen::Dynamic>& rGlobalDofsColumn);

    //! @brief calculates the coefficient matrix for the rTimeDerivative derivative in the differential equation
    //! for a mechanical problem, this corresponds to the stiffness, damping or mass matrix (rTimeDerivative=0,1,2)
    void ElementCoefficientMatrix(int rElementId,
    		                        int rTimeDerivative,
                                    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rResult,
                                    NuTo::FullVector<int,Eigen::Dynamic>& rGlobalDofsRow,
                                    NuTo::FullVector<int,Eigen::Dynamic>& rGlobalDofsColumn);

    //! @brief calculates the coefficient matrix for the 0-th derivative in the differential equation
    //! and compares it to the matrix using central differences
    //! for a mechanical problem, this corresponds to the stiffness matrix
    //! @param rElementId element
    //! @param rDelta  delta step for finite differences
    //! @return maximum difference between analytical and central difference method
    double ElementCoefficientMatrix_0_Check(int rElementId, double rDelta, NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rDifference);

#ifndef SWIG
    //! @brief calculates the coefficient matrix for the 0-th derivative in the differential equation using central differences
    //! @param rElementPtr element
    //! @param rDelta  delta step for finite differences
    //! @param stiffnessCDF  stiffness from central differences (return value, size should be allocated correctly before entering the routine)
    //! @return maximum difference between analytical and central difference method
    void ElementCoefficientMatrix_0_Resforce(ElementBase* rElementPtr, double rDelta, NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& stiffnessCDF);
#endif //SWIG

    //! @brief calculates the coefficient matrix for the 0-th derivative in the differential equation
    //! @param rElementId elementId
    //! @param rDelta  delta step for finite differences
    //! @param stiffnessCDF  stiffness from central differences (return value, size should be allocated correctly before entering the routine)
    //! @return maximum difference between analytical and central difference method
    void ElementCoefficientMatrix_0_Resforce(int rElementId, double rDelta, NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& stiffnessCDF);

    //! @brief calculates the coefficient matrix for the 0-th derivative in the differential equation
    //! and compares it to the matrix using central differences
    //! for a mechanical problem, this corresponds to the stiffness matrix
    //! @param rDelta  delta step for finite differences
    //! @return element with maximum error
    int ElementTotalCoefficientMatrix_0_Check(double rDelta, NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rDifference);

    //! @brief Checks the global CoefficientMatrix_0 and the internal forces vector by the internal energy
    //! @return false, if stiffness is not correct
    virtual bool CheckStiffness();


    //! @brief Compares the global CoefficientMatrix_0 to the matrix using central differences
    //! @param rDelta delta step for finite differences
    //! @param rPrintResult 'true' prints the result
    bool CheckCoefficientMatrix_0(double rDelta, bool rPrintResult);

    //! @brief Compares the element CoefficientMatrix_0 to the matrix using central differences
    //! for all elements
    //! @param rDelta delta step for finite differences
    bool ElementCheckCoefficientMatrix_0(double rDelta);

    //! @brief Compares the element CoefficientMatrix_0 to the matrix using central differences
    //! for a specific element
    //! @param rDelta delta step for finite differences
    //! @param rElementId element id
    //! @param rDifference difference to CDF solution
    //! @param rPrintResult 'true' prints the result
    bool ElementCheckCoefficientMatrix_0(double rDelta,
            int rElementId, NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rDifference, bool rPrintResult);

    //! @brief calculates the gradient of the internal potential
    //! for a mechanical problem, this corresponds to the internal force vector
    void ElementGradientInternalPotential(int rElementId,
                                          NuTo::FullVector<double, Eigen::Dynamic>& rResult,
                                          NuTo::FullVector<int, Eigen::Dynamic>& rGlobalDofsRow);

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
#endif //SWIG




    //! @brief calculates the engineering strain
    //! @param rElemIdent  element number
    //! @param rEngineerungStrain engineering strain (return value, always 6xnumIp matrix)
    void ElementGetEngineeringStrain(int rElementId, NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rEngineeringStrain);

    //! @brief calculates the engineering plastic strain
    //! @param rElemIdent  element number
    //! @param rEngineerungStrain engineering plastic strain (return value, always 6xnumIp matrix)
    void ElementGetEngineeringPlasticStrain(int rElementId, NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rEngineeringPlasticStrain);

    //! @brief calculates the engineering stress
    //! @param rElemIdent  element number
    //! @param rEingineeringStress Engineering Stress (return value, always 6xnumIp matrix)
    void ElementGetEngineeringStress(int rElementId, NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rEngineeringStress);

    //! @brief calculates the damage
    //! @param rElemIdent  identifier for the element
    //! @param rDamage (return value, always 1xnumIp matrix)
    void ElementGetDamage(int rElementId, NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rDamage);

    //! @brief calculates the global integration point coordinates
    //! @param rElemIdent  identifier for the element
    //! @param rCoordinates integration point coordinates (return value, always 3xnumIp matrix)
    void ElementGetIntegrationPointCoordinates(int rElementId, NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rCoordinates);

    //! @brief calculates the maximum damage in all elements
    //! @param rElemIdent  identifier for the element
    //! @return max damage value
    double ElementTotalGetMaxDamage();

    //! @brief updates the history data of a all elements
    NuTo::Error::eError ElementTotalUpdateStaticData();

    //! @brief updates the temprory static data of a all elements
    //! its is a const function, since only mutuable data (instead of const) is updated (kind of temporary data)
    NuTo::Error::eError ElementTotalUpdateTmpStaticData();

    //! @brief saves static data of a all elements
    NuTo::Error::eError ElementFatigueSaveStaticData();

    //! @brief restores static data of a all elements
    NuTo::Error::eError ElementFatigueRestoreStaticData();

    //! @brief extrapolates static data of a all elements
    NuTo::Error::eError ElementFatigueExtrapolateStaticData();

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

    //! @brief writes the Lagrange multiplier and Slack variables (inequalities) of a constraint to the prescribed matrix
    //! @param ConstraintId constraint id
    //! @param rMultiplier Lagrange multiplier (first col Lagrange, evtl. second col Slackvariables)
    void ConstraintLagrangeGetMultiplier(int ConstraintId, NuTo::FullVector<double,Eigen::Dynamic>& rMultiplier)const;

    //! @brief sets the penalty stiffness of the augmented Lagragian to the prescribed value
    //! @param ConstraintId constraint id
    //! @param rPenalty penalty parameter
    void ConstraintLagrangeSetPenaltyStiffness(int ConstraintId, double rPenalty);

    //! @brief adds a displacement constraint equation for a node group solved using Lagrange multiplier
    //! @param rGroupId group id
    //! @param rDirection direction of the constraint (in 2D a point with 2 entries, in 3D 3 entries, in 1D not used)
    //! @param rValue prescribed value (e.g. zero to fix a displacement to zero)
    //! @return integer id to delete or modify the constraint
    int ConstraintLagrangeSetDisplacementNodeGroup(int rGroupId, const NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rDirection, const std::string& rSign, double rValue);

#ifndef SWIG
    //! @brief adds a displacement constraint equation for a node group solved using Lagrange multiplier
    //! @param rGroup group pointer
    //! @param rDirection direction of the constraint (in 2D a point with 2 entries, in 3D 3 entries, in 1D not used)
    //! @param rValue prescribed value (e.g. zero to fix a displacement to zero)
    //! @return integer id to delete or modify the constraint
    int ConstraintLagrangeSetDisplacementNodeGroup(Group<NodeBase>* rGroup, const NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rDirection, NuTo::Constraint::eEquationSign rEquationSign, double rValue);

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
    //! @param rAttribute displacements, rotations, temperatures
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

    //! @brief returns the number of constraint equations
    //! @return number of constraints
    int ConstraintGetNumLinearConstraints()const;

    //! @brief calculates the constraint matrix that builds relations between the nodal dagrees of freedom (before gauss elimination)
    //! @param rConstraintMatrix constraint matrix
    void ConstraintGetConstraintMatrixBeforeGaussElimination(NuTo::SparseMatrixCSRGeneral<double>& rConstraintMatrix);

    //! @brief returns the constraint matrix  (after gauss elimination)
    //! @param rConstraintMatrix constraint matrix
    void ConstraintGetConstraintMatrixAfterGaussElimination(NuTo::SparseMatrixCSRGeneral<double>& rConstraintMatrix);

    //! @brief returns the constraint vector after gauss elimination
    //! rConstraintMatrix*DOFS = RHS
    //! @param rRHS rhs
    void ConstraintGetRHSAfterGaussElimination(NuTo::FullVector<double,Eigen::Dynamic>& rRHS);

    //! @brief returns the constraint vector after gauss elimination
    //! rConstraintMatrix*DOFS = RHS
    //! @param rConstraintMatrix constraint matrix
    void ConstraintGetRHSBeforeGaussElimination(NuTo::FullVector<double,Eigen::Dynamic>& rhsBeforeGaussElimination);

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

#ifndef SWIG
    //!@brief number the free DOFS in the constraints (Lagrange multipliers)
    //!@param rDOF current maximum DOF number, increased in the number
    void ConstraintNumberGlobalDofs(int& rDOF);

    //! @brief renumber the dofs of the Lagrange multipliers according to predefined ordering
    //! @param rMappingInitialToNewOrdering ... mapping from initial ordering to the new ordering
    void ConstraintRenumberGlobalDofs(const std::vector<int>& mappingInitialToNewOrdering);

    //! @brief extract dof values from the node (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    void ConstraintExtractGlobalDofValues(NuTo::FullVector<double,Eigen::Dynamic>& rActiveDofValues, NuTo::FullVector<double,Eigen::Dynamic>& rDependentDofValues)const;

    //! @brief write dof values to the Lagrange multipliers (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    void ConstraintMergeGlobalDofValues(const NuTo::FullVector<double,Eigen::Dynamic>& rActiveDofValues, const NuTo::FullVector<double,Eigen::Dynamic>& dependentDofValues);

#endif


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
    void ConstraintLinearEquationCreate(int rConstraint, int rNode, NuTo::Node::eAttributes rDofType, int rDofComponent, double rCoefficient, double rRHS = 0);
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
    void ConstraintLinearEquationAddTerm(int rConstraint, int rNode, NuTo::Node::eAttributes rDofType, int rDofComponent, double rCoefficient);

    //! @brief ... add the contribution of Lagrange multipliers to the global system of equations
    //! @param rMatrixJJ ... matrix jj
    //! @param rMatrixJK ... matrix jk
    void ConstraintsBuildGlobalCoefficientSubMatrices0General(NuTo::SparseMatrix<double>& rMatrixJJ, NuTo::SparseMatrix<double>& rMatrixJK)const;

    //! @brief ... add the contribution of Lagrange multipliers to the global system of equations
    //! @param rMatrixJJ ... matrix jj
    //! @param rMatrixJK ... matrix jk
    //! @param rMatrixKJ ... matrix kj
    //! @param rMatrixKK ... matrix kk
    void ConstraintBuildGlobalCoefficientSubMatrices0General(NuTo::SparseMatrix<double>& rMatrixJJ, NuTo::SparseMatrix<double>& rMatrixJK, NuTo::SparseMatrix<double>& rMatrixKJ, NuTo::SparseMatrix<double>& rMatrixKK)const;

    //! @brief ... add the contribution of Lagrange multipliers to the global system of equations
    //! @param rMatrixJJ ... matrix jj
    //! @param rMatrixJK ... matrix jk
    void ConstraintBuildGlobalCoefficientSubMatrices0Symmetric(NuTo::SparseMatrix<double>& rMatrixJJ, NuTo::SparseMatrix<double>& rMatrixJK)const;

    //! @brief ... add the contribution of Lagrange multipliers to the global system of equations
    //! @param rMatrixJJ ... matrix jj
    //! @param rMatrixJK ... matrix jk
    //! @param rMatrixKK ... matrix kk
    void ConstraintBuildGlobalCoefficientSubMatrices0Symmetric(NuTo::SparseMatrix<double>& rMatrixJJ, NuTo::SparseMatrix<double>& rMatrixJK, NuTo::SparseMatrix<double>& rMatrixKK) const;

    //! @brief ... add the contribution of Lagrange multipliers to the global system of equations
    //! @param rActiveDofGradientVector ... gradient of active dofs
    //! @param rDependentDofGradientVector ... gradient of dependent dofs
    void ConstraintBuildGlobalGradientInternalPotentialSubVectors(NuTo::FullVector<double,Eigen::Dynamic>& rActiveDofGradientVector, NuTo::FullVector<double,Eigen::Dynamic>& rDependentDofGradientVector) const;

    //! @brief ... calculates the potential due to the addition of Lagrange terms
    double ConstraintTotalGetTotalEnergy()const;

    //! @brief info about the nodes in the Structure
    void ConstraintInfo(int mVerboseLevel)const;

#endif

private:
    //! @brief ... convert input string to dof type and dof component
    //! @brief rDof ... input string
    //! @param rDofType ... type of dof (e.g DISPLACEMENTS, ROTATIONS, TEMPERATURES)
    //! @param rDofComponent ... dof component (0, 1, 2)
    void ConstraintEquationGetDofInformationFromString(const std::string& rDof, NuTo::Node::eAttributes& rDofType, int& rDofComponent);

public:
    //*************************************************
    //************ Load routines        ***************
    //***  defined in structures/StructureLoad.cpp  ***
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
    bool ConstitutiveLawGetVariableBool(int rIdent, const std::string& rIdentifier) const;

    //! @brief ... sets a variable of the constitutive law which is selected by a string
    //! @param rIdent ... constitutive law identifier
    //! @param rIdentifier ... String to identify the requested variable
    //! @param rValue ... new value for requested variable
    void ConstitutiveLawSetVariableBool(int rIdent, const std::string& rIdentifier, bool rValue);

    //! @brief ... gets a variable of the constitutive law which is selected by a string
    //! @param rIdent ... constitutive law identifier
    //! @param rIdentifier ... String to identify the requested variable
    //! @return ... value of the requested variable
    double ConstitutiveLawGetVariableDouble(int rIdent, const std::string& rIdentifier) const;

    //! @brief ... sets a variable of the constitutive law which is selected by a string
    //! @param rIdent ... constitutive law identifier
    //! @param rIdentifier ... String to identify the requested variable
    //! @param rValue ... new value for requested variable
    void ConstitutiveLawSetVariableDouble(int rIdent, const std::string& rIdentifier, double rValue);

    //! @brief ... gets a variable of the constitutive law which is selected by a string
    //! @param rIdent ... constitutive law identifier
    //! @param rIdentifier ... String to identify the requested variable
    //! @return ... value of the requested variable
    NuTo::FullVector<double, Eigen::Dynamic> ConstitutiveLawGetVariableFullVectorDouble(int rIdent, const std::string& rIdentifier) const;

    //! @brief ... sets a variable of the constitutive law which is selected by a string
    //! @param rIdent ... constitutive law identifier
    //! @param rIdentifier ... String to identify the requested variable
    //! @param rValue ... new value for requested variable
    void ConstitutiveLawSetVariableFullVectorDouble(int rIdent, const std::string& rIdentifier, NuTo::FullVector<double, Eigen::Dynamic>  rValue);


#ifndef SWIG

    //! @brief ... gets a variable of the constitutive law which is selected by an enum
    //! @param rIdent ... constitutive law identifier
    //! @param rIdentifier ... Enum to identify the requested variable
    //! @return ... value of the requested variable
    bool ConstitutiveLawGetVariableBool(int rIdent, Constitutive::eConstitutiveVariable rIdentifier) const;

    //! @brief ... sets a variable of the constitutive law which is selected by an enum
    //! @param rIdent ... constitutive law identifier
    //! @param rIdentifier ... Enum to identify the requested variable
    //! @param rValue ... new value for requested variable
    void ConstitutiveLawSetVariableBool(int rIdent, Constitutive::eConstitutiveVariable rIdentifier, bool rValue);

    //! @brief ... gets a variable of the constitutive law which is selected by an enum
    //! @param rIdent ... constitutive law identifier
    //! @param rIdentifier ... Enum to identify the requested variable
    //! @return ... value of the requested variable
    double ConstitutiveLawGetVariableDouble(int rIdent, Constitutive::eConstitutiveVariable rIdentifier) const;

    //! @brief ... sets a variable of the constitutive law which is selected by an enum
    //! @param rIdent ... constitutive law identifier
    //! @param rIdentifier ... Enum to identify the requested variable
    //! @param rValue ... new value for requested variable
    void ConstitutiveLawSetVariableDouble(int rIdent, Constitutive::eConstitutiveVariable rIdentifier, double rValue);

    //! @brief ... gets a variable of the constitutive law which is selected by an enum
    //! @param rIdent ... constitutive law identifier
    //! @param rIdentifier ... Enum to identify the requested variable
    //! @return ... value of the requested variable
    NuTo::FullVector<double, Eigen::Dynamic> ConstitutiveLawGetVariableFullVectorDouble(int rIdent, Constitutive::eConstitutiveVariable rIdentifier) const;

    //! @brief ... sets a variable of the constitutive law which is selected by an enum
    //! @param rIdent ... constitutive law identifier
    //! @param rIdentifier ... Enum to identify the requested variable
    //! @param rValue ... new value for requested variable
    void ConstitutiveLawSetVariableFullVectorDouble(int rIdent, Constitutive::eConstitutiveVariable rIdentifier, NuTo::FullVector<double, Eigen::Dynamic>  rValue);
#endif // SWIG

    //! @brief ... set density
    //! @param rIdent ... constitutive law identifier
    //! @param rRho ... density
    void ConstitutiveLawSetDensity(int rIdent, double rRho);

    //! @brief ... get densuty
    //! @param rIdent ... constitutive law identifier
    //! @return ... density
    double ConstitutiveLawGetDensity(int rIdent) const;

    //! @brief ... set Young's modulus
    //! @param rIdent ... constitutive law identifier
    //! @param rE ... Young's modulus
    void ConstitutiveLawSetYoungsModulus(int rIdent, double rE);

    //! @brief ... get Young's modulus
    //! @param rIdent ... constitutive law identifier
    //! @return ... Young's modulus
    double ConstitutiveLawGetYoungsModulus(int rIdent) const;

    //! @brief ... set Poisson's ratio
    //! @param rIdent ... constitutive law identifier
    //! @param rNu ... Poisson's ratio
    void ConstitutiveLawSetPoissonsRatio(int rIdent, double rNu);

    //! @brief ... get Poisson's ratio
    //! @param rIdent ... constitutive law identifier
    //! @return ... Poisson's ratio
    double ConstitutiveLawGetPoissonsRatio(int rIdent) const;

    //! @brief ... get initial yield strength
    //! @return ... yield strength
    double ConstitutiveLawGetInitialYieldStrength(int rIdent) const;

    //! @brief ... set initial yield strength
    //! @param rSigma ...  yield strength
    void ConstitutiveLawSetInitialYieldStrength(int rIdent, double rSigma);

    //! @brief ... get yield strength for multilinear response
    //! @return ... first column: equivalent plastic strain
    //! @return ... second column: corresponding yield strength
    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> ConstitutiveLawGetYieldStrength(int rIdent) const;

    //! @brief ... add yield strength
    //! @param rEpsilon ...  equivalent plastic strain
    //! @param rSigma ...  yield strength
    void ConstitutiveLawAddYieldStrength(int rIdent, double rEpsilon, double rSigma);

    //! @brief ... get initial hardening modulus
    //! @return ... hardening modulus
    double ConstitutiveLawGetInitialHardeningModulus(int rIdent) const;

    //! @brief ... set hardening modulus
    //! @param rH ...  hardening modulus
    void ConstitutiveLawSetInitialHardeningModulus(int rIdent, double rH);

    //! @brief ... get hardening value
    //! @return ... hardening value
    double ConstitutiveLawGetHardeningValue(int rIdent) const;

    //! @brief ... set hardening value
    //! @return ... hardening value
    void ConstitutiveLawSetHardeningValue(int rIdent, double rHardening);

    //! @brief ... get hardening exponent
    //! @return ... hardening exponent
    double ConstitutiveLawGetHardeningExponent(int rIdent) const;

    //! @brief ... set hardening exponent
    //! @return ... hardening exponent
    void ConstitutiveLawSetHardeningExponent(int rIdent, double rHardeningExponent);

    //! @brief ... get hardening modulus for multilinear response
    //! @return ... first column: equivalent plastic strain
    //! @return ... second column: corresponding hardening modulus
    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> ConstitutiveLawGetHardeningModulus(int rIdent) const;

    //! @brief ... add hardening modulus
    //! @param rEpsilon ...  equivalent plastic strain
    //! @param rSigma ...  hardening modulus
    void ConstitutiveLawAddHardeningModulus(int rIdent, double rEpsilon, double rH);

    //! @brief ... get nonlocal radius
    //! @return ... nonlocal radius
    double ConstitutiveLawGetNonlocalRadius(int rIdent) const;

    //! @brief ... set nonlocal radius
    //! @param rRadius ...  nonlocal radius
    void ConstitutiveLawSetNonlocalRadius(int rIdent, double rRadius);

    //! @brief ... get nonlocal radius parameter
    //! @return ... nonlocal radius parameter
    double ConstitutiveLawGetNonlocalRadiusParameter(int rIdent) const;

    //! @brief ... set nonlocal radius parameter
    //! @param rRadius ...  nonlocal radius parameter
    void ConstitutiveLawSetNonlocalRadiusParameter(int rIdent, double rRadiusParameter);


    //! @brief ... get tensile strength
    //! @param rTensileStrength ...  tensile strength
    double ConstitutiveLawGetTensileStrength(int rIdent);

    //! @brief ... set tensile strength
    //! @param rTensileStrength ...  tensile strength
    void ConstitutiveLawSetTensileStrength(int rIdent, double rTensileStrength);

    //! @brief ... get shear strength
    //! @param rIdent ...  constitutive model
    double ConstitutiveLawGetShearStrength(int rIdent);

    //! @brief ... set shear strength
    //! @param rShearStrength ...  shear strength
    void ConstitutiveLawSetShearStrength(int rIdent, double rShearStrength);

    //! @brief ... get compressive strength
    //! @param rCompressiveStrength ...  compressive strength
    double ConstitutiveLawGetCompressiveStrength(int rIdent);

    //! @brief ... set compressive strength
    //! @param rCompressiveStrength ...  compressive strength
    void ConstitutiveLawSetCompressiveStrength(int rIdent, double rCompressiveStrength);

    //! @brief ... get biaxial compressive strength
    //! @param rBiaxialCompressiveStrength ...  biaxial compressive strength
    double ConstitutiveLawGetBiaxialCompressiveStrength(int rIdent);

    //! @brief ... set biaxial compressive strength
    //! @param rBiaxialCompressiveStrength ...  biaxial compressive strength
    void ConstitutiveLawSetBiaxialCompressiveStrength(int rIdent, double rBiaxialCompressiveStrength);

    //! @brief ... get fracture energy
    //! @param rFractureEnergy ...  fracture energy
    double ConstitutiveLawGetFractureEnergy(int rIdent);

    //! @brief ... set fracture energy
    //! @param rFractureEnergy ...  fracture energy
    void ConstitutiveLawSetFractureEnergy(int rIdent, double rFractureEnergy);

    //! @brief ... get friction coefficient
    //! @param rIdent ...  constitutive model
    double ConstitutiveLawGetFrictionCoefficient(int rIdent);

    //! @brief ... set friction coefficient
    //! @param rFrictionCoefficient ...  friction coefficient
    void ConstitutiveLawSetFrictionCoefficient(int rIdent, double rFrictionCoefficient);

    //! @brief ... set heat cpacity
    void ConstitutiveLawSetHeatCapacity(int rIdent, double rHeatCapacity);

    //! @brief ... get heat cpacity
    double ConstitutiveLawGetheatCapacity(int rIdent) const;

    //! @brief ... set thermal conductivity
    void ConstitutiveLawSetThermalConductivity(int rIdent, double rThermalConductivity);

    //! @brief ... get thermal conductivity
    double ConstitutiveLawGetThermalConductivity(int rIdent) const;

    //! @brief ... get viscosity
    //! @return ... viscosity
    double ConstitutiveLawGetViscosity(int rIdent) const;

    //! @brief ... set viscosity
    //! @param ... viscosity
    void ConstitutiveLawSetViscosity(int rIdent, double rViscosity);

    //! @brief ... get viscosity exponent
    //! @return ... viscosity exponent
    double ConstitutiveLawGetViscosityExponent(int rIdent) const;

    //! @brief ... set fatigue extrapolation flag
    //! @param ... fatigue extrapolation flag
    void ConstitutiveLawSetFatigueExtrapolation(int rIdent, bool rFatigueExtrapolation);

    //! @brief ... get fatigue extrapolation flag
    //! @param ... fatigue extrapolation flag
    bool ConstitutiveLawGetFatigueExtrapolation(int rIdent) const;

    //! @brief ... set viscosity exponent
    //! @param ... viscosity exponent
    void ConstitutiveLawSetViscosityExponent(int rIdent, double rViscosityExponent);

    //! @brief ... get damage distribution (determines the portion of damage via viscoplasticity and plasticity)
    //! @return ... damage distribution
    double ConstitutiveLawGetDamageDistribution(int rIdent) const;

    //! @brief ... set damage distribution (determines the portion of damage via viscoplasticity and plasticity)
    //! @param ... damage distribution
    void ConstitutiveLawSetDamageDistribution(int rIdent, double rDamageDistribution);


    //! @brief ... get damage law
    //! @return ... damage law
    NuTo::FullVector<double, Eigen::Dynamic> ConstitutiveLawGetDamageLaw(int rIdent) const;

    //! @brief ... set damage law
    //! @param rDamageLaw ... damage law <BR>
    //! ============================================================================================<BR>
    //! rDamageLaw[0] = Constitutive::Constitutive::eDamageLawType::ISOTROPIC_NO_SOFTENING <BR>
    //! w(k) = 1 - e_0/k <BR>
    //! rDamageLaw[1] = e_0 // strain at elastic limit <BR>
    //! ============================================================================================<BR>
    //! rDamageLaw[0] = Constitutive::Constitutive::eDamageLawType::ISOTROPIC_LINEAR_SOFTENING <BR>
    //! w(k) = e_c/k * (k-e_0) / (e_c-e_0) <BR>
    //! rDamageLaw[1] = e_0 // strain at elastic limit <BR>
    //! rDamageLaw[1] = e_c // strain at full damage <BR>
    //! ============================================================================================<BR>
    //! rDamageLaw[0] = Constitutive::Constitutive::eDamageLawType::ISOTROPIC_EXPONENTIAL_SOFTENING <BR>
    //! w(k) = 1 - e_0/k exp{ (e_0 - k) / e_f } <BR>
    //! rDamageLaw[1] = e_0 // strain at elastic limit <BR>
    //! rDamageLaw[1] = e_f // post-peak slope parameter <BR>
    void ConstitutiveLawSetDamageLaw(int rIdent, const NuTo::FullVector<double, Eigen::Dynamic> rDamageLaw);

    //! @brief ... get viscoplastic yield surface offset with respect to the plastic yield surface
    //! @return ... viscoplastic yield surface offset
    double ConstitutiveLawGetViscoplasticYieldSurfaceOffset(int rIdent) const;

    //! @brief ... set viscoplastic yield surface offset with respect to the plastic yield surface
    //! @param ... viscoplastic yield surface offset
    void ConstitutiveLawSetViscoplasticYieldSurfaceOffset(int rIdent, double rViscoplasticYieldSurfaceOffset);

    //! @brief ... gets the equilibrium water volume fraction depend on the relative humidity
    //! @param rIdent ... constitutive law identifier
    //! @param rRelativeHumidity ... relative humidity
    //! @param rCoeffs ... polynomial coefficients of the sorption curve
    //! @return ... equilibrium water volume fraction
    double ConstitutiveLawGetEquilibriumWaterVolumeFraction(int rIdent, double rRelativeHumidity, NuTo::FullVector<double,Eigen::Dynamic> rCoeffs) const;

#ifndef SWIG

    //! @brief ... create a new section
    //! @param rIdent ... section identifier
    //! @param rType ... section type
    void ConstitutiveLawCreate(int rIdent, Constitutive::eConstitutiveType rType);

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
    //! @param rThickness ... cross-section thickness
    void SectionSetThickness(int rId, double rThickness);

    //! @brief ... get section thickness
    //! @param rIdent ... section identifier
    //! @return section thickness
    double SectionGetThickness(int rId) const;

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
    int SectionCreate(Section::eSectionType rType);

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
    int GroupCreate(NuTo::Groups::eGroupId rEnumType);

    //! @brief ... Creates a group for the structure
    //! @param ... rIdent identifier for the group
    //! @param ... rType  type of the group
    void GroupCreate(int id, NuTo::Groups::eGroupId rEnumType);
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
    void GroupAddNodeCoordinateRange(int rIdentGroup, int rDirection, double rMin, double rMax);

    //! @brief ... Adds an element to an element group
    //! @param ... rIdentGroup identifier for the group
    //! @param ... rIdentNode  identifier for the element
    void GroupAddElement(int rIdentGroup, int rIdElement);

    //! @brief ... Adds all elements to a group whose nodes are in the given node group
    //! @param ... rElementGroupId identifier for the element group
    //! @param ... rNodeGroupId idenbtifier for the node group
    //! @param ... rHaveAllNodes if set to true, the element is only selected when all element nodes are in the node group, if set
    //! to false, the element is select if at least one node is in the node group
    void GroupAddElementsFromNodes(int rElementGroupId, int rNodeGroupId, bool rHaveAllNodes);

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
    NuTo::IntegrationTypeBase* GetPtrIntegrationType(NuTo::IntegrationType::eIntegrationType rIdentIntegrationType);
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
    void SetNumExtrapolatedCycles(int rNumber);

    //! @brief get the number of cycles to be extrapolated in the cycle jump routine
    int GetNumExtrapolatedCycles() const;

    //! @brief absolute tolerance for entries of the global stiffness matrix (coefficientMatrix0)
    //! values smaller than that one will not be added to the global matrix
    void SetToleranceStiffnessEntries(double rToleranceStiffnessEntries);

    //! @brief absolute tolerance for entries of the global stiffness matrix (coefficientMatrix0)
    //! values smaller than that one will not be added to the global matrix
    double GetToleranceStiffnessEntries()const;

    //! @brief returns the number of degrees of freedom
    //! @return ... number of degrees of freedom
    int GetNumDofs()const;

    //! @brief returns the number of active degrees of freedom
    //! @return ... number of active degrees of freedom
    int GetNumActiveDofs()const;

    //! @brief returns the a reference to the constraint matrix
    const NuTo::SparseMatrixCSRGeneral<double>& GetConstraintMatrix()const;

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

    //! @brief do a postprocessing step after each converged load step (for Newton Raphson iteration) overload this function to use Newton Raphson
    //virtual void PostProcessDataAfterUpdate(int rLoadStep, int rNumNewtonIterations, double rLoadFactor, double rDeltaLoadFactor, double rResidual);

    //! @brief do a postprocessing step after each converged load step (for Newton Raphson iteration) overload this function to use Newton Raphson
    //virtual void PostProcessDataAfterConvergence(int rLoadStep, int rNumNewtonIterations, double rLoadFactor, double rDeltaLoadFactor, double rResidual);

    //! @brief do a postprocessing step after each line search within the load step(for Newton Raphson iteration) overload this function to use Newton Raphson
    //virtual void PostProcessDataAfterLineSearch(int rLoadStep, int rNewtonIteration, double rLineSearchFactor, double rLoadFactor, double rResidual, const NuTo::FullVector<double,Eigen::Dynamic>& rResidualVector);

    //! @brief do a postprocessing step after each line search within the load step(for Newton Raphson iteration) overload this function to use Newton Raphson
    //virtual void PostProcessDataInLineSearch(int rLoadStep, int rNewtonIteration, double rLineSearchFactor, double rLoadFactor, double rResidual, double rPrevResidual);

    //! @brief initialize some stuff before a new load step (e.g. output directories for visualization, if required)
    //virtual void InitBeforeNewLoadStep(int rLoadStep);


    void SetUpdateTmpStaticDataRequired()
    {
    	mUpdateTmpStaticDataRequired = true;
    }

    //! @brief is only true for structure used as multiscale (structure in a structure)
    //! @parameters rTypeOfSpecimen 0 box, 1 dogbone
    //! @parameters rBoundingBox box for the spheres (3*2 matrix)
    //! @parameters rSeed seed for the random number generator
    //! @parameters rRadiusBoundaryParticles radius particles simulated on the boundary
    //! @parameters rDistanceBoundaryParticles distance of the boundary particles
    //! @parameters rTypeOfSpecimen 0 box, 1 dogbone
    //! @return ... matrix with spheres (coordinates x y z and radius)
    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> CreateSpheresOnSpecimenBoundary(int rTypeOfSpecimen, NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rBoundingBox, int rSeed,
    		double rRadiusBoundaryParticles, double rDistanceBoundaryParticles);


    //! @brief is only true for structure used as multiscale (structure in a structure)
    //! @parameters rTypeOfSpecimen 0 box, 1 dogbone
    //! @parameters rBoundingBox box for the spheres (3*2 matrix)
    //! @parameters rRelParticleVolume percentage of particle volume inside the specimen
    //! @parameters rGradingCurve matrix with each line min_diameter, max_diameter, volume percentage of that sieve size
    //! @parameters relativeDistance scaling factor to increase the diameter when inserting the sphere to ensure a minimum distance
    //! @parameters absoluteDistance distance to increase the diameter when inserting the sphere to ensure a minimum distance
    //! @parameters rSeed seed for the random number generator
    //! @parameters rSpheresBoundary particles simulated on the boundary e.g. created with CreateSpheresOnBoxBoundary (they do not contribute to the grading curve)
    //! @return ... matrix with spheres (coordinates x y z and radius)
    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> CreateSpheresInSpecimen(int rTypeOfSpecimen, NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rBoundingBox, double rRelParticleVolume, NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rGradingCurve,
    		double relativeDistance, double absoluteDistance, int rSeed, NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rSpheresBoundary);

    //! @brief cut spheres at a given z-coordinate to create circles (in 2D)
    //! @parameters rSpheres matrix with the spheres (x,y,z,r)
    //! @parameters rZCoord z coordinate (where to cut)
    //! @parameters rMinRadius minimal radius of the circle
    //! @return ... matrix with the circles (x,y,r)
    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> CutSpheresZ(NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rSpheres, double rZCoord, double rMinRadius);

    //! @brief sets the Hessian to be constant or variable
    //! @parameters rTimeDerivative (0 = stiffness, 1 damping, 2 mass)
    //! @parameters rValue (true = const false=variable)
    void SetHessianConstant(int rTimeDerivative, bool rValue);

	//! @brief sets the Hessian to be constant or variable
	//! @parameters rTimeDerivative (0 = stiffness, 1 damping, 2 mass)
    //! @return (true = const false=variable)
	bool GetHessianConstant(int rTimeDerivative)const;


protected:
    //! @brief ... number of time derivatives (0 : static, 1: velocities, 2: accelerations)
	int mNumTimeDerivatives;

	//! @brief ... storing the beginning of the time increment
	double mPrevTime;

    //! @brief ... storing the end of the time increment (current time)
	double mTime;

    int mDimension;

    //! @brief ... number of cycles applied for extrapolation in the cycle jump.
    int mNumExtrapolatedCycles;

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
#ifdef ENABLE_VISUALIZE
    //! @brief ... map storing the components (displacements, strains, nonlocal weights etc) to be included in the output (VTK) file
    boost::ptr_list<NuTo::VisualizeComponentBase> mVisualizeComponents;
#endif //ENABLE_VISUALIZE
    //! @brief ... total number of degrees of freedom of the structure
    int mNumDofs;

    //! @brief ... active dofs
    int mNumActiveDofs;

    //!brief ... renumbering of nodal DOFs required or not
    bool mNodeNumberingRequired;

    //! @brief constraint matrix relating the prescibed nodal unknowns to the free parameters
    SparseMatrixCSRGeneral<double> mConstraintMatrix;

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
    SparseMatrixCSRGeneral<double> mConstraintMappingRHS;

    //! @brief right hand side of the constraint equations
    FullVector<double,Eigen::Dynamic> mConstraintRHS;

    //! @brief is set to true, if at least one constitutive model requires an update of tmpStaticData before stress and stiffness routines are called
    bool mHaveTmpStaticData;

    //! @brief is set to false, if the structure is changed (nodes, elements) or (DOFs at the nodes)
    bool mUpdateTmpStaticDataRequired;

    //! @brief absolute tolerance for entries of the global stiffness matrix (coefficientMatrix0)
    //! values smaller than that one will not be added to the global matrix
    double mToleranceStiffnessEntries;

    //! @brief parameters of the time integration scheme indicating, if the hessian is constant (0 stiffness, 1 damping, 2 mass)
    //! note that if a matrix is constant, the corresponding term is no longer considered in the gradient calculation
    //! because the term can be obtained from the global solution procedure
    bool mHessianConstant[3];

#ifdef _OPENMP
    //@brief maximum independent sets used for parallel assembly of the stiffness resforce etc.
    mutable std::vector<std::vector<ElementBase*> > mMIS;
    //set to true to use MIS (enables parallel adding of element matrices to the global matrix - faster for elements with almost identical effort)
    //set to false to assemble the global stiffness matrix using a barrier in OMP - faster for elements with very uneven effort)
    bool mUseMIS;
    //@brief number of processors used in an openmp simulation
    int mNumProcessors;
#endif

    //! @brief logger class to redirect the output to some file or the console (or both), can be changed even for const routines
    mutable Logger mLogger;

    //! @brief ... standard constructor just for the serialization routine
    StructureBase()
    {}

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

    //! @brief ... check for dof numbering and build of tmpStaticData
    void BuildGlobalCoefficientMatrixCheck();


    //! @brief ... inserts a particle into subboxes to increase efficiency when performing overlap checks
    void InsertParticleIntoBox(NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rParticles, int rTheParticle, std::vector<std::vector<int > >& rSubBox, std::array<int,3>& rNSubBox,std::array<double,3>& rLSubBox, FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rBoundingBox);
};
} //namespace NuTo
#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
#include <boost/serialization/assume_abstract.hpp>
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::StructureBase)
#endif // SWIG
#endif  // ENABLE_SERIALIZATION
#endif // STRUCTUREBASE_H
