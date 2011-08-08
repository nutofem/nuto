// $Id$

#ifndef STRUCTUREBASE_H
#define STRUCTUREBASE_H

#include <ctime>

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

#include <boost/ptr_container/ptr_map.hpp>
#include <boost/ptr_container/ptr_list.hpp>

#include <string>
#include "nuto/base/Logger.h"
#include "nuto/base/NuToObject.h"
#include "nuto/math/SparseMatrixCSRGeneral_Def.h"
#include "nuto/mechanics/constitutive/ConstitutiveBase.h"
#include "nuto/mechanics/constraints/ConstraintBase.h"
#include "nuto/mechanics/groups/GroupBase.h"
#include "nuto/mechanics/elements/IpDataEnum.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeBase.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeEnum.h"
#include "nuto/mechanics/loads/LoadBase.h"
#include "nuto/mechanics/nodes/NodeEnum.h"
#include "nuto/mechanics/sections/SectionBase.h"
#include "nuto/mechanics/sections/SectionEnum.h"

#ifdef ENABLE_VISUALIZE
#include "nuto/visualize/VisualizeBase.h"
#include "nuto/visualize/VisualizeComponentBase.h"
#endif // ENABLE_VISUALIZE

namespace NuTo
{
class ElementBase;
template <class T> class FullMatrix;
class NodeBase;
template<class T> class SparseMatrixCSRSymmetric;
template<class T> class SparseMatrixCSRVector2General;
template<class T> class SparseMatrixCSRVector2Symmetric;
class EngineeringStrain2D;
class NewtonRaphsonAuxRoutinesBase;
class CrackBase;
class ConstitutiveStaticDataMultiscale2DPlaneStrain;
class VisualizeUnstructuredGrid;

class VisualizeComponentBase;

//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... standard abstract class for all mechanical structures
class StructureBase : public NuToObject
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
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
    int GetDimension();

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

    //! @brief ... Add visualization of principal stresses to the internal list, which is finally exported via the ExportVtkDataFile command
    void AddVisualizationComponentPrincipalEngineeringStress();

    //! @brief ... Add nonlocal weights to the internal list, which is finally exported via the ExportVtkDataFile command
    //! @param rElementId ... Element id
    //! @param rIp ... local ip number
    void AddVisualizationComponentNonlocalWeights(int rElementId, int rIp);

    //! @brief ... clear all visualization components
    void ClearVisualizationComponents();

    //! @brief ... export the entire structure to Vtk data file
    //! @param rFileName ... file name
    //! @param rWhat ... string which describes what to plot
    void ExportVtkDataFile(const std::string& rFileName) const;

    //Visualize for all integration points the fine scale structure
    void ElementGroupVisualizeIpMultiscale(int rGroupIdent, const std::string& rFileName, bool rVisualizeDamage)const;

    //Visualize for all integration points the fine scale structure
    void ElementGroupVisualizeIpMultiscaleDamage(int rGroupIdent, const std::string& rFileName)const;

    //Visualize for all integration points the fine scale structure
    void ElementGroupVisualizeIpMultiscaleHomogeneous(int rGroupIdent, const std::string& rFileName)const;

    //! @brief ... export an element group to Vtk data file
    //! @param rGroupIdent ... group ident
    //! @param rFileName ... file name
    //! @param rWhat ... string which describes what to plot
    void ElementGroupExportVtkDataFile(int rGroupIdent, const std::string& rFileName) const;

#ifndef SWIG
    //! @brief ... define the data sets (scalar, vector etc for the visualize routine based on the mVisualizecomponents
    void DefineVisualizeData(VisualizeUnstructuredGrid& rVisualize, const boost::ptr_list<NuTo::VisualizeComponentBase>& rWhat)const;

    //! @brief ... adds all the elements in the vector to the data structure that is finally visualized
    void ElementVectorAddToVisualize(VisualizeUnstructuredGrid& rVisualize, const boost::ptr_list<NuTo::VisualizeComponentBase>& rWhat, const std::vector<const ElementBase*>& rElements) const;

    //! @brief ... adds all the elements in the vector to the data structure that is finally visualized
    void ElementTotalAddToVisualize(VisualizeUnstructuredGrid& rVisualize, const boost::ptr_list<NuTo::VisualizeComponentBase>& rWhat) const;

    //! @brief ... adds all the elements in a group to the data structure that is finally visualized
    void ElementGroupAddToVisualize(int rGroupId, VisualizeUnstructuredGrid& rVisualize, const boost::ptr_list<NuTo::VisualizeComponentBase>& rWhat) const;

    //! @brief ... adds all the elements in the vector to the data structure that is finally visualized
    void ElementVectorAddToVisualizeIpMultiscale(VisualizeUnstructuredGrid& rVisualize, const boost::ptr_list<NuTo::VisualizeComponentBase>& rWhat,
    		const std::vector<const ElementBase*>& rElements, bool rVisualizeDamage) const;

    //! @brief ... adds all the elements in the vector to the data structure that is finally visualized
    void ElementTotalAddToVisualizeIpMultiscale(VisualizeUnstructuredGrid& rVisualize,
    		const boost::ptr_list<NuTo::VisualizeComponentBase>& rWhat, bool rVisualizeDamage) const;

    //! @brief ... adds all the elements in a group to the data structure that is finally visualized
    void ElementGroupAddToVisualizeIpMultiscale(int rGroupId, VisualizeUnstructuredGrid& rVisualize,
    		const boost::ptr_list<NuTo::VisualizeComponentBase>& rWhat, bool rVisualizeDamage) const;
#endif //SWIG
#endif // ENABLE_VISUALIZE

    //@brief determines the maximum independent sets and stores it at the structure
    // is only relevant for openmp, otherwise the routine is just empty
    void CalculateMaximumIndependentSets();

    //! @brief ... build global coefficient matrix (e.g stiffness) for primary dofs (e.g displacements, rotations, temperature)
    //! @param rMatrix ... global coefficient matrix (nonsymmetric)
    //! @param rVector ... global equivalent load vector (e.g. due to prescribed displacements)
    void BuildGlobalCoefficientMatrix0(NuTo::SparseMatrixCSRGeneral<double>& rMatrix, NuTo::FullMatrix<double>& rVector);

    //! @brief ... build global coefficient matrix (e.g stiffness) for primary dofs (e.g displacements, rotations, temperature)
    //! @param rMatrix ... global coefficient matrix (symmetric)
    //! @param rVector ... global equivalent load vector (e.g. due to prescribed displacements)
    void BuildGlobalCoefficientMatrix0(NuTo::SparseMatrixCSRSymmetric<double>& rMatrix, NuTo::FullMatrix<double>& rVector);

    //! @brief ... build global coefficient matrix (e.g stiffness) for primary dofs (e.g displacements, rotations, temperature)
    //! @param rMatrix ... global coefficient matrix (nonsymmetric)
    //! @param rVector ... global equivalent load vector (e.g. due to prescribed displacements)
    void BuildGlobalCoefficientMatrix0(NuTo::SparseMatrixCSRVector2General<double>& rMatrix, NuTo::FullMatrix<double>& rVector);

    //! @brief ... build global coefficient matrix (e.g stiffness) for primary dofs (e.g displacements, rotations, temperature)
    //! @param rMatrix ... global coefficient matrix (nonsymmetric)
    //! @param rVector ... global equivalent load vector (e.g. due to prescribed displacements)
    void BuildGlobalCoefficientMatrix0(SparseMatrixCSRVector2Symmetric<double>& rMatrix, FullMatrix<double>& rVector);

    //! @brief ... build global external load vector
    //! @param rVector ... external load vector
    void BuildGlobalExternalLoadVector(NuTo::FullMatrix<double>& rVector);

    //! @brief ... build global gradient of the internal potential (e.g. the internal forces)
    //! @param rVector ... global gradient of the internal potential (e.g. internal force vector)
    void BuildGlobalGradientInternalPotentialVector(NuTo::FullMatrix<double>& rVector);

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

    //! @brief ... store all elements connected to this node in a vector
    //! @param rNodeId (Input) 			... node id
    //! @param rElementNumbers (Output) ... vector of element ids
    void NodeGetElements(const int rNodeId, NuTo::FullMatrix<int>& rElementNumbers);

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
    void NodeSetDisplacements(int rId,const NuTo::FullMatrix<double>& rDisplacements);

    //! @brief sets the displacements of a group of nodes
    //! @param rIdent node group identifier
    //! @param rDisplacements matrix (one column) with the displacements
    void NodeGroupSetDisplacements(int rIdent, const FullMatrix<double>& rDisplacements);

    //! @brief gets the coordinates of a node
    //! @param rNode node identifier
    //! @param rCoordinates matrix (one column) with the coordinates
    void NodeGetCoordinates(int rNode, NuTo::FullMatrix<double>& rCoordinates)const;

    //! @brief gets the coordinates of a group of nodes (be careful, the order of the nodes in a group might change between different runs)
    //! @param rNodeGroup node group identifier
    //! @param rCoordinates matrix (rows/nodes columns/coordinates)
    void NodeGroupGetCoordinates(int rNodeGroup, NuTo::FullMatrix<double>& rCoordinates);

    //! @brief gets the displacements of a node
    //! @param rNode node identifier
    //! @param rDisplacements matrix (one column) with the displacements
    void NodeGetDisplacements(int rNode, NuTo::FullMatrix<double>& rDisplacements)const;

    //! @brief gets the displacements of a group of nodes (be careful, the order of the nodes in a group might change between different runs)
    //! @param rNodeGroup node group identifier
    //! @param rDisplacements matrix (rows/nodes columns/rDisplacements)
    void NodeGroupGetDisplacements(int rNodeGroup, NuTo::FullMatrix<double>& rDisplacements);

    //! @brief extract dof values (e.g. displacements, temperatures to the nodes)
    //! @param rActiveDofValues ... vector of global active dof values (ordering according to global dofs, size is number of active dofs)
    //! @param rDependentDofValues ... vector of global dependent dof values (ordering according to (global dofs) - (number of active dofs), size is (total number of dofs) - (number of active dofs))
    virtual void NodeExtractDofValues(NuTo::FullMatrix<double>& rActiveDofValues, NuTo::FullMatrix<double>& rDependentDofValues) const = 0;

    //! @brief write dof values (e.g. displacements, temperatures to the nodes)
    //! @param rActiveDofValues ... vector of global dof values (ordering according to global dofs, size is number of active dofs)
    virtual void NodeMergeActiveDofValues(const NuTo::FullMatrix<double>& rActiveDofValues)=0;

    //! @brief calculate the internal force vector for a node
    //! @param rId ... node id
    //! @param rGradientInternalPotential ...vector for all the dofs the corresponding internal force (return value)
    void NodeInternalForce(int rId, NuTo::FullMatrix<double>& rNodeForce) const;

    //! @brief calculate the internal force vector for a node group of nodes
    //! @param rIdent ... group id
    //! @param rGradientInternalPotential ...vector for all the dofs the corresponding internal force (return value)
    void NodeGroupInternalForce(int rIdent, NuTo::FullMatrix<double>& rNodeForce) const;

#ifndef SWIG
    //! @brief calculate the internal force vector for a node
    //! @param rNodePtr  node for which this has to be calculated
    //! @param rGradientInternalPotential ...vector for all the dofs the corresponding internal force (return value)
    void NodeInternalForce(const NodeBase* rNodePtr, NuTo::FullMatrix<double>& rNodeForce) const;
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
    void ElementStiffness(int rElementId, NuTo::FullMatrix<double>& rResult ,
                          NuTo::FullMatrix<int>& rGlobalDofsRow,
                          NuTo::FullMatrix<int>& rGlobalDofsColumn)const;

    //! @brief calculates the coefficient matrix for the 0-th derivative in the differential equation
    //! for a mechanical problem, this corresponds to the stiffness matrix
    void ElementCoefficientMatrix_0(int rElementId,
                                    NuTo::FullMatrix<double>& rResult,
                                    NuTo::FullMatrix<int>& rGlobalDofsRow,
                                    NuTo::FullMatrix<int>& rGlobalDofsColumn)const;

    //! @brief calculates the coefficient matrix for the 0-th derivative in the differential equation
    //! and compares it to the matrix using central differences
    //! for a mechanical problem, this corresponds to the stiffness matrix
    //! @param rElementId element
    //! @param rDelta  delta step for finite differences
    //! @return maximum difference between analytical and central difference method
    double ElementCoefficientMatrix_0_Check(int rElementId, double rDelta, NuTo::FullMatrix<double>& rDifference);

#ifndef SWIG
    //! @brief calculates the coefficient matrix for the 0-th derivative in the differential equation using central differences
    //! @param rElementPtr element
    //! @param rDelta  delta step for finite differences
    //! @param stiffnessCDF  stiffness from central differences (return value, size should be allocated correctly before entering the routine)
    //! @return maximum difference between analytical and central difference method
    void ElementCoefficientMatrix_0_Resforce(ElementBase* rElementPtr, double rDelta, NuTo::FullMatrix<double>& stiffnessCDF);
#endif //SWIG

    //! @brief calculates the coefficient matrix for the 0-th derivative in the differential equation
    //! @param rElementId elementId
    //! @param rDelta  delta step for finite differences
    //! @param stiffnessCDF  stiffness from central differences (return value, size should be allocated correctly before entering the routine)
    //! @return maximum difference between analytical and central difference method
    void ElementCoefficientMatrix_0_Resforce(int rElementId, double rDelta, NuTo::FullMatrix<double>& stiffnessCDF);

    //! @brief calculates the coefficient matrix for the 0-th derivative in the differential equation
    //! and compares it to the matrix using central differences
    //! for a mechanical problem, this corresponds to the stiffness matrix
    //! @param rDelta  delta step for finite differences
    //! @return element with maximum error
    int ElementTotalCoefficientMatrix_0_Check(double rDelta, NuTo::FullMatrix<double>& rDifference);

    //! @brief similiar to above, but this time the global matrix is checked, not the element matrices
    //! @return false, if stiffness is not correct
    virtual bool CheckStiffness();

    //! @brief calculates the coefficient matrix for the 1-th derivative in the differential equation
    //! for a mechanical problem, this corresponds to the damping matrix
    void ElementCoefficientMatrix_1(int rElementId,
                                    NuTo::FullMatrix<double>& rResult,
                                    NuTo::FullMatrix<int>& rGlobalDofsRow,
                                    NuTo::FullMatrix<int>& rGlobalDofsColumn)const;

    //! @brief calculates the coefficient matrix for the 2-th derivative in the differential equation
    //! for a mechanical problem, this corresponds to the Mass matrix
    void ElementCoefficientMatrix_2(int rElementId,
                                    NuTo::FullMatrix<double>& rResult,
                                    NuTo::FullMatrix<int>& rGlobalDofsRow,
                                    NuTo::FullMatrix<int>& rGlobalDofsColumn)const;

    //! @brief calculates the gradient of the internal potential
    //! for a mechanical problem, this corresponds to the internal force vector
    void ElementGradientInternalPotential(int rElementId,
                                          NuTo::FullMatrix<double>& rResult,
                                          NuTo::FullMatrix<int>& rGlobalDofsRow)const;

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
#endif //SWIG

     //! @brief modifies the fine scale model at an ip for a multiscale approach
     //! @param rElementIdent identifier for the element
     //! @param rIp integration point
     //! @param rFileName binary file to be deserialize the structure from
     void ElementIpSetFineScaleModel(int rElementId, int rIp, std::string rFileName);

     //! @brief modifies the fine scale model for a group of elements for a multiscale approach
     //! @param rGroupIdent identifier for the group of elements
     //! @param rFileName binary file to be deserialize the structure from
     void ElementGroupSetFineScaleModel(int rGroupIdent, std::string rFileName);

     //  @brief modifies the fine scale model for all element ips for a multiscale approach
     //! @param rFileName binary file to be deserialize the structure from
     void ElementTotalSetFineScaleModel(std::string rFileName);

 #ifndef SWIG
     //! @brief modifies the constitutive law of a single element
     //! @param rElement element pointer
     //! @param rIp integration point number
     //! @param rFileName binary file to be deserialize the structure from
     void ElementIpSetFineScaleModel(ElementBase* rElement, int rIp, std::string rFileName, double rLengthCoarseScale, std::string rIPName);
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

    //! @brief modifies the constitutive law of a single element
    //! @param rElement element pointer
    //! @param rConstitutive material pointer
    void ElementSetIntegrationType(ElementBase* rElement, const IntegrationTypeBase* rIntegrationType, NuTo::IpData::eIpDataType rIpDataType);
#endif //SWIG

    //! @brief modifies the section of a single element
    //! @param rElementIdent element number
    //! @param rSectionIdent identifier for the section
    void ElementSetIntegrationType(int rElementId, const std::string& rIntegrationTypeIdent, std::string rIpDataTypeStr);

    //! @brief modifies the section of a group of elements
    //! @param rGroupIdent identifier for the group of elements
    //! @param rSectionId identifier for the section
    void ElementGroupSetIntegrationType(int rGroupIdent, const std::string& rIntegrationTypeIdent, std::string rIpDataTypeStr);

    //! @brief modifies the section of a all elements
    //! @param rSectionIdent identifier for the section
    void ElementTotalSetIntegrationType(const std::string& rIntegrationTypeIdent, std::string rIpDataTypeStr);


#ifndef SWIG
    //! @brief modifies the section of a single element
    //! @param rElement element pointer
    //! @param rConstitutive section
    void ElementSetSection(ElementBase* rElement, SectionBase* rSection);
#endif //SWIG

    //! @brief calculates the engineering strain
    //! @param rElemIdent  element number
    //! @param rEngineerungStrain engineering strain (return value, always 6xnumIp matrix)
    void ElementGetEngineeringStrain(int rElementId, NuTo::FullMatrix<double>& rEngineeringStrain)const;

    //! @brief calculates the engineering plastic strain
    //! @param rElemIdent  element number
    //! @param rEngineerungStrain engineering plastic strain (return value, always 6xnumIp matrix)
    void ElementGetEngineeringPlasticStrain(int rElementId, NuTo::FullMatrix<double>& rEngineeringPlasticStrain)const;

    //! @brief calculates the engineering stress
    //! @param rElemIdent  element number
    //! @param rEingineeringStress Engineering Stress (return value, always 6xnumIp matrix)
    void ElementGetEngineeringStress(int rElementId, NuTo::FullMatrix<double>& rEngineeringStress)const;

    //! @brief calculates the damage
    //! @param rElemIdent  identifier for the element
    //! @param rDamage (return value, always 1xnumIp matrix)
    void ElementGetDamage(int rElementId, FullMatrix<double>& rDamage)const;

    //! @brief calculates the maximum damage in all elements
    //! @param rElemIdent  identifier for the element
    //! @return max damage value
    double ElementTotalGetMaxDamage()const;

    //! @brief updates the history data of a all elements
    void ElementTotalUpdateStaticData();

    //! @brief updates the temprory static data of a all elements
    //! its is a const function, since only mutuable data (instead of const) is updated (kind of temporary data)
    void ElementTotalUpdateTmpStaticData();

    //! @brief calculates the average stress
    //! @param rVolume  volume of the structure in 3D /area in 2D/ length in 1D
    //! this is a parameter of the model, since holes have to be considered (zero stress, but still nonzero area)
    //! @param rEngineeringStress  average stress (return value)
    void ElementTotalGetAverageStress(double rVolume, NuTo::FullMatrix<double>& rEngineeringStress)const;

    //! @brief calculates the average stress
    //! @param rGroupId  group number
    //! @param rVolume  volume of the structure in 3D /area in 2D/ length in 1D
    //! this is a parameter of the model, since holes have to be considered (zero stress, but still nonzero area)
    //! @param rEngineeringStress  average stress (return value)
    void ElementGroupGetAverageStress(int rGroupId, double rVolume, NuTo::FullMatrix<double>& rEngineeringStress)const;

    //! @brief calculates the average strain
    //! @param rVolume  volume of the structure in 3D /area in 2D/ length in 1D
    //! this is a parameter of the model, since holes have to be considered (zero stress, but still nonzero area)
    //! @param rEngineeringStrain  average strain (return value)
    void ElementTotalGetAverageStrain(double rVolume, NuTo::FullMatrix<double>& rEngineeringStrain)const;

    //! @brief calculates the average strain
    //! @param rGroupId  group number
    //! @param rVolume  volume of the structure in 3D /area in 2D/ length in 1D
    //! this is a parameter of the model, since holes have to be considered (zero stress, but still nonzero area)
    //! @param rEngineeringStrain  average strain (return value)
    void ElementGroupGetAverageStrain(int rGroupId, double rVolume, NuTo::FullMatrix<double>& rEngineeringStrain)const;

    //! @brief calculates the total energy of the system
    //! @return total energy
    virtual double ElementTotalGetTotalEnergy()const;

    //! @brief calculates the total energy of the system
    //! @return total energy
    virtual double ElementGroupGetTotalEnergy(int rGroupId)const;

    //! @brief calculates the elastic energy of the system
    //! @return elastic energy
    double ElementTotalGetElasticEnergy()const;

    //! @brief sets the parameters  of the finescale model (structure ip) for a group of elements
    //! @parameter rElement Element pointer
    //! @parameter rName name of the parameter, e.g. YoungsModulus
    //! @parameter rParameter value of the parameter
    void ElementGroupSetFineScaleParameter(int rGroupId, std::string rName, double rParameter);

    //! @brief sets the parameters  of the finescale model (structure ip) for a group of elements
    //! @parameter rElement Element pointer
    //! @parameter rName name of the parameter, e.g. YoungsModulus
    //! @parameter rParameter value of the parameter
    void ElementGroupSetFineScaleParameter(int rGroupId, std::string rName, std::string rParameter);

    //! @brief sets the parameters  of the finescale model (structure ip) for all elements
    //! @parameter rElement Element pointer
    //! @parameter rName name of the parameter, e.g. YoungsModulus
    //! @parameter rParameter value of the parameter
    void ElementTotalSetFineScaleParameter(std::string rName, double rParameter);

    //! @brief sets the parameters  of the finescale model (structure ip) for all elements
    //! @parameter rElement Element pointer
    //! @parameter rName name of the parameter, e.g. YoungsModulus
    //! @parameter rParameter value of the parameter
    void ElementTotalSetFineScaleParameter(std::string rName, std::string rParameter);

#ifndef SWIG
    //! @brief sets the parameters  of the finescale model (structure ip)
    //! @parameter rElement Element pointer
    //! @parameter rName name of the parameter, e.g. YoungsModulus
    //! @parameter rParameter value of the parameter
    void ElementSetFineScaleParameter(ElementBase* rElement, std::string rName, double rParameter);

    //! @brief sets the parameters  of the finescale model (structure ip)
    //! @parameter rElement Element pointer
    //! @parameter rName name of the parameter, e.g. YoungsModulus
    //! @parameter rParameter value of the parameter
    void ElementSetFineScaleParameter(ElementBase* rElement, std::string rName, std::string rParameter);

    //! @brief modifies the finescale model of a multiscale element from the initial elastic solution phase to the nonlinear phase if principal stress is larger than tensile strength
    //! @return true, if an adaptation has been performed, otherwise false
    void ElementTotalMultiscaleSwitchToNonlinear();
#endif

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
    void ConstraintLagrangeGetMultiplier(int ConstraintId, NuTo::FullMatrix<double>& rMultiplier)const;

    //! @brief sets the penalty stiffness of the augmented Lagragian to the prescribed value
    //! @param ConstraintId constraint id
    //! @param rPenalty penalty parameter
    void ConstraintLagrangeSetPenaltyStiffness(int ConstraintId, double rPenalty);

    //! @brief adds a displacement constraint equation for a node group solved using Lagrange multiplier
    //! @param rGroupId group id
    //! @param rDirection direction of the constraint (in 2D a point with 2 entries, in 3D 3 entries, in 1D not used)
    //! @param rValue prescribed value (e.g. zero to fix a displacement to zero)
    //! @return integer id to delete or modify the constraint
    int ConstraintLagrangeSetDisplacementNodeGroup(int rGroupId, const NuTo::FullMatrix<double>& rDirection, const std::string& rSign, double rValue);

#ifndef SWIG
    //! @brief adds a displacement constraint equation for a node group solved using Lagrange multiplier
    //! @param rGroup group pointer
    //! @param rDirection direction of the constraint (in 2D a point with 2 entries, in 3D 3 entries, in 1D not used)
    //! @param rValue prescribed value (e.g. zero to fix a displacement to zero)
    //! @return integer id to delete or modify the constraint
    int ConstraintLagrangeSetDisplacementNodeGroup(Group<NodeBase>* rGroup, const NuTo::FullMatrix<double>& rDirection, NuTo::Constraint::eEquationSign rEquationSign, double rValue);

    //! @brief adds a displacement constraint equation for a node
    //! @param rNode pointer to node
    //! @param rDirection direction of the constraint (in 2D a point with 2 entries, in 3D 3 entries, in 1D not used)
    //! @param rValue prescribed value (e.g. zero to fix a displacement to zero)
    //! @return integer id to delete or modify the constraint
    int ConstraintLinearSetDisplacementNode(NodeBase* rNode, const NuTo::FullMatrix<double>& rDirection, double rValue);
#endif

    //! @brief adds a displacement constraint equation for a node
    //! @param rNode identifier for node
    //! @param rComponent e.g. the first (count from zero) displacement component
    //! @param rValue prescribed value (e.g. zero to fix a displacement to zero)
    int  ConstraintLinearSetDisplacementNode(int rIdent, const NuTo::FullMatrix<double>& rDirection, double rValue);

    #ifndef SWIG
    //! @brief adds a fine scale displacement constraint equation for a node
    //! @param rNode pointer to node
    //! @param rDirection direction of the constraint (in 2D a point with 2 entries, in 3D 3 entries, in 1D not used)
    //! @param rValue prescribed value (e.g. zero to fix a displacement to zero)
    //! @return integer id to delete or modify the constraint
    int ConstraintLinearSetFineScaleDisplacementNode(NodeBase* rNode, const NuTo::FullMatrix<double>& rDirection, double rValue);
#endif

    //! @brief adds a displacement constraint equation for a node
    //! @param rNode identifier for node
    //! @param rComponent e.g. the first (count from zero) displacement component
    //! @param rValue prescribed value (e.g. zero to fix a displacement to zero)
    int  ConstraintLinearSetFineScaleDisplacementNode(int rIdent, const NuTo::FullMatrix<double>& rDirection, double rValue);

#ifndef SWIG
    //! @brief adds a displacement constraint equation for a group of node
    //! @param rNode pointer to group of nodes
    //! @param rDirection direction of the constraint (in 2D a point with 2 entries, in 3D 3 entries, in 1D not used)
    //! @param rValue prescribed value (e.g. zero to fix a displacement to zero)
    //! @return integer id to delete or modify the constraint
    int ConstraintLinearSetDisplacementNodeGroup(Group<NodeBase>* rGroup, const NuTo::FullMatrix<double>& rDirection, double rValue);
#endif

    //! @brief adds a constraint equation for a group of nodes
    //! @param rGroupIdent identifier for group of nodes
    //! @param rAttribute displacements, rotations, temperatures
    //! @param rComponent e.g. the first (count from zero) displacement component
    //! @param rValue prescribed value (e.g. zero to fix a displacement to zero)
    int ConstraintLinearSetDisplacementNodeGroup(int rGroupIdent, const NuTo::FullMatrix<double>& rDirection, double rValue);

    #ifndef SWIG
    //! @brief adds a fine scale displacement constraint equation for a group of node
    //! @param rNode pointer to group of nodes
    //! @param rDirection direction of the constraint (in 2D a point with 2 entries, in 3D 3 entries, in 1D not used)
    //! @param rValue prescribed value (e.g. zero to fix a displacement to zero)
    //! @return integer id to delete or modify the constraint
    int ConstraintLinearSetFineScaleDisplacementNodeGroup(Group<NodeBase>* rGroup, const NuTo::FullMatrix<double>& rDirection, double rValue);
#endif

    //! @brief adds a constraint equation for a group of nodes
    //! @param rGroupIdent identifier for group of nodes
    int ConstraintLinearSetFineScaleDisplacementNodeGroup(int rGroupIdent, const NuTo::FullMatrix<double>& rDirection, double rValue);

    //! @brief returns the number of constraint equations
    //! @return number of constraints
    int ConstraintGetNumLinearConstraints()const;

    //! @brief calculates the constraint matrix that builds relations between the nodal dagrees of freedom
    //! rConstraintMatrix*DOFS = RHS
    //! @param rConstraintMatrix constraint matrix
    //! @param rRHS right hand side
    void ConstraintGetConstraintMatrix(NuTo::SparseMatrixCSRGeneral<double>& rConstraintMatrix, NuTo::FullMatrix<double>& rRHS);

    //!@brief sets/modifies the right hand side of the constraint equations
    //!@param rRHS new right hand side
    //!@param rRHS new right hand side
    void ConstraintSetRHS(int rConstraintEquation, double rRHS);

    //!@brief sets/modifies the strain of a constraint equation (works only for periodic bc)
    //!@param rConstraintEquation id of the constraint equation
    //!@param rStrain new strain
    void ConstraintPeriodicSetStrain(int rConstraintEquation, NuTo::FullMatrix<double> rStrain);

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
    void ConstraintExtractGlobalDofValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues)const;

    //! @brief write dof values to the Lagrange multipliers (based on global dof number)
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    void ConstraintMergeGlobalDofValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& dependentDofValues);

    //!@brief sets/modifies the strain of a constraint equation (works only for periodic bc)
    //!@param rConstraintEquation id of the constraint equation
    //!@param rStrain new strain
    void ConstraintPeriodicSetStrain2D(int rConstraintEquation, const NuTo::EngineeringStrain2D& rStrain);
#endif

    //!@brief sets/modifies the crack opening of a constraint equation (works only for periodic bc)
    //!@param rConstraintEquation id of the constraint equation
    //!@param rCrackOpening new crack opening (x,y)
    void ConstraintPeriodicSetCrackOpening(int rConstraintEquation,
            NuTo::FullMatrix<double> rCrackOpeningdouble);

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
    int ConstraintLinearDisplacementsSetPeriodic2D(double angle, NuTo::FullMatrix<double> rStrain,
            NuTo::FullMatrix<double> rCrackOpening, double rRadiusToCrackWithoutConstraints,
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
    void ConstraintsBuildGlobalCoefficientSubMatrices0General(SparseMatrix<double>& rMatrixJJ, SparseMatrix<double>& rMatrixJK)const;

    //! @brief ... add the contribution of Lagrange multipliers to the global system of equations
    //! @param rMatrixJJ ... matrix jj
    //! @param rMatrixJK ... matrix jk
    //! @param rMatrixKJ ... matrix kj
    //! @param rMatrixKK ... matrix kk
    void ConstraintBuildGlobalCoefficientSubMatrices0General(SparseMatrix<double>& rMatrixJJ, SparseMatrix<double>& rMatrixJK, SparseMatrix<double>& rMatrixKJ, SparseMatrix<double>& rMatrixKK)const;

    //! @brief ... add the contribution of Lagrange multipliers to the global system of equations
    //! @param rMatrixJJ ... matrix jj
    //! @param rMatrixJK ... matrix jk
    void ConstraintBuildGlobalCoefficientSubMatrices0Symmetric(SparseMatrix<double>& rMatrixJJ, SparseMatrix<double>& rMatrixJK)const;

    //! @brief ... add the contribution of Lagrange multipliers to the global system of equations
    //! @param rMatrixJJ ... matrix jj
    //! @param rMatrixJK ... matrix jk
    //! @param rMatrixKK ... matrix kk
    void ConstraintBuildGlobalCoefficientSubMatrices0Symmetric(SparseMatrix<double>& rMatrixJJ, SparseMatrix<double>& rMatrixJK, SparseMatrix<double>& rMatrixKK) const;

    //! @brief ... add the contribution of Lagrange multipliers to the global system of equations
    //! @param rActiveDofGradientVector ... gradient of active dofs
    //! @param rDependentDofGradientVector ... gradient of dependent dofs
    void ConstraintBuildGlobalGradientInternalPotentialSubVectors(NuTo::FullMatrix<double>& rActiveDofGradientVector, NuTo::FullMatrix<double>& rDependentDofGradientVector) const;

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
    int LoadCreateNodeForce(int rNodeIdent, const NuTo::FullMatrix<double>& rDirection, double rValue);

    //! @brief adds a force for a node group
    //! @param rGroupIdent ... identifier for node group
    //! @param rDirection ... direction of the force
    //! @param rValue ... force
    //! @return integer id to delete or modify the load
    int LoadCreateNodeGroupForce(int rGroupIdent, const NuTo::FullMatrix<double>& rDirection, double rValue);

    //! @brief delete load
    //! @param rIdent ... load identifier
    void LoadDelete(int rIdent);

#ifndef SWIG
    //! @brief adds a force for a node
    //! @param rNode ... pointer to node
    //! @param rDirection ... direction of the force
    //! @param rValue ... force
    //! @return integer id to delete or modify the load
    int LoadCreateNodeForce(const NodeBase* rNode, const NuTo::FullMatrix<double>& rDirection, double rValue);

    //! @brief adds a force for a node grpup
    //! @param rNodeGroup ... pointer to node group
    //! @param rDirection ... direction of the force
    //! @param rValue ... force
    //! @return integer id to delete or modify the load
    int LoadCreateNodeGroupForce(const Group<NodeBase>* rNodeGroup, const NuTo::FullMatrix<double>& rDirection, double rValue);

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
    NuTo::FullMatrix<double> ConstitutiveLawGetYieldStrength(int rIdent) const;

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

    //! @brief ... get hardening modulus for multilinear response
    //! @return ... first column: equivalent plastic strain
    //! @return ... second column: corresponding hardening modulus
    NuTo::FullMatrix<double> ConstitutiveLawGetHardeningModulus(int rIdent) const;

    //! @brief ... add hardening modulus
    //! @param rEpsilon ...  equivalent plastic strain
    //! @param rSigma ...  hardening modulus
    void ConstitutiveLawAddHardeningModulus(int rIdent, double rEpsilon, double rH);

    //! @brief ... get nonlocal radius
    //! @return ... nonlocal radius
    double ConstitutiveLawGetNonlocalRadius(int rIdent) const;

    //! @brief ... set nonlocal radius
    //! @param rH ...  nonlocal radius
    void ConstitutiveLawSetNonlocalRadius(int rIdent, double rRadius);

    //! @brief ... get tensile strength
    //! @param rTensileStrength ...  tensile strength
    double ConstitutiveLawGetTensileStrength(int rIdent);

    //! @brief ... set tensile strength
    //! @param rTensileStrength ...  tensile strength
    void ConstitutiveLawSetTensileStrength(int rIdent, double rTensileStrength);

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

    //! @brief ... get elastic stiffness
    //! @return ...  elastic stiffness
    NuTo::FullMatrix<double> ConstitutiveLawGetElasticStiffness(int rIdent);

    //! @brief ... set fracture energy
    //! @param rElasticStiffness ...  fracture energy
    void ConstitutiveLawSetElasticStiffness(int rIdent, NuTo::FullMatrix<double> rElasticStiffness);

    //! @brief ... get elastic stiffness
    //! @param rIdent ...  fracture energy
    std::string ConstitutiveLawGetMultiscaleFile(int rIdent);

    //! @brief ... set fracture energy
    //! @param rFractureEnergy ...  file name
    void ConstitutiveLawSetMultiscaleFile(int rIdent, std::string rFileName);

    //! @brief ... get crack transition radius
    //! @param rIdent ...  identifier
    double ConstitutiveLawGetCrackTransitionRadius(int rIdent);

    //! @brief ... set crack transition radius
    //! @param rCrackTransitionRadius ...  fracture energy
    void ConstitutiveLawSetCrackTransitionRadius(int rIdent, double rCrackTransitionRadius);

    //! @brief ... get PenaltyStiffnessCrackAngle
    //! @param rIdent ...  identifier
    double ConstitutiveLawGetPenaltyStiffnessCrackAngle(int rIdent);

    //! @brief ... set PenaltyStiffnessCrackAngle
    //! @param rPenaltyStiffnessCrackAngle ...  PenaltyStiffnessCrackAngle
    void ConstitutiveLawSetPenaltyStiffnessCrackAngle(int rIdent, double rPenaltyStiffnessCrackAngle);

    //! @brief ... get scaling factor for the crack angle
    //! @param rIdent ...  identifier
    double ConstitutiveLawGetScalingFactorCrackAngle(int rIdent);

    //! @brief ... set scaling factor for the crack angle
    //! @param rScalingFactorCrackAngle ...  scaling factor
    void ConstitutiveLawSetScalingFactorCrackAngle(int rIdent, double rScalingFactorCrackAngle);

    //! @brief ... get scaling factor for the crack opening
    //! @param rIdent ...  identifier
    double ConstitutiveLawGetScalingFactorCrackOpening(int rIdent);

    //! @brief ... set scaling factor for the crack opening
    //! @param rScalingFactorCrackAngle ...  scaling factor
    void ConstitutiveLawSetScalingFactorCrackOpening(int rIdent, double rScalingFactorCrackOpening);

    //! @brief ... get scaling factor for the total strain
    //! @param rIdent ...  identifier
    double ConstitutiveLawGetScalingFactorEpsilon(int rIdent);

    //! @brief ... set PenaltyStiffnessCrackAngle
    //! @param rScalingFactorCrackAngle ...  scaling factor
    void ConstitutiveLawSetScalingFactorEpsilon(int rIdent, double rScalingFactorEpsilon);

    //! @brief ... get result directory for fine scale models in multiscale simulation
    //! @param rIdent ...  identifier
    std::string ConstitutiveLawGetResultDirectory(int rIdent);

    //! @brief ... set ResultDirectory for fine scale models in multiscale simulation
    //! @param rResultDirectory ...  ResultDirectory
    void ConstitutiveLawSetResultDirectory(int rIdent, std::string rResultDirectory);

    //! @brief ... get load step macro for fine scale models in multiscale simulation
    //! @param rIdent ...  identifier
    int ConstitutiveLawGetLoadStepMacro(int rIdent);

    //! @brief ... set LoadStepMacro for fine scale models in multiscale simulation
    //! @param rLoadStepMacro ...  LoadStepMacro
    void ConstitutiveLawSetLoadStepMacro(int rIdent, int rLoadStepMacro);

    //! @brief ... get if the fine scale model is to be used with the linear elastic periodic boundary shape functions
    //! @return rUseAdditionalPeriodicShapeFunctions
    bool ConstitutiveLawGetUseAdditionalPeriodicShapeFunctions(int rIdent)const;

    //! @brief ... set if the fine scale model is to be used with the linear elastic periodic boundary shape functions
    //! @param rUseAdditionalPeriodicShapeFunctions ...  true or false
    void ConstitutiveLawSetUseAdditionalPeriodicShapeFunctions(int rIdent, bool rUseAdditionalPeriodicShapeFunctions);

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
#endif

    //! @brief ... Creates a group for the structure
    //! @param ... rIdent identifier for the group
    //! @param ... rType  type of the group, e.g. "NODES" or "ELEMENTS"
    int GroupCreate(const std::string& rType);

    //! @brief ... Creates a group for the structure
    //! @param ... rIdent identifier for the group
    //! @param ... rType  type of the group
    void GroupCreate(int id, NuTo::Groups::eGroupId rEnumType);

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

    //! @brief ... Adds all nodes to a group whose coordinates are in the specified range
    //! @param ... rIdentGroup identifier for the group
    //! @param ... rCenter center of the selection circle
    //! @param ... rMin ... minimum radius
    //! @param ... rMax ... maximum radius
    void GroupAddNodeRadiusRange(int rIdentGroup, NuTo::FullMatrix<double> rCenter, double rMin, double rMax);

    //! @brief ... Adds an element to an element group
    //! @param ... rIdentGroup identifier for the group
    //! @param ... rIdentElement  identifier for the element
    void GroupAddElement(int rIdentGroup, int rIdElement);

    //! @brief ... Returns the number of members in a group
    //! @param ... rIdentGroup identifier for the group
    //! @return ... number of members
    int GroupGetNumMembers(int rIdentGroup)const;

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
    //! @brief ... Info routine that prints general information about the object (detail according to verbose level)
    virtual void Info()const;

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

    //! @brief set parameters for the  Newton Raphson iteration
    //! @parameters rToleranceResidualForce  convergence criterion for the norm of the residual force vector
    void SetNewtonRaphsonToleranceResidualForce(double rToleranceResidualForce)
    {
    	if (rToleranceResidualForce<=0)
    	{
    		std::cout << "tolerance residual force " << rToleranceResidualForce << std::endl;
    		throw MechanicsException("[StructureBase::SetNewtonRaphsonToleranceResidualForce] tolerance has to be positive.");
    	}
    	mToleranceResidualForce = rToleranceResidualForce;
    }

    //! @brief set parameters for the  Newton Raphson iteration
    //! @parameters rAutomaticLoadstepControl true, if the step length should be adapted
    void SetNewtonRaphsonAutomaticLoadStepControl(bool rAutomaticLoadstepControl)
    {
    	mAutomaticLoadstepControl = rAutomaticLoadstepControl;
    }

    //! @brief set parameters for the  Newton Raphson iteration
    //! @parameters maximum of the delta_step length
    void SetNewtonRaphsonMaxDeltaLoadFactor(double rMaxDeltaLoadFactor)
    {
    	if (rMaxDeltaLoadFactor<=0 || rMaxDeltaLoadFactor>1)
    		throw MechanicsException("[StructureBase::SetNewtonRaphsonMaxDeltaLoadFactor] factor has to be in the internal (0,1].");
    	mMaxDeltaLoadFactor = rMaxDeltaLoadFactor;
    }

    //! @brief set parameters for the  Newton Raphson iteration
    //! @parameters rMaxNumNewtonIterations maximum number of iterations per Newton step
    void SetNewtonRaphsonMaxNumNewtonIterations(double rMaxNumNewtonIterations)
    {
    	if (rMaxNumNewtonIterations<=0)
    		throw MechanicsException("[StructureBase::SetNewtonRaphsonMaxDeltaLoadFactor] maximum number of newton iterations has to be positive.");
    	mMaxNumNewtonIterations = rMaxNumNewtonIterations;
    }

    //! @brief set parameters for the  Newton Raphson iteration
    //! @parameters rDecreaseFactor decrease the load factor in case of no convergence with the prescribed number of Newton iterations
    void SetNewtonRaphsonDecreaseFactor(double rDecreaseFactor)
    {
    	if (rDecreaseFactor<=0 || rDecreaseFactor>=1)
    		throw MechanicsException("[StructureBase::SetNewtonRaphsonDecreaseFactor] factor has to be in the internal (0,1).");
    	mDecreaseFactor= rDecreaseFactor;
    }

    //! @brief set parameters for the  Newton Raphson iteration
    //! @parameters
    void SetNewtonRaphsonMinNumNewtonIterations(double rMinNumNewtonIterations)
    {
    	if (rMinNumNewtonIterations<=0 )
    		throw MechanicsException("[StructureBase::SetNewtonRaphsonMinNumNewtonIterations] minimum number of newton iterations has to be positive.");
    	mMinNumNewtonIterations = rMinNumNewtonIterations;
    }

    //! @brief set parameters for the  Newton Raphson iteration
    //! @parameters rIncreaseFactor if convergence is achieved in less than rMinNumNewtonIterations, the step length is increased
    void SetNewtonRaphsonIncreaseFactor(double rIncreaseFactor)
    {
    	if (rIncreaseFactor<=1)
    		throw MechanicsException("[StructureBase::SetNewtonRaphsonIncreaseFactor] factor has to be greater than 1.");
    	mIncreaseFactor= rIncreaseFactor;
    }

    //! @brief set parameters for the  Newton Raphson iteration
    //! @parameters rMinDeltaLoadFactor if the load factor is smaller the procedure is assumed to diverge (throwing an exception)
    void SetNewtonRaphsonMinDeltaLoadFactor(double rMinDeltaLoadFactor)
    {
    	if (rMinDeltaLoadFactor<=0 || rMinDeltaLoadFactor>=1)
    		throw MechanicsException("[StructureBase::SetNewtonRaphsonMinDeltaLoadFactor] factor has to be in the interval (0,1).");
    	mMinDeltaLoadFactor= rMinDeltaLoadFactor;
    }

    //! @brief set parameters for the  Newton Raphson iteration
    //! @parameters rMinLineSearchFactor smallest line search factor, if normResidual does not fulfill the linesearch criterion with d(t+1)=d(t)+alpha*delta(t), step size is reduced
    void SetNewtonRaphsonMinLineSearchFactor(double rMinLineSearchFactor)
    {
    	if (rMinLineSearchFactor<=0 || rMinLineSearchFactor>=1)
    		throw MechanicsException("[StructureBase::SetNewtonRaphsonMinLineSearchFactor] factor has to be in the interval (0,1).");
    	mMinLineSearchFactor= rMinLineSearchFactor;
    }

    //! @brief performs a Newton Raphson iteration (displacement and/or load control) no save for updates in between time steps
    void NewtonRaphson();

    //! @brief performs a Newton Raphson iteration (displacement and/or load control)
    //! @parameters rSaveStructureBeforeUpdate if set to true, save the structure (done in a separate routine to be implemented by the user) before an update is performed
    //! @parameters rSaveStringStream stringstream the routine is saved to
    //! @parameters rIsSaved return parameter describing, if the routine actually had to to a substepping and consequently had to store the initial state into rSaveStringStream
    virtual void NewtonRaphson(bool rSaveStructureBeforeUpdate,
            std::stringstream& rSaveStringStream,
            bool& rIsSaved);

    //! @brief performs an adaption of the model
    //! @returns true, if the model has actually be changed, or false if no change has been made
    virtual bool AdaptModel()
    {
    	return false;
    }

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

    //! @brief set the load factor (load or displacement control) overload this function to use Newton Raphson
    //! @param load factor
    virtual void SetLoadFactor(double rLoadFactor);

    //! @brief do a postprocessing step after each converged load step (for Newton Raphson iteration) overload this function to use Newton Raphson
    virtual void PostProcessDataAfterUpdate(int rLoadStep, int rNumNewtonIterations, double rLoadFactor, double rDeltaLoadFactor, double rResidual)const;

    //! @brief do a postprocessing step after each converged load step (for Newton Raphson iteration) overload this function to use Newton Raphson
    virtual void PostProcessDataAfterConvergence(int rLoadStep, int rNumNewtonIterations, double rLoadFactor, double rDeltaLoadFactor, double rResidual)const;

    //! @brief do a postprocessing step after each line search within the load step(for Newton Raphson iteration) overload this function to use Newton Raphson
    virtual void PostProcessDataAfterLineSearch(int rLoadStep, int rNewtonIteration, double rLineSearchFactor, double rLoadFactor, double rResidual, const NuTo::FullMatrix<double>& rResidualVector)const;

    //! @brief do a postprocessing step after each line search within the load step(for Newton Raphson iteration) overload this function to use Newton Raphson
    virtual void PostProcessDataInLineSearch(int rLoadStep, int rNewtonIteration, double rLineSearchFactor, double rLoadFactor, double rResidual, double rPrevResidual)const;

    //! @brief initialize some stuff before a new load step (e.g. output directories for visualization, if required)
    virtual void InitBeforeNewLoadStep(int rLoadStep);

    //! @brief only for debugging, info at some stages of the Newton Raphson iteration
    virtual void NewtonRaphsonInfo(int rVerboseLevel)const
    {}

    //! @brief is only true for structure used as multiscale (structure in a structure)
    virtual bool IsMultiscaleStructure()const
    {
    	return false;
    }

    //! @brief is only true for structure used as multiscale (structure in a structure)
    virtual void ScaleCoordinates(double rCoordinates[3])const
    {
    	throw MechanicsException("[NuTo::StructureBase::ScaleCoordinates] only implemented for multiscale structures.");
    }

    void SetUpdateTmpStaticDataRequired()
    {
    	mUpdateTmpStaticDataRequired = true;
    }


protected:
    int mDimension;

    //! @brief ... map storing the name and the pointer to the constitutive law
    //! @sa ConstitutiveBase
    boost::ptr_map<int,ConstitutiveBase> mConstitutiveLawMap;

    //! @brief ... map storing the constraints
    //! @sa ConstraintBase
    boost::ptr_map<int,ConstraintBase> mConstraintMap;

    //! @brief ... map storing node loads
    //! @sa LoadBase
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

    //! @brief ... a mapping from the enums of the predefined integration types to their corresponding string name
    std::vector<std::string> mMappingIntEnum2String;

    //! @brief ... map storing the components (displacements, strains, nonlocal weights etc) to be included in the output (VTK) file
    boost::ptr_list<NuTo::VisualizeComponentBase> mVisualizeComponents;

    //! @brief ... total number of degrees of freedom of the structure
    int mNumDofs;

    //! @brief ... active dofs
    int mNumActiveDofs;

    //!brief ... renumbering of nodal DOFs required or not
    bool mNodeNumberingRequired;

    //! @brief constraint matrix relating the prescibed nodal unknowns to the free parameters
    SparseMatrixCSRGeneral<double> mConstraintMatrix;

    //! @brief right hand side of the constraint equations
    FullMatrix<double> mConstraintRHS;

    //! @brief is set to true, if at least one constitutive model requires an update of tmpStaticData before stress and stiffness routines are called
    bool mHaveTmpStaticData;

    //! @brief is set to false, if the structure is changed (nodes, elements) or (DOFs at the nodes)
    bool mUpdateTmpStaticDataRequired;

    //! @brief absolute tolerance for entries of the global stiffness matrix (coefficientMatrix0)
    //! values smaller than that one will not be added to the global matrix
    double mToleranceStiffnessEntries;

    //*********************************************
    //parameters of the Newton Raphson iteration
    //! @brief convergence criterion for the force norm (or the maximum value of the residual force, both are checked)
    double mToleranceResidualForce;
    //! @brief use (true) automatic load step control (increase, decrease load step)
    bool mAutomaticLoadstepControl;
    //! @brief maximum delta load factor for automatic load step control
    double mMaxDeltaLoadFactor;
    //! @brief maximum number of Newton iterations before the loadstep is decreases (automatic load control) or error of no convergence
    int mMaxNumNewtonIterations;
    //! @brief decrease factor of automatic load control, if no convergence is achieved
    double mDecreaseFactor;
    //! @brief if number of Newton iterations is smaller than this value, the deltaLoadFactor is multiplied by mIncreaseFactor for the next iteration
    int mMinNumNewtonIterations;
    //! @brief see mMinNumNewtonIterations
    double mIncreaseFactor;
    //! @brief smallest load increment, if deltaLoadFactor<mMinDeltaLoadFactor -> no convergence
    double mMinDeltaLoadFactor;
    //! @brief smallest line search factor, if normResidual does not fulfill the linesearch criterion with d(t+1)=d(t)+alpha*delta(t), step size is reduced
    double mMinLineSearchFactor;
    //*********************************************

#ifdef _OPENMP
    //@brief maximum independent sets used for parallel assembly of the stiffness resforce etc.
    mutable std::vector<std::vector<ElementBase*> > mMIS;
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

    //! @brief ... store all elements of a group in a vector
    //! @param rElementGroup ... element group
    //! @param rElements ... vector of element pointer
    void GetElementsByGroup(const Group<ElementBase>* rElementGroup, std::vector<const ElementBase*>& rElements) const;

    //! @brief ... check for dof numbering and build of tmpStaticData
    void BuildGlobalCoefficientMatrixCheck();

    //! @brief ... based on the global dofs build submatrices of the global coefficent matrix0
    //! @param rMatrixJJ ... submatrix jj (number of active dof x number of active dof)
    //! @param rMatrixJK ... submatrix jk (number of active dof x number of dependent dof)
    virtual void BuildGlobalCoefficientSubMatrices0General(NuTo::SparseMatrix<double>& rMatrixJJ, NuTo::SparseMatrix<double>& rMatrixJK) const = 0;

    //! @brief ... based on the global dofs build submatrices of the global coefficent matrix0
    //! @param rMatrixJJ ... submatrix jj (number of active dof x number of active dof)
    //! @param rMatrixJK ... submatrix jk (number of active dof x number of dependent dof)
    //! @param rMatrixKJ ... submatrix kj (number of dependent dof x number of active dof)
    //! @param rMatrixKK ... submatrix kk (number of dependent dof x number of dependent dof)
    virtual void BuildGlobalCoefficientSubMatrices0General(NuTo::SparseMatrix<double>& rMatrixJJ, NuTo::SparseMatrix<double>& rMatrixJK, NuTo::SparseMatrix<double>& rMatrixKJ, NuTo::SparseMatrix<double>& rMatrixKK) const = 0;

    //! @brief ... based on the global dofs build submatrices of the global coefficent matrix0
    //! @param rMatrixJJ ... submatrix jj (number of active dof x number of active dof)
    //! @param rMatrixJK ... submatrix jk (number of active dof x number of dependent dof)
    virtual void BuildGlobalCoefficientSubMatrices0Symmetric(NuTo::SparseMatrix<double>& rMatrixJJ, NuTo::SparseMatrix<double>& rMatrixJK) const = 0;

    //! @brief ... based on the global dofs build submatrices of the global coefficent matrix0
    //! @param rMatrixJJ ... submatrix jj (number of active dof x number of active dof)
    //! @param rMatrixJK ... submatrix jk (number of active dof x number of dependent dof)
    //! @param rMatrixKK ... submatrix kk (number of dependent dof x number of dependent dof)
    virtual void BuildGlobalCoefficientSubMatrices0Symmetric(NuTo::SparseMatrix<double>& rMatrixJJ, NuTo::SparseMatrix<double>& rMatrixJK, NuTo::SparseMatrix<double>& rMatrixKK) const = 0;

    //! @brief ... based on the global dofs build sub-vectors of the global internal potential gradient
    //! @param rActiveDofGradientVector ... global internal potential gradient which corresponds to the active dofs
    //! @param rDependentDofGradientVector ... global internal potential gradient which corresponds to the dependent dofs
    virtual void BuildGlobalGradientInternalPotentialSubVectors(NuTo::FullMatrix<double>& rActiveDofGradientVector, NuTo::FullMatrix<double>& rDependentDofGradientVector) const = 0;

};
} //namespace NuTo
#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
#include <boost/serialization/assume_abstract.hpp>
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::StructureBase)
#endif // SWIG
#endif  // ENABLE_SERIALIZATION
#endif // STRUCTUREBASE_H
