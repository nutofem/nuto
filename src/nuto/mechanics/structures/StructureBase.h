// $Id$

#ifndef STRUCTUREBASE_H
#define STRUCTUREBASE_H

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#endif // ENABLE_SERIALIZATION

#include <boost/ptr_container/ptr_map.hpp>
#include <boost/ptr_container/ptr_list.hpp>
#include <string>
#include "nuto/base/NuToObject.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"
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

    //! @brief ... Add engineering strains to the internal list, which is finally exported via the ExportVtkDataFile command
    void AddVisualizationComponentEngineeringStrain();

    //! @brief ... Add engineering plastic strains to the internal list, which is finally exported via the ExportVtkDataFile command
    void AddVisualizationComponentEngineeringPlasticStrain();

    //! @brief ... Add engineering stress to the internal list, which is finally exported via the ExportVtkDataFile command
    void AddVisualizationComponentEngineeringStress();

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

    //! @brief ... export an element group to Vtk data file
    //! @param rGroupIdent ... group ident
    //! @param rFileName ... file name
    //! @param rWhat ... string which describes what to plot
    void ElementGroupExportVtkDataFile(const std::string& rGroupIdent, const std::string& rFileName) const;

#ifndef SWIG
    //! @brief ... export elements to Vtk data file
    //! @param rElements ... vector of elements
    //! @param rFileName ... file name
    void ExportVtkDataFile(const std::vector<const ElementBase*>& rElements, const std::string& rFileName) const;
#endif //SWIG
#endif // ENABLE_VISUALIZE

    //! @brief ... build global coefficient matrix (e.g stiffness) for primary dofs (e.g displacements, rotations, temperature)
    //! @param rMatrix ... global coefficient matrix (nonsymmeetric)
    //! @param rVector ... global equivalent load vector (e.g. due to prescribed displacements)
    void BuildGlobalCoefficientMatrix0(NuTo::SparseMatrixCSRGeneral<double>& rMatrix, NuTo::FullMatrix<double>& rVector);

    //! @brief ... build global coefficient matrix (e.g stiffness) for primary dofs (e.g displacements, rotations, temperature)
    //! @param rMatrix ... global coefficient matrix (symmetric)
    //! @param rVector ... global equivalent load vector (e.g. due to prescribed displacements)
    void BuildGlobalCoefficientMatrix0(NuTo::SparseMatrixCSRSymmetric<double>& rMatrix, NuTo::FullMatrix<double>& rVector);

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

    //! @brief gives the identifier of a node
    //! @param reference to a node
    //! @return identifier
    virtual int NodeGetId(const NodeBase* rNode)const=0;
#endif //SWIG

    //! @brief delete node
    //! @param rIdent ... node identifier
    virtual void NodeDelete(const int rIdent)=0;

    //! @brief info about the nodes in the Structure
    virtual void NodeInfo(int mVerboseLevel)const=0;

    //! @brief numbers the dofs in the structure
    virtual void NodeBuildGlobalDofs()=0;

    //! @brief returns the number of global dofs
    //! @return number of global dofs
    int NodeGetNumberGlobalDofs()const;

    //! @brief returns the number of active dofs
    //! @return number of active dofs
    int NodeGetNumberActiveDofs()const;

    //! @brief sets the displacements of a node
    //! @param rIdent node identifier
    //! @param rDisplacements matrix (one column) with the displacements
    void NodeSetDisplacements(int rId,const NuTo::FullMatrix<double>& rDisplacements);

    //! @brief gets the displacements of a node
    //! @param rIdent node identifier
    //! @param rDisplacements matrix (one column) with the displacements
    void NodeGetDisplacements(int rId, NuTo::FullMatrix<double>& rDisplacements)const;

    //! @brief extract dof values (e.g. displacements, temperatures to the nodes)
    //! @param rActiveDofValues ... vector of global active dof values (ordering according to global dofs, size is number of active dofs)
    //! @param rDependentDofValues ... vector of global dependent dof values (ordering according to (global dofs) - (number of active dofs), size is (total number of dofs) - (number of active dofs))
    virtual void NodeExtractDofValues(NuTo::FullMatrix<double>& rActiveDofValues, NuTo::FullMatrix<double>& rDependentDofValues) const = 0;

//*************************************************
//************ Element routines     ***************
//**  defined in structures/StructureElement.cpp **
//*************************************************
    //! @brief returns the number of nodes
    //! @return number of nodes
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
#endif //SWIG

    //! @brief delete element
    //! @param rIdent ... element identifier
    virtual void ElementDelete(const int rIdent)=0;

    //! @brief info about the elements in the Structure
    virtual void ElementInfo(int mVerboseLevel)const=0;

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
    void ElementGroupSetConstitutiveLaw(const std::string& rGroupIdent, int rConstitutiveLawIdent);

    //! @brief modifies the constitutive law of a all elements
    //! @param rConstitutiveLawIdent identifier for the material
    void ElementTotalSetConstitutiveLaw(int rConstitutiveLawIdent);

#ifndef SWIG
    //! @brief modifies the constitutive law of a single element
    //! @param rElement element pointer
    //! @param rConstitutive material pointer
    void ElementSetConstitutiveLaw(ElementBase* rElement, ConstitutiveBase* rConstitutive);
#endif //SWIG

    //! @brief modifies the section of a single element
    //! @param rElementIdent element number
    //! @param rSectionIdent identifier for the section
    void ElementSetSection(int rElementId, const std::string& rSectionIdent);

    //! @brief modifies the section of a group of elements
    //! @param rGroupIdent identifier for the group of elements
    //! @param rSectionId identifier for the section
    void ElementGroupSetSection(const std::string& rGroupIdent, const std::string& rSectionIdent);

    //! @brief modifies the section of a all elements
    //! @param rSectionIdent identifier for the section
    void ElementTotalSetSection(const std::string& rSectionIdent);

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
    void ElementGroupSetIntegrationType(const std::string& rGroupIdent, const std::string& rIntegrationTypeIdent, std::string rIpDataTypeStr);

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

    //! @brief updates the history data of a all elements
    void ElementTotalUpdateStaticData();

    //! @brief updates the temprory static data of a all elements
    //! its is a const function, since only mutuable data (instead of const) is updated (kind of temporary data)
    void ElementTotalUpdateTmpStaticData();

    //*************************************************
    //************ Constraint routines     ***************
    //**  defined in StructureBaseConstraints.cpp **
    //*************************************************
#ifndef SWIG
    //! @brief adds a displacement constraint equation for a node
    //! @param rNode pointer to node
    //! @param rDirection direction of the constraint (in 2D a point with 2 entries, in 3D 3 entries, in 1D not used)
    //! @param rValue prescribed value (e.g. zero to fix a displacement to zero)
    //! @return integer id to delete or modify the constraint
    int ConstraintSetDisplacementNode(NodeBase* rNode, const NuTo::FullMatrix<double>& rDirection, double rValue);
#endif

    //! @brief adds a displacement constraint equation for a node
    //! @param rNode identifier for node
    //! @param rComponent e.g. the first (count from zero) displacement component
    //! @param rValue prescribed value (e.g. zero to fix a displacement to zero)
    int  ConstraintSetDisplacementNode(int rIdent, const NuTo::FullMatrix<double>& rDirection, double rValue);

#ifndef SWIG
    //! @brief adds a displacement constraint equation for a group of node
    //! @param rNode pointer to group of nodes
    //! @param rDirection direction of the constraint (in 2D a point with 2 entries, in 3D 3 entries, in 1D not used)
    //! @param rValue prescribed value (e.g. zero to fix a displacement to zero)
    //! @return integer id to delete or modify the constraint
    int ConstraintSetDisplacementNodeGroup(Group<NodeBase>* rGroup, const NuTo::FullMatrix<double>& rDirection, double rValue);
#endif

    //! @brief adds a constraint equation for a group of nodes
    //! @param rGroupIdent identifier for group of nodes
    //! @param rAttribute displacements, rotations, temperatures
    //! @param rComponent e.g. the first (count from zero) displacement component
    //! @param rValue prescribed value (e.g. zero to fix a displacement to zero)
    int ConstraintSetDisplacementNodeGroup(std::string rGroupIdent, const NuTo::FullMatrix<double>& rDirection, double rValue);

    //! @brief returns the number of constraint equations
    //! @return number of constraints
    int ConstraintGetNumConstraintEquations()const;

    //! @brief calculates the constraint matrix that builds relations between the nodal dagrees of freedom
    //! rConstraintMatrix*DOFS = RHS
    //! @param rConstraintMatrix constraint matrix
    //! @param rRHS right hand side
    void ConstraintGetConstraintMatrix(NuTo::SparseMatrixCSRGeneral<double>& rConstraintMatrix, NuTo::FullMatrix<double>& rRHS);

    //!@brief sets/modifies the right hand side of the constraint equations
    //!@param rRHS new right hand side
    //!@param rRHS new right hand side
    void ConstraintSetRHS(int rConstraintEquation, double rRHS);

    //! @brief ... create a constraint equation
    //! @param rNode ... node id in the first constraint equation term
    //! @param rDof ... dof in the first constraint equation term (e.g "X_DISPLACEMENT", "Z_Rotation", "Temperature")
    //! @param rCoefficient ... weight factor of this term
    //! @param rRHS ... prescribed right hand side value
    //! @return integer id of the constraint
    int ConstraintEquationCreate(int rNode, const std::string& rDof, double rCoefficient, double rRHS = 0);

    //! @brief ... create a constraint equation
    //! @param rConstraint ... constraint id
    //! @param rNode ... node id in the first constraint equation term
    //! @param rDof ... dof in the first constraint equation term (e.g "X_DISPLACEMENT", "Z_Rotation", "Temperature")
    //! @param rCoefficient ... weight factor of this term
    //! @param rRHS ... prescribed right hand side value
    void ConstraintEquationCreate(int rConstraint, int rNode, const std::string& rDof, double rCoefficient, double rRHS = 0);

#ifndef SWIG
    //! @brief ... create a constraint equation
    //! @param rConstraint ... constraint id
    //! @param rNode ... node id in the first constraint equation term
    //! @param rDofType ... type of dof in the first constraint equation term (e.g DISPLACEMENTS, ROTATIONS, TEMPERATURES)
    //! @param rDofComponent ... dof component (0, 1, 2)
    //! @param rCoefficient ... weight factor of this term
    //! @param rRHS ... prescribed right hand side value
    void ConstraintEquationCreate(int rConstraint, int rNode, NuTo::Node::eAttributes rDofType, int rDofComponent, double rCoefficient, double rRHS = 0);
#endif

    //! @brief ... add a term to a constraint equation
    //! @param rConstraint ... constraint id
    //! @param rNode ... node id
    //! @param rDof ... dof (e.g "X_DISPLACEMENT", "Z_Rotation", "Temperature")
    //! @param rCoefficient ... weight factor of this term
    void ConstraintEquationAddTerm(int rConstraint, int rNode, const std::string& rDof, double rCoefficient);

#ifndef SWIG
    //! @brief ... add a term to a constraint equation
    //! @param rConstraint ... constraint id
    //! @param rNode ... node id
    //! @param rDofType ... type of dof (e.g DISPLACEMENTS, ROTATIONS, TEMPERATURES)
    //! @param rDofComponent ... dof component (0, 1, 2)
    //! @param rCoefficient ... weight factor of this term
    void ConstraintEquationAddTerm(int rConstraint, int rNode, NuTo::Node::eAttributes rDofType, int rDofComponent, double rCoefficient);
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
    int LoadCreateNodeGroupForce(std::string rGroupIdent, const NuTo::FullMatrix<double>& rDirection, double rValue);

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
    void SectionCreate(const std::string& rIdent, const std::string& rType);

    //! @brief ... delete an existing section
    //! @param rIdent ... section identifier
    void SectionDelete(const std::string& rIdent);

    //! @brief ... set section cross-section area
    //! @param rIdent ... section identifier
    //! @param rArea ... cross-section area
    void SectionSetArea(const std::string& rIdent, double rArea);

    //! @brief ... get section cross-section area
    //! @param rIdent ... section identifier
    //! @return section cross-section area
    double SectionGetArea(const std::string& rIdent) const;

    //! @brief ... set section thickness
    //! @param rIdent ... section identifier
    //! @param rThickness ... cross-section thickness
    void SectionSetThickness(const std::string& rIdent, double rThickness);

    //! @brief ... get section thickness
    //! @param rIdent ... section identifier
    //! @return section thickness
    double SectionGetThickness(const std::string& rIdent) const;

    //! @brief ... print information about all sections
    //! @param rVerboseLevel ... controls the verbosity of the information
    void SectionInfo(unsigned short rVerboseLevel) const;

    //! @brief ... print information of a single section
    //! @param rIdent ... section identifier
    //! @param rVerboseLevel ... controls the verbosity of the information
    void SectionInfo(const std::string& rIdent, unsigned short rVerboseLevel) const;

#ifndef SWIG
    //! @brief ... create a new section
    //! @param rIdent ... section identifier
    //! @param rType ... section type
    void SectionCreate(const std::string& rIdent, Section::eSectionType rType);

    //! @brief ... get the pointer to a section from the section identifier
    //! @param rIdent ... section identifier
    //! @return ... pointer to the section
    SectionBase* SectionGetSectionPtr(const std::string& rIdent);

    //! @brief ... get the pointer to a section from the section identifier
    //! @param rIdent ... section identifier
    //! @return ... pointer to the section
    const SectionBase* SectionGetSectionPtr(const std::string& rIdent) const;

    //! @brief ... get the identifier of an section from the pointer to section object
    //! @param rSectionPtr ... pointer to section object
    //! @return ... section identifier
    std::string SectionGetId(const SectionBase* rSectionPtr) const;
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
    virtual std::string GroupGetId(GroupBase* rGroup)const;

    //! @brief ... get the pointer to a group from the group identifier
    //! @param rIdent ... group identifier
    //! @return ... pointer to the group
    GroupBase* GroupGetGroupPtr(const std::string& rIdent);

    //! @brief ... get the pointer to a group from the group identifier
    //! @param rIdent ... group identifier
    //! @return ... pointer to the group
    const GroupBase* GroupGetGroupPtr(const std::string& rIdent) const;
#endif

    //! @brief ... Creates a group for the structure
    //! @param ... rIdent identifier for the group
    //! @param ... rType  type of the group, e.g. "NODES" or "ELEMENTS"
    void GroupCreate(const std::string& rIdent, const std::string& rType);

    //! @brief ... Deletes a group from the structure
    //! @param ... rIdent identifier for the group
    void GroupDelete(const std::string& rIdent);

    //! @brief ... Unites two groups and stores the result in a new group
    //! @param ... rIdentGroup1 identifier for the first group
    //! @param ... rIdentGroup2 identifier for the second group
    //! @result ... rIdentGroupResult identifier for the created result group
    void GroupUnion(const std::string& rIdentGroup1, const std::string& rIdentGroup2, const std::string& rIdentGroupResult);
    //! @brief ... Difference between two groups and stores the result in a new group
    //! @param ... rIdentGroup1 identifier for the first group
    //! @param ... rIdentGroup2 identifier for the second group
    //! @result ... rIdentGroupResult identifier for the created result group
    void GroupDifference(const std::string& rIdentGroup1, const std::string& rIdentGroup2, const std::string& rIdentGroupResult);

    //! @brief ... Calculates the intersection between two groups and stores the result in a new group
    //! @param ... rIdentGroup1 identifier for the first group
    //! @param ... rIdentGroup2 identifier for the second group
    //! @result ... rIdentGroupResult identifier for the created result group
    void GroupIntersection(const std::string& rIdentGroup1, const std::string& rIdentGroup2, const std::string& rIdentGroupResult);

    //! @brief ... Calculates the symmetric difference between two groups and stores the result in a new group
    //! @param ... rIdentGroup1 identifier for the first group
    //! @param ... rIdentGroup2 identifier for the second group
    //! @result ... rIdentGroupResult identifier for the created result group
    void GroupSymmetricDifference(const std::string& rIdentGroup1, const std::string& rIdentGroup2, const std::string& rIdentGroupResult);

    //! @brief ... Adds a node to a node group
    //! @param ... rIdentGroup identifier for the group
    //! @param ... rIdentNode  identifier for the node
    void GroupAddNode(const std::string& rIdentGroup, int rIdNode);

    //! @brief ... Adds an element to an element group
    //! @param ... rIdentGroup identifier for the group
    //! @param ... rIdentElement  identifier for the element
    void GroupAddElement(const std::string& rIdentGroup, int rIdElement);

    //! @brief ... Returns the number of members in a group
    //! @param ... rIdentGroup identifier for the group
    //! @return ... number of members
    int GroupGetNumMembers(const std::string& rIdentGroup)const;

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
    const NuTo::IntegrationTypeBase* GetPtrIntegrationType(const std::string& rIdentIntegrationType);

    //! @brief ... Returns a pointer to an integration type
    //! if the integration type does not exist (in the map), the integration type is created
    //! @param identIntegrationType Identifier for an integration type
    const NuTo::IntegrationTypeBase* GetPtrIntegrationType(NuTo::IntegrationType::eIntegrationType rIdentIntegrationType);
#endif //SWIG

    //*************************************************
    //************ Basic routines     ***************
    //**  defined in structures/StructureBase.cpp **
    //*************************************************
    //! @brief ... Info routine that prints general information about the object (detail according to verbose level)
    virtual void Info()const;


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
    boost::ptr_map<std::string,GroupBase> mGroupMap;

    //! @brief ... map storing the name and the pointer to the integration types
    //! @sa IntegrationTypeBase
    boost::ptr_map<std::string,IntegrationTypeBase> mIntegrationTypeMap;

    //! @brief ... map storing the section name and the pointer to the section object
    //! @sa SectionBase
    boost::ptr_map<std::string,SectionBase> mSectionMap;

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

    //! @brief ... store all elements of a structure in a vector
    //! @param rElements ... vector of element pointer
    virtual void GetElementsTotal(std::vector<const ElementBase*>& rElements) const = 0;

    //! @brief ... store all elements of a structure in a vector
    //! @param rElements ... vector of element pointer
    virtual void GetElementsTotal(std::vector<ElementBase*>& rElements) = 0;

    //! @brief ... store all elements of a group in a vector
    //! @param rElementGroup ... element group
    //! @param rElements ... vector of element pointer
    void GetElementsByGroup(const Group<ElementBase>* rElementGroup, std::vector<const ElementBase*>& rElements) const;

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
