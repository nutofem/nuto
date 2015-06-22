/*
 * Element2D.h
 *
 *  Created on: 9 Apr 2015
 *      Author: ttitsche
 */

#ifndef ELEMENT2D_H_
#define ELEMENT2D_H_

#include "nuto/mechanics/elements/ElementBase.h"

namespace NuTo
{


class StructureBase;
class DeformationGradient2D;
class TemperatureGradient2D;
class LocalEqStrain;
class NonlocalEqStrain;
class ConstitutiveTangentLocal3x3;
class EngineeringStress2D;
class HeatFlux2D;
template <int TNumRows, int TNumColumns> class ConstitutiveTangentLocal;
template <int TNumRows, int TNumColumns> class ConstitutiveTangentNonlocal;

class Element2D: public ElementBase
{

public:
    Element2D(const NuTo::StructureBase* rStructure,  const std::vector<NuTo::NodeBase* >& rNodes,
            ElementData::eElementDataType rElementDataType,IpData::eIpDataType rIpDataType, InterpolationType* rInterpolationType);

    virtual ~Element2D() {};

    //! @brief calculates output data for the element
    //! @param eOutput ... coefficient matrix 0 1 or 2  (mass, damping and stiffness) and internal force (which includes inertia terms)
    //!                    @param updateStaticData (with DummyOutput), IPData, globalrow/column dofs etc.
    Error::eError Evaluate(boost::ptr_multimap<NuTo::Element::eOutput, NuTo::ElementOutputBase>& rElementOutput) override;


    //! @brief returns the enum (type of the element)
    //! @return enum
    NuTo::Element::eElementType GetEnumType()const override;

    //! @brief returns the global dimension of the element
    //! this is required to check, if an element can be used in a 1d, 2D or 3D Structure
    //! there is also a routine GetLocalDimension, which is e.g. 2 for plane elements and 1 for truss elements
    //! @return global dimension
    int GetGlobalDimension()const override;

    //! @brief returns a pointer to the i-th node of the element
    //! @param local node number
    //! @return pointer to the node
    NodeBase* GetNode(int rLocalNodeNumber) override;

    //! @brief returns a pointer to the i-th node of the element
    //! @param local node number
    //! @return pointer to the node
    const NodeBase* GetNode(int rLocalNodeNumber)const override;

    //! @brief returns a pointer to the i-th node of the element
    //! @param local node number
    //! @brief rDofType dof type
    //! @return pointer to the node
    NodeBase* GetNode(int rLocalNodeNumber, Node::eAttributes rDofType) override;

    //! @brief returns a pointer to the i-th node of the element
    //! @param local node number
    //! @brief rDofType dof type
    //! @return pointer to the node
    const NodeBase* GetNode(int rLocalNodeNumber, Node::eAttributes rDofType)const override;

    //! @brief sets the rLocalNodeNumber-th node of the element
    //! @param local node number
    //! @param pointer to the node
    void SetNode(int rLocalNodeNumber, NodeBase* rNode) override;

    //! @brief resizes the node vector
    //! @param rNewNumNodes new number of nodes
    void ResizeNodes(int rNewNumNodes);

    //! brief exchanges the node ptr in the full data set (elements, groups, loads, constraints etc.)
    //! this routine is used, if e.g. the data type of a node has changed, but the restraints, elements etc. are still identical
    void ExchangeNodePtr(NodeBase* rOldPtr, NodeBase* rNewPtr) override;

    //! @brief sets the section of an element
    //! @param rSection pointer to section
    void SetSection(const SectionBase* rSection) override;

    //! @brief returns a pointer to the section of an element
    //! @return pointer to section
    const SectionBase* GetSection()const override;

    //! @brief Allocates static data for an integration point of an element
    //! @param rConstitutiveLaw constitutive law, which is called to allocate the static data object
    ConstitutiveStaticDataBase* AllocateStaticData(const ConstitutiveBase* rConstitutiveLaw)const override;

    //! @brief calculates the volume of an integration point (weight * detJac)
    //! @return rVolume  vector for storage of the ip volumes (area in 2D, length in 1D)
    const Eigen::VectorXd GetIntegrationPointVolume() const override;

    const Eigen::MatrixXd ExtractNodeValues(int rTimeDerivative, Node::eAttributes) const override;

    void ExtractNodeValues(Eigen::MatrixXd& rNodeValues, int rTimeDerivative, Node::eAttributes rDofType) const;

    //! @brief Calculates the the inverse of the Jacobian and its determinant
    //! @param rDerivativeShapeFunctions Derivatives of the shape functions (dN1dx, dN1dy, dN1dz, dN2dx, ..
    //! @param rNodeCoordinates Node coordinates (X1,Y1,Z1,X2,Y2,Z2,...
    //! @param rInvJacobian inverse Jacobian matrix (return value)
    //! @param rDetJac determinant of the Jacobian (return value)
    void CalculateJacobian(const Eigen::MatrixXd& rDerivativeShapeFunctions,
                           const Eigen::MatrixXd& rNodeCoordinates,
                           Eigen::Matrix2d& rInvJacobian,
                           double& rDetJac)const;

    double CalculateArea() const override;

    //! @brief ... extract global dofs from nodes (mapping of local row ordering of the element matrices to the global dof ordering)
    const Eigen::VectorXi CalculateGlobalRowDofs() const;

    //! @brief ... extract global dofs from nodes (mapping of local column ordering of the element matrices to the global dof ordering)
    const Eigen::VectorXi CalculateGlobalColumnDofs() const;

    //! @brief calculates the deformation gradient in 2D
    //! @param rDerivativeShapeFunctionsLocal derivatives of the shape functions with respect to local coordinates (plane world coordinates)
    //! @param rNodalDisplacements local displacements
    //! @return rDeformationGradient (return value)
    const DeformationGradient2D CalculateDeformationGradient(const Eigen::MatrixXd& rDerivativeShapeFunctionsLocal,
                                      const Eigen::MatrixXd& rNodalDisplacements)const;

    //! @brief adds to a matrix the product B^tCBnonlocal, where B contains the derivatives of the shape functions and C is the constitutive tangent and Bnonlocal is the nonlocal B matrix
    //! eventually include also area/width of an element
    //! @param rLocalDerivativeShapeFunctions derivatives of the shape functions with respect to global coordinates
    //! @param rNonlocalDerivativeShapeFunctions derivatives of the shape functions with respect to global coordinates
    //! @param ConstitutiveTangentBase constitutive tangent matrix
    //! @param rFactor factor including determinant of Jacobian and IP weight
    //! @param rCoefficientMatrix to be added to
    //! &param rFirstCol first column of the coefficient matrix to be modified (corresponding to the current nonlocal element)
    void AddDetJBtCB(const Eigen::MatrixXd& rLocalDerivativeShapeFunctionsLocal,const Eigen::MatrixXd& rNonlocalDerivativeShapeFunctionsLocal,
                                  const ConstitutiveTangentLocal<3,3>& rConstitutiveTangent, double rFactor,
                                  FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rCoefficientMatrix, int rFirstRow, int rFirstCol)const;

    //! @brief adds to a matrix the product B^tCB, where B contains the derivatives of the shape functions and C is the constitutive tangent
    //! eventually include also area/width of an element (that's the thermal solution)
    //! @param rDerivativeShapeFunctions derivatives of the shape functions with respect to global coordinates
    //! @param ConstitutiveTangentBase constitutive tangent matrix
    //! @param rFactor factor including determinant of Jacobian and IP weight
    //! @param rRow row, where to start to add the submatrix
    //! @param rCoefficientMatrix to be added to
    void AddDetJBtCB(const Eigen::MatrixXd& rDerivativeShapeFunctionsGlobal,
                                  const ConstitutiveTangentLocal<2,2>& rConstitutiveTangent, double rFactor,
                                  int rRow, int rCol,
                                  FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rCoefficientMatrix)const;

    //! @brief adds up the internal force vector
    //! @param rDerivativeShapeFunctions derivatives of the shape functions with respect to global coordinates
    //! @param rEngineeringStress stress
    //! @param factor factor including det Jacobian area and integration point weight
    //! @param rRow start row (in case of a multifield problem)
    //! @param rResult resforce vector
    void AddDetJBtSigma(const Eigen::MatrixXd& rDerivativeShapeFunctionsGlobal,
                        const EngineeringStress2D& rEngineeringStress,
                        double factor,
                        int rRow,
                        FullVector<double,Eigen::Dynamic>& rResult)const;

    //! @brief calculates the Kee matrix
    //! @param rShapeFunctions of the ip for all shape functions
    //! @param rDerivativeShapeFunctions of the ip for all shape functions
    //! @param nonlocal gradient radius xi
    //! @param rFactor multiplication factor (detJ area..)
    //! @param Kee return matrix with detJ * (Nt N + cBtB)
    void CalculateKee(Eigen::VectorXd rShapeFunctions,
                        const Eigen::MatrixXd& rDerivativeShapeFunctions,
                        ConstitutiveTangentLocal<1, 1>& rNonlocalParameter,
                        double rFactor,
                        FullMatrix<double,
                        Eigen::Dynamic,
                        Eigen::Dynamic>& rKee) const;

    //! @brief add Kee*nonlocalEqStrain-detJ*N.T*localEqStrain (detJ is already included in Kee)
    //! @param rShapeFunctions of the ip for all shape functions
    //! @param rLocalEqStrain local eq. strain values
    //! @param rKee stiffness matrix Kee
    //! @param rNodeNonlocalEqStrain nodal nonlocal eq strain values
    //! @param rFactor factor including detJ and area
    //! @param rResult result
    void AddDetJRnonlocalEqStrain(const Eigen::VectorXd& rShapeFunctions,
                        LocalEqStrain& rLocalEqStrain,
                        FullMatrix<double,
                        Eigen::Dynamic, Eigen::Dynamic>& rKee,
                        Eigen::Matrix<double, Eigen::Dynamic,
                        Eigen::Dynamic>& rNodeNonlocalEqStrain,
                        double rFactor,
                        int startIndexNonlocalEqStrain,
                        FullVector<double, Eigen::Dynamic>& rResult) const;

    //! @brief add detJ B.T dSigma/dnonlocalEqStrain N
    //! @param derivativeShapeFunctions of the ip for all shape functions
    //! @param tangentStressNonlocalEqStrain derivative of the stress with respect to the nonlocal eq strain
    //! @param rShapeFunctions of the ip for all shape functions
    //! @param rFactor factor including detJ and area
    //! @param rResult result
    void AddDetJBtdSigmadNonlocalEqStrainN(const Eigen::MatrixXd& rDerivativeShapeFunctions,
            ConstitutiveTangentLocal<3, 1>& rTangentStressNonlocalEqStrain,
            Eigen::VectorXd rShapeFunctions,
            double rFactor,
            int rRow,
            int rCol,
            FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rResult) const;

    //! @brief add detJ N_transpose dEqStrain/dEpsilon B
    //! @param rShapeFunctions of the ip for the nonlocal eq strain dofs
    //! @param rTangentLocalEqStrainStrain derivative of the local eq strains with respect to the strain
    //! @param rderivativeShapeFunctions of the ip for the displacement dofs
    //! @param rFactor factor including detJ and area
    //! @param rResult result
    void AddDetJNtdLocalEqStraindEpsilonB(Eigen::VectorXd rShapeFunctions, ConstitutiveTangentLocal<3, 1>& rTangentLocalEqStrainStrain, const Eigen::MatrixXd& rDerivativeShapeFunctions,
            double rFactor, int rRow, int rCol, FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rResult) const;

    //! @brief cast the base pointer to an Element2D, otherwise throws an exception
    const NuTo::Element2D* AsElement2D()const override
    {
        return this;
    }

    //! @brief cast the base pointer to an Element2D, otherwise throws an exception
    NuTo::Element2D* AsElement2D()
    {
        return this;
    }

protected:
    //! @brief ... reorder nodes such that the sign of the length/area/volume of the element changes
    void ReorderNodes() override;

    //! @brief ... check if the element is properly defined (check node dofs, nodes are reordered if the element length/area/volum is negative)
    void CheckElement() override;




private:
    std::vector<NodeBase*> mNodes;

    const SectionBase *mSection;
};

} /* namespace NuTo */

#endif /* ELEMENT2D_H_ */
