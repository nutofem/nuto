/*
 * Element1DInXD.h
 *
 *  Created on: 06 August 2015
 *      Author: phuschke
 */

#pragma once

#include "nuto/mechanics/elements/Element1D.h"

namespace NuTo
{

class StructureBase;
class DeformationGradient1D;
class TemperatureGradient1D;
class LocalEqStrain;
class NonlocalEqStrain;
class ConstitutiveTangentLocal1x1;
class EngineeringStress1D;
class HeatFlux1D;
template<int TNumRows, int TNumColumns> class ConstitutiveTangentLocal;
template<int TNumRows, int TNumColumns> class ConstitutiveTangentNonlocal;

class Element1DInXD: public Element1D
{

#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:
    Element1DInXD(const NuTo::StructureBase* rStructure, const std::vector<NuTo::NodeBase*>& rNodes, ElementData::eElementDataType rElementDataType, IpData::eIpDataType rIpDataType, InterpolationType* rInterpolationType);

    virtual ~Element1DInXD()
    {
    }
    ;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif  // ENABLE_SERIALIZATION

    //! @brief returns the enum (type of the element)
    //! @return enum
    NuTo::Element::eElementType GetEnumType() const override;

    const Eigen::MatrixXd ExtractNodeValues(int rTimeDerivative, Node::eDof) const override;
    const Eigen::MatrixXd ExtractGlobalNodeValues(int rTimeDerivative, Node::eDof rDofType) const;

    const Eigen::VectorXd InterpolateDofGlobal(const Eigen::VectorXd& rNaturalCoordinates, Node::eDof rDofType) const override;
    const Eigen::VectorXd InterpolateDofGlobal(int rTimeDerivative, const Eigen::VectorXd& rNaturalCoordinates, Node::eDof rDofType) const override;

    //! @brief adds to a matrix the product B^tCB, where B contains the derivatives of the shape functions and C is the constitutive tangent
    //! eventually include also area/width of an element (that's the thermal solution)
    //! @param rDerivativeShapeFunctions derivatives of the shape functions with respect to global coordinates
    //! @param ConstitutiveTangentBase constitutive tangent matrix
    //! @param rFactor factor including determinant of Jacobian and IP weight
    //! @param rRow row, where to start to add the submatrix
    //! @param rCoefficientMatrix to be added to
    void AddDetJBtCB(const Eigen::MatrixXd& rDerivativeShapeFunctionsGlobal, const ConstitutiveTangentLocal<1, 1>& rConstitutiveTangent, double rFactor, int rRow, int rCol, FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rCoefficientMatrix) const override;

    //! @brief adds up the internal force vector
    //! @param rDerivativeShapeFunctions derivatives of the shape functions with respect to global coordinates
    //! @param rEngineeringStress stress
    //! @param factor factor including det Jacobian area and integration point weight
    //! @param rRow start row (in case of a multifield problem)
    //! @param rResult resforce vector
    void AddDetJBtSigma(const Eigen::MatrixXd& rDerivativeShapeFunctionsGlobal, const EngineeringStress1D& rEngineeringStress, double factor, int rRow, FullVector<double, Eigen::Dynamic>& rResult) const override;

    //! @brief ... extract global dofs from nodes (mapping of local row ordering of the element matrices to the global dof ordering)
    const Eigen::VectorXi CalculateGlobalRowDofs() const override;

    //! @brief cast the base pointer to an Element1D, otherwise throws an exception
    const NuTo::Element1DInXD* AsElement1DInXD() const
    {
        return this;
    }

    //! @brief cast the base pointer to an Element1D, otherwise throws an exception
    NuTo::Element1DInXD* AsElement1DInXD()
    {
        return this;
    }

protected:

    //! @brief ... check if the element is properly defined (check node dofs, nodes are reordered if the element length/area/volum is negative)
    void CheckElement() override;

private:

    //! @brief constructor for serialization
    Element1DInXD() = default;

    //! @brief Calculates the rotation matrix. 2x2 in 2D, 3x3 in 3D
    //! @return mRotationMatrix
    Eigen::MatrixXd CalculateRotationMatrix();

    //! @brief Calculates the transformation matrix. The size of the matrix depends on the number of nodes and the global dimension of the element.
    //! @param rGlobalDimension: global dimension of the element, i.e. 2 or 3
    //! @param rNumberOfNodes: number of nodes of the element
    //! @return  transformationMatrix
    Eigen::MatrixXd CalculateTransformationMatrix(unsigned int rGlobalDimension, unsigned int rNumberOfNodes) const;

    //! @brief Returns the number of dofs for each node depending on the dof type, i.e. scalar or vector quantity
    //! @return number of dofs per node
    int GetNumDofsPerNode(Node::eDof rDofType) const;

    Eigen::MatrixXd mRotationMatrix;
};

} /* namespace NuTo */

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::Element1DInXD)
#endif
