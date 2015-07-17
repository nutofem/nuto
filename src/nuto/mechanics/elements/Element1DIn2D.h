/*
 * Element1DIn2D.h
 *
 *  Created on: 13 July 2015
 *      Author: phuschke
 */

#ifndef ELEMENT1DIN2D_H_
#define ELEMENT1DIN2D_H_

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

class Element1DIn2D: public Element1D
{

public:
    Element1DIn2D(const NuTo::StructureBase* rStructure, const std::vector<NuTo::NodeBase*>& rNodes, ElementData::eElementDataType rElementDataType, IpData::eIpDataType rIpDataType, InterpolationType* rInterpolationType);

    virtual ~Element1DIn2D()
    {
    }
    ;

    //! @brief returns the enum (type of the element)
    //! @return enum
    NuTo::Element::eElementType GetEnumType() const override;

    //! @brief returns the global dimension of the element
    //! this is required to check, if an element can be used in a 1d, 2D or 3D Structure
    //! there is also a routine GetLocalDimension, which is e.g. 2 for plane elements and 1 for truss elements
    //! @return global dimension
    int GetGlobalDimension() const override;

    //! @brief returns the local dimension of the element
    //! this is required to check, if an element can be used in a 1d, 2D or 3D Structure
    //! @return local dimension
    int GetLocalDimension() const override;

    const Eigen::MatrixXd ExtractNodeValues(int rTimeDerivative, Node::eAttributes) const override;
    const Eigen::MatrixXd ExtractGlobalNodeValues(int rTimeDerivative, Node::eAttributes rDofType) const;

    const Eigen::VectorXd InterpolateDofGlobal(const Eigen::VectorXd& rNaturalCoordinates, Node::eAttributes rDofType) const override;
    const Eigen::VectorXd InterpolateDofGlobal(int rTimeDerivative, const Eigen::VectorXd& rNaturalCoordinates, Node::eAttributes rDofType) const override;

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
    const NuTo::Element1DIn2D* AsElement1DIn2D() const
    {
        return this;
    }

    //! @brief cast the base pointer to an Element1D, otherwise throws an exception
    NuTo::Element1DIn2D* AsElement1DIn2D()
    {
        return this;
    }

    //! @brief gets element length
    //! @return element length
    double GetElementLength() const;

    //! @brief sets element length
    //! @param rElementLength element length
    void SetElementLength(const double rElementLength);

protected:

    //! @brief ... check if the element is properly defined (check node dofs, nodes are reordered if the element length/area/volum is negative)
    void CheckElement() override;

private:
    //! @brief calculates the element length
    //! @return mElementLength
    double CalculateElementLength();

    //! @brief calculates the rotation matrix
    //! @return mRotationMatrix
    Eigen::Matrix2d CalculateRotationMatrix();

    double mElementLength;
    Eigen::Matrix2d mRotationMatrix;
};

} /* namespace NuTo */

#endif /* ELEMENT1DIN2D_H_ */
