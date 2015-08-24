/*
 * Element1DSpring.h
 *
 *  Created on: 18 August 2015
 *      Author: phuschke
 */

#pragma once

#include "nuto/mechanics/elements/Element1D.h"

namespace NuTo
{

class Element1DSpring: public Element1D
{

public:
    Element1DSpring(const NuTo::StructureBase* rStructure, const std::vector<NuTo::NodeBase*>& rNodes, ElementData::eElementDataType rElementDataType, IpData::eIpDataType rIpDataType, InterpolationType* rInterpolationType);

    virtual ~Element1DSpring()
    {
    }
    ;

    //! @brief returns the enum (type of the element)
    //! @return enum
    NuTo::Element::eElementType GetEnumType() const override;

    const Eigen::MatrixXd ExtractNodeValues(int rTimeDerivative, Node::eAttributes) const override;

    //! @brief calculates output data for the element
    Error::eError Evaluate(boost::ptr_multimap<NuTo::Element::eOutput, NuTo::ElementOutputBase>& rElementOutput) override;

    void AddStiffnessMatrix(double rSpringStiffness, NuTo::FullVector<double, Eigen::Dynamic> rSpringDirection, unsigned int rRow, unsigned int rCol, FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rCoefficientMatrix) const;

    const Eigen::VectorXi CalculateGlobalRowDofs() const override;


protected:

    //! @brief ... check if the element is properly defined
    void CheckElement() override;

};

} /* namespace NuTo */

