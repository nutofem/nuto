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

    const Eigen::MatrixXd ExtractNodeValues(int rTimeDerivative, Node::eDof) const override;

    //! @brief calculates output data for the element
    Error::eError Evaluate(boost::ptr_multimap<NuTo::Element::eOutput, NuTo::ElementOutputBase>& rElementOutput) override;


    // VHIRTHAM TODO TEMP BEGIN

    //! @brief calculates output data for the element
    //! @param eOutput ... coefficient matrix 0 1 or 2  (mass, damping and stiffness) and internal force (which includes inertia terms)
    //!                    @param updateStaticData (with DummyOutput), IPData, globalrow/column dofs etc.
    Error::eError Evaluate(std::map<Element::eOutput, std::shared_ptr<ElementOutputBase>>& rElementOutput) override
    {
        throw NuTo::MechanicsException(std::string("[")+__PRETTY_FUNCTION__ + "] --- temp implementation missing");
    }

    // VHIRTHAM TODO TEMP END





    void AddStiffnessMatrix(double rSpringStiffness, NuTo::FullVector<double, Eigen::Dynamic> rSpringDirection, unsigned int rRow, unsigned int rCol, FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rCoefficientMatrix) const;

    const Eigen::VectorXi CalculateGlobalRowDofs() const override;


protected:

    //! @brief ... check if the element is properly defined
    void CheckElement() override;

};

} /* namespace NuTo */

