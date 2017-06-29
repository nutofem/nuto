#pragma once

#include "mechanics/constitutive/inputoutput/ConstitutiveVector.h"

namespace NuTo
{

class ConstitutiveScalar: public ConstitutiveVector<1>
{
public:
    ConstitutiveScalar()                                        = default;
    ConstitutiveScalar(const ConstitutiveScalar& )              = default;
    ConstitutiveScalar(      ConstitutiveScalar&&)              = default;

    virtual std::unique_ptr<ConstitutiveIOBase> clone() override
    {
        return std::make_unique<ConstitutiveScalar>(*this);
    }

    virtual ~ConstitutiveScalar()                               = default;

    ConstitutiveScalar& operator=(const ConstitutiveScalar& )   = default;
    ConstitutiveScalar& operator=(      ConstitutiveScalar&&)   = default;

    virtual Eigen::Matrix<double, 1, 1>& AsVector() override
    {
        throw Exception(std::string("[")+__PRETTY_FUNCTION__+"] You are calling a vector method on a scalar.");
    }

    Eigen::Matrix<double, 1, 1>& AsScalar()
    {
        return *this;
    }

};

} /* namespace NuTo */
