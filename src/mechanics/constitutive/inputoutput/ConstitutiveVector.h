#pragma once

#include "mechanics/constitutive/inputoutput/ConstitutiveMatrix.h"

namespace NuTo
{

template <int TRows>
class ConstitutiveVector : public ConstitutiveMatrix<TRows, 1>
{
public:
    ConstitutiveVector() = default;
    ConstitutiveVector(const ConstitutiveVector&) = default;
    ConstitutiveVector(ConstitutiveVector&&) = default;

    virtual ~ConstitutiveVector() = default;

    virtual std::unique_ptr<ConstitutiveIOBase> clone() override
    {
        return std::make_unique<ConstitutiveVector<TRows>>(*this);
    }

    ConstitutiveVector& operator=(const ConstitutiveVector&) = default;
    ConstitutiveVector& operator=(ConstitutiveVector&&) = default;


    double& operator[](int rRow) override
    {
        return this->Eigen::Matrix<double, TRows, 1>::operator[](rRow);
    }


    double operator[](int rRow) const override
    {
        return this->Eigen::Matrix<double, TRows, 1>::operator[](rRow);
    }

    virtual Eigen::Matrix<double, TRows, 1>& AsVector()
    {
        return *this;
    }


    //! @brief Convert vector of arbitrary dimension to 3d
    //! This sets the rest to zero. Useful for visualization.
    ConstitutiveVector<3> ConvertTo3DVector()
    {
        ConstitutiveVector<3> result;
        result.setZero();
        result.topRows(TRows) = this->AsVector();
        return result;
    }
};

} /* namespace NuTo */
