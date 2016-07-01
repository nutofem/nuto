#pragma once

#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveIOBase.h"
#include <eigen3/Eigen/Core>

namespace NuTo
{

template <int TRows, int TCols>
class ConstitutiveMatrix : public ConstitutiveIOBase, public Eigen::Matrix<double, TRows, TCols>
{
public:

    ConstitutiveMatrix()                                        = default;
    ConstitutiveMatrix(const ConstitutiveMatrix& )              = default;
    ConstitutiveMatrix(      ConstitutiveMatrix&&)              = default;

    virtual ~ConstitutiveMatrix()                               = default;

    virtual std::unique_ptr<ConstitutiveIOBase> clone() override
    {
        return std::make_unique<ConstitutiveMatrix<TRows, TCols>>(*this);
    }

    ConstitutiveMatrix& operator=(const ConstitutiveMatrix& )   = default;
    ConstitutiveMatrix& operator=(      ConstitutiveMatrix&&)   = default;


    double& operator ()(int rRow, int rCol) override
    {
        return this->Eigen::Matrix<double, TRows, TCols>::operator () (rRow, rCol);
    }


    double operator ()(int rRow, int rCol) const override
    {
        return this->Eigen::Matrix<double, TRows, TCols>::operator () (rRow, rCol);
    }

    void SetZero() override
    {
        this->setZero();
    }

    int GetNumRows() const override
    {
        if (TRows == Eigen::Dynamic)
            return this->rows();
        return TRows;
    }

    int GetNumColumns() const override
    {
        return TCols;
    }

    virtual Eigen::Matrix<double, TRows, TCols>& AsMatrix()
    {
        return *this;
    }
};

} /* namespace NuTo */
