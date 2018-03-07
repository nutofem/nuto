#pragma once

#include "mechanics/dofs/DofVector.h"

namespace NuTo
{

class GlobalDofVector
{
public:
    NuTo::DofVector<double> J;
    NuTo::DofVector<double> K;

    double& operator()(DofType dof, int globalDofNumber)
    {
        const int numIndependent = J[dof].size();
        if (globalDofNumber < numIndependent)
            return J[dof][globalDofNumber];
        else
            return K[dof][globalDofNumber - numIndependent];
    }

    double operator()(DofType dof, int globalDofNumber) const
    {
        const int numIndependent = J[dof].size();
        if (globalDofNumber < numIndependent)
            return J[dof][globalDofNumber];
        else
            return K[dof][globalDofNumber - numIndependent];
    }

    std::vector<double> operator()(DofType dof, std::vector<int> globalDofNumbers) const
    {
        std::vector<double> v;
        v.reserve(globalDofNumbers.size());
        for (int globalDofNumber : globalDofNumbers)
            v.push_back((*this)(dof, globalDofNumber));
        return v;
    }

    GlobalDofVector& operator+=(const GlobalDofVector& rhs)
    {
        this->J += rhs.J;
        this->K += rhs.K;
        return *this;
    }

    GlobalDofVector& operator-=(const GlobalDofVector& rhs)
    {
        this->J.AddScaled(rhs.J, -1.);
        this->K.AddScaled(rhs.K, -1.);
        return *this;
    }
};

} /* NuTo */
