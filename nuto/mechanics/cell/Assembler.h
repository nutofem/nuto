#pragma once
#include "nuto/mechanics/dofs/DofInfo.h"
#include "nuto/mechanics/dofs/DofContainer.h"
#include "nuto/mechanics/dofs/DofVector.h"
#include "nuto/mechanics/dofs/DofMatrix.h"
#include "nuto/mechanics/dofs/DofMatrixSparse.h"

namespace NuTo
{

std::vector<DofType> AvailableDofTypes(const DofVector<int> n, std::vector<DofType> dofTypes)
{
    if (dofTypes.empty())
        return n.DofTypes();

    std::vector<DofType> intersection;
    for (const auto& d1 : n.DofTypes())
        for (const auto& d2 : dofTypes)
            if (d1.Id() == d2.Id())
            {
                intersection.push_back(d1);
                continue;
            }
    return intersection;
}

class VectorAssembler
{
public:
    VectorAssembler(DofContainer<int> sizes = {})
    {
        Resize(sizes);
    }

    void Resize(DofContainer<int> sizes)
    {
        for (auto dofSize : sizes)
            mVector[dofSize.first].setZero(dofSize.second);
    }

    void Add(const DofVector<double>& v, const DofVector<int>& numbering, std::vector<DofType> dofTypes = {})
    {
        for (const auto& dofType : AvailableDofTypes(numbering, dofTypes))
        {
            for (unsigned localDofNumber = 0; localDofNumber < numbering[dofType].size(); ++localDofNumber)
            {
                const int globalDofNumber = numbering[dofType][localDofNumber];
                mVector[dofType][globalDofNumber] += v[dofType][localDofNumber];
            }
        }
    }

    void Reset()
    {
        for (const auto& dofType : mVector.DofTypes())
            mVector[dofType].setZero();
    }

    const DofVector<double>& Get() const
    {
        return mVector;
    }

private:
    DofVector<double> mVector;
};

class MatrixAssembler
{
public:
private:
    DofVector<double> mVector;
};


} /* NuTo */
