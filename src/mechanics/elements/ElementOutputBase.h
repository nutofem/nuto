// $Id $
#pragma once

#include <vector>

namespace NuTo
{
class ElementOutputIpData;


template <typename T>
class BlockFullMatrix;
template <typename T>
class BlockFullVector;

//! @author Joerg F. Unger
//! @date Apr 29, 2010
//! @brief ...
class ElementOutputBase
{
public:
    ElementOutputBase();

    virtual ~ElementOutputBase();

    virtual BlockFullMatrix<double>& GetBlockFullMatrixDouble();

    virtual BlockFullVector<double>& GetBlockFullVectorDouble();

    virtual BlockFullVector<int>& GetBlockFullVectorInt();

    virtual std::vector<int>& GetVectorInt();

    virtual ElementOutputIpData& GetIpData();

    virtual void SetSymmetry(bool rSymmetric);

    virtual bool GetSymmetry() const;

    virtual void SetConstant(bool rConstant);

    virtual bool GetConstant() const;

    virtual ElementOutputBase* Clone() const = 0;
};
}
