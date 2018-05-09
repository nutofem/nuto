#pragma once
#include <vector>

namespace NuTo
{
namespace Math
{

class Term
{
public:
    Term(size_t id, double value)
        : mData(id, value)
    {
    }
    size_t& Id()
    {
        return mData.first;
    }

    double& Coefficient()
    {
        return mData.second;
    }

private:
    std::pair<size_t, double> mData;
};

class Equation
{
public:
    std::vector<Term>& Terms()
    {
        return mData.first;
    }

    double& Value()
    {
        return mData.second;
    }

private:
    std::pair<std::vector<Term>, double> mData;
};

} // namespace Math
} // namespace NuTo
