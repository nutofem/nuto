#pragma once


namespace NuTo
{
namespace Node
{

//! @brief provides a hash for the eDof enum - used e.g. in unordered maps
struct eDofHash
{
    template <typename T>
    std::size_t operator()(T t) const
    {
        return static_cast<std::size_t>(t);
    }
};

//! @brief provides a hash for the eDof enum - used e.g. in unordered maps
struct eDofPairHash
{
    template <typename T>
    std::size_t operator()(std::pair<T, T> t) const
    {
        return CombineDofs(t.first, t.second);
    }
};

}// namespace Node
}// namespace NuTo
