#pragma once
#include <eigen3/Eigen/Core>
#include <vector>
#include <array>
#include <map>
#include <unordered_map>
#include <memory>


namespace NuTo
{

template <typename TVector>
struct CompareVector
{
    bool operator()(const TVector& l, const TVector& r) const
    {
        assert(l.size() == r.size());
        //int size = std::min(l.size(), r.size());
        int size = l.size(); 
        int iDim = 0;
        for (; iDim < size - 1; ++iDim)
            if (l[iDim] != r[iDim])
                return l[iDim] < r[iDim];
        // iDim = size -1 = last entry;
        return l[iDim] < r[iDim];
    }
};

// https://en.wikipedia.org/wiki/Memoization: Not to be confused with Memorization.
template <typename TResult, typename TNaturalCoords, typename TCompare = CompareVector<TNaturalCoords>>
class NaturalCoordinateMemoizerMap
{
public:
    //! @brief ctor
    //! @param function ... function for memoization
    NaturalCoordinateMemoizerMap(std::function<TResult(TNaturalCoords)> function)
        : mFunction(function)
    {
    }

    //! @brief returns the value of the function for the given arguments.
    //!        Repeated calls of Get() with the same arguments will return the memoized
    //!        result.
    //! @param v ... argument
    //! @return reference to the 'memoized' result
    const TResult& Get(const TNaturalCoords& v) const
    {
        auto it = mCache.find(v);
        if (it != mCache.end())
            return it->second;

        auto emplace_result = mCache.emplace(v, mFunction(v));
        return emplace_result.first->second;
    }

    void ClearCache() const
    {
        mCache.clear();
    }

private:
    //! @brief mutable to mark Get() and ClearCache() const
    mutable std::map<TNaturalCoords, TResult, TCompare> mCache;

    std::function<TResult(TNaturalCoords)> mFunction;
};


} /* NuTo */
