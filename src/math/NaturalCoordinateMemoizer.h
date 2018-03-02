#pragma once
#include <Eigen/Core>
#include <map>
namespace NuTo
{

template <typename TVector>
struct CompareVector
{

    //! @brief defines a compare operator, following the compare concept
    //! http://en.cppreference.com/w/cpp/concept/Compare
    //! @remark Looks wierd, is fast. std::tie performs elementwise comparisons.
    bool operator()(const TVector& l, const TVector& r) const
    {
        assert(l.rows() == r.rows());
        switch (l.rows())
        {
        case 1:
            return l[0] < r[0];
        case 2:
            return std::tie(l[0], l[1]) < std::tie(r[0], r[1]);
        case 3:
            return std::tie(l[0], l[1], l[2]) < std::tie(r[0], r[1], r[2]);
        default:
            throw;
        }
    }
};

// https://en.wikipedia.org/wiki/Memoization: Not to be confused with Memorization.
template <typename TResult, typename TNaturalCoords, typename TCompare = CompareVector<TNaturalCoords>>
class NaturalCoordinateMemoizerMap
{
public:
    //! @brief ctor
    //! @param function function for memoization
    NaturalCoordinateMemoizerMap(std::function<TResult(TNaturalCoords)> function)
        : mFunction(function)
    {
    }

    //! @brief returns the value of the function for the given arguments.
    //!        Repeated calls of Get() with the same arguments will return the memoized
    //!        result.
    //! @param v argument
    //! @return reference to the 'memoized' result
    const TResult& Get(const TNaturalCoords& v) const
    {
        auto it = mCache.find(v);
        if (it == mCache.end())
        {
// The memoizer is shared between all threads in the parallel assembly.
// If it is still empty, all threads want to calculate and emplace the values
// into the cache. This will cause problems like infinite recursions, segfaults
// during map::rebalancing and so on. Thus: omp critical.
// The 'hot' path however, this the code above, which is hopefully not
// affected by the performance loss of omp critical.
#ifdef _OPENMP
#pragma omp critical
#endif
            it = mCache.emplace_hint(it, v, mFunction(v));
        }

        return it->second;
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
