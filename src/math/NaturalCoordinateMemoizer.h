#pragma once
#include <eigen3/Eigen/Core>
#include <vector>
#include <array>
#include <map>
#include <unordered_map>
#include <memory>


namespace NuTo
{

//! @brief transforms a vector of natural node coordinates to a reasonable id
template <int TRaster>
struct NaturalCoordianteToId
{
    //! @brief 'rasterize' the natural coordinates (always [-1, 1]) into MaxId() equidistant blocks
    template <typename TNaturalCoords>
    size_t operator()(const TNaturalCoords& v) const
    {
        int id = 0;
        constexpr float eps = 1.e-2;
        constexpr float factor = 0.5 * TRaster - 2 * eps;
        // reverse the loop to get
        // id = x + 10*y + 100*z
        // instead of (unreversed)
        // id = 100*x + 10*y + z
        // IMO the first version is better to keep the data in 1D/2D closer together
        for (int i = v.size() - 1; i >= 0; --i)
        {
            assert(v[i] > -1. - eps);
            assert(v[i] < 1. + eps);
            id *= TRaster;
            id += factor * (v[i] + 1 - eps); // transformation from float[-1, 1] --> int[0, TRaster]
        }
        return id;
    }

    static constexpr size_t MaxId()
    {
        return TRaster * TRaster * TRaster;
    }
};


//! @brief provides a memoization of the a std::function<TResult(TNaturalCoords)>, mainly used for
//!        element shape functions and their derivatives
//! @param TResult ... return type of the function
//! @param TNaturalCoords ... vector type of natural coordinates
//! @param TIdHash ... function object that transforms TNaturalCoords into a vector id
template <typename TResult, typename TNaturalCoords, typename TIdHash = NaturalCoordianteToId<32>>
class NaturalCoordinateMemoizer
{
public:
    //! @brief ctor
    //! @param function ... function for memoization
    NaturalCoordinateMemoizer(std::function<TResult(TNaturalCoords)> function)
        : mFunction(function)
    {
        mCache.resize(TIdHash::MaxId());
#ifndef NDEBUG
        mCoordinates.resize(TIdHash::MaxId());
#endif
    }

    //! @brief returns the value of the function for the given arguments.
    //!        Repeated calls of Get() with the same arguments will return the memoized
    //!        result.
    //! @param v ... argument
    //! @return reference to the 'memoized' result
    const TResult& Get(const TNaturalCoords& v) const
    {
        size_t id = TIdHash()(v);
        if (mCache[id] == nullptr)
        {
            mCache[id] = std::make_unique<TResult>(mFunction(v));
#ifndef NDEBUG
            mCoordinates[id] = v;
#endif
        }

#ifndef NDEBUG
        if ((v - mCoordinates[id]).norm() > 1.e-15)
            throw;
#endif

        return *mCache[id];
    }

    void ClearCache() const
    {
        for (auto& ptr : mCache)
            ptr = nullptr;
#ifndef NDEBUG
        mCoordinates.clear();
#endif
    }

private:
    //! @brief mutable to mark the Get() method const
    mutable std::vector<std::unique_ptr<TResult>> mCache;

#ifndef NDEBUG
    mutable std::vector<TNaturalCoords> mCoordinates;
#endif

    std::function<TResult(TNaturalCoords)> mFunction;
};


template <typename TVector>
struct CompareVector
{
    bool operator()(const TVector& l, const TVector& r) const
    {
        assert(l.size() == r.size());
        int iDim = 0;
        for (; iDim < l.size() - 1; ++iDim)
            if (l[iDim] != r[iDim])
                return l[iDim] < r[iDim];
        return l[iDim] < r[iDim];
    }
};


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
    //! @brief mutable to mark the Get() method const
    mutable std::map<TNaturalCoords, TResult, TCompare> mCache;
    std::function<TResult(TNaturalCoords)> mFunction;
};


template <typename TVector>
struct VectorHash
{
    size_t operator()(const TVector& v) const
    {
        union {
            double d = 0;
            unsigned long ul;
        } c[3];
        for (int i = 0; i < v.size(); ++i)
            c[i].d = v[i];
        return static_cast<size_t>((3 * c[0].ul) ^ (5 * c[1].ul) ^ (7 * c[2].ul));
    }
};

template <typename TVector>
struct VectorEqual
{
    bool operator()(const TVector& lhs, const TVector& rhs) const
    {
        assert(lhs.size() == rhs.size());
        for (int i = 0; i < lhs.size(); ++i)
            if (lhs[i] != rhs[i])
                return false;
        return true;
    }
};


template <typename TResult, typename TNaturalCoords, typename THash = VectorHash<TNaturalCoords>,
          typename TEqual = VectorEqual<TNaturalCoords>>
class NaturalCoordinateMemoizerUnorderedMap
{
public:
    //! @brief ctor
    //! @param function ... function for memoization
    NaturalCoordinateMemoizerUnorderedMap(std::function<TResult(TNaturalCoords)> function)
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
    //! @brief mutable to mark the Get() method const
    mutable std::unordered_map<TNaturalCoords, TResult, THash, TEqual> mCache;
    std::function<TResult(TNaturalCoords)> mFunction;
};


} /* NuTo */
