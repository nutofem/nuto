#pragma once
#include <eigen3/Eigen/Core>
#include <vector>
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
        constexpr float factor = 0.5 * TRaster;
        constexpr float eps = 1.e-6;
        // reverse the loop to get
        // id = x + 10*y + 100*z
        // instead of (unreversed)
        // id = 100*x + 10*y + z
        // IMO the first version is better to keep the data in 1D/2D closer together
        for (int i = v.size() - 1; i >= 0; --i)
        {
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
template <typename TResult, typename TNaturalCoords, typename TIdHash = NaturalCoordianteToId<16>>
class NaturalCoordinateMemoizer
{
public:

    //! @brief ctor
    //! @param function ... function for memoization
    NaturalCoordinateMemoizer(std::function<TResult(TNaturalCoords)> function)
        : mFunction(function)
    {
        mCache.resize(TIdHash::MaxId());
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
            mCache[id] = std::make_unique<TResult>(mFunction(v));
        return *mCache[id];
    }

private:
    //! @brief mutable to mark the Get() method const
    mutable std::vector<std::unique_ptr<TResult>> mCache;

    std::function<TResult(TNaturalCoords)> mFunction;
};
} /* NuTo */
