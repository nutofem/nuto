#include <benchmark/benchmark.h>

#include "mechanics/elements/ElementShapeFunctions.h"
#include "math/NaturalCoordinateMemoizer.h"

#include <memory>
#include <vector>
#include <unordered_map>

/*
 * Benchmarks various implementations of the ShapeFunctionMemoization based on
 *  - vector
 *      Obviously, vector is fastest but causes problems based on its 'rasterized' implementation. IMO this requires
 *      that user know exactly how this is implemented to properly identify bugs.
 *  - unordered_map
 *      Requires a not obvious hash function floating point coordinates and is not faster that the much simpler map ...
 *  - map
 *
 */

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
    }

    //! @brief returns the value of the function for the given arguments.
    //!        Repeated calls of Get() with the same arguments will return the memoized
    //!        Eigen::MatrixXd.
    //! @param v ... argument
    //! @return reference to the 'memoized' Eigen::MatrixXd
    const TResult& Get(const TNaturalCoords& v) const
    {
        size_t id = TIdHash()(v);
        if (mCache[id] == nullptr)
            mCache[id] = std::make_unique<TResult>(mFunction(v));
        return *mCache[id];
    }

    void ClearCache() const
    {
        for (auto& ptr : mCache)
            ptr = nullptr;
    }

private:
    //! @brief mutable to mark the Get() method const
    mutable std::vector<std::unique_ptr<TResult>> mCache;
    std::function<TResult(TNaturalCoords)> mFunction;
};


struct Lobatto
{
    std::vector<Eigen::Vector3d> ips;
    Lobatto()
    {
        std::vector<double> ip1D = {-1., -0.654653670707977087, 0., +0.654653670707977087, +1.};
        int num = ip1D.size();
        ips.reserve(num * num * num);
        for (int i = 0; i < num; i++)
            for (int j = 0; j < num; j++)
                for (int k = 0; k < num; k++)
                    ips.push_back({ip1D[i], ip1D[j], ip1D[k]});
    }
};

auto testFunction = NuTo::ShapeFunctions3D::DerivativeShapeFunctionsTetrahedronOrder2;

template <typename TMemoizer>
void Run(benchmark::State& state)
{
    TMemoizer memo(testFunction);
    Lobatto l;
    for (auto _ : state)
        for (const auto& ip : l.ips)
            memo.Get(ip);
}

template <typename TVector>
struct VectorHash
{
    size_t operator()(const TVector& v) const
    {
        // taken from http://www.beosil.com/download/CollisionDetectionHashing_VMV03.pdf
        static std::array<size_t, 3> magicNumbers{73856093, 19349663, 83492791};

        union Binary {
            double d;
            unsigned long ul;
        };

        size_t vectorHash = 0;
        for (int i = 0; i < v.size(); ++i)
            vectorHash ^= magicNumbers[i] * Binary({v[i]}).ul;

        return vectorHash;
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
    //!        Eigen::MatrixXd.
    //! @param v ... argument
    //! @return reference to the 'memoized' Eigen::MatrixXd
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

BENCHMARK_TEMPLATE(Run, NaturalCoordinateMemoizer<Eigen::MatrixXd, Eigen::Vector3d>);
BENCHMARK_TEMPLATE(Run, NuTo::NaturalCoordinateMemoizerMap<Eigen::MatrixXd, Eigen::Vector3d>);
BENCHMARK_TEMPLATE(Run, NaturalCoordinateMemoizerUnorderedMap<Eigen::MatrixXd, Eigen::Vector3d>);

BENCHMARK_MAIN();
