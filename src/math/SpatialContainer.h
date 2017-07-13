#pragma once

#include <ANN/ANN.h>
#include <eigen3/Eigen/Core>

#include <set>

#include <vector>
#include <memory>

namespace NuTo
{

//! @brief Container that builds a KD-tree data structure using ANN for various neighbor searches
//! @tparam T object type
//! @tparam TDim struct that defines
//!              Eigen::VectorXd operator() (const T&)
//!              to extract the coordinate of T
template <typename T, typename Coordinate>
class SpatialContainer
{
public:
    //! @brief ctor. saves rValues and builds the KD-tree for each coordinate in rValues
    //! @param rValues objects
    SpatialContainer(const std::vector<T>& rValues)
        : mData(rValues)
    {
        assert(not mData.empty());

        const int dim = Coordinate()(mData[0]).rows();
        const int size = static_cast<int>(mData.size());


        mPoints = annAllocPts(size, dim);
        for (int iRow = 0; iRow < size; ++iRow)
        {
            Eigen::VectorXd coordinate = Coordinate()(mData[iRow]);
            for (int iDim = 0; iDim < dim; ++iDim)
                mPoints[iRow][iDim] = coordinate[iDim];
        }
        mTree = std::make_unique<ANNkd_tree>(mPoints, size, dim);
    }

    ~SpatialContainer()
    {
        annDeallocPts(mPoints);
        annClose();
    }


    //! @brief checks if there are entries at rCoordinate in the radius rRadius
    //! @param rCoordinate querry point
    //! @param rRadius nearest neighbor search radius
    //! @return true (at least one coordinate found) or false (nothing found)
    bool HasEntryAtCoordinate(Eigen::VectorXd rCoordinate, double rRadius) const
    {
        ANNpoint querryPoint = rCoordinate.data();
        ANNdist radiusSquared = rRadius * rRadius;
        return mTree->annkFRSearch(querryPoint, radiusSquared, 0) > 0;
    }

    //! @brief finds and returns values that have the same coordinates within the radius rRadius
    //! @param rRadius nearest neighbor search radius
    //! @return outer vector has the size of unique coordinates
    //!         inner vector contains all objects with the same coordinates
    std::vector<std::vector<T>> GetAllDuplicateValues(double rRadius) const
    {
        auto duplicateIds = GetAllDuplicateIDs(rRadius);
        std::vector<std::vector<T>> values;
        values.reserve(duplicateIds.size());

        for (auto& duplicateIdVector : duplicateIds)
        {
            std::vector<T> equalValues;
            equalValues.reserve(duplicateIdVector.size());

            for (int id : duplicateIdVector)
                equalValues.push_back(mData[id]);

            values.push_back(equalValues);
        }
        return values;
    }

    //! @brief completely similar to public function @GetAllDuplicateValues but returns IDs instead of values
    //! @param rRadius nearest neighbor search radius
    //! @return see @GetAllDuplicateValues (with IDs)
    std::vector<std::vector<int>> GetAllDuplicateIDs(double rRadius) const
    {
        // define a set of searchIDs. All points in the neighborhood of a previously searched ID
        // are removed from this set. (to prevent duplicate duplicates. u see?)
        std::set<int> searchIds;
        for (int i = 0; i < static_cast<int>(mData.size()); ++i)
            searchIds.insert(i);

        std::vector<std::vector<int>> duplicates;
        while (not searchIds.empty())
        {
            auto duplicatesAtId = FindIDsWithinRadius(*searchIds.begin(), rRadius);
            duplicates.push_back(duplicatesAtId);

            for (int duplicate : duplicatesAtId)
                searchIds.erase(duplicate); // do not search duplicate ids again
        }
        return duplicates;
    }

    //! @brief finds all ids in the neighborhood of point rIndex
    //! @param rRadius nearest neighbor search radius
    //! @return all ids in the neighborhood of point rIndex
    std::vector<int> FindIDsWithinRadius(int rIndex, double rRadius) const
    {
        ANNpoint querryPoint = mPoints[rIndex];
        ANNdist radiusSquared = rRadius * rRadius;

        // the return value of this function is ALWAYS the number of points in the radiusSquared
        int numPointsInDistance = mTree->annkFRSearch(querryPoint, radiusSquared, 0);

        // tiny hack: ANNidx == int. so allocate std::vector<int> directly
        std::vector<ANNidx> nearestNeighbourIds(numPointsInDistance);
        mTree->annkFRSearch(querryPoint, radiusSquared, numPointsInDistance, nearestNeighbourIds.data());
        return nearestNeighbourIds;
    }

private:
    //! @brief data reference
    const std::vector<T>& mData;

    //! @brief coordinates
    ANNpointArray mPoints;

    //! @brief KD-tree
    std::unique_ptr<ANNkd_tree> mTree;
};

} // namespace NuTo
