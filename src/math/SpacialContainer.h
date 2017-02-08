#pragma once

#include <ANN/ANN.h>
#include <eigen3/Eigen/Dense>

#include <set>

#include <vector>


namespace NuTo
{

//! @brief
//! @tparam T
//! @tparam TDim
template <typename T, int TDim>
class SpacialContainer
{
public:
    SpacialContainer(const std::vector<Eigen::Matrix<double, TDim, 1>>& rCoordinates, const std::vector<T>& rValues)
        : mData(rValues), mSize(static_cast<int>(mData.size()))
    {
        assert(rCoordinates.size() == rValues.size());
        mPoints = annAllocPts(mSize, TDim);
        for (int iRow = 0; iRow < mSize; ++iRow)
            for (int iDim = 0; iDim < TDim; ++iDim)
                mPoints[iRow][iDim] = rCoordinates[iRow][iDim];

        mTree = new ANNkd_tree(mPoints, mSize, TDim);
    }

    ~SpacialContainer()
    {
        delete mTree;
        annDeallocPts(mPoints);
        annClose();
    }

    std::vector<std::vector<T>> GetAllDuplicateValues(double rTolerance)
    {
        auto duplicateIds = GetAllDuplicateIDs(rTolerance);
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


    std::vector<std::vector<int>> GetAllDuplicateIDs(double rTolerance)
    {
        std::vector<std::vector<int>> duplicates;

        std::set<int> searchIds;
        for (int i = 0; i < mSize; ++i)
            searchIds.insert(i);

        while (not searchIds.empty())
        {
            auto duplicatesAtId = GetDuplicateIDs(*searchIds.begin(), rTolerance);
            duplicates.push_back(duplicatesAtId);

            for (int duplicate : duplicatesAtId)
                searchIds.erase(duplicate); // do not search duplicate ids again
        }
        return duplicates;
    }

    std::vector<int> GetDuplicateIDs(int rIndex, double rTolerance)
    {
        ANNpoint querryPoint = mPoints[rIndex];
        ANNdist distanceSquared = rTolerance * rTolerance;

        // the return value of this function is ALWAYS the number of points in the distanceSquared
        int numPointsInDistance = mTree->annkFRSearch(querryPoint, distanceSquared, 0);

        // tiny hack: ANNidx == int. so allocate std::vector<int> directly
        std::vector<ANNidx> nearestNeighbourIds(numPointsInDistance);
        mTree->annkFRSearch(querryPoint, distanceSquared, numPointsInDistance, nearestNeighbourIds.data());
        return nearestNeighbourIds;
    }



private:

    const std::vector<T>& mData;

    ANNpointArray mPoints;

    ANNkd_tree* mTree;

    int mSize; // annoying "unsigned comparison...
};

} // namespace NuTo
