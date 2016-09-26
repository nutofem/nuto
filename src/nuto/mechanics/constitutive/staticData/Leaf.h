#pragma once

#include <vector>
#include <typeinfo>
#include "nuto/mechanics/constitutive/staticData/Component.h"
#include "nuto/mechanics/MechanicsException.h"

namespace NuTo
{
namespace Constitutive
{
namespace StaticData
{

template<typename T>
class Leaf : public Component
{
public:
    static Leaf<T>* Create(const T& rInitialValue)
    {
        return new Leaf<T>(rInitialValue);
    }

    virtual Leaf<T>* Clone() const override
    {
        Leaf<T>* newLeaf = new Leaf<T>();
        newLeaf->mData = this->mData;
        return newLeaf;
    }

    bool operator==(const Leaf<T>& rhs) const
    {
        return this->mData == rhs.mData;
    }

    bool operator!=(const Leaf<T>& rhs) const
    {
        return !this->operator==(rhs); 
    }
     
    friend std::ostream& operator<<(std::ostream& os, const Leaf<T>& leaf)
    {
        os << "Static data of type " << typeid(T).name() << ":\n";
        unsigned int i = 0;
        for (auto value : leaf.mData)
        {
            os << "Timestep " << i << ": " << value << "\n";
            i++;
        }
        return os;
    }

    //! @brief Set a new value for the current static data.
    //! @param newData New value for current static data.
    void SetData(T newData)
    {
        if (mData.empty())
            mData.push_back(newData);
        else
            mData[0] = newData;
    }

    //! @brief Get the data at `timestep`.
    //! @param timeStep Timestep at which to retrieve the data. The timestep is optional and defaults to 0, that is,
    //!                 the current time step.
    T& GetData(int timeStep = 0)
    {
        if (timeStep > GetNumData() - 1)
            throw MechanicsException(__PRETTY_FUNCTION__, "You requested time step " + std::to_string(timeStep) + ". Number of allocated time steps: " + std::to_string(GetNumData()));
        return mData.at(timeStep);
    }

    //! @brief Copies the "current" ([0]) static data `numAdditionalData` times.
    void AllocateAdditionalData(int numAdditionalData)
    {
        if (mData.empty())
            throw MechanicsException(__PRETTY_FUNCTION__, "No static data allocated yet.");

        for (int i = 0; i < numAdditionalData; ++i)
        {
            mData.push_back(mData[0]);
        }
    }

    //! @brief Puts current static data to previous static data, previous to pre-previous, etc.
    //! The current data is copied to the previous data, all others are moved.
    void ShiftToPast() override
    {
        if (mData.size() < 2)
            throw MechanicsException(__PRETTY_FUNCTION__, "There need to be at least two time steps allocated.");

        mData.pop_back();
        mData.insert(mData.begin(), mData[0]);
    }

    //! @brief Puts previous static data to current static data, pre-previous to previous, etc.
    void ShiftToFuture() override
    {
        if (mData.size() < 2)
            throw MechanicsException(__PRETTY_FUNCTION__, "There need to be at least two time steps allocated.");
        
        mData.erase(mData.begin());
        mData.push_back(mData.back());
    }

    //! @brief Returns the total number of static data sets.
    int GetNumData() const
    {
        return mData.size();
    }

private:
    Leaf() = default;

    Leaf(const T& rInitialValue)
    {
        mData.clear();
        mData.push_back(rInitialValue);
    }

    std::vector<T> mData;
};

} // namespace StaticData
} // namespace Constitutive
} // namespace NuTo

