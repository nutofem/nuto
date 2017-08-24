
#pragma once

#include <map>
#include <vector>

namespace NuTo
{

template <class T>
// \brief Transforms a std::map<key,value> into a std::map<value,std::vector<key>> where key and value are of type T
class ReverseMap
{
public:
    using vector_t = std::vector<T>;
    using map_t = std::map<T, T>;
    using reverseMap_t = std::map<T, vector_t>;

    ReverseMap() = default;

    ReverseMap(const std::map<T, T>& map)
    {
        addMap(map);
    };

    void addMap(const std::map<T, T>& map)
    {
        for (const auto& pair : map)
        {
            auto ret = mMap.emplace(pair.second, vector_t({pair.first}));

            if (not ret.second)
                mMap[pair.second].push_back(pair.first);
        }
    }

    typename std::map<T, vector_t>::const_iterator cbegin()
    {
        return mMap.begin();
    }

    typename std::map<T, vector_t>::const_iterator cend()
    {
        return mMap.end();
    }


    typename std::map<T, vector_t>::iterator begin()
    {
        return mMap.begin();
    }

    typename std::map<T, vector_t>::iterator end()
    {
        return mMap.end();
    }

    std::vector<T> operator[](T key)
    {
        return mMap[key];
    }

    size_t size()
    {
        return mMap.size();
    }

    auto find(T key)
    {
        return mMap.find(key);
    }

private:
    std::map<T, vector_t> mMap;
};

} // namespace NuTo
