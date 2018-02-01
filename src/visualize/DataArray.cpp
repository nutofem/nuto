#include "visualize/DataArray.h"
#include <iostream>
#include <fstream>
#include <typeinfo>

using namespace NuTo::Visualize;

template <typename TDataType>
DataArray<TDataType>::DataArray(std::string name, int numComponents, std::vector<TDataType>&& data)
    : mName(name)
    , mNumComponents(numComponents)
    , mData(std::move(data))
{
}

namespace NuTo
{
namespace Visualize
{
template <>
std::string DataArray<double>::VtkTypeString() const
{
    static_assert(sizeof(double) == 8, "unexpected type size");
    return "Float64";
}

template <>
std::string DataArray<float>::VtkTypeString() const
{
    static_assert(sizeof(float) == 4, "unexpected type size");
    return "Float32";
}

template <>
std::string DataArray<unsigned>::VtkTypeString() const
{
    static_assert(sizeof(unsigned) == 4, "unexpected type size");
    return "UInt32";
}

template <>
std::string DataArray<uint8_t>::VtkTypeString() const
{
    static_assert(sizeof(uint8_t) == 1, "unexpected type size");
    return "UInt8";
}

} // namespace Visualize
} // namespace NuTo

template <typename TDataType>
int DataArray<TDataType>::Offset() const
{
    return 8 + sizeof(TDataType) * mData.size();
}

template <typename TDataType>
void DataArray<TDataType>::WriteCommonHeader(std::ostream& file) const
{
    file << "<DataArray ";
    file << "type=\"" << VtkTypeString() << "\" ";
    file << "Name=\"" << mName << "\" ";
    if (mNumComponents != 0)
        file << " NumberOfComponents=\"" << mNumComponents << "\" ";
}

template <typename TDataType>
void DataArray<TDataType>::WriteAscii(std::ostream& file) const
{
    WriteCommonHeader(file);
    file << "format=\""
         << "ascii"
         << "\" ";
    file << ">\n";

    for (auto data : mData)
        // https://stackoverflow.com/questions/14644716/how-to-output-a-character-as-an-integer-through-cout
        if (typeid(TDataType) == typeid(uint8_t))
            file << " " << static_cast<int>(data);
        else
            file << " " << data;

    file << "\n</DataArray>";
}

template <typename TDataType>
void DataArray<TDataType>::WriteBinaryHeader(std::ostream& file, int* offset) const
{
    WriteCommonHeader(file);
    file << "format=\""
         << "appended"
         << "\" ";
    file << "offset=\"" << *offset << "\" ";
    file << "/>\n";
    *offset += Offset();
}

template <typename TDataType>
void DataArray<TDataType>::WriteBinaryData(std::ofstream& file) const
{
    const uint64_t size = Offset() - 8;
    file.write(reinterpret_cast<const char*>(&size), sizeof(uint64_t));
    file.write(reinterpret_cast<const char*>(mData.data()), mData.size() * sizeof(TDataType));
}

template class NuTo::Visualize::DataArray<double>;
template class NuTo::Visualize::DataArray<float>;
template class NuTo::Visualize::DataArray<unsigned>;
template class NuTo::Visualize::DataArray<uint8_t>;
