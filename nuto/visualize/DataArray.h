#include <string>
#include <vector>
#include <iosfwd>

namespace NuTo
{
namespace Visualize
{
template <typename TDataType>
//! @brief collects data for a vtk <DataArray> and writes it in ascii or raw binary appended format
class DataArray
{
public:
    //! @param name "name" attribute
    //! @param numComponents  "NumberOfComponents" attribute
    //! @param data data as a xvalue. You'll have to std::move() the data in here.
    DataArray(std::string name, int numComponents, std::vector<TDataType>&& data);

    //! @return matching vtk type string to <typename TDataType>
    std::string VtkTypeString() const;

    //! @return calculates the byte offset in the raw binary appended format via
    //! offset = headerSize (usually 8) + #data * sizeof(data)
    int Offset() const;

    //! @brief writes the header for both ascii and binary format including "name", "NumberOfComponents" and "type"
    //! @param file outstream to write into
    void WriteCommonHeader(std::ostream& file) const;

    //! @brief writes to ascii including the common header, "type" and the actual data
    //! @param file outstream to write into
    void WriteAscii(std::ostream& file) const;

    //! @brief writes the binary header including the common header, "type", "offset" and calculates the new offset
    //! @param file outstream to write into
    //! @param offset "offset" attribute that gets updated with the current offset
    //! @remark By using the * syntax it is clear (from the callers side) that offset is changed
    //!         [called with &int]
    void WriteBinaryHeader(std::ostream& file, int* offset) const;

    //! @brief writes the data in uncompressed raw binary
    //! @param file file to write the raw binary
    void WriteBinaryData(std::ofstream& file) const;

private:
    std::string mName;
    int mNumComponents;
    std::vector<TDataType> mData;
};
} /* Visualize */
} /* NuTo */
