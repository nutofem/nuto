#include <fstream>
#include <iomanip>

#include "visualize/XMLWriter.h"
#include "visualize/UnstructuredGrid.h"
#include "visualize/VisualizeException.h"

#include "visualize/base64.h"

int ToVtkCellType(NuTo::eCellTypes type)
{
    switch (type)
    {
    case NuTo::eCellTypes::VERTEX:
        return 1;
    case NuTo::eCellTypes::HEXAHEDRON:
        return 12;
    case NuTo::eCellTypes::LINE:
        return 3;
    case NuTo::eCellTypes::QUAD:
        return 9;
    case NuTo::eCellTypes::TETRAEDER:
        return 10;
    case NuTo::eCellTypes::TRIANGLE:
        return 5;
    case NuTo::eCellTypes::POLYGON:
        return 7;
    }
}

Eigen::VectorXd TransformData(Eigen::VectorXd data)
{
    if (data.rows() == 6)
    {
        //                     0  1  2  3  4  5
        // NuTo voigt format: xx yy zz yz xz xy
        std::swap(data[3], data[5]);
        // swap to:           xx yy zz xy xz yz
        std::swap(data[4], data[5]);
        // VTK voigt format:  xx yy zz xy yz xz
    }
    return data;
}

#ifdef HAVE_ZLIB
#include <zlib.h>
extern "C" {
int compress(Bytef* dest, uLongf* destLen, const Bytef* source, uLong sourceLen);
}

template <typename T>
std::vector<unsigned char> compress_data(const std::vector<T>& data)
{
    // Compute length of uncompressed data
    const unsigned long uncompressed_size = data.size() * sizeof(T);

    // Compute maximum length of compressed data
    unsigned long compressed_size = (uncompressed_size + (((uncompressed_size) / 1000) + 1) + 12);

    // Allocate space for compressed data
    std::vector<unsigned char> compressed_data(compressed_size);

    // Compress data
    if (compress((Bytef*)compressed_data.data(), &compressed_size, (const Bytef*)data.data(), uncompressed_size) !=
        Z_OK)
    {
        NuTo::VisualizeException(__PRETTY_FUNCTION__, "Zlib error while compressing data");
    }

    compressed_data.resize(compressed_size);

    return compressed_data;
}
#endif

template <typename T>
std::string encode_base64(const std::vector<T>& data)
{
#ifdef HAVE_ZLIB
    std::uint32_t header[4];
    header[0] = 1;
    header[1] = data.size() * sizeof(T);
    header[2] = 0;

    std::vector<unsigned char> compressed_data = compress_data(data);

    header[3] = compressed_data.size();

    auto prefix = base64_encode((const unsigned char*)&header[0], 4 * sizeof(std::uint32_t));
    auto dataString = base64_encode(compressed_data.data(), compressed_data.size());
#else
    const std::uint32_t size = data.size() * sizeof(T);
    auto prefix = base64_encode((const unsigned char*)&size, sizeof(std::uint32_t));
    auto dataString = base64_encode((const unsigned char*)&data[0], size);
#endif
    return prefix + dataString;
}


void WriteDataLine(std::ofstream& file, const Eigen::VectorXd& vec)
{
    for (int i = 0; i < vec.rows(); ++i)
        file << " " << vec[i];
    file << '\n';
}


void NuTo::Visualize::XMLWriter::Export(std::string filename, const UnstructuredGrid& unstructuredGrid, bool binary)
{
    const auto& points = unstructuredGrid.mPoints;
    const auto& cells = unstructuredGrid.mCells;

    const auto format = binary ? std::quoted("binary") : std::quoted("ascii");
    std::ofstream file(filename);
    if (!file.is_open())
        throw NuTo::VisualizeException(__PRETTY_FUNCTION__, "Error opening file " + filename);

    // header /////////////////////////////////////////////////////////////////
    file << "<VTKFile type=" << std::quoted("UnstructuredGrid") << " version=" << std::quoted("0.1")
         << " byte_order=" << std::quoted("LittleEndian");
#ifdef HAVE_ZLIB
    file << " compressor=" << std::quoted("vtkZLibDataCompressor");
#endif
    file << ">\n";
    file << "  <UnstructuredGrid>\n";
    file << "    <Piece NumberOfPoints=\"" << points.size() << "\" NumberOfCells=\"" << cells.size() << "\">\n";
    ///////////////////////////////////////////////////////////////////////////


    // point data /////////////////////////////////////////////////////////////////
    file << "      <PointData>\n";

    for (auto pointDataName : unstructuredGrid.mPointDataNames)
    {
        const int index = unstructuredGrid.GetPointDataIndex(pointDataName);
        const int numComponents = points[0].GetData(index).rows();
        file << "        <DataArray type=\"Float64\" Name=\"" << pointDataName << "\" NumberOfComponents=\""
             << numComponents << "\" format=" << format << ">\n";
        if (binary)
        {
            std::vector<double> allPointData;
            allPointData.reserve(points.size() * numComponents);
            for (const auto& point : points)
            {
                const Eigen::VectorXd& pointData = TransformData(point.GetData(index));
                assert(pointData.rows() == numComponents);
                for (int i = 0; i < pointData.size(); ++i)
                    allPointData.push_back(pointData[i]);
            }
            file << encode_base64(allPointData) << "\n";
        }
        else
        {
            for (const auto& point : points)
            {
                const Eigen::VectorXd& pointData = TransformData(point.GetData(index));
                assert(pointData.rows() == numComponents);
                WriteDataLine(file, pointData);
            }
        }
        file << "        </DataArray>\n";
    }
    file << "      </PointData>\n";

    // cell data /////////////////////////////////////////////////////////////////
    file << "      <CellData>\n";

    for (auto cellDataName : unstructuredGrid.mCellDataNames)
    {
        const int index = unstructuredGrid.GetCellDataIndex(cellDataName);
        const int numComponents = cells[0].GetData(index).rows();
        file << "        <DataArray type=\"Float64\" Name=\"" << cellDataName << "\" NumberOfComponents=\""
             << numComponents << "\" format=" << format << ">\n";
        if (binary)
        {
            std::vector<double> allCellData;
            allCellData.reserve(cells.size() * numComponents);
            for (const auto& cell : cells)
            {
                const Eigen::VectorXd& cellData = TransformData(cell.GetData(index));
                assert(cellData.rows() == numComponents);
                for (int i = 0; i < cellData.size(); ++i)
                    allCellData.push_back(cellData[i]);
            }
            file << encode_base64(allCellData) << "\n";
        }
        else
        {
            for (const auto& cell : cells)
            {
                const Eigen::VectorXd& cellData = TransformData(cell.GetData(index));
                assert(cellData.rows() == numComponents);
                WriteDataLine(file, cellData);
            }
        }
        file << "        </DataArray>\n";
    }
    file << "      </CellData>\n";

    // points /////////////////////////////////////////////////////////////////
    file << "      <Points>\n";
    file << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=" << format << ">\n";
    if (binary)
    {
        std::vector<double> allCoordinates;
        allCoordinates.reserve(points.size() * points[0].GetCoordinates().size());
        for (const auto& point : points)
        {
            for (int i = 0; i < point.GetCoordinates().size(); ++i)
                allCoordinates.push_back(point.GetCoordinates()[i]);
        }
        file << encode_base64(allCoordinates) << "\n";
    }
    else
    {
        for (const auto& point : points)
        {
            WriteDataLine(file, point.GetCoordinates());
        }
    }
    file << "        </DataArray>\n";
    file << "      </Points>\n";
    // end points /////////////////////////////////////////////////////////////////////

    // Cells - connectivity ///////////////////////////////////////////////////////////
    file << "      <Cells>\n";
    file << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=" << format << ">\n";
    if (binary)
    {
        std::vector<int> allCellPointIds;
        for (const auto& cell : cells)
        {
            for (int id : cell.GetPointIds())
                allCellPointIds.push_back(id);
        }
        file << encode_base64(allCellPointIds) << "\n";
    }
    else
    {
        for (const auto& cell : cells)
        {
            for (int id : cell.GetPointIds())
                file << " " << id;
            file << '\n';
        }
    }
    file << "        </DataArray>\n";

    // Cells - offsets //////////////////////////////////////////////////////////////
    file << "        <DataArray type=\"Int32\" Name=\"offsets\" format=" << format << ">\n";
    if (binary)
    {
        std::vector<int> offsets;
        int offset = 0;
        for (const auto& cell : cells)
        {
            offset += cell.GetNumPoints();
            offsets.push_back(offset);
        }
        file << encode_base64(offsets) << "\n";
    }
    else
    {
        int offset = 0;
        for (const auto& cell : cells)
            file << " " << (offset += cell.GetNumPoints()) << "\n";
    }
    file << "        </DataArray>\n";

    // Cells - types //////////////////////////////////////////////////////////////
    file << "        <DataArray type=\"Int32\" Name=\"types\" format=" << format << ">\n";
    if (binary)
    {
        std::vector<int> cellTypes;
        cellTypes.reserve(cells.size());
        for (const auto& cell : cells)
            cellTypes.push_back(ToVtkCellType(cell.GetCellType()));
        file << encode_base64(cellTypes) << "\n";
    }
    else
    {
        for (const auto& cell : cells)
            file << ToVtkCellType(cell.GetCellType()) << '\n';
    }
    file << "        </DataArray>\n";
    file << "      </Cells>\n";
    // end cells /////////////////////////////////////////////////////////////////////////

    // end header /////////////////////////////////////////////////////////////////
    file << "    </Piece>\n";
    file << "  </UnstructuredGrid>\n";
    file << "</VTKFile>\n";

    file.close();
}
