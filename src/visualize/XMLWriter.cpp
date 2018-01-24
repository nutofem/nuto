#include "visualize/XMLWriter.h"

#include <fstream>
#include <iomanip>

#include "visualize/UnstructuredGrid.h"
#include "visualize/DataArray.h"
#include "base/Exception.h"

using namespace NuTo;

int ToVtkCellType(eCellTypes type)
{
    switch (type)
    {
    case eCellTypes::VERTEX:
        return 1;
    case eCellTypes::HEXAHEDRON:
        return 12;
    case eCellTypes::LINE:
        return 3;
    case eCellTypes::QUAD:
        return 9;
    case eCellTypes::TETRAEDER:
        return 10;
    case eCellTypes::TRIANGLE:
        return 5;
    case eCellTypes::POLYGON:
        return 7;
    case eCellTypes::WEDGE:
        return 13;
    case eCellTypes::PYRAMID:
        return 14;
    case eCellTypes::LINE2NDORDER:
        return 21;
    case eCellTypes::TRIANGLE2NDORDER:
        return 22;
    case eCellTypes::QUAD2NDORDER:
        return 23;
    case eCellTypes::TETRAEDER2NDORDER:
        return 24;
    case eCellTypes::HEXAHEDRON2NDORDER:
        return 25;
    }
    throw Exception(__PRETTY_FUNCTION__, "Unknown cell type");
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

using DataType = float;

//! @brief extracts point data and cell data into a DataArray object
template <typename TContainer>
Visualize::DataArray<DataType> ExtractData(const TContainer& container, std::string name, int dataIndex)
{
    const int numComponents = container.begin()->GetData(dataIndex).rows();

    std::vector<DataType> data;
    data.reserve(container.size() * numComponents);
    for (const auto& item : container)
    {
        auto dataItem = TransformData(item.GetData(dataIndex));
        assert(numComponents == dataItem.size());
        for (int i = 0; i < dataItem.size(); ++i)
            data.push_back(dataItem[i]);
    }
    return Visualize::DataArray<DataType>(name, numComponents, std::move(data));
}

Visualize::DataArray<DataType> ExtractCoordinates(const std::vector<Visualize::Point>& points)
{
    std::vector<DataType> data;
    data.reserve(points.size() * 3);
    for (const auto& point : points)
    {
        const Eigen::Vector3d coordinates = point.GetCoordinates();
        data.push_back(coordinates[0]);
        data.push_back(coordinates[1]);
        data.push_back(coordinates[2]);
    }
    return Visualize::DataArray<DataType>("points", 3, std::move(data));
}

//! @brief collection of DataArrays for the return type of ExtractCellInfos()
struct CellInfos
{
    Visualize::DataArray<unsigned> connectivity;
    Visualize::DataArray<unsigned> offsets;
    Visualize::DataArray<uint8_t> types;
};

CellInfos ExtractCellInfos(const std::vector<Visualize::Cell>& cells)
{
    std::vector<unsigned> connectivity;
    std::vector<unsigned> offsets;
    std::vector<uint8_t> vtkCellTypes;

    unsigned offset = 0;
    for (const auto& cell : cells)
    {
        offset += cell.GetNumPoints();
        const auto& pointIds = cell.GetPointIds();
        int cellType = ToVtkCellType(cell.GetCellType());

        connectivity.insert(connectivity.end(), pointIds.begin(), pointIds.end());
        offsets.push_back(offset);
        vtkCellTypes.push_back(cellType);
    }

    Visualize::DataArray<unsigned> dataConnectivity("connectivity", 0, std::move(connectivity));
    Visualize::DataArray<unsigned> dataOffsets("offsets", 0, std::move(offsets));
    Visualize::DataArray<uint8_t> dataCellTypes("types", 0, std::move(vtkCellTypes));
    return CellInfos({std::move(dataConnectivity), std::move(dataOffsets), std::move(dataCellTypes)});
}

template <typename TDataTye>
void WriteDataArray(std::ofstream& file, const Visualize::DataArray<TDataTye>& dataArray, bool binary, int* offset)
{
    if (binary)
        dataArray.WriteBinaryHeader(file, offset);
    else
        dataArray.WriteAscii(file);
}

void Visualize::XMLWriter::Export(std::string filename, const UnstructuredGrid& unstructuredGrid, bool binary)
{
    const auto& points = unstructuredGrid.mPoints;
    const auto& cells = unstructuredGrid.mCells;

    std::vector<DataArray<DataType>> pointData;
    for (auto name : unstructuredGrid.mPointDataNames)
        pointData.push_back(ExtractData(points, name, unstructuredGrid.GetPointDataIndex(name)));

    std::vector<DataArray<DataType>> cellData;
    for (auto name : unstructuredGrid.mCellDataNames)
        cellData.push_back(ExtractData(cells, name, unstructuredGrid.GetCellDataIndex(name)));

    DataArray<DataType> pointCoordinates = ExtractCoordinates(points);
    CellInfos cellInfos = ExtractCellInfos(cells);

    std::ofstream file(filename);
    if (!file.is_open())
        throw Exception(__PRETTY_FUNCTION__, "Error opening file " + filename);
    // header #########################################################################################################
    file << R"(<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian")";
    if (binary)
        file << R"( header_type="UInt64")";
    file << ">\n <UnstructuredGrid>\n";
    file << "  <Piece NumberOfPoints=\"" << points.size() << "\" NumberOfCells=\"" << cells.size() << "\">\n";
    int offset = 0;
    // point data #####################################################################################################
    file << "   <PointData>\n";
    for (const auto& dataArray : pointData)
        WriteDataArray(file, dataArray, binary, &offset);
    file << "   </PointData>\n";

    // cell data ######################################################################################################
    file << "   <CellData>\n";
    for (const auto& dataArray : cellData)
        WriteDataArray(file, dataArray, binary, &offset);
    file << "   </CellData>\n";

    // points #########################################################################################################
    file << "   <Points>\n";
    WriteDataArray(file, pointCoordinates, binary, &offset);
    file << "   </Points>\n";

    // Cells - ########################################################################################################
    file << "   <Cells>\n";
    WriteDataArray(file, cellInfos.connectivity, binary, &offset);
    WriteDataArray(file, cellInfos.offsets, binary, &offset);
    WriteDataArray(file, cellInfos.types, binary, &offset);
    file << "   </Cells>\n";
    file << "  </Piece>\n";
    file << " </UnstructuredGrid>\n";

    if (binary)
    {
        file << " <AppendedData encoding=\"raw\">\n_";

        for (const auto& dataArray : pointData)
            dataArray.WriteBinaryData(file);
        for (const auto& dataArray : cellData)
            dataArray.WriteBinaryData(file);
        pointCoordinates.WriteBinaryData(file);
        cellInfos.connectivity.WriteBinaryData(file);
        cellInfos.offsets.WriteBinaryData(file);
        cellInfos.types.WriteBinaryData(file);

        file << "\n </AppendedData>\n";
    }
    file << "</VTKFile>\n";
    file.close();
}
