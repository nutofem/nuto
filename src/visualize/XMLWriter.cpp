#include "visualize/XMLWriter.h"

#include <fstream>
#include <iomanip>

#include "visualize/UnstructuredGrid.h"
#include "visualize/DataArray.h"
#include "visualize/VisualizeException.h"


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

using DataType = float;

//! @brief extracts point data and cell data into a DataArray object
template <typename TContainer>
NuTo::Visualize::DataArray<DataType> ExtractData(const TContainer& container, std::string name, int dataIndex)
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
    return NuTo::Visualize::DataArray<DataType>(name, numComponents, std::move(data));
}

NuTo::Visualize::DataArray<DataType> ExtractCoordinates(const std::vector<NuTo::Visualize::Point>& points)
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
    return NuTo::Visualize::DataArray<DataType>("points", 3, std::move(data));
}

//! @brief collection of DataArrays for the return type of ExtractCellInfos()
struct CellInfos
{
    NuTo::Visualize::DataArray<unsigned> connectivity;
    NuTo::Visualize::DataArray<unsigned> offsets;
    NuTo::Visualize::DataArray<uint8_t> types;
};

CellInfos ExtractCellInfos(const std::vector<NuTo::Visualize::Cell>& cells)
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

    NuTo::Visualize::DataArray<unsigned> dataConnectivity("connectivity", 0, std::move(connectivity));
    NuTo::Visualize::DataArray<unsigned> dataOffsets("offsets", 0, std::move(offsets));
    NuTo::Visualize::DataArray<uint8_t> dataCellTypes("types", 0, std::move(vtkCellTypes));
    return CellInfos({std::move(dataConnectivity), std::move(dataOffsets), std::move(dataCellTypes)});
}

template <typename TDataTye>
void WriteDataArray(std::ofstream& file, const NuTo::Visualize::DataArray<TDataTye>& dataArray, bool binary,
                    int* offset)
{
    if (binary)
        dataArray.WriteBinaryHeader(file, offset);
    else
        dataArray.WriteAscii(file);
}

void NuTo::Visualize::XMLWriter::Export(std::string filename, const UnstructuredGrid& unstructuredGrid, bool binary)
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
        throw NuTo::VisualizeException(__PRETTY_FUNCTION__, "Error opening file " + filename);
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
