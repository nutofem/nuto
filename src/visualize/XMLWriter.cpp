#include "visualize/XMLWriter.h"
#include "visualize/UnstructuredGrid.h"
#include "visualize/VisualizeException.h"

#include <fstream>
#include <iomanip>

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

void WriteDataLine(std::ofstream& file, const Eigen::VectorXd& data)
{
    for (int i = 0; i < data.rows(); ++i)
        file << " " << data[i];
    file << '\n';
}

void NuTo::Visualize::XMLWriter::Export(std::string filename, const UnstructuredGrid& unstructuredGrid, bool asBinary)
{
    const auto& points = unstructuredGrid.mPoints;
    const auto& cells = unstructuredGrid.mCells;

    const auto format = asBinary? std::quoted("ascii") : std::quoted("ascii");

    std::ofstream file(filename);
    if (!file.is_open())
    	throw NuTo::VisualizeException(__PRETTY_FUNCTION__, "Error opening file " + filename);

    // header /////////////////////////////////////////////////////////////////
    file << "<VTKFile type=" << std::quoted("UnstructuredGrid") 
         << " version=" << std::quoted("0.1") 
         << " byte_order=" << std::quoted("LittleEndian") << ">\n";
    file << "  <UnstructuredGrid>\n";
    file << "    <Piece NumberOfPoints=\"" << points.size() << "\" NumberOfCells=\"" << cells.size() << "\">\n";
    ///////////////////////////////////////////////////////////////////////////


    // point data /////////////////////////////////////////////////////////////////
    file << "      <PointData>\n";

    for (auto pointDataName : unstructuredGrid.mPointDataNames)
    {
        const int index = unstructuredGrid.GetPointDataIndex(pointDataName);
        const int numComponents = points[0].GetData(index).rows(); 
        file << "        <DataArray type=\"Float32\" Name=\"" << pointDataName << "\" NumberOfComponents=\"" << numComponents << "\" Format=" << format << ">\n";
        for (const auto& point : points)
        {
            const Eigen::VectorXd& pointData = TransformData(point.GetData(index));
            assert(pointData.rows() == numComponents);
            WriteDataLine(file, pointData);
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
        file << "        <DataArray type=\"Float32\" Name=\"" << cellDataName << "\" NumberOfComponents=\"" << numComponents << "\" Format=\"ascii\">\n";
        for (const auto& cell : cells)
        {
            const Eigen::VectorXd& cellData = TransformData(cell.GetData(index));
            assert(cellData.rows() == numComponents);
            WriteDataLine(file, cellData);
        }
        file << "        </DataArray>\n";
    }
    file << "      </CellData>\n";

    // points /////////////////////////////////////////////////////////////////
    file << "      <Points>\n";
    file << "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" Format=\"ascii\">\n";
    for (const auto& point : points)
        WriteDataLine(file, point.GetCoordinates());
    file << "        </DataArray>\n";
    file << "      </Points>\n";
    // end points /////////////////////////////////////////////////////////////////////

    // Cells - connectivity ///////////////////////////////////////////////////////////
    file << "      <Cells>\n";
    file << "        <DataArray type=\"Int32\" Name=\"connectivity\" Format=\"ascii\">\n";
    for (const auto& cell : cells) 
    {
        for (int id : cell.GetPointIds())
            file << " " << id;
        file << '\n';
    }
    file << "        </DataArray>\n";

    // Cells - offsets //////////////////////////////////////////////////////////////
    file << "        <DataArray type=\"Int32\" Name=\"offsets\" Format=\"ascii\">\n";
    int offset = 0;
    for (const auto& cell : cells)
        file << " " << (offset += cell.GetNumPoints()) << "\n"; 
    file << "        </DataArray>\n";

    // Cells - types //////////////////////////////////////////////////////////////
    file << "        <DataArray type=\"Int32\" Name=\"types\" Format=\"ascii\">\n";
    for (const auto& cell : cells)
        file << ToVtkCellType(cell.GetCellType()) << '\n';
    file << "        </DataArray>\n";
    file << "      </Cells>\n";
    //end cells /////////////////////////////////////////////////////////////////////////

    // end header /////////////////////////////////////////////////////////////////
    file << "    </Piece>\n";
    file << "  </UnstructuredGrid>\n";
    file << "</VTKFile>\n";

    file.close();
}
