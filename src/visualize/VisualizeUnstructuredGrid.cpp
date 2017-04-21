// $Id$
#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <boost/tokenizer.hpp>

#include "visualize/CellHexahedron.h"
#include "visualize/CellVertex.h"
#include "visualize/CellLine.h"
#include "visualize/CellPyramid.h"
#include "visualize/CellQuad.h"
#include "visualize/CellTetra.h"
#include "visualize/CellTriangle.h"
#include "visualize/Point.h"
#include "visualize/VisualizeDataBase.h"
#include "visualize/VisualizeDataType.h"
#include "visualize/VisualizeException.h"
#include "visualize/VisualizeUnstructuredGrid.h"



NuTo::VisualizeUnstructuredGrid::VisualizeUnstructuredGrid()
{}

NuTo::VisualizeUnstructuredGrid::~VisualizeUnstructuredGrid()
{}

// add point
unsigned int NuTo::VisualizeUnstructuredGrid::AddPoint(const double* rCoordinates)
{
    this->mPoints.push_back(new Point(rCoordinates));
    unsigned int PointId = this->mPoints.size() - 1;
    for (unsigned int PointDataCount = 0; PointDataCount < this->mPointData.size(); PointDataCount++)
    {
        switch (this->mPointData[PointDataCount].GetDataType())
        {
        case NuTo::eVisualizeDataType::SCALAR:
            this->mPoints[PointId].AddDataScalar(PointDataCount);
            break;
        case NuTo::eVisualizeDataType::VECTOR:
            this->mPoints[PointId].AddDataVector(PointDataCount);
            break;
        case NuTo::eVisualizeDataType::TENSOR:
            this->mPoints[PointId].AddDataTensor(PointDataCount);
            break;
        case NuTo::eVisualizeDataType::FIELD:
            this->mPoints[PointId].AddDataField(PointDataCount, this->mPointData[PointDataCount].GetNumData());
            break;
        default:
            throw NuTo::VisualizeException("[NuTo::VisualizeUnstructuredGrid::AddPoint] Unsupported data type.");
        }
    }
    return PointId;
}

// add vertex cell
unsigned int NuTo::VisualizeUnstructuredGrid::AddVertexCell(const unsigned int* rPoints)
{
    this->CheckPoints(1,rPoints);
    this->mCells.push_back(new CellVertex(rPoints,this->mCellData));
    return this->mCells.size() - 1;
}

// add line cell
unsigned int NuTo::VisualizeUnstructuredGrid::AddLineCell(const unsigned int* rPoints)
{
    this->CheckPoints(2,rPoints);
    this->mCells.push_back(new CellLine(rPoints,this->mCellData));
    return this->mCells.size() - 1;
}

// add triangle cell
unsigned int NuTo::VisualizeUnstructuredGrid::AddTriangleCell(const unsigned int* rPoints)
{
    this->CheckPoints(3,rPoints);
    this->mCells.push_back(new CellTriangle(rPoints,this->mCellData));
    return this->mCells.size() - 1;
}

// add quadrilateral cell
unsigned int NuTo::VisualizeUnstructuredGrid::AddQuadCell(const unsigned int* rPoints)
{
    this->CheckPoints(4,rPoints);
    this->mCells.push_back(new CellQuad(rPoints,this->mCellData));
    return this->mCells.size() - 1;
}

// add tetraeder cell
unsigned int NuTo::VisualizeUnstructuredGrid::AddTetraCell(const unsigned int* rPoints)
{
    this->CheckPoints(4,rPoints);
    this->mCells.push_back(new CellTetra(rPoints,this->mCellData));
    return this->mCells.size() - 1;
}

// add pyramid cell
unsigned int NuTo::VisualizeUnstructuredGrid::AddPyramidCell(const unsigned int* rPoints)
{
    this->CheckPoints(5,rPoints);
    this->mCells.push_back(new CellPyramid(rPoints,this->mCellData));
    return this->mCells.size() - 1;
}

// add hexahedron cell
unsigned int NuTo::VisualizeUnstructuredGrid::AddHexahedronCell(const unsigned int* rPoints)
{
    this->CheckPoints(8,rPoints);
    this->mCells.push_back(new CellHexahedron(rPoints,this->mCellData));
    return this->mCells.size() - 1;
}

// export to Vtu Datafile (XML based file format)
void NuTo::VisualizeUnstructuredGrid::ExportVtuDataFile(const std::string& rFilename) const
{
    std::ofstream file(rFilename.c_str());
    if (!file.is_open())
    {
    	throw NuTo::VisualizeException(std::string("[NuTo::VisualizeUnstructuredGrid::ExportVtuDatafile] Error opening file ")+rFilename.c_str());
    }
    // header /////////////////////////////////////////////////////////////////
    file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    file << "  <UnstructuredGrid>\n";
    file << "    <Piece NumberOfPoints=\"" << this->mPoints.size() << "\" NumberOfCells=\"" << this->mCells.size() << "\">\n";
    ///////////////////////////////////////////////////////////////////////////


    // point data /////////////////////////////////////////////////////////////////
    file << "      <PointData>\n";

    for (unsigned int PointDataCount = 0; PointDataCount < this->mPointData.size(); PointDataCount++)
    {
        file << "        <DataArray type=\"Float32\" Name=\"" << this->mPointData[PointDataCount].GetIdent() << "\" NumberOfComponents=\"" << this->mPointData[PointDataCount].GetNumData() << "\" Format=\"ascii\">\n";
        boost::ptr_vector<Point>::const_iterator PointIter = this->mPoints.begin();
        while (PointIter != this->mPoints.end())
        {
            const VisualizeDataBase* tmpData = PointIter->GetData(PointDataCount);
            if (tmpData->GetDataType() != this->mPointData[PointDataCount].GetDataType())
            {
                throw NuTo::VisualizeException("[NuTo::VisualizeUnstructuredGrid::ExportVtkDatafile] mismatch in point data.");
            }
            file <<"          " << *tmpData << '\n';
            PointIter++;
        }
        file << "        </DataArray>\n";
    }

    file << "      </PointData>\n";

    // cell data /////////////////////////////////////////////////////////////////
    file << "      <CellData>\n";

    for (unsigned int CellDataCount = 0; CellDataCount < this->mCellData.size(); CellDataCount++)
    {
        file << "        <DataArray type=\"Float32\" Name=\"" << this->mCellData[CellDataCount].GetIdent() << "\" NumberOfComponents=\"" << this->mCellData[CellDataCount].GetNumData() << "\" Format=\"ascii\">\n";
        boost::ptr_vector<CellBase>::const_iterator CellIter = this->mCells.begin();
        while (CellIter != this->mCells.end())
        {
            const VisualizeDataBase* tmpData = CellIter->GetData(CellDataCount);
            if (tmpData->GetDataType() != this->mCellData[CellDataCount].GetDataType())
            {
                throw NuTo::VisualizeException("[NuTo::VisualizeUnstructuredGrid::ExportVtkDatafile] mismatch in Cell data.");
            }
            file <<"          " << *tmpData << '\n';
            CellIter++;
        }
        file << "        </DataArray>\n";
    }

    file << "      </CellData>\n";

    // points /////////////////////////////////////////////////////////////////
    // modify format for point coordinate output
    std::ios_base::fmtflags OriginalFlags = file.flags(); // store original format
    file.setf(std::ios_base::uppercase | std::ios_base::scientific | std::ios_base::showpos);

    file << "      <Points>\n";
    file << "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" Format=\"ascii\">\n";
    // write point coordinates to file
    boost::ptr_vector<Point>::const_iterator PointIter = this->mPoints.begin();
    while (PointIter != this->mPoints.end())
    {
        const double* PointCoordinates = PointIter->GetCoordinates();
        file <<"          " << PointCoordinates[0] << " " << PointCoordinates[1] << " " << PointCoordinates[2] << " \n";
        PointIter++;
    }
    // reset original format
    file.flags(OriginalFlags);
    file << "        </DataArray>\n";
    file << "      </Points>\n";
    // end points /////////////////////////////////////////////////////////////////////


    // Cells //////////////////////////////////////////////////////////////////
    file << "      <Cells>\n";
    file << "        <DataArray type=\"Int32\" Name=\"connectivity\" Format=\"ascii\">\n";
    boost::ptr_vector<CellBase>::const_iterator CellIter = this->mCells.begin();
    while (CellIter != this->mCells.end())
    {
    	unsigned int NumPoints = CellIter->GetNumPoints();
    	const unsigned int *Points = CellIter->GetPoints();
    	file <<"          ";
        for (unsigned int PointCount = 0; PointCount < NumPoints; PointCount++)
        {
            file << " " << Points[PointCount];
        }
        file << '\n';
        CellIter++;
    }
    file << "        </DataArray>\n";
    file << "        <DataArray type=\"Int32\" Name=\"offsets\" Format=\"ascii\">\n";
    CellIter = this->mCells.begin();
    int offset(0);
    while (CellIter != this->mCells.end())
    {
        offset+=CellIter->GetNumPoints();
        file << " " << offset << '\n';
        CellIter++;
    }
    file << "        </DataArray>\n";
    file << "        <DataArray type=\"Int32\" Name=\"types\" Format=\"ascii\">\n";
    CellIter = this->mCells.begin();
    while (CellIter != this->mCells.end())
    {
        file << CellIter->GetVtkCellType() << '\n';
        CellIter++;
    }
    file << "        </DataArray>\n";
    file << "      </Cells>\n";
    //end cells /////////////////////////////////////////////////////////////////////////

    // end header /////////////////////////////////////////////////////////////////
    file << "    </Piece>\n";
    file << "  </UnstructuredGrid>\n";
    file << "</VTKFile>\n";

    file.close();
}


// check points
void NuTo::VisualizeUnstructuredGrid::CheckPoints(const unsigned int rNumPoints, const unsigned int *rPoints) const
{
    for (unsigned int PointCount = 0; PointCount < rNumPoints; PointCount++)
    {
        if (rPoints[PointCount] >= this->mPoints.size())
        {
            throw NuTo::VisualizeException("[NuTo::VisualizeUnstructuredGrid::CheckPoints] Point id is out of range.");
        }
    }
}

// define scalar point data
void NuTo::VisualizeUnstructuredGrid::DefinePointDataScalar(const std::string& rIdent)
{
    this->CheckPointDataIdent(rIdent);
    this->mPointData.push_back(NuTo::VisualizeDataType(rIdent, NuTo::eVisualizeDataType::SCALAR));
    unsigned int DataIndex = this->mPointData.size() - 1;
    // add data to all points
    boost::ptr_vector<Point>::iterator PointIter = this->mPoints.begin();
    while (PointIter != this->mPoints.end())
    {
        PointIter->AddDataScalar(DataIndex);
        PointIter++;
    }
}

// define scalar cell data
void NuTo::VisualizeUnstructuredGrid::DefineCellDataScalar(const std::string& rIdent)
{
    this->CheckCellDataIdent(rIdent);
    this->mCellData.push_back(NuTo::VisualizeDataType(rIdent, NuTo::eVisualizeDataType::SCALAR));
    unsigned int DataIndex = this->mCellData.size() - 1;
    // add data to all cells
    boost::ptr_vector<CellBase>::iterator CellIter = this->mCells.begin();
    while (CellIter != this->mCells.end())
    {
        CellIter->AddDataScalar(DataIndex);
        CellIter++;
    }
}

// define vector point data
void NuTo::VisualizeUnstructuredGrid::DefinePointDataVector(const std::string& rIdent)
{
    this->CheckPointDataIdent(rIdent);
    this->mPointData.push_back(NuTo::VisualizeDataType(rIdent, NuTo::eVisualizeDataType::VECTOR));
    unsigned int DataIndex = this->mPointData.size() - 1;
    // add data to all points
    boost::ptr_vector<Point>::iterator PointIter = this->mPoints.begin();
    while (PointIter != this->mPoints.end())
    {
        PointIter->AddDataVector(DataIndex);
        PointIter++;
    }
}

// define vector cell data
void NuTo::VisualizeUnstructuredGrid::DefineCellDataVector(const std::string& rIdent)
{
    this->CheckCellDataIdent(rIdent);
    this->mCellData.push_back(NuTo::VisualizeDataType(rIdent, NuTo::eVisualizeDataType::VECTOR));
    unsigned int DataIndex = this->mCellData.size() - 1;
    // add data to all cells
    boost::ptr_vector<CellBase>::iterator CellIter = this->mCells.begin();
    while (CellIter != this->mCells.end())
    {
        CellIter->AddDataVector(DataIndex);
        CellIter++;
    }
}

// define tensor point data
void NuTo::VisualizeUnstructuredGrid::DefinePointDataTensor(const std::string& rIdent)
{
    this->CheckPointDataIdent(rIdent);
    this->mPointData.push_back(NuTo::VisualizeDataType(rIdent, NuTo::eVisualizeDataType::TENSOR));
    unsigned int DataIndex = this->mPointData.size() - 1;
    // add data to all points
    boost::ptr_vector<Point>::iterator PointIter = this->mPoints.begin();
    while (PointIter != this->mPoints.end())
    {
        PointIter->AddDataTensor(DataIndex);
        PointIter++;
    }
}

// define tensor cell data
void NuTo::VisualizeUnstructuredGrid::DefineCellDataTensor(const std::string& rIdent)
{
    this->CheckCellDataIdent(rIdent);
    this->mCellData.push_back(NuTo::VisualizeDataType(rIdent, NuTo::eVisualizeDataType::TENSOR));
    unsigned int DataIndex = this->mCellData.size() - 1;
    // add data to all cells
    boost::ptr_vector<CellBase>::iterator CellIter = this->mCells.begin();
    while (CellIter != this->mCells.end())
    {
        CellIter->AddDataTensor(DataIndex);
        CellIter++;
    }
}

// define field point data
void NuTo::VisualizeUnstructuredGrid::DefinePointDataField(const std::string& rIdent, unsigned int rNumData)
{
    this->CheckPointDataIdent(rIdent);
    this->mPointData.push_back(NuTo::VisualizeDataType(rIdent, NuTo::eVisualizeDataType::FIELD));
    this->mPointData[this->mPointData.size() - 1].SetNumData(rNumData);
    unsigned int DataIndex = this->mPointData.size() - 1;
    // add data to all points
    boost::ptr_vector<Point>::iterator PointIter = this->mPoints.begin();
    while (PointIter != this->mPoints.end())
    {
        PointIter->AddDataField(DataIndex, rNumData);
        PointIter++;
    }
}

// define field cell data
void NuTo::VisualizeUnstructuredGrid::DefineCellDataField(const std::string& rIdent, unsigned int rNumData)
{
    this->CheckCellDataIdent(rIdent);
    this->mCellData.push_back(NuTo::VisualizeDataType(rIdent, NuTo::eVisualizeDataType::FIELD));
    this->mCellData[this->mCellData.size() - 1].SetNumData(rNumData);
    unsigned int DataIndex = this->mCellData.size() - 1;
    // add data to all points
    boost::ptr_vector<CellBase>::iterator CellIter = this->mCells.begin();
    while (CellIter != this->mCells.end())
    {
        CellIter->AddDataField(DataIndex, rNumData);
        CellIter++;
    }
}

// check point data identifier
void NuTo::VisualizeUnstructuredGrid::CheckPointDataIdent(const std::string& rIdent) const
{
    this->CheckDataIdent(rIdent);
    std::vector<VisualizeDataType>::const_iterator iter = mPointData.begin();
    while (iter != mPointData.end())
    {
        if (iter->IsIdent(rIdent))
        {
            break;
        }
        iter++;
    }
    if (iter != mPointData.end())
    {
        throw NuTo::VisualizeException("[NuTo::VisualizeUnstructuredGrid::CheckPointDataIdent] data identifier already exist for point data.");
    }
}

// check cell data identifier
void NuTo::VisualizeUnstructuredGrid::CheckCellDataIdent(const std::string& rIdent) const
{
    this->CheckDataIdent(rIdent);
    std::vector<VisualizeDataType>::const_iterator iter = mCellData.begin();
    while (iter != mCellData.end())
    {
        if (iter->IsIdent(rIdent))
        {
            break;
        }
        iter++;
    }
    if (iter != mCellData.end())
    {
        throw NuTo::VisualizeException("[NuTo::VisualizeUnstructuredGrid::CheckCellDataIdent] data identifier already exist for cell data.");
    }
}

// check strings
void NuTo::VisualizeUnstructuredGrid::CheckDataIdent(const std::string& rIdent) const
{
    std::vector<std::string> tokens;
    std::istringstream iss(rIdent);
    std::copy(std::istream_iterator<std::string>(iss), std::istream_iterator<std::string>(), std::back_inserter<std::vector<std::string> >(tokens));
    if (tokens.size() != 1)
    {
        //std::cout << "data ident : " << rIdent << '\n';
    	throw NuTo::VisualizeException("[NuTo::VisualizeUnstructuredGrid::CheckDataIdent] data identifier must be a single word.");
    }
}

// set scalar point data
void NuTo::VisualizeUnstructuredGrid::SetPointDataScalar(unsigned int rPointIndex, const std::string& rDataIdent, double rData)
{
    if (rPointIndex >= this->mPoints.size())
    {
        throw NuTo::VisualizeException("[NuTo::VisualizeUnstructuredGrid::SetPointDataScalar] invalid point index.");
    }
    unsigned int PointDataIndex = this->GetPointDataIndex(rDataIdent);
    this->mPoints[rPointIndex].SetDataScalar(PointDataIndex, rData);
}

// set vector point data
void NuTo::VisualizeUnstructuredGrid::SetPointDataVector(unsigned int rPointIndex, const std::string& rDataIdent, double rData[3])
{
    if (rPointIndex >= this->mPoints.size())
    {
        throw NuTo::VisualizeException("[NuTo::VisualizeUnstructuredGrid::SetPointDataVector] invalid point index.");
    }
    unsigned int PointDataIndex = this->GetPointDataIndex(rDataIdent);
    this->mPoints[rPointIndex].SetDataVector(PointDataIndex, rData);
}

// set scalar cell data
void NuTo::VisualizeUnstructuredGrid::SetCellDataScalar(unsigned int rCellIndex, const std::string& rDataIdent, double rData)
{
    if (rCellIndex >= this->mCells.size())
    {
        throw NuTo::VisualizeException("[NuTo::VisualizeUnstructuredGrid::SetCellDataScalar] invalid cell index.");
    }
    unsigned int CellDataIndex = this->GetCellDataIndex(rDataIdent);
    this->mCells[rCellIndex].SetDataScalar(CellDataIndex, rData);
}

// set vector cell data
void NuTo::VisualizeUnstructuredGrid::SetCellDataVector(unsigned int rCellIndex, const std::string& rDataIdent, double rData[3])
{
    if (rCellIndex >= this->mCells.size())
    {
        throw NuTo::VisualizeException("[NuTo::VisualizeUnstructuredGrid::SetCellDataTensor] invalid cell index.");
    }
    unsigned int CellDataIndex = this->GetCellDataIndex(rDataIdent);
    this->mCells[rCellIndex].SetDataVector(CellDataIndex, rData);
}

// set tensor cell data
void NuTo::VisualizeUnstructuredGrid::SetCellDataTensor(unsigned int rCellIndex, const std::string& rDataIdent, double rData[9])
{
    if (rCellIndex >= this->mCells.size())
    {
        throw NuTo::VisualizeException("[NuTo::VisualizeUnstructuredGrid::SetCellDataTensor] invalid cell index.");
    }
    unsigned int CellDataIndex = this->GetCellDataIndex(rDataIdent);
    this->mCells[rCellIndex].SetDataTensor(CellDataIndex, rData);
}

// get point data index
unsigned int NuTo::VisualizeUnstructuredGrid::GetPointDataIndex(const std::string& rIdent) const
{
    for (unsigned int PointDataCount = 0; PointDataCount < this->mPointData.size(); PointDataCount++)
    {
        if (this->mPointData[PointDataCount].IsIdent(rIdent))
        {
            return PointDataCount;
        }
    }
    std::cout << rIdent << '\n';
    throw NuTo::VisualizeException("[NuTo::VisualizeUnstructuredGrid::GetPointDataIndex] data identifier not found.");
}

// get cell data index
unsigned int NuTo::VisualizeUnstructuredGrid::GetCellDataIndex(const std::string& rIdent) const
{
    for (unsigned int CellDataCount = 0; CellDataCount < this->mCellData.size(); CellDataCount++)
    {
        if (this->mCellData[CellDataCount].IsIdent(rIdent))
        {
            return CellDataCount;
        }
    }
    std::cout << rIdent << '\n';
    throw NuTo::VisualizeException("[NuTo::VisualizeUnstructuredGrid::GetCellDataIndex] data identifier not found.");
}
