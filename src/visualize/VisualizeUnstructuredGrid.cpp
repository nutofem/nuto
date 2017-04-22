// $Id$
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>

#include "visualize/Point.h"
#include "visualize/VisualizeException.h"
#include "visualize/VisualizeUnstructuredGrid.h"


int NuTo::VisualizeUnstructuredGrid::AddPoint(Eigen::Vector3d coordinates)
{
    mPoints.push_back(Point(coordinates, mPointDataNames.size()));
    return mPoints.size() - 1;
}

int NuTo::VisualizeUnstructuredGrid::AddCell(CellBase cell)
{
    mCells.push_back(cell);
    return mCells.size() - 1;
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

void NuTo::VisualizeUnstructuredGrid::DefinePointData(std::string name)
{
    if (std::find(mPointDataNames.begin(), mPointDataNames.end(), name) != mPointDataNamesToId.end())
        throw VisualizeException(__PRETTY_FUNCTION__, "data identifier already exist for point data.");

    if (not mPoints.empty())
        throw VisualizeException(__PRETTY_FUNCTION__, "define all data fields _before_ adding points");

    mPointDataNames.push_back(name);
}

void NuTo::VisualizeUnstructuredGrid::DefineCellData(std::string name)
{
    if (std::find(mCellDataNames.begin(), mCellDataNames.end(), name) != mCellDataNamesToId.end())
        throw NuTo::VisualizeException(__PRETTY_FUNCTION__, "data identifier already exist for point data.");
    
    if (not mCells.empty())
        throw VisualizeException(__PRETTY_FUNCTION__, "define all data fields _before_ adding cells");

    mCellDataNames.push_back(name);
}

void NuTo::VisualizeUnstructuredGrid::SetPointData(int pointIndex, const std::string& name, double data)
{
    SetPointData(pointIndex, name, Eigen::Matrix<double, 1, 1>::Constant(data));
}

void NuTo::VisualizeUnstructuredGrid::SetPointData(int pointIndex, const std::string& name, Eigen::VectorXd data)
{
    mPoints[pointIndex].SetData(GetPointDataIndex(name), data); 
}

void NuTo::VisualizeUnstructuredGrid::SetCellData(int cellIndex, const std::string& name, double data)
{
    SetCellData(cellIndex, name, Eigen::Matrix<double, 1, 1>::Constant(data));
}

void NuTo::VisualizeUnstructuredGrid::SetCellData(int cellIndex, const std::string& name, Eigen::VectorXd data)
{
    mCells[cellIndex].SetData(GetCellDataIndex(name), data); 
}

int NuTo::VisualizeUnstructuredGrid::GetPointDataIndex(const std::string& name) const
{
    auto it = std::find(mPointDataNames.begin(), mPointDataNames.end(), name);
    if (it == mPointDataNamesToId.end())
        throw NuTo::VisualizeException(__PRETTY_FUNCTION__, "data " + name + " not defined");
    return std::distance(mPointDataNames.begin(), it);
}

int NuTo::VisualizeUnstructuredGrid::GetCellDataIndex(const std::string& name) const
{
    auto it = std::find(mCellDataNames.begin(), mCellDataNames.end(), name);
    if (it == mCellDataNamesToId.end())
        throw NuTo::VisualizeException(__PRETTY_FUNCTION__, "data " + name + " not defined");
    return std::distance(mCellDataNames.begin(), it);
}
