/**\file
 * Result viewer result data.
 */
/*
 * Written 2011 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#include "common.h"
#include "Data.h"

#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkPointData.h>

namespace nutogui
{
  
  ResultViewerImpl::Data::Data (vtkDataSet* dataset)
   : dataset (dataset)
  {
    CollectArrays ();
  }
  
  void ResultViewerImpl::Data::CollectArrays ()
  {
    {
      vtkCellData* cellData = dataset->GetCellData();
      for (int i = 0; i < cellData->GetNumberOfArrays(); i++)
      {
	DataArray data;
	data.name = wxString (cellData->GetArrayName (i), wxConvLibc);
	data.assoc = perCell;
	data.arrayIndex = i;
	vtkDataArray* array = cellData->GetArray (i);
	data.arrayComps = array->GetNumberOfComponents ();
	ComputeDisplayDataRange (data, array);
	arrays.push_back (data);
      }
    }
    
    {
      vtkPointData* pointData = dataset->GetPointData();
      for (int i = 0; i < pointData->GetNumberOfArrays(); i++)
      {
	DataArray data;
	data.name = wxString (pointData->GetArrayName (i), wxConvLibc);
	data.assoc = perPoint;
	data.arrayIndex = i;
	vtkDataArray* array = pointData->GetArray (i);
	data.arrayComps = array->GetNumberOfComponents ();
	ComputeDisplayDataRange (data, array);
	arrays.push_back (data);
      }
    }
  }
  
  void ResultViewerImpl::Data::ComputeDisplayDataRange (DataArray& data,
							vtkDataArray* arrayData)
  {
    double range[2];
    arrayData->GetRange (range, -1);
    data.magMin = range[0];
    data.magMax = range[1];
    
    arrayData->GetRange (range, 0);
    for (int comp = 1; comp < arrayData->GetNumberOfComponents(); comp++)
    {
      double compRange[2];
      arrayData->GetRange (compRange, comp);
      range[0] = std::min (range[0], compRange[0]);
      range[1] = std::max (range[1], compRange[1]);
    }
    data.compMin = range[0];
    data.compMax = range[1];
  }
  
  vtkDataSet* ResultViewerImpl::Data::GetDataSet() const
  {
    return dataset;
  }
  
  size_t ResultViewerImpl::Data::GetNumDataArrays () const
  {
    return arrays.size ();
  }
  
  const wxString& ResultViewerImpl::Data::GetDataArrayName (size_t i) const
  {
    return arrays[i].name;
  }
  
  ResultViewerImpl::Data::DataArrayAssociation
  ResultViewerImpl::Data::GetDataArrayAssociation (size_t i) const
  {
    return arrays[i].assoc;
  }
  
  int ResultViewerImpl::Data::GetDataArrayIndex (size_t i) const
  {
    return arrays[i].arrayIndex;
  }
  
  int ResultViewerImpl::Data::GetDataArrayComponents  (size_t i) const
  {
    return arrays[i].arrayComps;
  }

  void ResultViewerImpl::Data::GetDataArrayMagnitudeRange (size_t i, double* range) const
  {
    range[0] = arrays[i].magMin;
    range[1] = arrays[i].magMax;
  }
  
  void ResultViewerImpl::Data::GetDataArrayCompValueRange (size_t i, double* range) const
  {
    range[0] = arrays[i].compMin;
    range[1] = arrays[i].compMax;
  }
  
  bool ResultViewerImpl::Data::IsDataArrayDisplacement (size_t i) const
  {
    return (arrays[i].arrayComps == 3) && (arrays[i].name == wxT ("Displacements"));
  }
  
  vtkDataArray* ResultViewerImpl::Data::GetDataArrayRawData (size_t i) const
  {
    if (arrays[i].assoc == perPoint)
    {
      vtkPointData* pointData = dataset->GetPointData();
      return pointData->GetArray (arrays[i].arrayIndex);
    }
    else
    {
      vtkCellData* cellData = dataset->GetCellData();
      return cellData->GetArray (arrays[i].arrayIndex);
    }
  }
} // namespace nutogui
