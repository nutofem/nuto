/**\file
 * Result viewer result data.
 */
/*
 * Written 2011 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#ifndef __NUTOGUI_RESULTVIEWER_DATA_H__
#define __NUTOGUI_RESULTVIEWER_DATA_H__

#include "ResultViewerImpl.h"

#include <vector>

#include <vtkSmartPointer.h>

class vtkDataArray;
class vtkDataSet;

namespace nutogui
{
  class ResultViewerImpl::Data
  {
    vtkSmartPointer<vtkDataSet> dataset;
    
    struct DataArray;
    std::vector<DataArray> arrays;
    
    void CollectArrays ();
    void ComputeDisplayDataRange (DataArray& data, vtkDataArray* arrayData);
  public:
    Data (vtkDataSet* dataset);
    
    vtkDataSet* GetDataSet() const;
    
    size_t GetNumDataArrays () const;
    const wxString& GetDataArrayName (size_t i) const;
    enum DataArrayAssociation
    {
      perCell, perPoint
    };
    DataArrayAssociation GetDataArrayAssociation (size_t i) const;
    int GetDataArrayIndex (size_t i) const;
    int GetDataArrayComponents  (size_t i) const;

    void GetDataArrayMagnitudeRange (size_t i, double* range) const;
    void GetDataArrayCompValueRange (size_t i, double* range) const;
    
    bool IsDataArrayDisplacement (size_t i) const;
    vtkDataArray* GetDataArrayRawData (size_t i) const;
  private:
    struct DataArray
    {
      wxString name;
      DataArrayAssociation assoc;
      int arrayIndex;
      int arrayComps;
      
      double magMin, magMax;
      double compMin, compMax;
    };
  };
} // namespace nutogui

#endif // __NUTOGUI_RESULTVIEWER_DATA_H__
