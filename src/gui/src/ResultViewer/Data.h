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
    ResultDataSourceVTKPtr dataSource;
    
    struct DataArray;
    std::vector<DataArray> arrays;
    
    void CollectArrays ();
    void ComputeDisplayDataRange (DataArray& data, vtkDataArray* arrayData);
  public:
    Data (const ResultDataSourceVTKPtr& dataSource);
    
    size_t GetDataSetNum() const;
    const wxString& GetDataSetName (size_t index) const;
    vtkDataSet* GetDataSet (size_t index) const;
    
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
    vtkDataArray* GetDataArrayRawData (size_t setIndex, size_t i) const;
    
    /// Special component IDs
    enum
    {
      compMagnitude = -1
    };
    /// Get display name (for human consumption) for an array component
    wxString GetDataArrayComponentDisplayName (size_t i, int comp) const;
    
    /**
     * Query the original cell ID for a cell from the dataset.
     * The dataset must have some “ancestry” in a dataset a Data object.
     * Returns an invalid cell ID otherwise.
     */
    vtkIdType GetOriginalCell (vtkDataSet* dataset, vtkIdType cellID) const;
  private:
    struct DataArray
    {
      wxString name;
      DataArrayAssociation assoc;
      int arrayIndex;
      int arrayComps;
      
      double magMin, magMax;
      double compMin, compMax;
      
      DataArray();
    };
    
    /**
     * Add a cell data array that contains the cell IDs.
     * Useful when the original IDs have to be determined after the cells
     * changed (e.g. after clipping).
     */
    void AddOriginalCellIDs (vtkDataSet* dataset);
  };
} // namespace nutogui

#endif // __NUTOGUI_RESULTVIEWER_DATA_H__
