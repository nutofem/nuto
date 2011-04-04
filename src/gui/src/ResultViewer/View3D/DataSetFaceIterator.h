/**\file
 * Class to iterate over all faces of an VTK data set.
 */
/*
 * Written 2011 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#ifndef __NUTOGUI_DATASETFACEITERATOR_H__
#define __NUTOGUI_DATASETFACEITERATOR_H__

#include <vtkSmartPointer.h>
#include <vtkType.h>

class vtkCell;
class vtkDataSet;
class vtkGenericCell;

namespace nutogui
{
  /**
   * Iterate over all faces of a vtkDataSet.
   */
  class DataSetFaceIterator
  {
    vtkSmartPointer<vtkDataSet> dataset;
    
    /// ID of current cell in data set.
    vtkIdType currentCellId;
    /// Current cell in data set.
    vtkSmartPointer<vtkCell> currentCell;
    
    /// ID of current face in current cell.
    int currentFaceCellId;
    /// Next cell returned
    vtkSmartPointer<vtkCell> nextFaceCell;
    
    /// Advance to next point
    void Advance ();
    
    /**
     * Return whether cells of the given cell type have faces.
     * (Cells without faces are returned straight away.)
     */
    static bool HasCellTypeFaces (int cellType);
    
    /**
     * Cells to receive cell data from data set.
     * Some data sets always return the same vtkCell* pointer,
     * but we need different ones each time.
     * So cycle through a number of preallocated cells that are filled.
     */
    vtkSmartPointer<vtkGenericCell> cellStorage[2];
    /// Next member of cellStorage to use
    unsigned int nextStorage;
    
    /// Get a cell from data set. Uses cellStorage to save cell data.
    vtkCell* GetCell (vtkIdType index);
  public:
    /// Construct iterator.
    DataSetFaceIterator (vtkDataSet* dataset);
    
    /// Return whether there is a 'next' face available.
    bool HasNext() const;
    /// Return next face
    vtkCell* Next ();
    /// Return next face + originating cell
    vtkCell* Next (vtkIdType& origCell);
  };
} // namespace nutogui

#endif // __NUTOGUI_DATASETFACEITERATOR_H__
