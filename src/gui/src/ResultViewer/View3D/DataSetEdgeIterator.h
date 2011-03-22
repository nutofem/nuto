/**\file
 * Class to iterate over all edges of an VTK data set.
 */
/*
 * Written 2011 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#ifndef __NUTOGUI_DATASETEDGEITERATOR_H__
#define __NUTOGUI_DATASETEDGEITERATOR_H__

#include <vtkSmartPointer.h>
#include <vtkType.h>

#include <utility>

class vtkCell;
class vtkDataSet;

namespace nutogui
{
  /**
   * Iterate over all edges of a vtkDataSet.
   */
  class DataSetEdgeIterator
  {
    vtkSmartPointer<vtkDataSet> dataset;
    
    /// ID of current cell in data set.
    vtkIdType currentCellId;
    /// Current cell in data set.
    vtkSmartPointer<vtkCell> currentCell;
    
    /// ID of current edge in current cell.
    int currentEdgeCellId;
    /// Current edge in current cell.
    vtkSmartPointer<vtkCell> currentEdgeCell;
    
    /// Current point in edge. (Note: first point is skipped)
    vtkIdType currentEdgeCellPointId;
    
    /// The last point fetched (used as first point in the pair of points that is an edge)
    vtkIdType lastPoint;
    
    /// Advance to next point
    void Advance ();
  public:
    /// Construct iterator.
    DataSetEdgeIterator (vtkDataSet* dataset);
    
    /// Edge type.
    typedef std::pair<vtkIdType, vtkIdType> Edge;
    /// Return whether there is a 'next' edge available.
    bool HasNext() const;
    /// Return next edge
    Edge Next ();
    /// Return next edge + originating cell
    Edge Next (vtkIdType& origCell);
  };
} // namespace nutogui

#endif // __NUTOGUI_DATASETEDGEITERATOR_H__
