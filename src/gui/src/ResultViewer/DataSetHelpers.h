/**\file
 * Helpers to deal with data sets.
 */
/*
 * Written 2011 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#ifndef __NUTOGUI_DATASETHELPERS_H__
#define __NUTOGUI_DATASETHELPERS_H__

class vtkDataSet;
class vtkPointSet;

namespace nutogui
{
  struct DataSetHelpers
  {
    /**
     * Copy point data from some data set onto a point set.
     */
    static void CopyPoints (vtkPointSet* dest, vtkDataSet* src);
  };
} // namespace nutogui

#endif // __NUTOGUI_DATASETHELPERS_H__
