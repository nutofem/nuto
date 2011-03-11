/**\file
 * Interface to access result (i.e. visualization) data from NuTo as a VTK data set
 */
/*
 * Written 2011 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#ifndef __NUTOGUI_RESULTDATASOURCEVTK_H__
#define __NUTOGUI_RESULTDATASOURCEVTK_H__

class vtkDataSet;

#include <boost/shared_ptr.hpp>

namespace nutogui
{
  /// Interface to access result (i.e. visualization) data from NuTo as a VTK data set
  struct ResultDataSourceVTK
  {
    /// Get data as a VTK DataSet
    virtual vtkDataSet* QueryDataSet () = 0;
  };
  typedef boost::shared_ptr<ResultDataSourceVTK> ResultDataSourceVTKPtr;
} // namespace nutogui

#endif // __NUTOGUI_RESULTDATASOURCEVTK_H__

