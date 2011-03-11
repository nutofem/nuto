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
class wxString;

#include <boost/shared_ptr.hpp>

namespace nutogui
{
  /// Interface to access result (i.e. visualization) data from NuTo as a VTK data set
  struct ResultDataSourceVTK
  {
    /// Get number of data sets
    virtual size_t GetNumDataSets () const = 0;
    /// Get data as a VTK DataSet
    virtual vtkDataSet* GetDataSet (size_t index) = 0;
    /// Query name of a data set
    virtual const wxString& GetDataSetName (size_t index) const = 0;
  };
  typedef boost::shared_ptr<ResultDataSourceVTK> ResultDataSourceVTKPtr;
} // namespace nutogui

#endif // __NUTOGUI_RESULTDATASOURCEVTK_H__

