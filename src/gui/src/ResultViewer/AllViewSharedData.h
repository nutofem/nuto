/**\file
 * Data shared by pretty much all views.
 */
/*
 * Written 2011 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#ifndef __NUTOGUI_RESULTVIEWER_ALLVIEWSHAREDDATA_H__
#define __NUTOGUI_RESULTVIEWER_ALLVIEWSHAREDDATA_H__

#include "ResultViewerImpl.h"

#include "ViewPanel.h"

#include <boost/unordered_set.hpp>

#include <vtkType.h>

namespace nutogui
{
  struct ResultViewerImpl::AllViewSharedData : public ViewPanel::SharedDataBase
  {
    typedef boost::unordered_set<vtkIdType> SelectedCellsSet;
    /// Set of selected cells
    SelectedCellsSet selectedCellIDs;
    
    static AllViewSharedDataPtr Setup (ViewPanel* panel);
  };
} // namespace nutogui

#endif // __NUTOGUI_RESULTVIEWER_ALLVIEWSHAREDDATA_H__


