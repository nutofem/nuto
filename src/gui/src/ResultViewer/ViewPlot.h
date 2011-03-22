/**\file
 * Result viewer panel for plots.
 */
/*
 * Written 2011 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#ifndef __NUTOGUI_RESULTVIEWER_VIEWPLOT_H__
#define __NUTOGUI_RESULTVIEWER_VIEWPLOT_H__

#include "ResultViewerImpl.h"
#include "ViewPanelContentVTK.h"

#include <vtkSmartPointer.h>

class vtkRenderer;

namespace nutogui
{
  class ResultViewerImpl::ViewPlot : public ViewPanelContentVTK
  {
  protected:
    vtkSmartPointer<vtkRenderer> renderer;
    void SetupVTKRenderer ();
    void SetupRenderer ();
  public:
    ViewPlot (ViewPanel* parent);

    void SetData (const DataConstPtr& data);
    wxWindow* CreateTopTools (wxWindow* parentWindow);
  };
} // namespace nutogui

#endif // __NUTOGUI_RESULTVIEWER_VIEWPLOT_H__
