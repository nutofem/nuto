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

class vtkChartXY;
class vtkColorSeries;
class vtkContextActor;
class vtkRenderer;
class vtkTextActor;

namespace nutogui
{
  class ResultViewerImpl::ViewPlot : public ViewPanelContentVTK
  {
  protected:
    DataConstPtr data;
    AllViewSharedDataPtr sharedAllData;
    
    vtkSmartPointer<vtkColorSeries> plotColors;
    vtkSmartPointer<vtkRenderer> renderer;
    vtkSmartPointer<vtkContextActor> contextActor;
    vtkSmartPointer<vtkChartXY> chart;
    vtkSmartPointer<vtkTextActor> noDataMessage;
    void SetupVTKRenderer ();
    void SetupRenderer ();
    
    void SetupChart ();
    
    void OnSelectedCellsChanged (wxCommandEvent& event);
  public:
    ViewPlot (ViewPanel* parent);

    void SetData (const DataConstPtr& data);
    wxWindow* CreateTopTools (wxWindow* parentWindow);
			      
    DECLARE_EVENT_TABLE()
  };
} // namespace nutogui

#endif // __NUTOGUI_RESULTVIEWER_VIEWPLOT_H__
