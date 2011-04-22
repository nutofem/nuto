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

#include "uicommon/ControlWithItemsClientDataWrapper.h"
#include <vector>

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
    AllViewSharedDataPtr sharedAllData;
    
    vtkSmartPointer<vtkColorSeries> plotColors;
    vtkSmartPointer<vtkRenderer> renderer;
    vtkSmartPointer<vtkContextActor> contextActor;
    vtkSmartPointer<vtkChartXY> chart;
    vtkSmartPointer<vtkTextActor> noDataMessage;
    void SetupVTKRenderer ();
    void SetupRenderer ();
    
    wxAuiToolBar* toolbar;
    uicommon::ControlWithItemsClientDataWrapper<size_t> displayDataChoice;
    wxChoice* visChoiceCtrl;
    void SetData (const DataConstPtr& data);
    
    bool SetupChart (size_t dataArray, int viscomp);
    
    void OnSelectedCellsChanged (wxCommandEvent& event);
    
    size_t currentDataArray;
    int currentVisComp;
    struct VisOpt
    {
      int comp;
      
      VisOpt() : comp (-1) {}
    };
    std::vector<VisOpt> lastVisOpt;
    
    void OnDisplayDataChanged (wxCommandEvent& event);
    void OnVisOptionChanged (wxCommandEvent& event);
    void UpdateVisOptionChoice (size_t dataIndex, int initialSel);
  public:
    ViewPlot (ViewPanel* parent);

    wxWindow* CreateTopTools (wxWindow* parentWindow);
    void DestroyTopTools (wxWindow* tools);
			      
    DECLARE_EVENT_TABLE()
  };
} // namespace nutogui

#endif // __NUTOGUI_RESULTVIEWER_VIEWPLOT_H__
