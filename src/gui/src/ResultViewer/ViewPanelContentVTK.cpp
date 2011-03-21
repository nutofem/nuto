/**\file
 * Result viewer view panel: class for child with content displayed using VTK
 */
/*
 * Written 2011 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#include "common.h"

#include "ViewPanelContentVTK.h"

#include "vtkwx.h"

#include <vtkRenderWindow.h>

namespace nutogui
{
  BEGIN_DECLARE_EVENT_TYPES()
    DECLARE_EVENT_TYPE(EVENT_RENDER_WIDGET_REALIZED, 0)
  END_DECLARE_EVENT_TYPES()
  DEFINE_EVENT_TYPE(EVENT_RENDER_WIDGET_REALIZED)
  
  class ResultViewerImpl::ViewPanelContentVTK::RenderWidget : public vtkwx::RenderWidget
  {
    bool painted;
    
    void OnPaint (wxPaintEvent& event);
  public:
    RenderWidget (wxWindow *parent, wxWindowID id = -1,
      const wxPoint& pos = wxDefaultPosition,
      const wxSize& size = wxDefaultSize,
      long style = 0, const wxString& name = wxT("RenderWidget"))
     : vtkwx::RenderWidget (parent, id, pos, size, style, name), painted (false) {}
     
    DECLARE_EVENT_TABLE()
  };
  
  BEGIN_EVENT_TABLE(ResultViewerImpl::ViewPanelContentVTK::RenderWidget, vtkwx::RenderWidget)
    EVT_PAINT(RenderWidget::OnPaint)
  END_EVENT_TABLE()
  
  void ResultViewerImpl::ViewPanelContentVTK::RenderWidget::OnPaint (wxPaintEvent& event)
  {
    if (!painted)
    {
      painted = true;
      // First paint: post event to parent so it can set up the VTK renderer etc.
      wxCommandEvent event (EVENT_RENDER_WIDGET_REALIZED);
      event.SetEventObject (this);
      wxPostEvent (this, event);
    }
    
    event.Skip();
  }

  //-------------------------------------------------------------------
  
  BEGIN_EVENT_TABLE(ResultViewerImpl::ViewPanelContentVTK, ViewPanel::Content)
    EVT_COMMAND(wxID_ANY, EVENT_RENDER_WIDGET_REALIZED, ResultViewerImpl::ViewPanelContentVTK::OnRenderWidgetRealized)
  END_EVENT_TABLE()

  ResultViewerImpl::ViewPanelContentVTK::ViewPanelContentVTK (ViewPanel* parent)
    : ViewPanel::Content (parent), renderWidget (nullptr) {}
    

  vtkwx::RenderWidget* ResultViewerImpl::ViewPanelContentVTK::CreateRenderWidget (wxWindow* parent)
  {
    renderWidget = new RenderWidget (parent);
    return renderWidget;
  }
  
  vtkwx::RenderWidget* ResultViewerImpl::ViewPanelContentVTK::GetRenderWidget ()
  {
    return renderWidget;
  }
  
  void ResultViewerImpl::ViewPanelContentVTK::Render ()
  {
    renderWidget->GetRenderWindow()->Render ();
  }
  
  void ResultViewerImpl::ViewPanelContentVTK::OnRenderWidgetRealized (wxCommandEvent& event)
  {
    SetupVTKRenderer();
  }
} // namespace nutogui
