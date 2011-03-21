/**\file
 * Result viewer view panel: class for child with content displayed using VTK
 */
/*
 * Written 2011 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#ifndef __RESULTVIEWERIMPL_VIEWPANELCONTENTVTK_H__
#define __RESULTVIEWERIMPL_VIEWPANELCONTENTVTK_H__

#include "ViewPanelContent.h"

namespace vtkwx
{
  class RenderWidget;
} // namespace vtkwx

namespace nutogui
{
  class ResultViewerImpl::ViewPanelContentVTK : public ViewPanel::Content
  {
  public:
    ViewPanelContentVTK (ViewPanel* parent);
    
    DECLARE_EVENT_TABLE()
  protected:
    /// Creates the actual render widget
    vtkwx::RenderWidget* CreateRenderWidget (wxWindow* parent);
    /// Set up actual VTK rendering. Called after VTK display control was created
    virtual void SetupVTKRenderer () = 0;
    
    /// Get render widget
    vtkwx::RenderWidget* GetRenderWidget ();
    /// Convenience method: re-render current view
    void Render ();
  private:
    class RenderWidget;
    /// Control with VTK display
    RenderWidget* renderWidget;
    
    void OnRenderWidgetRealized (wxCommandEvent& event);
  };
} // namespace nutogui

#endif // __RESULTVIEWERIMPL_VIEWPANELCONTENTVTK_H__
