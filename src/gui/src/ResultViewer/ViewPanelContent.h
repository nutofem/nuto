/**\file
 * Result viewer view panel: class for child with content
 */
/*
 * Written 2011 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#ifndef __RESULTVIEWERIMPL_VIEWPANELCONTENT_H__
#define __RESULTVIEWERIMPL_VIEWPANELCONTENT_H__

#include "ViewPanel.h"

namespace nutogui
{
  class ResultViewerImpl::ViewPanel::Content : public wxPanel
  {
  public:
    Content (ViewPanel* parent);
    
    /// Create window to be displayed in top tool bar (between panel splitting controls)
    virtual wxWindow* CreateTopTools (wxWindow* parentWindow) = 0;
    /**
     * Destroy any data associated with the window created by CreateTopTools.
     * The window itself must not be destroyed/deleted!
     */
    virtual void DestroyTopTools (wxWindow* tools) {}
  protected:
    /// Post an event to all other view panel content children
    void PostToOthers (wxEvent& event);
    
    /**
     * Helper function: update minimal size of a control in a toolbar
     */
    void UpdateToolbarControlMinSize (wxWindow* control,
				      wxAuiToolBar* toolbar,
				      int forceHeight = 0);
  };
} // namespace nutogui

#endif // __RESULTVIEWERIMPL_VIEWPANELCONTENT_H__
