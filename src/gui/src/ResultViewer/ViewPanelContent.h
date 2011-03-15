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
    
    virtual void SetData (const DataConstPtr& data) = 0;
    virtual Content* Clone (ViewPanel* parent) = 0;
    virtual wxWindow* CreateTopTools (wxWindow* parentWindow) = 0;
  protected:
    /// Post an event to all other view panel content children
    void PostToOthers (wxEvent& event);
  };
} // namespace nutogui

#endif // __RESULTVIEWERIMPL_VIEWPANELCONTENT_H__
