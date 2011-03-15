/**\file
 * Result viewer view panel: class for child with content
 */
/*
 * Written 2011 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#include "common.h"
#include "ViewPanelContent.h"

namespace nutogui
{
  ResultViewerImpl::ViewPanel::Content::Content (ViewPanel* parent) : wxPanel (parent)
  {
  }
  
  void ResultViewerImpl::ViewPanel::Content::PostToOthers (wxEvent& event)
  {
    static_cast<ViewPanel*> (GetParent())->PostToOthers (event);
  }
} // namespace nutogui
