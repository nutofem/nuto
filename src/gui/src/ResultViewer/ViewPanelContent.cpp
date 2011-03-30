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

#include <wx/aui/auibar.h>

namespace nutogui
{
  ResultViewerImpl::ViewPanel::Content::Content (ViewPanel* parent) : wxPanel (parent)
  {
  }
  
  void ResultViewerImpl::ViewPanel::Content::PostToOthers (wxEvent& event)
  {
    static_cast<ViewPanel*> (GetParent())->PostToOthers (event);
  }
  
  void ResultViewerImpl::ViewPanel::Content::UpdateToolbarControlMinSize (wxWindow* control,
									  wxAuiToolBar* toolbar,
									  int forceHeight)
  {
    wxAuiToolBarItem* ctrl_item = nullptr;
    for (size_t i = 0; i < toolbar->GetToolCount(); i++)
    {
      wxAuiToolBarItem* item = toolbar->FindToolByIndex (i);
      if (item->GetWindow() == control)
      {
	ctrl_item = item;
	break;
      }
    }
    if (ctrl_item)
    {
      wxSize newMinSize (control->GetBestSize());
      if (forceHeight > 0) newMinSize.SetHeight (forceHeight);
      ctrl_item->SetMinSize (newMinSize);
      toolbar->Realize();
    }
  }

} // namespace nutogui
