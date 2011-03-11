/**\file
 * Slider control with a "buddy" (another control to display+enter data)
 */
/*
 * Written 2010 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#include "common.h"
#include "uicommon/BuddySlider.h"

#include <wx/wx.h>
#include <wx/tooltip.h>

namespace uicommon
{
    
  BEGIN_EVENT_TABLE(BuddySlider, wxWindow)
    EVT_SIZE(BuddySlider::OnResize)
    EVT_COMMAND_SCROLL(wxID_ANY, BuddySlider::OnSliderScroll)
    EVT_SYS_COLOUR_CHANGED(BuddySlider::OnSysColourChanged)
  END_EVENT_TABLE()
  
  BuddySlider::BuddySlider (wxWindow* parent, wxWindowID id,
			    const wxPoint& pos, const wxSize& size)
   : wxWindow (parent, id, pos, size)
  {
    sliderCtrl = new wxSlider (this, wxID_ANY, 0, 0, 1);
    sliderCtrl->SetEventHandler (this);
    
    // Create sizer ...
    wxBoxSizer* sizer = new wxBoxSizer (wxHORIZONTAL);
    sizer->Add (sliderCtrl, 1, wxEXPAND);
    // Dummy spacer, later replaced by buddy
    sizer->AddSpacer (0);
    
    SetSizer (sizer);
  }

  void BuddySlider::SetToolTip (const wxString& tip)
  {
    wxWindow::SetToolTip (tip);
    sliderCtrl->SetToolTip (tip);
    wxWindow* buddy = GetBuddy();
    if (buddy) buddy->SetToolTip (tip);
  }
  
  void BuddySlider::SetBuddy (wxWindow* buddy)
  {
    // Put buddy into sizer ...
    GetSizer()->Detach (1);
    GetSizer()->Add (buddy);
    if (GetToolTip())
      buddy->SetToolTip (GetToolTip()->GetTip());
    
    UpdateSliderRange ();
    SliderValueChanged (sliderCtrl->GetValue ());
  }
  
  void BuddySlider::UpdateSliderRange ()
  {
    sliderCtrl->SetRange (GetRangeMin(), GetRangeMax());
    SliderValueChanged (sliderCtrl->GetValue ());
  }
  
  void BuddySlider::UpdateSliderValue (int value)
  {
    sliderCtrl->SetValue (value);
  }
  
  void BuddySlider::OnResize (wxSizeEvent& event)
  {
    // FIXME: shouldn't that happen, like, automatically?
    if ((event.GetSize().x == 0) || (event.GetSize().y == 0)) return;
    // Don't want resize events from children ...
    if (event.GetEventObject() != this)
    {
    #if wxABI_VERSION < 20900
      event.Skip ();
    #else
      // WX 2.9: triggers infinite loop
    #endif
      return;
    }
    GetSizer()->SetDimension (0, 0, event.GetSize().x, event.GetSize().y);
    GetSizer()->Layout();
  }
  
  void BuddySlider::OnSliderScroll (wxScrollEvent& event)
  {
    // Event will be propagated in any case
    event.Skip ();
    
    // Don't handle our own events
    if (event.GetEventObject() == this) return;
    
    // Update buddy control ...
    SliderValueChanged (event.GetPosition());
    // ...and pretend it's from 'this' window
    event.SetEventObject (this);
    event.SetId (GetId ());
  }
  
  void BuddySlider::OnSysColourChanged(wxSysColourChangedEvent& event)
  {
    /* wxGTK propagates SysColourChanged events to it's children.
       In our special case though this gets back to this window -
       so the default wxWindow SysColourChanged event handling creates
       an infinite loop which we break manually.

       The 0 check is actually a workaround for an apparent bug where
       aforementioned propagation code doesn't set the event object.
     */
    wxWindow* evWindow = static_cast<wxWindow*> (event.GetEventObject());
    if ((evWindow == 0) || (evWindow->GetParent() == this)) return;
    event.Skip();
  }
  
} // namespace uicommon
