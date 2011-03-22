/**\file
 * Slider control with an (integer) text control "buddy"
 */
/*
 * Written 2010 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#include "common.h"
#include "uicommon/TextCtrlBuddySlider.h"

#include <wx/wx.h>

namespace uicommon
{
    
  BEGIN_EVENT_TABLE(TextCtrlBuddySlider, BuddySlider)
    EVT_TEXT(wxID_ANY, TextCtrlBuddySlider::OnTextChanged)
    EVT_KEY_DOWN(TextCtrlBuddySlider::OnKeyDown)
  END_EVENT_TABLE()
  
  // Right-aligning edit seems to be broken on MSW
  #ifndef __WXMSW__
  #define EDIT_STYLE  wxTE_RIGHT
  #else
  #define EDIT_STYLE  0
  #endif

  TextCtrlBuddySlider::TextCtrlBuddySlider (wxWindow* parent, wxWindowID id, 
					    int value, int minValue, int maxValue,
					    const wxPoint& pos, const wxSize& size,
					    long sliderStyle)
   : BuddySlider (parent, id, pos, size, sliderStyle), buddy (0), rangeMin (minValue), rangeMax (maxValue)
  {
    buddy = new wxTextCtrl (this, wxID_ANY, wxEmptyString,
			    wxDefaultPosition, wxDefaultSize,
			    EDIT_STYLE);
    buddy->SetValidator (wxTextValidator (wxFILTER_NUMERIC));
    oldBuddyEvtHandler = buddy->GetEventHandler();
    buddy->SetEventHandler (this);
    
    /* Make text control font a bit smaller (text controls are usually quite
       a bit higher than sliders) */
    {
      wxFont buddyFont (buddy->GetFont ());
      buddyFont.SetPointSize (buddyFont.GetPointSize() * 0.9);
      buddy->SetFont (buddyFont);
    }
    
    UpdateTextCtrlMinSize ();
    SetBuddy (buddy);
    SetValue (value);
  }
  
  TextCtrlBuddySlider::~TextCtrlBuddySlider ()
  {
    buddy->SetEventHandler (oldBuddyEvtHandler);
  }
  
  void TextCtrlBuddySlider::SliderValueChanged (int newValue)
  {
    currentValue = newValue;
    buddy->ChangeValue (wxString::Format (wxT("%d"), newValue));
  }

  wxWindow* TextCtrlBuddySlider::GetBuddy() const
  {
    return buddy;
  }

  void TextCtrlBuddySlider::UpdateTextCtrlMinSize ()
  {
    wxSize buddyMinSize (buddy->GetMinSize ());
    {
      wxScreenDC tempDC;
      tempDC.SetFont (buddy->GetFont ());
      // Measure space needed for range min value + some slack ...
      wxSize numbersExtent1 (
	tempDC.GetTextExtent (wxString::Format (wxT("%d00"), rangeMin)));
      // ... as well as range max value
      wxSize numbersExtent2 (
	tempDC.GetTextExtent (wxString::Format (wxT("%d00"), rangeMax)));
      buddyMinSize.SetWidth (std:: max (numbersExtent1.GetWidth (),
					numbersExtent2.GetWidth ()));
    }
    buddy->SetMinSize (buddyMinSize);
  }
  
  int TextCtrlBuddySlider::GetValue() const
  {
    return currentValue;
  }
  
  void TextCtrlBuddySlider::SetValue (int value)
  {
    if (currentValue != value)
    {
      UpdateSliderValue (value);
      // slider value may be != 'value' due range limiting
      currentValue = GetSlider()->GetValue();
      buddy->ChangeValue (wxString::Format (wxT("%d"), currentValue));
    }
  }
		     
  void TextCtrlBuddySlider::SetRange (int minValue, int maxValue)
  {
    rangeMin = minValue;
    rangeMax = maxValue;
    UpdateSliderRange();
    UpdateTextCtrlMinSize ();
    Layout ();
  }
  
  void TextCtrlBuddySlider::OnTextChanged (wxCommandEvent& event)
  {
    // Seems this event is triggered during text control creation, when buddy is not yet valid
    if (!buddy) return;

    long value;
    if (buddy->GetValue ().ToLong (&value))
    {
      UpdateSliderValue (value);
      // slider value may be != 'value' due range limiting
      currentValue = GetSlider()->GetValue();
      if (currentValue != value)
	buddy->ChangeValue (wxString::Format (wxT("%d"), currentValue));
      
      wxScrollEvent ev (wxEVT_SCROLL_CHANGED, GetId(), currentValue, wxHORIZONTAL);
      ev.SetEventObject (this);
      GetEventHandler()->ProcessEvent (ev);   
    }
  }

  void TextCtrlBuddySlider::OnKeyDown (wxKeyEvent& event)
  {
    int dir = 0;
    if (event.GetKeyCode() == WXK_UP)
      dir = -1;
    else if (event.GetKeyCode() == WXK_DOWN)
      dir = 1;
    
    if (dir == 0)
    {
      // Skip all other keys (so the text control still works ... ;P)
      event.Skip ();
      return;
    }
    // Emulate what an actual slider does (change by line size on up/down)
    int newValue = GetSlider()->GetValue() + GetSlider()->GetLineSize()*dir;
    SetValue (newValue);
    
    wxScrollEvent ev (wxEVT_SCROLL_CHANGED, GetId(), currentValue, wxHORIZONTAL);
    ev.SetEventObject (this);
    GetEventHandler()->ProcessEvent (ev);   
  }

} // namespace uicommon
