/**\file
 * Slider control with an (integer) text control "buddy"
 */
/*
 * Written 2010 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#ifndef __UICOMMON_TEXTCTRLBUDDYSLIDER_H__
#define __UICOMMON_TEXTCTRLBUDDYSLIDER_H__

#include "export.h"
#include "BuddySlider.h"

class wxTextCtrl;

namespace uicommon
{
    
  class UICOMMON_EXPORTED TextCtrlBuddySlider : public BuddySlider
  {
  protected:
    wxTextCtrl* buddy;
    int rangeMin, rangeMax;
    int currentValue;
    
    virtual int GetRangeMin() { return rangeMin; }
    virtual int GetRangeMax() { return rangeMax; }
    virtual void SliderValueChanged (int newValue);
    virtual wxWindow* GetBuddy() const;
    
    /**
     * Adjusts the minimum size of the text control to hold all digits for
     * the minimum and maximum values plus one digit
     */
    void UpdateTextCtrlMinSize ();
    
    void OnTextChanged (wxCommandEvent& event);
    void OnKeyDown (wxKeyEvent& event);
  public:
    TextCtrlBuddySlider (wxWindow* parent, wxWindowID id, 
		         int value, int minValue, int maxValue,
		         const wxPoint& pos = wxDefaultPosition,
		         const wxSize& size = wxDefaultSize,
			 long sliderStyle = wxSL_HORIZONTAL);
    
    int GetValue() const;
    void SetValue (int value);
    
    int GetMin() const { return rangeMin; }
    int GetMax() const { return rangeMax; }
    void SetRange (int minValue, int maxValue);

    DECLARE_EVENT_TABLE()
  };
  
} // namespace uicommon

#endif // __UICOMMON_TEXTCTRLBUDDYSLIDER_H__
