/**\file
 * Slider control with a "buddy" (another control to display+enter data)
 */
/*
 * Written 2010 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#ifndef __UICOMMON_BUDDYSLIDER_H__
#define __UICOMMON_BUDDYSLIDER_H__

#include "export.h"
#include <wx/slider.h>
#include <wx/window.h>

class wxSlider;

namespace uicommon
{
    
  class UICOMMON_EXPORTED BuddySlider : public wxWindow
  {
  private:
    wxSlider* sliderCtrl;
    
  protected:
    /// Get slider minimum value from derived implementations
    virtual int GetRangeMin() = 0;
    /// Get slider maximum value from derived implementations
    virtual int GetRangeMax() = 0;
    /// Tell derived implementations the slider value has changed
    virtual void SliderValueChanged (int newValue) = 0;
    /// Query buddy control from derived implementation
    virtual wxWindow* GetBuddy() const = 0;
    
    /// Set the buddy control (called by derived implementations)
    void SetBuddy (wxWindow* buddy);
    /// Update the range of the slider control (called by derived implementations)
    void UpdateSliderRange ();
    /**
     * Update the value of the slider control, usually after the buddy control
     * was changed (called by derived implementations)
     */
    void UpdateSliderValue (int value);
    
    wxSlider* GetSlider() const { return sliderCtrl; }
    
    void OnResize (wxSizeEvent& event);
    void OnSliderScroll (wxScrollEvent& event);
    void OnSysColourChanged(wxSysColourChangedEvent& event);
  public:
    BuddySlider (wxWindow* parent, wxWindowID id, const wxPoint& pos = wxDefaultPosition,
		 const wxSize& size = wxDefaultSize, long sliderStyle = wxSL_HORIZONTAL);

    // Shadow wxWindow method, forward tip to children
    void SetToolTip (const wxString& tip);

    DECLARE_EVENT_TABLE()
  };
  
} // namespace uicommon

#endif // __UI_BUDDYSLIDER_H__
