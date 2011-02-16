/**\file
 * Text control to facilitate visual feedback by color change
 */
/*
 * Written 2011 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#ifndef __UICOMMON_COLORFEEDBACKTEXTCTRL_H__
#define __UICOMMON_COLORFEEDBACKTEXTCTRL_H__

#include "export.h"

#include <wx/wx.h>

namespace uicommon
{
  /**
   * Text control subclass to facilitate visual feedback on input changes by changing the control colors.
   * This control has, in addition to the default foreground/background colors,
   * so-called "signal" foreground/background colors. How much of these colors is shown
   * depends on the "mix" factor. A "mix" factor of 0 means the original colors are shown;
   * a factor of 1 means the signal colors are at their maximum contribution.
   * (They never replace the original colors entirely to avoid crass visuals.)
   *
   * The default signal colors are red for the background and a light yellow for the foreground.
   */
  class UICOMMON_EXPORTED ColorFeedbackTextCtrl : public wxTextCtrl
  {
    DECLARE_DYNAMIC_CLASS(ColorFeedbackTextCtrl)
    
    wxColour originalBackgroundColor;
    wxColour originalForegroundColor;
    
    wxColour signalBackgroundColor;
    wxColour signalForegroundColor;
    float currentMix;
    
    void UpdateActualColors (bool repaint);
    void SetMix (float mix);
  public:
    ColorFeedbackTextCtrl ();
    ColorFeedbackTextCtrl (wxWindow* parent, wxWindowID id,
			   const wxString& value = wxEmptyString,
			   const wxPoint& pos = wxDefaultPosition,
			   const wxSize& size = wxDefaultSize,
			   long style = 0,
			   const wxValidator& validator = wxDefaultValidator,
			   const wxString& name = wxTextCtrlNameStr);

    /// Reset mix factor to show only original colors.
    void ResetFeedbackColorFactor () { SetMix (0); }
    /// Set the amount of signal colors to mix into foreground and background.
    void SetFeedbackColorFactor (float factor);
    
    /*@{*/
    /// Get/set original colors.
    bool SetBackgroundColour (const wxColour& colour);
    bool SetForegroundColour (const wxColour& colour);
    wxColour GetBackgroundColour() { return originalBackgroundColor; }
    wxColour GetForegroundColour() { return originalForegroundColor; }
    /*@}*/
    
    /// Set signal background color.
    void SetSignalBackgroundColour (const wxColour& colour);
    /// Set signal foreground color.
    void SetSignalForegroundColour (const wxColour& colour);
    /// Get signal background color.
    wxColour GetSignalBackgroundColour() const { return signalBackgroundColor; }
    /// Get signal foreground color.
    wxColour GetSignalForegroundColour() const { return signalForegroundColor; }
  };
} // namespace uicommon

#endif // __UICOMMON_COLORFEEDBACKTEXTCTRL_H__
