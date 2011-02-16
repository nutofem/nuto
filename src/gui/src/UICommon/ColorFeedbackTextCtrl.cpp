/**\file
 * Text control to facilitate visual feedback by color change
 */
/*
 * Written 2011 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#include "common.h"

#include "uicommon/ColorFeedbackTextCtrl.h"

#include <float.h>

namespace uicommon
{
  IMPLEMENT_DYNAMIC_CLASS (ColorFeedbackTextCtrl, wxTextCtrl)
  
  // Maximum fraction of signal color to apply.
  static const float maxMix = 0.6f;
  
  ColorFeedbackTextCtrl::ColorFeedbackTextCtrl ()
   : signalBackgroundColor (255, 0, 0), signalForegroundColor (255, 255, 208), currentMix (0)
  {
  }

  ColorFeedbackTextCtrl::ColorFeedbackTextCtrl (wxWindow* parent, wxWindowID id,
						const wxString& value,
						const wxPoint& pos,
						const wxSize& size,
						long style,
						const wxValidator& validator,
						const wxString& name)
   : wxTextCtrl (parent, id, value, pos, size, style, validator, name),
     signalBackgroundColor (255, 0, 0), signalForegroundColor (255, 255, 208), currentMix (0)
  {
    originalBackgroundColor = GetBackgroundColour ();
    originalForegroundColor = GetForegroundColour ();
  }
  
  static inline wxColour ColourLerp (const wxColour& a, const wxColour& b, float factor)
  {
    const float f1 = 1-factor;
    const float f2 = factor;
    return wxColour (a.Red()*f1   + b.Red()*f2,
		     a.Green()*f1 + b.Green()*f2,
		     a.Blue()*f1  + b.Blue()*f2,
		     a.Alpha()*f1 + b.Alpha()*f2);
  }

  void ColorFeedbackTextCtrl::UpdateActualColors (bool repaint)
  {
    wxTextCtrl::SetBackgroundColour (ColourLerp (originalBackgroundColor, signalBackgroundColor, currentMix));
    wxTextCtrl::SetForegroundColour (ColourLerp (originalForegroundColor, signalForegroundColor, currentMix));
    if (repaint) Refresh (true);
  }

  void ColorFeedbackTextCtrl::SetMix (float mix)
  {
    if (fabsf (mix - currentMix) > FLT_EPSILON)
    {
      currentMix = mix;
      UpdateActualColors (true);
    }
  }

  void ColorFeedbackTextCtrl::SetFeedbackColorFactor (float factor)
  {
    SetMix (factor * maxMix);
  }

  bool ColorFeedbackTextCtrl::SetBackgroundColour (const wxColour& colour)
  {
    originalBackgroundColor = colour;
    return wxTextCtrl::SetBackgroundColour (ColourLerp (originalBackgroundColor, signalBackgroundColor, currentMix));
  }

  bool ColorFeedbackTextCtrl::SetForegroundColour (const wxColour& colour)
  {
    originalForegroundColor = colour;
    return wxTextCtrl::SetForegroundColour (ColourLerp (originalForegroundColor, signalForegroundColor, currentMix));
  }

  void ColorFeedbackTextCtrl::SetSignalBackgroundColour (const wxColour& colour)
  {
    signalBackgroundColor = colour;
    wxTextCtrl::SetBackgroundColour (ColourLerp (originalBackgroundColor, signalBackgroundColor, currentMix));
  }

  void ColorFeedbackTextCtrl::SetSignalForegroundColour (const wxColour& colour)
  {
    signalForegroundColor = colour;
    wxTextCtrl::SetForegroundColour (ColourLerp (originalForegroundColor, signalForegroundColor, currentMix));
  }

} // namespace uicommon
