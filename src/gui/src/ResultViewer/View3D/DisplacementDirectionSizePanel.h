/**\file
 * Panel for controlling the scaling of displacement direction displaying.
 */
/*
 * Written 2011 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#ifndef __NUTOGUI_DISPLACEMENTDIRECTIONSIZEPANEL_H__
#define __NUTOGUI_DISPLACEMENTDIRECTIONSIZEPANEL_H__

#include <wx/wx.h>

#include "DirScaleDragThing.h"

#include <locale>

namespace uicommon
{
  class ColorFeedbackTextCtrl;
} // namespace uicommon

namespace nutogui
{
  BEGIN_DECLARE_EVENT_TYPES()
    DECLARE_EVENT_TYPE(EVENT_DIRECTION_SCALE_CHANGED, 0)
  END_DECLARE_EVENT_TYPES()

  #define EVT_DIRECTION_SCALE_CHANGED(fn) \
      DECLARE_EVENT_TABLE_ENTRY( EVENT_DIRECTION_SCALE_CHANGED, wxID_ANY, -1, \
      (wxObjectEventFunction) \
      wxStaticCastEvent( wxEventFunction, & fn ), (wxObject *) NULL ),
  
  class DisplacementDirectionSizePanel : public wxPanel
  {
    class ScaleChangedEvent;
    
    uicommon::ColorFeedbackTextCtrl* scaleInput;
    
    std::locale locale; // Used for numeric input
    float scaleValue;
    
    void OnScaleInputEnter (wxCommandEvent& event);
    void OnScaleInputChanged (wxCommandEvent& event);
    
    float scaleDragStartScale;
    void OnScaleDragBegin (wxCommandEvent& event);
    void OnScaleDragChanged (DirScaleDragThing::ChangedEvent& event);
    void SetRelativeScale (float scale);
    
    long lastWheelTimeStamp;
    double wheelDist;
    void OnMouseWheel (wxMouseEvent& event);
  public:
    DisplacementDirectionSizePanel (wxWindow* parent, wxPoint position);
    
    float GetDisplacementScale() const { return scaleValue; }
    void SetDisplacementScale (float value);
    
    DECLARE_EVENT_TABLE()
  };
} // namespace nutogui

#endif // __NUTOGUI_DISPLACEMENTDIRECTIONSIZEPANEL_H__
