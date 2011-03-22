/**\file
 * Control to change the scale of displacement direction display by dragging the mouse.
 */
/*
 * Written 2011 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#ifndef __NUTOGUI_DIRSCALEDRAGTHING_H__
#define __NUTOGUI_DIRSCALEDRAGTHING_H__

#include <wx/wx.h>

namespace nutogui
{
  BEGIN_DECLARE_EVENT_TYPES()
    DECLARE_EVENT_TYPE(EVENT_SCALEDRAG_BEGIN, 0)
    DECLARE_EVENT_TYPE(EVENT_SCALEDRAG_CHANGED, 0)
    DECLARE_EVENT_TYPE(EVENT_SCALEDRAG_END, 0)
  END_DECLARE_EVENT_TYPES()
  
  /**
   * Control to change the scale of displacement direction display by dragging the mouse.
   *
   * How it works (API-wise):
   * - When the user starts dragging, a EVENT_SCALEDRAG_BEGIN event is emitted.
   * - If the user drags the mouse, a EVENT_SCALEDRAG_CHANGED event is emitted.
   *   It is of type DirScaleDragThing::ChangedEvent and carries a factor that
   *   gives the scale computed from drag distance. That scale is normalized
   *   and must be multiplied by the scale _at the time the drag started_
   *   by the event handler!
   * - When dragging ends (mouse button up...), a EVENT_SCALEDRAG_END is emitted.
   */
  class DirScaleDragThing : public wxControl
  {
    wxBitmap image;
    wxBitmap hotImage;
    bool isHot;
    
    void OnPaint (wxPaintEvent& event);
    
    void OnMouseLeftDown (wxMouseEvent& event);
    void OnMouseLeftUp (wxMouseEvent& event);
    void OnMouseMove (wxMouseEvent& event);
    void OnMouseCaptureLost (wxMouseCaptureLostEvent& event);
    void OnMouseEnter (wxMouseEvent& event);
    void OnMouseLeave (wxMouseEvent& event);
    
    wxPoint dragStart;
    void DragEnd ();
    
    wxTimer hotFadeTimer;
    wxStopWatch hotFadeStopWatch;
    void OnHotFadeTimer (wxTimerEvent& event);
    
    void PaintImage (wxDC& dc, float fade);
    void SetHot (bool hotState);
    void DrawBitmapsMixed (wxDC& dc, wxCoord x, wxCoord y,
			   const wxBitmap& image1, const wxBitmap& image2,
			   float factor);
  public:
    DirScaleDragThing (wxWindow* parent, wxWindowID id,
		       wxPoint position = wxDefaultPosition,
		       wxSize size = wxDefaultSize);

    DECLARE_EVENT_TABLE()

  public:
    /**
     * Event sent when the scale factor was changed by dragging.
     */
    class ChangedEvent : public wxCommandEvent
    {
      float totalScale;
    public:
      ChangedEvent (wxWindowID id, float totalScale);
      
      wxEvent* Clone() const;
      
      /**
       * Get the factor to be applied to the scale value, relative to it's value
       * when scaling was started.
       */
      float GetTotalScale() const { return totalScale; }
    };
    typedef void (wxEvtHandler::*ChangedEventFunction)(ChangedEvent&);
  };

  #define EVT_SCALEDRAG_CHANGED(id, fn) \
      DECLARE_EVENT_TABLE_ENTRY( EVENT_SCALEDRAG_CHANGED, id, -1, \
        (wxObjectEventFunction) (wxEventFunction) (wxCommandEventFunction) \
          wxStaticCastEvent( DirScaleDragThing::ChangedEventFunction, & fn ), \
        (wxObject *)nullptr ),
} // namespace nutogui

#endif // __NUTOGUI_DIRSCALEDRAGTHING_H__
