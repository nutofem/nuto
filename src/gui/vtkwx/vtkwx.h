/**\file
 * vtkwx: VTK rendering widget for WX
 */
/*
 * Written 2009 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#ifndef __VTKWX_H__
#define __VTKWX_H__

#include <vtkRenderWindowInteractor.h>

class vtkRenderWindow;

namespace vtkwx
{
  class Interactor;

  /**
   * A WX widget that encapsulates a VTK rendering window in a widget.
   */
  class RenderWidget : public wxWindow
  {
  public:
    RenderWidget (wxWindow *parent, wxWindowID id = -1,
      const wxPoint& pos = wxDefaultPosition,
      const wxSize& size = wxDefaultSize,
      long style = 0, const wxString& name = wxT("vtkwx__RenderWidget"));
    ~RenderWidget();

    /**
     * Get the VTK render window associated with this widget.
     * Creates the VTK render window, if necessary.
     * \remarks
     *  - You should not change the interactor of the window.
     */
    vtkRenderWindow* GetRenderWindow();
    /// Check whether the VTK render window exists yet.
    bool RenderWindowCreated() const { return renderWindow != 0; }

    /**
     * Get the 'ID' of this widget in a format that can be passed to the
     * VTK RendererWindow.
     */
    void* GetVTKWindowId () { return GetVTKWindowId (glContainer); }
    
    /// Enable/disable keyboard event handling. If enabled, keyboard events are sent to VTK.
    void EnableKeyboardHandling (bool flag) { passKeyboardEvents = flag; }
    /// Query state of keyboard event handling.
    bool IsKeyboardHandlingEnabled() const { return passKeyboardEvents; }
    
    DECLARE_EVENT_TABLE()
  protected:
    /// Need a redraw?
    bool needUpdate;
    /// VTK render window object
    vtkRenderWindow* renderWindow;
    /// Window that (supposedly) contains GL view
    wxWindow* glContainer;
    /// Our specific interactor
    Interactor* interactor;
    /// Window remap needed? (e.g. on GTK after reparent...)
    bool needRemap;
    /// Pass keyboard events to VTK?
    bool passKeyboardEvents;

    void OnPaint (wxPaintEvent& event);
    void OnSize (wxSizeEvent& event);
    void OnEraseBackground (wxEraseEvent& event);
    void OnIdle (wxIdleEvent& event);
    void OnMouse (wxMouseEvent& event);
    void OnKey (wxKeyEvent& event);

    /**
     * Cause the display to be rendered (potentially delayed to an opportune
     * time)
     */
    void UpdateDisplay();
    /**
     * Cause the display to be rendered immediately
     */
    void UpdateDisplayNow();

    /**
     * Get the 'ID' of a widget in a format that can be passed to the
     * VTK RenderWindow.
     */
    void* GetVTKWindowId (wxWindow* w);
    
    /// VTK RenderWindow as well as native window present?
    bool HaveWindow ();
    
  #if defined(__WXGTK__) && defined (__VTKWX_BUILD)
    /* __VTKWX_BUILD is a hack to avoid this methos to appear in
       headers where GTK isn't included */
    static gboolean signal_unrealize (GtkWidget *widget,
				      gpointer   user_data);
  #endif
  };

  /**
   * Custom render window interactor for use with RenderWidget.
   */
  class Interactor : public wxEvtHandler,
		     public vtkRenderWindowInteractor
  {
  public:
    static Interactor *New();

    DECLARE_EVENT_TABLE()
  protected:
    friend class RenderWidget;

    Interactor ();
    ~Interactor ();

    virtual int InternalCreateTimer (int timerId, int timerType, unsigned long duration);
    virtual int InternalDestroyTimer (int platformTimerId);

    void OnTimer(wxTimerEvent& event);

    /* In VTK, timers are identified by integers. WX timers are encapsulated
       in wxTimer objects. Thus we need to map between the VTK integer ID and
       the represented wxTimer object - the job of this helper class.
     */
    class TimerIdMapper
    {
      WX_DECLARE_HASH_MAP (int, wxTimer*, wxIntegerHash, wxIntegerEqual, IdToTimerHash);
      IdToTimerHash idToTimer;
    public:
      void SetIdForTimer (wxTimer* const timer, int id);
      wxTimer* const GetTimerForId (int id);

      void Delete (int id);
    };
    TimerIdMapper timerIDs;
  private:
    Interactor(const Interactor&);  // Not implemented.
    void operator=(const Interactor&);  // Not implemented.
  };

} // namespace vtkwx

#endif // __VTKWX_H__
