/**\file
 * vtkwx: VTK rendering widget for WX
 */
/*
 * Written 2009 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#include <wx/wx.h>

#if defined(__WXGTK__)
#include <gtk/gtkwidget.h>
#include <gdk/gdkx.h>
#endif

#define __VTKWX_BUILD
#include "vtkwx.h"

#include <vtkCommand.h>
#include <vtkRenderWindow.h>

namespace vtkwx
{
  #if defined(__WXGTK__)
  static GdkWindow* GetGdkWindow (wxWindow* w)
  {
    // On GTK GetHandle() returns a GtkWidget
    GtkWidget* gwidget = static_cast<GtkWidget*> (w->GetHandle());
    GdkWindow* gwindow;
  #ifdef HAS_gtk_widget_get_window
    // gtk_widget_get_window() only exists in newer GTK+ version
    gwindow = gtk_widget_get_window (gwidget);
  #else
    // For older versions, directly access the member
    gwindow = gwidget->window;
  #endif
    return gwindow;
  }
  #endif

  BEGIN_EVENT_TABLE(RenderWidget, wxWindow)
    EVT_PAINT(RenderWidget::OnPaint)
    EVT_SIZE(RenderWidget::OnSize)
    EVT_ERASE_BACKGROUND(RenderWidget::OnEraseBackground)
    EVT_IDLE(RenderWidget::OnIdle)
    EVT_MOUSE_EVENTS(RenderWidget::OnMouse)
    EVT_KEY_DOWN(RenderWidget::OnKey)
    EVT_KEY_UP(RenderWidget::OnKey)
    EVT_CHAR(RenderWidget::OnKey)
  END_EVENT_TABLE()

  RenderWidget::RenderWidget (wxWindow *parent, wxWindowID id,
    const wxPoint& pos,
    const wxSize& size,
    long style, const wxString& name)
   : wxWindow (parent, id, pos, size, style | wxWANTS_CHARS, name),
     needUpdate (true), renderWindow (0), needRemap (false)
  {
  #if defined(__WXGTK__) && (wxABI_VERSION < 20900)
    /* @@@ Odd (and I don't understand it yet) what's happening on GTK:
	When hijacking the GL container for the VTK render window all painting
        seems to go to the _parent_ window, not the 'intended' window!
       */
    // 2010-11-09: doesn't work (indeed, breaks things) on WX 2.9
    glContainer = new wxPanel (this, wxID_ANY, wxPoint (-10, -10), wxSize (0, 0));
  #else
    glContainer = this;
  #endif
  
  #if defined(__WXGTK__) 
    GtkWidget* gwidget = static_cast<GtkWidget*> (glContainer->GetHandle());
    g_signal_connect (gwidget, "unrealize",
		      G_CALLBACK (signal_unrealize),
		      this);
  #endif
    
    // Use a special interactor that wraps WX timers
    interactor = Interactor::New ();
  }

  RenderWidget::~RenderWidget()
  {
    if (renderWindow != 0) renderWindow->Delete();
    renderWindow = 0;
    interactor->Delete();
  }

  vtkRenderWindow* RenderWidget::GetRenderWindow()
  {
    // If no render window was created yet, create one
    if ((renderWindow == 0) || needRemap)
    {
      int w, h;
      GetSize (&w, &h);
      // Create & set size to current widget dimensions
      if (renderWindow == 0)
	renderWindow = vtkRenderWindow::New ();

      renderWindow->SetParentId (GetVTKWindowId (GetParent()));
      renderWindow->SetNextWindowId (GetVTKWindowId ());
  #if defined(__WXGTK__)
      /* Seems that occasionally a BadWindow error appears if we don't sync
         before letting VTK remap the window */
      GdkWindow* gwindow = GetGdkWindow (this);
      gdk_display_sync (gdk_drawable_get_display  (reinterpret_cast<GdkDrawable*> (gwindow)));
  #endif
      renderWindow->WindowRemap ();
      needRemap = false;

      renderWindow->SetSize (w, h);
      // Set up our own interactor.
      renderWindow->SetInteractor (interactor);
      if (!interactor->GetInitialized ())
	interactor->Initialize();
    }
    return renderWindow;
  }

  void RenderWidget::OnPaint (wxPaintEvent& event)
  {
    // Needed by WX
    wxPaintDC dc (this);

    UpdateDisplay();
  }

  void RenderWidget::OnSize (wxSizeEvent& event)
  {
    // Filter 0 dimensions, X11 doesn't like these
    if ((event.GetSize().x == 0) || (event.GetSize().y == 0)) return;

    if (HaveWindow ())
    {
      // Update render window and interactor sizes.
      renderWindow->SetSize (event.GetSize().x, event.GetSize().y);
      renderWindow->GetInteractor()->SetSize (event.GetSize().x, event.GetSize().y);

      UpdateDisplay();
    }
  }

  void RenderWidget::OnEraseBackground (wxEraseEvent& event)
  {
    if (!renderWindow) event.Skip (true);
  }

  void RenderWidget::OnIdle (wxIdleEvent& event)
  {
    /* Oddity #2: When calling Render() directly in reaction to changes
       to the window (resize, paint), the contents seem to 'lag' one frame.
       Delaying the rendering to the next idle event fixes this. */
    if (!renderWindow || !needUpdate) return;
    UpdateDisplayNow();
  }

  static wxWindow* GetTopLevelParent (wxWindow* w)
  {
    while (w)
    {
      if (w->IsTopLevel())
	break;
      w = w->GetParent();
    }
    return w;
  }

  void RenderWidget::OnMouse (wxMouseEvent& event)
  {
    if (renderWindow == 0) return;

    /* To make keyboard input work we need to grab focus.
       But that can be annoying, especially when the focused top-level
       window changes. To reduce annoyance, only change focus if we
       have the same top-level parent as the currently focused window. */
    if (passKeyboardEvents
	&& (GetTopLevelParent (FindFocus ()) == GetTopLevelParent (this)))
    {
      SetFocus ();
    }

    // Forward mouse events to VTK

    if (event.IsButton())
    {
      // Mouse button event
      interactor->SetEventInformationFlipY (event.GetX(), event.GetY(),
	event.ControlDown(), event.ShiftDown(), 0, event.ButtonDClick());
      interactor->SetAltKey (event.AltDown());
      if (event.ButtonDown())
      {
	switch (event.GetButton())
	{
	  case wxMOUSE_BTN_LEFT:
	    interactor->InvokeEvent(vtkCommand::LeftButtonPressEvent, 0);
	    break;
	  case wxMOUSE_BTN_RIGHT:
	    interactor->InvokeEvent(vtkCommand::RightButtonPressEvent, 0);
	    break;
	  case wxMOUSE_BTN_MIDDLE:
	    interactor->InvokeEvent(vtkCommand::MiddleButtonPressEvent, 0);
	    break;
	}
      }
      else
      {
	switch (event.GetButton())
	{
	  case wxMOUSE_BTN_LEFT:
	    interactor->InvokeEvent(vtkCommand::LeftButtonReleaseEvent, 0);
	    break;
	  case wxMOUSE_BTN_RIGHT:
	    interactor->InvokeEvent(vtkCommand::RightButtonReleaseEvent, 0);
	    break;
	  case wxMOUSE_BTN_MIDDLE:
	    interactor->InvokeEvent(vtkCommand::MiddleButtonReleaseEvent, 0);
	    break;
	}
      }
    }
    else if (event.GetEventType() == wxEVT_MOTION)
    {
      // Mouse move event
      interactor->SetEventInformationFlipY (event.GetX(), event.GetY(),
	event.ControlDown(), event.ShiftDown());
      interactor->SetAltKey (event.AltDown());
      interactor->InvokeEvent(vtkCommand::MouseMoveEvent, 0);
    }
    else if (event.GetEventType() == wxEVT_MOUSEWHEEL)
    {
      // Mouse wheel event
      if (event.GetWheelRotation() > 0) // FIXME ignores delta
	interactor->InvokeEvent (vtkCommand::MouseWheelForwardEvent, 0);
      else
	interactor->InvokeEvent (vtkCommand::MouseWheelBackwardEvent, 0);
    }
    else if (event.Entering())
    {
      // Mouse enters widget
      interactor->InvokeEvent (vtkCommand::EnterEvent, 0);
    }
    else if (event.Leaving())
    {
      // Mouse leaves widget
      interactor->InvokeEvent (vtkCommand::LeaveEvent, 0);
    }
  }
  
  static const char* WXKeyToKeySym (int wxkey)
  {
    // Lifted from vtkWin32RenderWindowInteractor
    static const char * const AsciiToKeySymTable[] = {
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      "space", "exclam", "quotedbl", "numbersign",
      "dollar", "percent", "ampersand", "quoteright",
      "parenleft", "parenright", "asterisk", "plus",
      "comma", "minus", "period", "slash",
      "0", "1", "2", "3", "4", "5", "6", "7",
      "8", "9", "colon", "semicolon", "less", "equal", "greater", "question",
      "at", "A", "B", "C", "D", "E", "F", "G",
      "H", "I", "J", "K", "L", "M", "N", "O",
      "P", "Q", "R", "S", "T", "U", "V", "W",
      "X", "Y", "Z", "bracketleft",
      "backslash", "bracketright", "asciicircum", "underscore",
      "quoteleft", "a", "b", "c", "d", "e", "f", "g",
      "h", "i", "j", "k", "l", "m", "n", "o",
      "p", "q", "r", "s", "t", "u", "v", "w",
      "x", "y", "z", "braceleft", "bar", "braceright", "asciitilde", "Delete"
    };
    static const int numAsciiToKeySyms = sizeof(AsciiToKeySymTable)/sizeof(const char*);

    const char* sym = 0;
    if (wxkey < numAsciiToKeySyms) sym = AsciiToKeySymTable[wxkey];
    
    if (sym == 0)
    {
      switch (wxkey)
      {
#define WX_TO_SYM(wxKey, keysym)  case WXK_ ## wxKey: return #keysym
	
	WX_TO_SYM(BACK, BackSpace);
	WX_TO_SYM(TAB, Tab);
	WX_TO_SYM(RETURN, Return);
	WX_TO_SYM(ESCAPE, Escape);
	WX_TO_SYM(DELETE, Delete);

	WX_TO_SYM(CANCEL, Cancel);
	WX_TO_SYM(CLEAR, Clear);
	WX_TO_SYM(SHIFT, Shift_L);
	WX_TO_SYM(ALT, Alt_L);
	WX_TO_SYM(CONTROL, Control_L);
	//WXK_MENU
	WX_TO_SYM(PAUSE, Pause);
	WX_TO_SYM(CAPITAL, Caps_Lock);
	WX_TO_SYM(END, End);
	WX_TO_SYM(HOME, Begin);
	WX_TO_SYM(LEFT, Left);
        WX_TO_SYM(UP, Up);
	WX_TO_SYM(RIGHT, Right);
	WX_TO_SYM(DOWN, Down);
	WX_TO_SYM(SELECT, Select);
	WX_TO_SYM(PRINT, Print);
	WX_TO_SYM(EXECUTE, Execute);
	//WXK_SNAPSHOT
	WX_TO_SYM(INSERT, Insert);
	WX_TO_SYM(HELP, Help);
	WX_TO_SYM(NUMPAD0, KP_0);
	WX_TO_SYM(NUMPAD1, KP_1);
	WX_TO_SYM(NUMPAD2, KP_2);
	WX_TO_SYM(NUMPAD3, KP_3);
	WX_TO_SYM(NUMPAD4, KP_4);
	WX_TO_SYM(NUMPAD5, KP_5);
	WX_TO_SYM(NUMPAD6, KP_6);
	WX_TO_SYM(NUMPAD7, KP_7);
	WX_TO_SYM(NUMPAD8, KP_8);
	WX_TO_SYM(NUMPAD9, KP_9);
	WX_TO_SYM(MULTIPLY, asterisk);
	WX_TO_SYM(ADD, plus);
	WX_TO_SYM(SEPARATOR, comma);
	WX_TO_SYM(SUBTRACT, minus);
	WX_TO_SYM(DECIMAL, period);
	WX_TO_SYM(DIVIDE, slash);
	WX_TO_SYM(F1, F1);
	WX_TO_SYM(F2, F2);
	WX_TO_SYM(F3, F3);
	WX_TO_SYM(F4, F4);
	WX_TO_SYM(F5, F5);
	WX_TO_SYM(F6, F6);
	WX_TO_SYM(F7, F7);
	WX_TO_SYM(F8, F8);
	WX_TO_SYM(F9, F9);
	WX_TO_SYM(F10, F10);
	WX_TO_SYM(F11, F11);
	WX_TO_SYM(F12, F12);
	WX_TO_SYM(F13, F13);
	WX_TO_SYM(F14, F14);
	WX_TO_SYM(F15, F15);
	WX_TO_SYM(F16, F16);
	WX_TO_SYM(F17, F17);
	WX_TO_SYM(F18, F18);
	WX_TO_SYM(F19, F19);
	WX_TO_SYM(F20, F20);
	WX_TO_SYM(F21, F21);
	WX_TO_SYM(F22, F22);
	WX_TO_SYM(F23, F23);
	WX_TO_SYM(F24, F24);
	WX_TO_SYM(NUMLOCK, Num_Lock);
	WX_TO_SYM(SCROLL, Scroll_Lock);
	WX_TO_SYM(PAGEUP, Page_Up);
	WX_TO_SYM(PAGEDOWN, Page_Down);

	WX_TO_SYM(NUMPAD_SPACE, KP_Space);
	WX_TO_SYM(NUMPAD_TAB, KP_Tab);
	WX_TO_SYM(NUMPAD_ENTER, KP_Enter);
	WX_TO_SYM(NUMPAD_F1, KP_F1);
	WX_TO_SYM(NUMPAD_F2, KP_F2);
	WX_TO_SYM(NUMPAD_F3, KP_F3);
	WX_TO_SYM(NUMPAD_F4, KP_F4);
	WX_TO_SYM(NUMPAD_HOME, KP_Home);
	WX_TO_SYM(NUMPAD_LEFT, KP_Left);
	WX_TO_SYM(NUMPAD_UP, KP_Up);
	WX_TO_SYM(NUMPAD_RIGHT, KP_Right);
	WX_TO_SYM(NUMPAD_DOWN, KP_Down);
	WX_TO_SYM(NUMPAD_PAGEUP, KP_Page_Up);
	WX_TO_SYM(NUMPAD_PAGEDOWN, KP_Page_Down);
	WX_TO_SYM(NUMPAD_END, KP_End);
	WX_TO_SYM(NUMPAD_BEGIN, KP_Begin);
	WX_TO_SYM(NUMPAD_INSERT, KP_Insert);
	WX_TO_SYM(NUMPAD_DELETE, KP_Delete);
	WX_TO_SYM(NUMPAD_EQUAL, KP_Equal);
	WX_TO_SYM(NUMPAD_MULTIPLY, KP_Multiply);
	WX_TO_SYM(NUMPAD_ADD, KP_Add);
	WX_TO_SYM(NUMPAD_SEPARATOR, KP_Separator);
	WX_TO_SYM(NUMPAD_SUBTRACT, KP_Subtract);
	WX_TO_SYM(NUMPAD_DECIMAL, KP_Decimal);
	WX_TO_SYM(NUMPAD_DIVIDE, KP_Divide);
    
	WX_TO_SYM(WINDOWS_LEFT, Win_L);
	WX_TO_SYM(WINDOWS_RIGHT, Win_R);
	WX_TO_SYM(MENU, App);
     	
      #undef WX_TO_SYM
      }
    }
    
    return sym;
  }
  
  void RenderWidget::OnKey (wxKeyEvent& event)
  {
    if (renderWindow == 0) return;

    if (passKeyboardEvents)
    {
      int keyChar;
      if (event.GetEventType() == wxEVT_CHAR)
      {
      #if defined(wxHAS_UNICODE) && wxHAS_UNICODE
	keyChar = event.GetUnicodeKey();
      #else
	keyChar = event.GetKeyCode();
	if (keyChar >= 128) keyChar = 0;
      #endif
      }
      else
	keyChar = 0;
      
      interactor->SetEventInformation (event.GetX(), event.GetY(), 
				      event.ControlDown(), event.ShiftDown(),
				      keyChar, 1,
				      WXKeyToKeySym (event.GetKeyCode()));
      interactor->SetAltKey (event.AltDown());
      if (event.GetEventType() == wxEVT_KEY_DOWN)
      {
	interactor->InvokeEvent (vtkCommand::KeyPressEvent, 0);
      }
      else if (event.GetEventType() == wxEVT_KEY_UP)
      {
	interactor->InvokeEvent (vtkCommand::KeyReleaseEvent, 0);
      }
      else if (event.GetEventType() == wxEVT_CHAR)
      {
	interactor->InvokeEvent (vtkCommand::CharEvent, 0);
      }
    }
    event.Skip();
  }

  void RenderWidget::UpdateDisplay()
  {
    needUpdate = true;
  }
  
  void RenderWidget::UpdateDisplayNow()
  {
    GetRenderWindow();
    renderWindow->Render();
    needUpdate = false;
  }

  void* RenderWidget::GetVTKWindowId (wxWindow* w)
  {
  #if defined(__WXGTK__)
    GdkWindow* gwindow = GetGdkWindow (w);
    if (!gwindow) return 0;
    // Finally, VTK wants an X11 'Window' value for window IDs
    Window xWindow = GDK_WINDOW_XWINDOW (gwindow);
    return reinterpret_cast<void*> (xWindow);
  #elif defined(__WXMSW__)
    // On Win32 GetHandle() returns a HWND, same that VTK wants
    return static_cast<void*> (w->GetHandle());
  #else
  #error Window ID translation not supported for this platform!
  #endif
  }

  bool RenderWidget::HaveWindow ()
  {
    return (renderWindow != 0) && (glContainer->GetHandle () != 0)
      && (GetVTKWindowId() != 0);
  }

#if defined(__WXGTK__) 
  gboolean RenderWidget::signal_unrealize (GtkWidget *widget,
					     gpointer   user_data)
  {
    RenderWidget* this_ = reinterpret_cast<RenderWidget*> (user_data);
    this_->needRemap = true;
    if (this_->renderWindow)
      this_->renderWindow->Finalize ();
    
    return false;
  }
#endif

  //-------------------------------------------------------------------------

  BEGIN_EVENT_TABLE(Interactor, wxEvtHandler)
    EVT_TIMER(wxID_ANY, Interactor::OnTimer)
  END_EVENT_TABLE()

  Interactor* Interactor::New ()
  {
    return new Interactor;
  }

  Interactor::Interactor ()
  {
  }
  
  Interactor::~Interactor()
  {
  }

  int Interactor::InternalCreateTimer (int timerId, int timerType, unsigned long duration)
  {
    // Creates a new wxTimer
    wxTimer* newTimer = new wxTimer (this, timerId);
    if (!newTimer->Start (duration, timerType == OneShotTimer))
    {
      delete newTimer;
      return 0;
    }
    // Store the wxTimer object for each timer ID
    timerIDs.SetIdForTimer (newTimer, timerId);
    return timerId;
  }

  int Interactor::InternalDestroyTimer (int platformTimerId)
  {
    // Get the wxTimer object for the timer ID
    wxTimer* timer = timerIDs.GetTimerForId (platformTimerId);
    if (timer == 0) return 0;
    // Release timer ID and wxTimer object
    timerIDs.Delete (platformTimerId);
    delete timer;
    return 1;
  }

  void Interactor::OnTimer(wxTimerEvent& event)
  {
    // Forward timer event to WX.
    // TODO Check if one-shot timers are destroyed by VTK or if we have to
    // do that manually.
    int timerId = event.GetId();
    InvokeEvent(vtkCommand::TimerEvent, (void*)&timerId);
  }

  void Interactor::TimerIdMapper::SetIdForTimer (wxTimer* const timer, int id)
  {
    idToTimer[id] = timer;
  }

  wxTimer* const Interactor::TimerIdMapper::GetTimerForId (int id)
  {
    IdToTimerHash::iterator timer = idToTimer.find (id);
    if (timer == idToTimer.end()) return 0;
    return timer->second;
  }

  void Interactor::TimerIdMapper::Delete (int id)
  {
    idToTimer.erase (id);
  }

} // namespace vtkwx
