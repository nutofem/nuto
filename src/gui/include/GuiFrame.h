#ifndef __NUTOGUI_GUIFRAME_H__
#define __NUTOGUI_GUIFRAME_H__

#include <boost/shared_ptr.hpp>
#include <wx/string.h>

class wxAcceleratorEntry;
class wxWindow;

namespace nutogui
{
  /**
   * Interface to control a NuToGUI frame window.
   */
  struct GuiFrame
  {
    /**
     * Tab callback - given to tabs so they can control some aspects of the
     * tab in the frame.
     */
    struct TabCallback
    {
      virtual ~TabCallback() {}

      /// Set caption of the tab
      virtual void SetCaption (const wxString& caption) = 0;
      
      /// Set closeable state of tab
      virtual void SetCloseable (bool flag) = 0;
      
      /// Activate tab
      virtual void Activate () = 0;
    };
    typedef boost::shared_ptr<TabCallback> TabCallbackPtr;
    typedef boost::weak_ptr<TabCallback> TabCallbackWeakPtr;

    /**
     * Interface for a tab in the frame.
     * Implement this to show something in the GUI frame.
     */
    struct Tab
    {
      virtual ~Tab() {}

      /**
       * Create tab contents.
       * Called by frame.
       */
      virtual wxWindow* CreateContents (const TabCallbackWeakPtr& callback,
					wxWindow* parentWindow) = 0;
					
      /// Tab is being activated
      virtual void TabActivate () = 0;
      /// Tab is being deactivated
      virtual void TabDeactivate () = 0;
    };
    typedef boost::shared_ptr<Tab> TabPtr;
    
    /**
     * Pane callback - given to tabs so they can control some aspects of the
     * pane in the frame.
     */
    struct PaneCallback
    {
      virtual ~PaneCallback() {}

      /// Set caption of the pane
      virtual void SetCaption (const wxString& caption) = 0;
      
      /// Add a keyboard accelerator
      virtual void AddAccelerator (const wxAcceleratorEntry& accel) = 0;
    };
    typedef boost::shared_ptr<PaneCallback> PaneCallbackPtr;
    typedef boost::weak_ptr<PaneCallback> PaneCallbackWeakPtr;

    /**
     * Interface for a pane in the frame.
     * Unlike tabs, all panes are always visible. Thus, they should only be
     * small (toolbars etc).
     */
    struct Pane
    {
      /**
       * Create pane contents.
       * Called by frame.
       */
      virtual wxWindow* CreateContents (const PaneCallbackWeakPtr& callback,
					wxWindow* parentWindow) = 0;
				
      /// Return whether pane is actually a toolbar.
      virtual bool IsToolbar() const = 0;
    };
    typedef boost::shared_ptr<Pane> PanePtr;

    virtual ~GuiFrame() {}

    /// Add a tab to the frame
    virtual void AddTab (const TabPtr& tab, bool activate = false) = 0;
    /// Activate a tab
    virtual void ActivateTab (Tab* tab) = 0;
    
    /// Possible default locations for a pane
    enum DefaultLocation
    {
      /// Left side
      Left,
      /// Right side
      Right,
      /// Top edge
      Top,
      /// Bottom edge
      Bottom
    };
    /// Add a pane to the frame
    virtual void AddPane (const wxString& name,
			  DefaultLocation defaultLoc,
			  const PanePtr& pane) = 0;
  };
  typedef boost::shared_ptr<GuiFrame> GuiFramePtr;

} // namespace nutogui

#endif // __NUTOGUI_GUIFRAME_H__
