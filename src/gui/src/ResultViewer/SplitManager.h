/**\file
 * Result viewer view split management.
 */
/*
 * Written 2010 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#ifndef __NUTOGUI_SPLITMANAGER_H__
#define __NUTOGUI_SPLITMANAGER_H__

#include "ResultViewerImpl.h"

#include <boost/weak_ptr.hpp>

namespace nutogui
{
  class ResultViewerImpl::SplitManager : public wxPanel
  {
    struct SplitNode;
    typedef boost::shared_ptr<SplitNode> SplitNodePtr;
    struct SplitNode
    {
      boost::weak_ptr<SplitNode> parent;
      /* Nodes: wxSplitterWindow
	 Leaves: View */
      wxWindow* window;
      SplitNodePtr children[2];
      
      SplitNode (const SplitNodePtr& parent, wxWindow* window)
       : parent (parent), window (window) {}
       
      bool IsLeaf() const { return !children[0] || !children[1]; }
      bool IsInnerNode() const { return children[0] && children[1]; }
    };
    SplitNodePtr splitTree;
    wxWindow* maximizedView;
    float oldSplitSashFrac;
    
    SplitNodePtr FindNode (const SplitNodePtr& tree, wxWindow* window);
    
    enum SplitDirection { splitV, splitH };
    void SplitCommon (wxWindow* view, SplitDirection direction,
		      wxWindow* newView);
    void RedistributeSize (const SplitNodePtr& node);
    std::pair<size_t, int> CountAndSizeInDirection (const SplitNodePtr& node,
						    int splitDir);
    int ComputeSizeInDirection (const SplitNodePtr& node, int splitDir, int viewSize);
    void SetSizeInDirection (const SplitNodePtr& node, int splitDir, int size);

    void PostToOthersUp (wxEvent& event, const SplitNodePtr& node);
    void PostToOthersDown (wxEvent& event, const SplitNodePtr& node);
  public:
    SplitManager (wxWindow* parent);
    
    void SetInitialView (wxWindow* firstView);
    
    bool CanSplit (wxWindow* view);
    void SplitVertically (wxWindow* view, wxWindow* newView);
    void SplitHorizontally (wxWindow* view, wxWindow* newView);
    bool CanUnsplit (wxWindow* view);
    void Unsplit (wxWindow* view);
    
    bool CanToggleMaximization (wxWindow* view);
    bool ToggleMaximization (wxWindow* view);
    
    void PostToOthers (wxEvent& event, wxWindow* source);
  };
} // namespace nutogui

#endif // __NUTOGUI_SPLITMANAGER_H__
