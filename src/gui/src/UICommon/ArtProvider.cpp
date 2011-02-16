/**\file
 * Custom WX art provider for art resources embedded in binary.
 */
/*
 * Written 2010 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#include "common.h"

#include <wx/wx.h>
#include <wx/filesys.h>
#include <wx/fs_mem.h>

#include "uicommon/ArtProvider.h"

namespace uicommon
{
  static const wxChar pathPrefix[] = wxT("/ArtProvider/");
  
  bool ArtProvider::handlersRegistered = false;
  
  void ArtProvider::RegisterHandlers ()
  {
    if (!handlersRegistered)
    {
      handlersRegistered = true;
      
      // PNG handler is used for image resources we embed
      wxImage::AddHandler (new wxPNGHandler);

      // Add file system handler used for embedding resources
      wxFileSystem::AddHandler(new wxMemoryFSHandler);
    }
  }
  
  ArtProvider::ArtProvider ()
    : failEarly (false)
  {
  }
  
  wxBitmap ArtProvider::CreateBitmap (const wxArtID& id, 
				      const wxArtClient& client,
				      const wxSize& size)
  {
    if (failEarly) return wxNullBitmap;
    
  #if defined(__WXGTK__)
    // GTK: prefer system art
    {
      failEarly = true;
      wxBitmap defaultArt (wxArtProvider::GetBitmap (id, client, size));
      failEarly = false;
      if (defaultArt.IsOk()) return defaultArt;
    }
  #endif
    
    wxString artName (id);
    if (artName.StartsWith (wxT ("wxART_"))) artName.erase (0, 6);
    artName.MakeLower ();
    
    wxString suffix;
    if (client == wxART_MENU)
      suffix = wxT("_menu");
    
    wxFileSystem fs;
    wxFSFile* artFile = 0;
    wxString memoryFSPrefix (wxT ("memory:"));
    memoryFSPrefix.Append (pathPrefix);
    if (!suffix.IsEmpty())
    {
      wxString artPath (memoryFSPrefix);
      artPath.Append (artName);
      artPath.Append (suffix);
      artPath.Append (wxT (".png"));
      artFile = fs.OpenFile (artPath);
    }
    if (!artFile)
    {
      wxString artPath (memoryFSPrefix);
      artPath.Append (artName);
      artPath.Append (wxT (".png"));
      artFile = fs.OpenFile (artPath);
    }
    
    if (!artFile)
      return wxNullBitmap;
    
    wxImage artImage (*(artFile->GetStream()));
    wxBitmap artBitmap (artImage);
    
    delete artFile;
    
    return artBitmap;
  }
  
  wxSize ArtProvider::DoGetSizeHint (const wxArtClient& client)
  {
    // Override toolbar, menu sizes with custom dimensions
    // @@@ FIXME: Maybe not that good for very high resolutions...
    if (client == wxART_TOOLBAR)
      return wxSize (24, 24);
    else if (client == wxART_MENU)
      return wxSize (16, 16);
    else
      return wxArtProvider::DoGetSizeHint (client);
  }
  
  void ArtProvider::RegisterArt (const wxChar* name, const char* data, size_t size,
				 const wxChar* mimeType)
  {
    RegisterHandlers();
    
    wxString path (pathPrefix);
    path.Append (name);
    wxMemoryFSHandler::AddFileWithMimeType (path, data, size, mimeType);
  }
  
} // namespace uicommon
