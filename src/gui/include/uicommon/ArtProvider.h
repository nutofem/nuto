/**\file
 * Custom WX art provider for art resources embedded in binary.
 */
/*
 * Written 2010 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#ifndef __UICOMMON_ARTPROVIDER_H__
#define __UICOMMON_ARTPROVIDER_H__

#include "export.h"

#include <wx/artprov.h>

namespace uicommon
{
  /**
   * Custom WX art provider for art resources embedded in binary.
   * Returns art provided with RegisterArt().
   */
  class UICOMMON_EXPORTED ArtProvider : public wxArtProvider
  {
    bool failEarly;
    
    static bool handlersRegistered;
    static void RegisterHandlers ();
  public:
    ArtProvider ();
    
    wxBitmap CreateBitmap (const wxArtID& id, 
			    const wxArtClient& client,
			    const wxSize& size);
    wxSize DoGetSizeHint(const wxArtClient& client);
    
    /**
     * Register an art resource.
     * \param name File name of art resource.
     * \param data Pointer to resource data.
     * \param size Size of resource data.
     * \param mimeType MIME type of resource data.
     */
    static void RegisterArt (const wxChar* name, const char* data, size_t size,
			     const wxChar* mimeType);
  };
} // namespace uicommon

#define ARTPROVIDER_REGISTER(Name, Data, Size, MIME)						\
  namespace {											\
    struct RegisterArt {									\
      RegisterArt() { uicommon::ArtProvider::RegisterArt (wxT (Name), Data, Size, wxT(MIME)); }	\
    };												\
    RegisterArt artRegistrator;									\
  }

#endif // __UICOMMON_ARTPROVIDER_H__
