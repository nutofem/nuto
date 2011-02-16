/**\file
 * Wrapper class to manage item client data in wxControlWithItems controls.
 */
/*
 * Written 2009 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#ifndef __UICOMMON_CONTROLWITHITEMSCLIENTDATAWRAPPER_H__
#define __UICOMMON_CONTROLWITHITEMSCLIENTDATAWRAPPER_H__

#include <wx/wx.h>

namespace uicommon
{
  /**
   * Wrapper class to manage item client data in wxControlWithItems controls.
   */
  template<typename T>
  class ControlWithItemsClientDataWrapper
  {
    wxControlWithItems* control;
    
    void* AllocClientData (const T& clientData)
    {
      void* newData;
      if (sizeof (T) > sizeof (void*))
      {
	newData = new T (clientData);
      }
      else
      {
	new (&newData) T (clientData);
      }
      return newData;
    }
    void FreeClientData (void* p)
    {
      if (sizeof (T) > sizeof (void*))
      {
	T* pT = reinterpret_cast<T*> (p);
	delete pT;
      }
      else
      {
	T* pT = reinterpret_cast<T*> (&p);
	pT->~T();
      }
    }
    T ConvertClientData (void* p) const
    {
      if (sizeof (T) > sizeof (void*))
      {
	T* pT = reinterpret_cast<T*> (p);
	return *pT;
      }
      else
      {
	T* pT = reinterpret_cast<T*> (&p);
	return *pT;
      }
    }
  public:
    ControlWithItemsClientDataWrapper () : control (0) {}
    ControlWithItemsClientDataWrapper (wxControlWithItems* control) : control (control) {}
    ~ControlWithItemsClientDataWrapper()
    {
      if (control != 0) Clear();
    }

    void SetControl (wxControlWithItems* control)
    {
      this->control = control;
    }
    wxControlWithItems* GetControl() const { return control; }
    bool HasControl() const { return control != 0; }

    void Clear()
    {
      for (unsigned int i = 0; i < control->GetCount(); i++)
      {
	FreeClientData (control->GetClientData (i));
      }
      control->Clear();
    }
    void Delete (unsigned int n)
    {
      FreeClientData (control->GetClientData (n));
      control->Delete (n);
    }

    int Append (const wxString& item, const T& clientData)
    {
      return control->Append (item, AllocClientData (clientData));
    }
    int Insert (const wxString& item, unsigned int pos, const T& clientData)
    {
      return control->Insert (item, pos, AllocClientData (clientData));
    }

    T GetClientData (unsigned int n) const
    {
      return ConvertClientData (control->GetClientData (n));
    }
    void SetClientData (unsigned int n, const T& data)
    {
      void* oldData = control->GetClientData (n);
      control->SetClientData (n, AllocClientData (data));
      FreeClientData (oldData);
    }
    bool HasClientData (unsigned int n) const
    {
      if (sizeof (T) > sizeof (void*))
      {
	return control->GetClientData (n) != 0;
      }
      else
	return true;
    }

    int GetSelection() const { return control->GetSelection(); }
    void SetSelection (int sel) { control->SetSelection (sel); }

    wxString GetString (int index) const { return control->GetString (index); }
    void SetString (int index, const wxString& str) { control->SetString (index, str); }
    wxString GetStringSelection () const { return control->GetStringSelection (); }

    unsigned int GetCount() const { return control->GetCount(); }
  };

} // namespace uicommon

#endif // __UICOMMON_CONTROLWITHITEMSCLIENTDATAWRAPPER_H__
