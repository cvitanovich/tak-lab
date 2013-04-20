#if !defined(__mtlrchad_h)              // Sentry, use file only if it's not already included.
#define __mtlrchad_h

/*  Project mtlrch
    
    Copyright © 1994. All Rights Reserved.

    SUBSYSTEM:    mtlrch.exe Application
    FILE:         mtlrchad.h
    AUTHOR:       


    OVERVIEW
    ========
    Class definition for mtlrchAboutDlg (TDialog).      
*/


#include <owl\owlpch.h>
#pragma hdrstop

#include "mtlrchap.rh"                  // Definition of all resources.


//{{TDialog = mtlrchAboutDlg}}
class mtlrchAboutDlg : public TDialog {
public:
    mtlrchAboutDlg (TWindow *parent, TResId resId = IDD_ABOUT, TModule *module = 0);
    virtual ~mtlrchAboutDlg ();

//{{mtlrchAboutDlgVIRTUAL_BEGIN}}
public:
    void SetupWindow ();
//{{mtlrchAboutDlgVIRTUAL_END}}
};    //{{mtlrchAboutDlg}}


// Reading the VERSIONINFO resource.
class ProjectRCVersion {
public:
    ProjectRCVersion (TModule *module);
    virtual ~ProjectRCVersion ();

    BOOL GetProductName (LPSTR &prodName);
    BOOL GetProductVersion (LPSTR &prodVersion);
    BOOL GetCopyright (LPSTR &copyright);
    BOOL GetDebug (LPSTR &debug);

protected:
    LPBYTE      TransBlock;
    void FAR    *FVData;

private:
    // Don't allow this object to be copied.
    ProjectRCVersion (const ProjectRCVersion &);
    ProjectRCVersion & operator =(const ProjectRCVersion &);
};


#endif                                      // __mtlrchad_h sentry.
