#if !defined(__mtlrchap_h)              // Sentry, use file only if it's not already included.
#define __mtlrchap_h

/*  Project mtlrch
    
    Copyright © 1994. All Rights Reserved.

    SUBSYSTEM:    mtlrch.exe Application
    FILE:         mtlrchap.h
    AUTHOR:       


    OVERVIEW
    ========
    Class definition for mtlrchApp (TApplication).      
*/


#include <owl\owlpch.h>
#pragma hdrstop

#include <owl\statusba.h>
#include <owl\controlb.h>
#include <owl\buttonga.h>
#include <owl\editview.h>
#include <owl\listview.h>
#include <owl\docmanag.h>
#include <owl\filedoc.h>
#include <owl\printer.h>

#include <classlib\bags.h>


#include "mtlrmdic.h"

#include "mtlrchap.rh"            // Definition of all resources.


// TFileDrop class Maintains information about a dropped file, its name, where it was dropped,
// and whether or not it was in the client area
class TFileDrop {
public:
    operator == (const TFileDrop& other) const {return this == &other;}

    char*   FileName;
    TPoint  Point;
    BOOL    InClientArea;

    TFileDrop (char*, TPoint&, BOOL, TModule*);
    ~TFileDrop ();

    const char* WhoAmI ();
private:
    //
    // hidden to prevent accidental copying or assignment
    //
    TFileDrop (const TFileDrop&);
    TFileDrop & operator = (const TFileDrop&);
};

typedef TIBagAsVector<TFileDrop> TFileList;
typedef TIBagAsVectorIterator<TFileDrop> TFileListIter;


//{{TApplication = mtlrchApp}}
class mtlrchApp : public TApplication {
private:

private:
    void SetupSpeedBar (TDecoratedMDIFrame *frame);
    void AddFiles (TFileList* files);

public:
    mtlrchApp ();
    virtual ~mtlrchApp ();

    mtlrchMDIClient  *mdiClient;

    // Public data members used by the print menu commands and Paint routine in MDIChild.
    TPrinter        *Printer;                           // Printer support.
    int             Printing;                           // Printing in progress.

//{{mtlrchAppVIRTUAL_BEGIN}}
public:
    virtual void InitMainWindow();
    virtual void InitInstance();
//{{mtlrchAppVIRTUAL_END}}

//{{mtlrchAppRSP_TBL_BEGIN}}
protected:
    void EvNewView (TView& view);
    void EvCloseView (TView& view);
    void CmHelpAbout ();
    void EvDropFiles (TDropInfo drop);
    void EvWinIniChange (char far* section);
//{{mtlrchAppRSP_TBL_END}}
  DECLARE_RESPONSE_TABLE(mtlrchApp);
};    //{{mtlrchApp}}


#endif                                      // __mtlrchap_h sentry.
