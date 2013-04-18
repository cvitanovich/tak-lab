#if !defined(__mtlrmdi1_h)              // Sentry, use file only if it's not already included.
#define __mtlrmdi1_h

/*  Project mtlrch
    
    Copyright © 1994. All Rights Reserved.

    SUBSYSTEM:    mtlrch.exe Application
    FILE:         mtlrmdi1.h
    AUTHOR:       


    OVERVIEW
    ========
    Class definition for mtlrchMDIChild (TMDIChild).      
*/


#include <owl\owlpch.h>
#pragma hdrstop

#include <owl\editfile.h>
#include <owl\listbox.h>

#include "mtlrchap.rh"            // Definition of all resources.


//{{TMDIChild = mtlrchMDIChild}}
class mtlrchMDIChild : public TMDIChild {
public:
    mtlrchMDIChild (TMDIClient &parent, const char far *title, TWindow *clientWnd, BOOL shrinkToClient = FALSE, TModule* module = 0);
    virtual ~mtlrchMDIChild ();

//{{mtlrchMDIChildVIRTUAL_BEGIN}}
public:
    virtual void Paint (TDC& dc, BOOL erase, TRect& rect);
//{{mtlrchMDIChildVIRTUAL_END}}
//{{mtlrchMDIChildRSP_TBL_BEGIN}}
protected:
    void EvGetMinMaxInfo (MINMAXINFO far& minmaxinfo);
//{{mtlrchMDIChildRSP_TBL_END}}
DECLARE_RESPONSE_TABLE(mtlrchMDIChild);
};    //{{mtlrchMDIChild}}


#endif                                      // __mtlrmdi1_h sentry.
