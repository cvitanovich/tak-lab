#if !defined(__mtlrmdic_h)              // Sentry, use file only if it's not already included.
#define __mtlrmdic_h

/*  Project mtlrch
    
    Copyright © 1994. All Rights Reserved.

    SUBSYSTEM:    mtlrch.exe Application
    FILE:         mtlrmdic.h
    AUTHOR:       


    OVERVIEW
    ========
    Class definition for mtlrchMDIClient (TMDIClient).      
*/


#include <owl\owlpch.h>
#pragma hdrstop

#include <owl\opensave.h>

#include "mtlrchap.rh"            // Definition of all resources.


//{{TMDIClient = mtlrchMDIClient}}
class mtlrchMDIClient : public TMDIClient {
public:
    int                     ChildCount;                 // Number of child window created. 

    mtlrchMDIClient ();
    virtual ~mtlrchMDIClient ();

    void OpenFile (const char *fileName = 0);

private:
    void LoadTextFile ();

//{{mtlrchMDIClientVIRTUAL_BEGIN}}
protected:
    virtual void SetupWindow ();
//{{mtlrchMDIClientVIRTUAL_END}}

//{{mtlrchMDIClientRSP_TBL_BEGIN}}
protected:
    void CmFilePrint ();
    void CmFilePrintSetup ();
    void CmFilePrintPreview ();
    void CmPrintEnable (TCommandEnabler &tce);
    void EvDropFiles (TDropInfo);
//{{mtlrchMDIClientRSP_TBL_END}}
DECLARE_RESPONSE_TABLE(mtlrchMDIClient);
};    //{{mtlrchMDIClient}}


#endif                                      // __mtlrmdic_h sentry.
