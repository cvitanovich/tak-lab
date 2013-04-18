/*  Project mtlrch
    
    Copyright © 1994. All Rights Reserved.

    SUBSYSTEM:    mtlrch.exe Application
    FILE:         mtlrmdic.cpp
    AUTHOR:       


    OVERVIEW
    ========
    Source file for implementation of mtlrchMDIClient (TMDIClient).      
*/


#include <owl\owlpch.h>
#pragma hdrstop

#include <dir.h>

#include "mtlrchap.h"
#include "mtlrmdic.h"
#include "mtlrmdi1.h"
#include "apxprint.h"
#include "apxprev.h"


//{{mtlrchMDIClient Implementation}}


//
// Build a response table for all messages/commands handled
// by mtlrchMDIClient derived from TMDIClient.
//
DEFINE_RESPONSE_TABLE1(mtlrchMDIClient, TMDIClient)
//{{mtlrchMDIClientRSP_TBL_BEGIN}}
    EV_COMMAND(CM_FILEPRINT, CmFilePrint),
    EV_COMMAND(CM_FILEPRINTERSETUP, CmFilePrintSetup),
    EV_COMMAND(CM_FILEPRINTPREVIEW, CmFilePrintPreview),
    EV_COMMAND_ENABLE(CM_FILEPRINT, CmPrintEnable),
    EV_COMMAND_ENABLE(CM_FILEPRINTERSETUP, CmPrintEnable),
    EV_COMMAND_ENABLE(CM_FILEPRINTPREVIEW, CmPrintEnable),
    EV_WM_DROPFILES,
//{{mtlrchMDIClientRSP_TBL_END}}
END_RESPONSE_TABLE;


//////////////////////////////////////////////////////////
// mtlrchMDIClient
// ===========
// Construction/Destruction handling.
 mtlrchMDIClient::mtlrchMDIClient ()
 : TMDIClient ()
{
    ChildCount = 0;

    // INSERT>> Your constructor code here.

}


 mtlrchMDIClient::~mtlrchMDIClient ()
{
    Destroy();

    // INSERT>> Your destructor code here.

}


//////////////////////////////////////////////////////////
// mtlrchMDIClient
// ===========
// MDIClient site initialization.
void mtlrchMDIClient::SetupWindow ()
{
    // Default SetUpWindow processing.
    TMDIClient::SetupWindow ();

    // Accept files via drag/drop in the client window.
    DragAcceptFiles(TRUE);
}


//////////////////////////////////////////////////////////
// mtlrchMDIClient
// ==========
// Menu File Print command
void mtlrchMDIClient::CmFilePrint ()
{
    //
    // Create Printer object if not already created.
    // 
    mtlrchApp *theApp = TYPESAFE_DOWNCAST(GetApplication(), mtlrchApp);
    if (theApp) {
        if (!theApp->Printer)
            theApp->Printer = new TPrinter;

        //
        // Create Printout window and set characteristics.
        //
        APXPrintOut printout(theApp->Printer, Title, GetActiveMDIChild(), TRUE);

        theApp->Printing++;

        //
        // Bring up the Print dialog and print the document.
        //
        theApp->Printer->Print(GetActiveMDIChild()->GetClientWindow(), printout, TRUE);

        theApp->Printing--;
    }
}


//////////////////////////////////////////////////////////
// mtlrchMDIClient
// ==========
// Menu File Print Setup command
void mtlrchMDIClient::CmFilePrintSetup ()
{
    mtlrchApp *theApp = TYPESAFE_DOWNCAST(GetApplication(), mtlrchApp);
    if (theApp) {
        if (!theApp->Printer)
            theApp->Printer = new TPrinter;

        //
        // Bring up the Print Setup dialog.
        //
        theApp->Printer->Setup(this);
    }
}


//////////////////////////////////////////////////////////
// mtlrchMDIClient
// ==========
// Menu File Print Preview command
void mtlrchMDIClient::CmFilePrintPreview ()
{
    mtlrchApp *theApp = TYPESAFE_DOWNCAST(GetApplication(), mtlrchApp);
    if (theApp) {
        if (!theApp->Printer)
            theApp->Printer = new TPrinter;

        theApp->Printing++;

        PreviewWindow *prevW = new PreviewWindow(Parent, theApp->Printer, GetActiveMDIChild(), "Print Preview", new TLayoutWindow(0));
        prevW->Create();

        GetApplication()->BeginModal(GetApplication()->GetMainWindow());

        // We must destroy the preview window explicitly.  Otherwise, the window will not be destroyed until
        // it's parent the MainWindow is destroyed.
        prevW->Destroy();
        delete prevW;

        theApp->Printing--;
    }
}


//////////////////////////////////////////////////////////
// mtlrchMDIClient
// ==========
// Menu enabler used by Print, Print Setup and Print Preview.
void mtlrchMDIClient::CmPrintEnable (TCommandEnabler &tce)
{
    if (GetActiveMDIChild()) {
        mtlrchApp *theApp = TYPESAFE_DOWNCAST(GetApplication(), mtlrchApp);
        if (theApp) {
            // If we have a Printer already created just test if all is okay.
            // Otherwise create a Printer object and make sure the printer
            // really exists and then delete the Printer object.
            if (!theApp->Printer) {
                theApp->Printer = new TPrinter;
                
                tce.Enable(theApp->Printer->GetSetup().Error == 0);
            } else
                tce.Enable(theApp->Printer->GetSetup().Error == 0);
        }
    } else
        tce.Enable(FALSE);
}


void mtlrchMDIClient::EvDropFiles (TDropInfo)
{
    Parent->ForwardMessage();
}
