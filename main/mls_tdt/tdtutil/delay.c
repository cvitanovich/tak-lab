#include <time.h>
#include "c:\matlabr11\extern\include\mex.h"


//void delay( clock_t wait );
/* Pauses for a specified number of milliseconds. */

void delay( clock_t wait )
{
   clock_t goal;
   goal = wait + clock();
   while( goal > clock() )
	   ;
}


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{ 
    int arg1; 
    /* Check for proper number of arguments */    
    if (nrhs != 1)
		mexErrMsgTxt("ONLY one input argument allowed (n_msecs)."); 

	arg1 = (int) mxGetScalar(prhs[0]);
   	delay(arg1);
	return;   
}