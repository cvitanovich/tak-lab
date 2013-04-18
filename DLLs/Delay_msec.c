#include <quicktime.h>
#include "c:\matlabr11\extern\include\mex.h"
#include <stdlib.h>
#include <stdio.h>

//void delay( clock_t wait );
/* Pauses for a specified number of milliseconds. */

void delay( clock_t wait )
{
   char *ticvar;
   clock_t goal;
   int i;

   printf("%d\n",clock());
   goal = wait + clock();
   while( goal > clock());
   printf("%d\n",clock());


   //printf("Start timing now \n\n");
   //for (i = 0; i < 100; i++)
//	printf("%d \n",clock());
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