/*

  DOS2UNIX_FLOAT.C	Sample .MEX file corresponding to YPRIME.M
		Solves simple 3 body orbit problem 

  The calling syntax is:

			[unix_float] = dos2unix_float(dos_float)


  Marc Ullman   June 22, 1987
  Copyright (c) 1987 The Mathworks, Inc.
  All Righs Reserved
*/

#include <math.h>
#include "mex.h"

/* Input Arguments */

#define	DOS_FLOAT	prhs[0]


/* Output Arguments */

#define	UNIX_FLOAT	plhs[0]


#define	max(A, B)	((A) > (B) ? (A) : (B))
#define	min(A, B)	((A) < (B) ? (A) : (B))


#define pi 3.14159265

static
void dos2unix_float(double unix_float[], double dos_float[], 
                    unsigned int m, unsigned int n)
{
        int j,i;
        unsigned long cur_dos_float= 0;
        float cur_unix_float;
        char *cpD, *cpU;
        
        for (j=0; j<m*n; j++) {

           cpU = (char *) &cur_unix_float; 
           cpD = (char *) &cur_dos_float;
           cpD += 3;
           if (dos_float[j] > 4.294967296e+09)
              printf("dos2unix_float error:  Input number exceeds 32 bytes.");
           cur_dos_float = (unsigned long) dos_float[j];
           for (i = 0; i < 4; i++) {
              *cpU = *cpD;
              cpU++;
              cpD--;
           }  
           unix_float[j] = (double) cur_unix_float;
	} /* j */

	return;
}


mexFunction(int nlhs, Matrix *plhs[], int nrhs, Matrix *prhs[])
{
	double	*unix_float;
	double	*dos_float;
	unsigned int	m,n;

	/* Check for proper number of arguments */

	if (nrhs != 1) {
		mexErrMsgTxt("DOS2UNIX_FLOAT requires one input argument.");
	} else if (nlhs > 1) {
		mexErrMsgTxt("DOS2UNIX_FLOAT requires one output argument.");
	}


	/* Check the dimensions of dos_float.  Dos_Float can be 1 X 1. */

	m = mxGetM(DOS_FLOAT);
	n = mxGetN(DOS_FLOAT);
	if (!mxIsNumeric(DOS_FLOAT) || mxIsComplex(DOS_FLOAT) || 
		!mxIsFull(DOS_FLOAT)  || !mxIsDouble(DOS_FLOAT)) {
		mexErrMsgTxt("DOS2UNIX_FLOAT: error in input argument format.");
	}


	/* Create a matrix for the return argument */

	UNIX_FLOAT = mxCreateFull(m, n, REAL);


	/* Assign pointers to the various parameters */

	unix_float = mxGetPr(UNIX_FLOAT);

	dos_float = mxGetPr(DOS_FLOAT);


	/* Do the actual computations in a subroutine */
 	dos2unix_float(unix_float,dos_float,m,n);
	return;
}


