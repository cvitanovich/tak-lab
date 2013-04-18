
 /*

  
  The calling syntax is:

			status = mtlwch(filename, channel, index)


 
  how to write C++ MEX code.  Changes to the C code are:

	* wrapped extern "C" around math.h and mex.h includes

	* removed #ifdef __STDC__ which accommodated both K&R and
	  ANSI code) - C++ code always has function prototches

	* wrapped mexFunction definition with extern "C" so that it
	  has "C linkage" for the benefit of MATLAB (within the body
	  of the mexFunction definition, full C++ syntax is allowed)

  To generate the MEX file, create the object file and then pass it
  through cmex:

	CC -o -I<directory containing mex.h> mtlrch.cc
	cmex chrime.o

  
*/

extern "C" {
#include <math.h>
#include "mex.h"
}

#include <fstream.h>

/* Input Arguments */

#define FILENAME        prhs[0]
#define CHANNEL prhs[1]
#define INDEX   prhs[2]

/* Output Arguments */



#define max(A, B)       ((A) > (B) ? (A) : (B))
#define min(A, B)       ((A) < (B) ? (A) : (B))



void swap_bytes(void *p,int size,int n)
//  void *p;   /* pointer to the first item to swap */
//  int size;  /* size of one item in bytes */
//  int n;     /* number of items */
{
  char *c0,*c1,d;
  int i;
  c0 = (char *)p;
  while(n-- > 0)
  {
    c1 = c0 + size -1;
    for(i = 0; i < size/2; i++)
      { d = *c0; *c0 = *c1; *c1 = d; c0++; c1--; }
    c0 += size/2 ;
  }
}





// matlab jumps here in

extern "C" {
void mexFunction(
	int             nlhs,
	Matrix  *plhs[],
	int             nrhs,
	Matrix  *prhs[]
	)
{
	double  *index;
	char    *filename;
	unsigned int    strlen,n,m;
	m=0;
	n=0;
	/* Check for proper number of arguments */

	if (nrhs != 3) {
		mexErrMsgTxt("mtlrch requires three input arguments.");
	} else if (nlhs > 1) {
		mexErrMsgTxt("mtlrch requires one output argument.");
	}


	
	if ( mxIsString(FILENAME)  != 1)
	 {
		mexErrMsgTxt("mtlrch requires that filename be string.");
	}
	    

	
	if (mxIsNumeric (INDEX) != 1  || mxGetM (INDEX)!=1 || mxGetN(INDEX)!=1)
	{
		mexErrMsgTxt("mtlrch requires that index be integer.");

	}
	
	if (mxIsNumeric (CHANNEL) != 1  )       
	{
		mexErrMsgTxt("mtlrch requires that index be integer.");

	}


	if (mxGetM (INDEX)!=1 && mxGetN(INDEX)!=1)
	{
		mexErrMsgTxt("mtlrch requires that one dimension is 1.");
	}

   
	strlen=mxGetN(FILENAME)+1;
	filename= new char[strlen];
	if (mxGetString(FILENAME,filename,strlen) < 0) 
		mexErrMsgTxt("Error using mxGetString.");


	index = mxGetPr (INDEX);

	/* Do the actual computations in a subroutine */

	int swapflag=0;
	unsigned char ant;
	short filetype, info_blocks, iscomplex, n_channels;
	int n_lines;
	double *pr;
	double *pi;


	fstream source (filename,ios::in|ios::out);
	if (!source) mexErrMsgTxt("unable to open file");

	source.seekg (79);
	source.read((unsigned char*) &ant,1);
	if (ant==26) swapflag=1;
	//      mexPrintf("ant %d \n", ant);

	source.seekg (80);
	source.read ((unsigned char *) &filetype,2);
	if (swapflag) swap_bytes(&filetype,2,1);
	//      mexPrintf("filetype %d \n", filetype);

	source.seekg (82);
	source.read ((unsigned char *) &info_blocks,2);
	if (swapflag) swap_bytes(&info_blocks,2,1);
	//      mexPrintf("info_blocks %d \n", info_blocks);

	source.seekg (84);
	source.read ((unsigned char *) &n_channels,2);
	if (swapflag) swap_bytes(&n_channels,2,1);
	//      mexPrintf("n_channels %d \n", n_channels);
	if (short(*index) > n_channels) mexErrMsgTxt("index to high");


	source.seekg (128);
	source.read ((unsigned char *) &n_lines,4);
	if (swapflag) swap_bytes(&n_lines,4,1);
	//      mexPrintf("n_lines %d \n", n_lines);

	switch (filetype)
	{
		case 2: iscomplex=0;
		break;
		case 3: iscomplex=1;
		break;
		default : mexErrMsgTxt("requested filetype is not supported\n");
	}
	//        mexPrintf("Test6\n");

	if (mxIsComplex(CHANNEL)!= iscomplex) mexErrMsgTxt("data type does not match file format.\n");

	pr=mxGetPr(CHANNEL);
	if (iscomplex) pi=mxGetPi(CHANNEL);

	source.seekp (info_blocks*256+n_lines*sizeof(float)*(iscomplex+1)*(int(*index)-1));

	float *ch_h= new float[n_lines*(iscomplex+1)];
	
	for (int i=0;i<n_lines; i++)
	{
		ch_h[i*(iscomplex+1)]=pr[i];
		//pr[i]=double (ch_h[i*(iscomplex+1)]);
		if (iscomplex) ch_h[i*(iscomplex+1)+1]=pi[i];
		//if (iscomplex) pi[i]=double (ch_h[i*(iscomplex+1)+1]);
	}       

	if (swapflag) swap_bytes(ch_h,4,n_lines*(iscomplex+1));

	source.write((unsigned char *) ch_h, n_lines*sizeof(float)*(iscomplex+1));

	return;
}
}
 /* end of extern "C" */

