/*	m110dx.c  -  mex file to interface MATLAB with the m110d.
*	based on m110d.c, with the exception that DATA
*	returns only spike times (no stimulus times or event markers)
*/

#include <dos.h>
#include <stdarg.h>
#include "e:\kip\code\c\include\mii.h"
#include "c:\matlabR11\extern\include\mex.h"
#include <time.h>
#include <string.h>

#define SHALLOW_BASE        0x218
#define DEEP_BASE           0x228

#define CNTRL_REG(x)   ( x )
#define CLCK_REG(x)    ( (x) + 0x01 )
#define STAT_REG(x)    ( (x) + 0x02 )
#define CONTROL(x)     ( (x) + 0x03 )
#define READ(x)        ( (x) + 0x04 )

union M110_STATUS m110_status;

static union
{
   struct
   {
      unsigned mode        : 1;
      unsigned en_data     : 1;
      unsigned event_pol   : 1;
      unsigned stim_pol    : 1;
      unsigned reset_mem   : 1;
      unsigned x1          : 1;
      unsigned en_half_int : 1;
      unsigned en_int      : 1;
   } bits;
   
   unsigned char byte;
} s_control_flags, d_control_flags;

static unsigned int s_control_reg = 0x88;
static unsigned int d_control_reg = 0x88;


/* Pauses for a specified number of milliseconds. */
void delay( clock_t wait )
{
   clock_t goal;
   goal = wait + clock();
   while( goal > clock() )
	   ;
}


int m110d( int op, ... )
{
   va_list ap;
   int count, parm, arg1, cnt;
   unsigned long temp, temp1, temp2, eventflag;
   double *buf;

   switch( op )
   {
      case C_INIT:
         _outp( CONTROL( DEEP_BASE ), 0x88 );
         d_control_reg = 0x88;
         _outp( CNTRL_REG( DEEP_BASE ), 0x10 );
         _outp( CNTRL_REG( DEEP_BASE ), 0x01 );
         delay( 100 );
         _outp( CNTRL_REG( DEEP_BASE ), 0x11 );
         d_control_flags.byte = 0x11;
         break;
          
      case C_MODE:
         va_start( ap, op );
         d_control_flags.bits.mode = va_arg( ap, int );
         va_end( ap );
         break;
      
      case C_POLARITY:
         va_start( ap, op );
         
         switch( va_arg( ap, int ) )
         {
            case M_EVENT:
               d_control_flags.bits.event_pol = va_arg( ap, int );
               break;
            
            case M_STIMULUS:
               d_control_flags.bits.stim_pol = va_arg( ap, int );
               break;
            
            default:
               va_end( ap );
               return( -1 );
         }
         
         va_end( ap );
         break;
      
      case C_INTERRUPT:
         va_start( ap, op );
         
         switch( va_arg( ap, int ) )
         {
            case M_HALFFULL:
               d_control_flags.bits.en_half_int = va_arg( ap, int );
               break;
            
            case M_DATA_AVAIL:
               d_control_flags.bits.en_int = va_arg( ap, int );
               break;
            
            default:
               va_end( ap );
               return( -1 );
         }
         
         va_end( ap );
         break;
      
      case C_START:
         d_control_flags.bits.en_data = 1;
         d_control_flags.bits.reset_mem = 1;
         _outp( CNTRL_REG( DEEP_BASE ), d_control_flags.byte );
         break;
         
      case C_STOP:
         d_control_flags.bits.en_data = 0;
         d_control_flags.bits.reset_mem = 1;
         _outp( CNTRL_REG( DEEP_BASE ), d_control_flags.byte );
         break;
      
      case C_CLOCK:
         va_start( ap, op );
         parm = va_arg( ap, int );
         va_end( ap );
         
         if( parm < 1 )     parm = 1;
         if( parm > 4095 )  parm = 4095;
                              
         _outp( CLCK_REG( DEEP_BASE ), parm );
         _outp( STAT_REG( DEEP_BASE ), parm >> 8 );
         
         _outp( CNTRL_REG( DEEP_BASE ), 0x31 );
         _outp( CNTRL_REG( DEEP_BASE ), d_control_flags.byte );
         break;
          
      case C_STATUS:
         return( _inp( STAT_REG( DEEP_BASE ) ) );
          
      case C_DATA:
						/* third arg is a ptr!!!! */
         va_start( ap, op );
         parm = va_arg( ap, int );
         buf = va_arg( ap, double *);
         
         count = 0;
          
         m110_status.whole = m110d( C_STATUS );
		 temp1 = 0;			// starting time
         while( count < parm && m110_status.flags.d_avail )
		{
            temp = _inp( READ( DEEP_BASE ) ) & 0xFF;
            temp2 = temp;
            temp = _inp( READ( DEEP_BASE ) + 2 ) & 0xFF;
            temp2 |= temp << 8;
            temp = _inp( READ( DEEP_BASE ) ) & 0xFF;
            temp2 |= temp<< 16;
            temp = _inp( READ( DEEP_BASE ) + 2 );
			eventflag = temp & 0xFF;
			temp = temp & 0x7F;
			if (eventflag >= 128)
			{
				temp2 |= temp << 24;
				buf[count] = (double) (temp2 - temp1);
				count++;
            }
			else
			{
				temp1 = temp << 24;
				for (cnt = 0; cnt < count; cnt++)
					buf[cnt] = 0;
				count = 0;
			}
						
            m110_status.whole = m110d( C_STATUS );
         }
		 
	     return;
         return( count );
      
      case C_8255:             /* 8255 control reg.
                                  arg1 = M_READ or M_SET
                                  arg2 = val, if arg1 == M_SET
                               */
         va_start( ap, op );
         
         switch( va_arg( ap, int ) )
         {
            case M_SET:
               d_control_reg = va_arg( ap, int );
               _outp( CONTROL( SHALLOW_BASE ), d_control_reg );
               break;
         
            case M_READ:
               va_end( ap );
               return( d_control_reg );
            
            default:
               va_end( ap );
               return( -1 );
         }
         
         va_end( ap );
         break;
 
      default:
         return( -1 );
   }
     
   return( 0 );
}

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
     
{ 
    int op, argout; 
    int arg1, arg2;
	long arg3; 
	double *buf_ptr;

    /* Check for proper number of arguments */    
    if (nrhs < 1) { 
	mexErrMsgTxt("At least one input argument required."); 
    } else if (nrhs >4) {
	mexErrMsgTxt("Too many input arguments."); 
	}
    
	op = (int) mxGetScalar(prhs[0]);

        
    /* The MII subroutine */
	switch (nrhs)
	{
	case 1:
	/* create matrix for output argument */
		plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);	
		argout = m110d(op);
		break;
	case 2:
	    arg1 = (int) mxGetScalar(prhs[1]);
		if (op == 8)			//C_.DATA
		{
	/* create matrix for output argument */
			plhs[0] = mxCreateDoubleMatrix(arg1, 1, mxREAL);	
			buf_ptr = mxGetPr(plhs[0]);
			m110d(op,arg1,buf_ptr);
		}
		else
		{
	/* create matrix for output argument */
			plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);	
			argout = m110d(op,arg1);
		}
		break;
	case 3:
	/* create matrix for output argument */
		plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);	
	    arg1 = (int) mxGetScalar(prhs[1]);
	    arg2 = (int) mxGetScalar(prhs[2]);
		argout = m110d(op,arg1,arg2);
		break;
	case 4:
		plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);	
	    arg1 = (int) mxGetScalar(prhs[1]);
	    arg2 = (int) mxGetScalar(prhs[2]);
	    arg3 = (long) mxGetScalar(prhs[3]);
		argout = m110d(op,arg1,arg2,arg3);
		break;
	}
	
	return;   
}