%M308x.m
%help file for M308x.dll
%
% returns -1 on error
%
% (defined in globals_mii.m as structures C_ and M_)
%
%	C_.INIT
%		M_.LEFT
%		or
%		M_.RIGHT
%
%	C_.DATA
%		M_.LEFT		or
%		M_.RIGHT
%			M_.SINE 
%				Freq (Hz)
%				SampPer (nominally 2 to 65000 usec)
%			or
%			M_.LOAD (loads array of values between 0 -> 4097)
%				MATLAB Array name
%				SampPer (nominally 2 to 65000 usec)
%
%	C_.CLOCK
%		M_.LEFT	
%		or
%		M_.RIGHT
%			M_.RAMP 
%				duration (in msec)
%			or
%			M_.FLAT
%				duration  (as (N+1)*100 usec)
%
%	C_.START (not implemented)
%
% example sequence sends ramped 100 msec sinewave of 1000 Hz:
%m308xx( C_.INIT, M_.LEFT );
%m308xx( C_.INIT, M_.RIGHT );
%m308xx( C_.DATA, M_.LEFT, M_.SINE, 1000, 20);
%m308xx( C_.DATA, M_.RIGHT, M_.SINE, 1000, 20);
%m308xx( C_.CLOCK, M_.LEFT, M_.RAMP, 5 );
%m308xx( C_.CLOCK, M_.RIGHT, M_.RAMP, 5 );
%m308xx( C_.CLOCK, M_.LEFT, M_.FLAT, 899 );
%m308xx( C_.CLOCK, M_.RIGHT, M_.FLAT, 899 );
