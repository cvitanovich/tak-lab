% M101xl.m
% help file for M101xl.dll
% Three main commands.
% (defined in globals_mii.m as strutures C_ and M_)
% returns -1 on error
%
% C_.INIT sets all bits as output
%
% C_.DATA is the main command
%		second arg specifies performance on byte (M_BYTE) or bit (M_BIT)
%		third arg specifies:
%			M_.READ input data from port
%			M_.SET output data to port
%			M_.CLEAR output all 0s to port
%			M_.PULSE pulses data to port then resets to original state
%				note if in M_.BIT mode this sets bit high then low
%
% details in MII Library functions
%
% modified to allow an additional (last) argument that sets the 
% approx pulse width in microsecs (accuracy may vary with machine speed)
%
% e.g. m100x( C_.INIT );
%		 m101xl( C_.DATA,M_.BIT,M_.PULSE,port,Nusecs)
