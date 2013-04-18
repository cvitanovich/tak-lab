% M210x.m
% help file for M210x.dll
%
% returns -1 on error
% accepts 6 commands
% (defined in globals_mii.m as strutures C_ and M_)
%
%	C_.INIT
%
%	C_.MODE sets operating mode
%		2nd arg one of: M_MEMORY (followed by M_.HOST or M_.EXTERNAL)
%			M_.BUFFER (followed by M_.CIRCULAR or M_.STOP)
%			or M_.ADDRESS sets number of values to read (3rd arg)
%
%	C_.DATA control data flow
%		2nd argone of: M_CLEAR (resets memory ptr to top)
%			M_.READ (reads value from buffer)
%			M_WRITE (writes value to buffer, value in arg 3)
%.
%	C_.INTERUPT followed by either M_.ENABLE or M_.DISABLE
%
%	C_.CONTROL allows access to control register
%		2nd arg either M_.SET (followed by value) or M_.READ
%
%	C_.C8255 allows access to 8255 DON'T DO THIS!!!!
%		2nd arg should be M_.READ!!!!

% details in MII Library functions
