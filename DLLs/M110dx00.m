% M110dx.m
% help file for M110dx.dll
%
% see also M110x.m
%
% returns -1 on error
% 10 main commands:
% (defined in globals_mii.m as strutures C_ and M_)
%	C_.INIT sets up 8255
%
%	C_.MODE sets to either M_.PST or M_.INTERVAL
%
%	C_.POLARITY sets polarities of BNC inputs
%		2nd arg either M_.EVENT or M_.STIMULUS
%		3rd arg either M_.POSITIVE or M_.NEGATIVE
%
%	C_.INTERUPT enables/disables interupt capability
%		2nd arg sets to interupt on either
%			M_.HALFFULL or M_.DATA_AVAIL
%		3rd arg specifies M_.ENABLE or M_.DISABLE
%
%	C_.START starts clock
%
%	C_.STOP stops clock and resets
%
%	C_.CLOCK sets clock speed in microsecs (1-4096)
%
%	C_.STATUS retruns header flags
%
%	C_.DATA reads memory buffer
%		2nd arg is integer number to read
%		3rd arg is pointer to where to put data (double)
%		returns # points read
%
%	C_.8255 allows direct access to 8255 DON'T DO THIS!
%		2nd arg must be M_.READ!!!

% details in MII Library functions
