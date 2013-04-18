% m202x.m
% help file for m202x.dll
%
% returns -1 on error
% accepts six commands
% (defined in globals_mii.m as strutures C_ and M_)
%	C_.INIT
%
%	C_.MODE sets memory, channel and trig modes
%		2nd arg should be one of:
%		M_.MEMORY (followed by M_.MANUAL or M_.AUTO)
%		M_.INPUT (followed by M_.AUTOSTEP or M_.NOSTEP, then 0-7)
%		M_.TRIGGER (followed by M_.SOFTWARE, M_.CLOCK, or M_.EXTERNAL)
%
%	C_.START
%
%	C_.DATA clears or sets to read the data (does not actually read)
%		2nd arg either M_C.LEAR or M_.READ
%
%	C_.CLOCK sets period between conversions
%		2nd arg is period
%
%	C_.C8255 allows access to 8255 DON'T DO THIS!!!!
%		2nd arg must be M_READ!!!!!!

% details in MII Library functions
