%M214x.m
%help file for M214x.dll
%
% returns -1 on error
%
% (defined in globals_mii.m as strutures C_ and M_)

%	C_.INIT initializes M214
%
%	C_.MODE six configurations:
%		M_.STOP
%		M_.GATED
%		M_.TRIGGERED
%		M_.TRIG_STOP
%		M_.CONVERT	(followed with either M_AUTO or M_MANUAL)
%		M_.SEQUENCE	(followed with either M_NORMAL or M_REVERSED)
%
%	C_.POLARITY
%		M_.PULSE
%		M_.TRIGOUT	(followed by either M_NEGATIVE or M_POSITIVE)
%
%	C_.INTERUPT
%		M_.ENABLE
%		M_.DISABLE
%
%	C_.CLOCK
%		M_.INTERVAL	(followed by value)
%		M_.DELAY		(followed by value)
%		M_.COUNT		(followed by value)
%
%	C_.CLEAR
%
%	C_.START
%
%	C_.STOP
%
%	C_.C8253 allows direct access to 8253 chip  DON'T DO THIS!!!!
