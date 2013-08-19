% M100x.m
% help file for M100x.dll
% routine accepts any of nine commands 
% (defined in globals_mii.m as strutures C_ and M_)
% a value of -1 is returned in the event of an error
%
% C_.C8253	allows direct access to control register for 8253
%		second arg should be value to write to register
%		usually used with C_CLOCK, M_RAWSET arguments
%
% C_.C8255 allows direct access to control register for 8255
%		second arg must be either M_.SET or M_.READ
%			if M_.SET: third arg is a new value to write
%
% C_.CLOCK reads or sets one of the timers
%		second argument spcifies exact operation to perform:
%			M_.SET sets a timer to mode 4(count in fourth arg)
%			M_.RAWSET sets a timer (w/o mode) (count in fourth arg)
%			M_.READ reads a timer
%		third arg denotes which timer (e.g. M_.TIMER0)
%
% C_.DATA reads or writes to the I/O port
%		second arg must be either M_.SET (third arg is new value)
%		or M_.READ
%
% C_.INIT initializes the M100
%
% C_.INTERUPT controls the timer 0 interupts:
%		M_.ENABLE enables timer 0 interupt
%		M_.DISABLE disables the timer 0 interupt
%		M_.CLEAR clears all interupts
%
% C_.START starts a clock
%		second arg should be e.g. M_.TIMER0
%
% C_.STATUS sets or reads the interupt flag
%		second arg is either M_.SET (third arg is value) or M_.READ
%
% C_.STOP stops a clock
%		second arg should be e.g. M_.TIMER0
%
% details in MII Library functions
