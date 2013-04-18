function [] = VCR(op)

%function [] = VCR(op)
%
% Helpfile for VCR.dll
% Allows control of the Panasonic AG-1960 VCR
% via the crimsonwave IR controller
% by placing calls to crimson.dll
%
% op:		1) PLAY
%		2) PAUSE/STILL
%		3) STOP
%		4) RECORD
%
% no other options are currently in place

% This DLL was written as ir2.cpp within
% d:\kip\code\c\crimson\samples\IRDevice
%
% considerable assistance was obtained from
% "Richard Wolf" <richard@graywolfsoftware.com>
